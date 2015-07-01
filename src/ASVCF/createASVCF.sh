#!/bin/bash

BAMLIST=$1
VCF=$2

ID=`cat /dev/urandom | tr -dc '[:alnum:]' | head -c 10`

TMP=`tempfile`
PASTE=pasteAS.tmp.sh
ALLASC=${VCF%.gz}.as.gz
NEWVCF=${VCF%.gz}.new.gz

# file search
for I in `cat $BAMLIST`
do
	if [ -f ${I%.bam}.$ID.gz ]
	then
		echo ${I%.bam}.$ID.gz exists.
		echo Need to move it first.
		exit
	fi
done
if [ -f $ALLASC ]
then
	echo $ALLASC exists.
	exit
fi
if [ -f $NEWVCF ]
then
	echo $NEWVCF exists.
	exit
fi
if [ -f $PASTE ]
then
	echo $PASTE exists.
	exit
fi



# count AS reads
for CHR in `tabix -l $VCF`
do
	tabix $VCF $CHR | cut -f 1-9 | gzip > $TMP.gz
	START=`zcat $TMP.gz | head -n 1 | cut -f 2`
	END=`zcat $TMP.gz | tail -n 1 | cut -f 2`
	echo $CHR:$START-$END
	for I in `cat $BAMLIST`
	do
		#echo $I
		samtools view -F 0x0100 $I $CHR:$START-$END | awk -v RASQUALDIR=$RASQUALDIR '$7=="="{cmd = RASQUALDIR"/src/ASVCF/parseCigar "$6; cmd | getline N; print $3"\t"$4"\t"$4+N-1"\t"$6"\t"$10; close(cmd);}' | \
			$RASQUALDIR/src/ASVCF/countAS $TMP.gz | awk '{print $5","$6}' | gzip >> ${I%.bam}.$ID.gz
	done
done
rm $TMP.gz


# paste all AS counts in one file
echo -n $RASQUALDIR/src/ASVCF/zpaste > $PASTE
for I in `cat $BAMLIST`
do
	echo -n " "${I%.bam}.$ID.gz >> $PASTE
done
sh $PASTE > $ALLASC
# clean up
rm $PASTE
for I in `cat $BAMLIST`
do
	rm ${I%.bam}.$ID.gz
done



# merge AS counts and VCF
$RASQUALDIR/src/ASVCF/pasteFiles $VCF $ALLASC | bgzip > $NEWVCF
# clean up
rm $ALLASC 


