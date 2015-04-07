# RASQUAL
RASQUAL (Robust Allele Specific Quantification and quality controL) maps QTLs using population and allele-specific signals for any sequenced based cellular tratis.

## How to build & install

Please make sure CLAPACK is installed in your environment.  Otherwise you can get the library at http://www.netlib.org/clapack/.  To build and install RASQUAL, firstly go to the _source_ directory (src), then set environment variables appropriately to point the CLAPACK library as follows.  Finally use "make" to build and install RASQUAL which will be installed in "$RASQUALDIR/bin".

	RASQUALDIR=/path/to/rasqualdir/
	cd $RASQUALDIR/src
	export CFLAGS="-I/nfs/users/nfs_n/nk5/team170/Natsuhiko/CLAPACK-3.1.1.1/INCLUDE -I/nfs/users/nfs_n/nk5/team170/Natsuhiko/CLAPACK-3.1.1.1/F2CLIBS"
	export LDFLAGS="-L/nfs/users/nfs_n/nk5/team170/Natsuhiko/CLAPACK-3.1.1.1 -L/nfs/users/nfs_n/nk5/team170/Natsuhiko/CLAPACK-3.1.1.1/F2CLIBS"
	make
	make install

## First eQTL mapping with RASQUAL

Mapping QTLs with RASQUAL is very straightforward.  You just need to prepare two data files (1) a fragment (read) count table and (2) a VCF file with imputed SNP genotypes and allele-specific counts at all SNP loci as standard input for RASQUAL.  RASQUAL also requires feature start(s) and end(s) (i.e., union of exons in this example) as inputs (-s and -e respectively).  To save memory usage, you are also required to count the number of testing SNPs (-l) and feature SNPs (-m) a priori.  Here are a couple of commands to map eQTLs (C11orf21 and TSPAN32 genes) to get the lead QTL SNP for each gene: 

        tabix data/chr11.gz 11 | bin/rasqual -y data/Y.bin -k data/K.bin -n 24 -j 1 -l 409 -m 63 \
                   -s 2316875,2320655,2321750,2321914,2324112 -e 2319151,2320937,2321843,2323290,2324279 -z -t -f C11orf21
        tabix data/chr11.gz 11 | bin/rasqual -y data/Y.bin -k data/K.bin -n 24 -j 2 -l 409 -m 61 \
                   -s 2323227,2323938,2324640,2325337,2328175,2329966,2330551,2331219,2334884,2335715,2338574,2339093 \
                   -e 2323452,2324188,2324711,2325434,2328220,2330040,2330740,2331248,2334985,2337897,2338755,2339430 \
                   -z -t -f TSPAN32

## Count data preparation

You can find an example expression data for C11orf21 and TSPAN32 genes in the _data_ directory.  There are two files: read count table (Y.txt) and sample specific offset data (K.txt), both has to be organized feature by sample format (i.e., row: gene; col: sample).  The following R script allows you to convert the count and offset files (text) into binary format used in RASQUAL software.  The script essentially converts a table data into a vector of values (double precision).

	cd $RASQUALDIR
	export RHOME=/software/R-3.0.0/bin/
	$RHOME/R --vanilla --quiet --args data/Y.txt data/K.txt < R/txt2bin.R > log

You also need to prepare custom VCF files containing the allele specific counts at all SNPs.  The files have to contain an additional subfield, "AS", located within the genotype field consisting of two integers, the reference and alternative allele counts, separated by a comma.  For example, sample **_i_** is heterozygous at a SNP and has 1 and 10 reads overlapping with the reference and alternative alleles respectively, the genotype field for the sample becomes

	... FORMAT ... Sample_i ...
	... GL:AS  ... 0|1:1,10 ...

You can find the example VCF (chr11.gz) in the _data_ directory.

## Genotype uncertainty

To maximize the ability of RASQUAL, we recommend to incorporate uncertainty in imputed genotypes.  There are 4 options: 

1. **Allelic probability** (AP) 

    Allelic probability can be obtained from the standard 2-step imputation scheme, where genotypes are first phased then imputed on haplotype-by-haplotype basis.  The custom AP field consists of the allelic probabilities (in Log 10 scale) of two haplotypes from one individual, separated by a comma:

        ... FORMAT ... Sample_i      ...
        ... AP:AS  ... 0.0,-5.0:1,10 ...

2. **Genotype likelihood** (GL)

    You may also use genotype likelihood from conventional genotype imputation in conjunction with phased genotype data:

        ... FORMAT    ... Sample_i                ...
        ... GT:GL:AS  ... 0|1:-4.0,-0.1,-0.6:1,10 ...

3. **Dosage** (DS)

    Genotype dosage can also be utilized to take account of genotype uncertainty:

        ... FORMAT    ... Sample_i     ...
        ... GT:DS:AS  ... 0|1:1.1:1,10 ...

4. **Imputation quality scoare** (R squre value; RSQ)

    Imputation methods often provid a quality score for each SNP locus that approximates squared correlation between true and observed genotypes (*e.g.*, *R*^2 from MaCH or Beagle; *I*^2 from IMPUTE2 ).  RASQUAL can convert the score into genotyping error rate to handle uncertainly:

        ... INFO            FORMAT ... Sample_i ...
        ... ...;RSQ=0.9;... GL:AS  ... 0|1:1,10 ...

We strongly recommend to use AP for QTL mapping.  If there are multiple subfield of AP, GL, DS, AP is prioritized than GL and DS and GL is prioritized than DS.  If you want to prioritize RSQ, you need to specify (-z) option for RASQUAL (see above example).

