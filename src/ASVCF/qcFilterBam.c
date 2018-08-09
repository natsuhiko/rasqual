//
//  qcFilterBam.c
//  

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "qcFilterBam.h"



int bestHitsN,gapOpensN,mismatchesN;
//unsigned short 
int minQual,minInsert,maxMismatch,maxGapOpen,maxBestHit,singleEnd=0,bowtie=1;//,hashSize,sortedByName=0,fragment=0;
//unsigned long 
int maxInsert;
char *oFormat=NULL, *skipMissing=NULL;

char* tag=NULL;
char* val=NULL;
char* strBuf=NULL;


void usage()
/* Explain usage and exit. */
{
    fprintf(stderr,
             "qcFilterBam - Fast reading and filtering of reads directly from bam file\n"
             "usage:\n"
             "   qcFilterBam [OPTIONS] <stdin> \n"
             "options:\n"
             "   -minQual XXX - minimum quality score (default 0)\n"
             "   -minInsert XXX - minimum insertSize score (default 0)\n"
             "   -maxInsert XXX - maximum insertSize score (default 1e7)\n"
             "   -maxMismatch XXX - maximum no. mismatches (default 3)\n"
             "   -maxGapOpen XXX - maximum no. gap opens (default 0)\n"
             "   -maxBestHit XXX - maximum no. best hits in genome (default 1)\n"
             "   -skipMissing (T|F) - Retain reads with missing tags (default is T - TRUE: skip reads with missing tags)\n"
             "   -singleEnd - Single end reads, ignores insert, read orientation and same chromosome flags\n"
             "   -bowtie - Alignments by bowtie, skips checks for X0 field\n"
             );
}

char* num=NULL;
char* op=NULL;

void parseTrueFalse (char *skipMissingArg) {
    if(!strcmp(skipMissingArg,"t"))
        strcpy(skipMissingArg,"T");
    if(!strcmp(skipMissingArg,"true"))
        strcpy(skipMissingArg,"T");
    if(!strcmp(skipMissingArg,"TRUE"))
        strcpy(skipMissingArg,"T");
    if(!strcmp(skipMissingArg,"f"))
        strcpy(skipMissingArg,"F");
    if(!strcmp(skipMissingArg,"false"))
        strcpy(skipMissingArg,"F");
    if(!strcmp(skipMissingArg,"FALSE"))
        strcpy(skipMissingArg,"F");
}

void checkClArgs() {
    parseTrueFalse(skipMissing);
    if((strcmp(skipMissing,"T"))&&(strcmp(skipMissing,"F"))){
        fprintf(stderr, "Invalid argument for option -skipMissing (must be t|f,T|F,true|false,TRUE|FALSE)");
        return;
    }
}


void copyRow(char *newRow[MAXROW], char *row[MAXROW], int wc) {
    int i;
    for(i=0;i<wc;i++) {
        if(strlen(row[i])>MAXFIELD){
	    row[i][MAXFIELD-1]='\0';
	}
        strcpy(newRow[i],row[i]);
    }
}

void printRow(char *newRow[MAXROW], int wc) {
    int i;
    for(i=0;i<wc;i++) {
        printf("%s\t",newRow[i]);
    }
    printf("\n");
}

int getLenFromCigarDG13( char *cigarString ) {
    int i,len=0;
    for(i=0;i<strlen(cigarString);i++) {
    }
    return len;
}

int getLenFromCigar(char* c1){
    int itr=0;
    //char op;
    long len;
    int offset=0;
    //char num[100];
    long res=0;
    int nchar;
    while((sscanf(c1+offset, "%ld%[MIDNSH=PX]", &len, op))>0){
        nchar = sprintf(num, "%ld", len);
        offset += nchar+1;
        itr++;
        if(op[0]=='D' || op[0]=='N' || op[0]=='P' || op[0]=='M' || op[0]=='=' || op[0]=='X'){res += len;};
    }
    return (int)res;
}

void trunkCigarSeq(int start, char* c1, char* seq1, int b){
    //char* op; // operation
    //op = (char*)calloc(100, sizeof(char));
    int len;// length in cigar
    int nseq=0;
    int nc1=0;
    //pc1 = c1;
    int end;
    int nchar;
    //fprintf(stderr, "%ld\n", at[0]);
    while((sscanf(c1+nc1, "%d%[MIDNSH=PX]", &len, op))>0){
        if(op[0]=='M' || op[0]=='=' || op[0]=='X'){
            end = start+len-1;
            if(start<b && b<=end){
                nchar = sprintf(c1+nc1, "%d", b-start);
                c1[nc1+nchar]=op[0];
                c1[nc1+nchar+1]='\0';
                seq1[nseq+b-start]='\0';
            }else if(start==b){
                c1[nc1]='\0';
                seq1[nseq]='\0';
            }
            start += len;
            nseq += len;
        }else if(op[0]=='D' || op[0]=='N' || op[0]=='P'){
            end = start+len-1;
            if(start<b && b<=end){
                nchar = sprintf(c1+nc1, "%d", b-start);
                c1[nc1+nchar]=op[0];
                c1[nc1+nchar+1]='\0';
                seq1[nseq]='\0';
            }else if(start==b){
                c1[nc1]='\0';
                seq1[nseq]='\0';
            }
            start += len;
        }else{// S, H and I
            nseq += len;
        }
        nchar = sprintf(num, "%d", len);
        nc1 += nchar+1;                       
    }
}

int findLeftRight(char *newRow[MAXROW], char *row[MAXROW], long unsigned int coords[2], char** cigar, char** seqs) {
    //int genomicLen=0;
    int iStart = atoi(newRow[3]);
    int jStart = atoi(row[3]);
    int isize = abs(atoi(newRow[8]));
    int jsize = abs(atoi(row[8]));
    int gap=0;
    if(iStart < jStart) {
        coords[0] = iStart;
        coords[1] = iStart + isize -1;
        cigar[0] = newRow[5];
        cigar[1] = row[5];
        seqs[0] = newRow[9];
        seqs[1] = row[9];
        gap=jStart-iStart-getLenFromCigar(cigar[0]);
	//HiC
	//if(isize==0){coords[1] = jStart + getLenFromCigar(cigar[1]) - 1;}
        if(gap<0){trunkCigarSeq(iStart, cigar[0], seqs[0], jStart);}
    } else {
        coords[0] = jStart;
        coords[1] = jStart + jsize -1;
        cigar[0] = row[5];
        cigar[1] = newRow[5];
        seqs[0] = row[9];
        seqs[1] = newRow[9];
        gap=iStart-jStart-getLenFromCigar(cigar[0]);
	//HiC
	//if(isize==0){coords[1] = iStart + getLenFromCigar(cigar[1]) - 1;}
        if(gap<0){trunkCigarSeq(jStart, cigar[0], seqs[0], iStart);}
    }
    return gap;
}

void getFieldTagAndValue(char *f, char *tag, char *val) {
  int i;
  char *ch = NULL, delim[] = ":";
  
  strcpy(strBuf,f);
  ch = strtok(strBuf,delim);
  strcpy(tag,ch);
  for(i=0;i<2;i++) {
    ch = strtok(NULL,delim);
  }
  strcpy(val,ch);
}


int checkAndConvertVal (char *tag, char *val) {
  int count = atoi(val);
  if(count>MAXCOUNT)
    fprintf(stderr, "Field value for tag %s found above default value (Max value %d, found %d)",tag,MAXCOUNT,count);
  return count;
}





int passesFilters(char *row[MAXROW], int wc) {
    
    short int bitFlag = atoi(row[1]); /* Bitwise flag */
    int Q = atoi(row[4]); /* Quality score */
    
    unsigned long int qst = strtoul(row[3],NULL,0); /* Query read start position */
    unsigned long int mst = strtoul(row[7],NULL,0); /* Mate read start position */
    
    //int qst = atoi(row[3]);
    //int mst = atoi(row[7]);
    
    //unsigned long int insertSize = strtoul(row[8],NULL,0);
    int insertSize = abs(atoi(row[8]));
    
    char qstr = (bitFlag & 16) == 0 ? '+' : '-'; /* Query strand */
    char mstr = (bitFlag & 32) == 0 ? '+' : '-'; /* Mate strand */
    //char pri = (bitFlag & 256) == 0 ? 'T' : 'F'; /* Mate strand */
    
    // proper pair
    //char pppr = (bitFlag & 2) == 0 ? '0' : '1';
    //if(pppr=='0'){return 0;}
    
    //pcr duplicates
    //char pcrd = (bitFlag & 1024) == 0 ? '0' : '1';
    //if(pcrd=='1'){return 0;}
    
    
    bestHitsN = mismatchesN = gapOpensN = MAXCOUNT; /* Global variables, set so that most of the time will exclude if not found */
    int i,found=0;
    for(i=11;i<wc;i++) {
        getFieldTagAndValue(row[i],tag,val);
        if(!strcmp(tag,"X0")) {
            bestHitsN = checkAndConvertVal(tag,val);
            found++;
        }
        if(!strcmp(tag,"XM")) {
            mismatchesN = checkAndConvertVal(tag,val);
            found++;
        }
        if(!strcmp(tag,"XO")) {
            gapOpensN = checkAndConvertVal(tag,val);
            found++;
        }
        if(found==3){break;}
    }
    if(!strcmp(skipMissing,"F")) { /* To include reads with missing tags we set their missing values to the thresholds, so read will pass QC */
        bestHitsN   = (bestHitsN == MAXCOUNT)   ? maxBestHit  : bestHitsN;
        gapOpensN   = (gapOpensN == MAXCOUNT)   ? maxGapOpen  : gapOpensN;
        mismatchesN = (mismatchesN == MAXCOUNT) ? maxMismatch : mismatchesN;
    }
    if(bowtie==1) { /* # best hits field is missing in bowtie alignments */
        bestHitsN = maxBestHit;
    } 
//fprintf(stderr, "%d ins=%d %d %d %d %c %c %ld %ld\n", Q, insertSize, mismatchesN, gapOpensN, bestHitsN, qstr, mstr, qst, mst);
    if(singleEnd==0){
        if(
	    // HiC 
	//(Q>=minQual) //&&((qstr != mstr)&&(((qstr == '+')&&(qst <= mst))||((mstr == '+')&&(mst <= qst))))){
            (Q>=minQual)&&
               (insertSize>=minInsert)&&
               (insertSize<=maxInsert)&&
               (strcmp(row[6],"=")==0)&&
               //(qstr != mstr)//&&
               ((qstr != mstr)&&(((qstr == '+')&&(qst <= mst))||((mstr == '+')&&(mst <= qst))))&&
               (mismatchesN<=maxMismatch)&&
               (gapOpensN<=maxGapOpen)&&
               (bestHitsN<=maxBestHit)

	    ){
            return 1;
        }
    } else {
        if((Q>=minQual)&&(mismatchesN<=maxMismatch)&&(gapOpensN<=maxGapOpen)&&(bestHitsN<=maxBestHit))
            return 1;
    }

    return 0;
}

int Min(int a, int b){return a >= b ? b : a; }

char* buffer;
int i0=0;
int ret=0;
int myLineFileChop(char** row){
	int wc=0, i=0, l=0;
	//int returnFlag=0;
	for(i=i0; i<ret; i++){
		if(buffer[i]=='\n'){
			if(wc<MAXROW)row[wc][l]='\0';
			wc++;
			i0=i+1;
			return Min(wc, MAXROW);
		}else if(buffer[i]=='\t'){
			if(wc<MAXROW)row[wc][l]='\0';
			l=0;
			wc++;
		}else{
			if(l<MAXFIELD-1 && wc<MAXROW)row[wc][l++]=buffer[i];
		}
	}
	while((ret=fread(buffer,sizeof(char),MAXFIELD*100,stdin))>0){
		buffer[ret]='\0';
		for(i=0;i<ret;i++){
			if(buffer[i]=='\n'){
				if(wc<MAXROW)row[wc][l]='\0';
                        	wc++;
                        	i0=i+1;
                        	return Min(wc, MAXROW);
                	}else if(buffer[i]=='\t'){
                        	if(wc<MAXROW)row[wc][l]='\0';
                        	l=0;
                        	wc++;
                	}else{
                        	if(l<MAXFIELD-1 && wc<MAXROW)row[wc][l++]=buffer[i];
                	}
		}
	}
	return wc;
}

void qcFilterBam(char *inFile) {
    /* qcFilterBam - Fast reading and filtering of reads directly from bam file. */
    unsigned int wc;
    long unsigned int goodReads=0,nReads=0;//,nHashKeys=0,missingTags=0
    char *row[MAXROW];
    int i;
    for(i=0;i<MAXROW;i++) row[i] = (char*)calloc(MAXFIELD, sizeof(char));
    char* cigar[2];
    char* seqs[2];
//fprintf(stderr, "%s\n", inFile);
    //struct lineFile *lf = lineFileOpen(inFile, TRUE);
    if(!singleEnd) {// paired end
        char mmstring='N';
        long ln=0;
        long unsigned int coords[2];
        char *prevRow[MAXROW];
        for(i=0;i<MAXROW;i++) prevRow[i] = (char*)calloc(MAXFIELD, sizeof(char));
        int pass1, pass2;
        //while ((wc = lineFileChop(lf, row)) != 0) {
        while ((wc = myLineFileChop(row)) != 0) {
//fprintf(stderr, "wc=%d\n", wc);
//fprintf(stderr, "%s\n", row[0]);
//fprintf(stderr, ">>>>%s\n\n", buffer+i0);
            pass1 = passesFilters(row, wc);
            ln++;
//fprintf(stderr, "%ld %s\n", ln, row[0]);
            nReads++;
            if(ln>1) {
                if(strcmp(prevRow[0],row[0])==0 && -atoi(prevRow[8])==atoi(row[8])) {// ID & ins size comparison
                    if((pass1+pass2)>0) {
                        /* printRow(row,wc); */
                        /* printRow(prevRow,wc); */
                        int gap = findLeftRight(row,prevRow,coords,cigar,seqs);
                        goodReads+= pass1 + pass2;
                        if(gap<=0){
                            printf("%s\t%ld\t%ld\t%s%s\t%s%s\n",row[2],coords[0],coords[1],cigar[0],cigar[1],seqs[0],seqs[1]);  
                        }else{
                            printf("%s\t%ld\t%ld\t%s%d%c%s\t%s%s\n",row[2],coords[0],coords[1],cigar[0],gap,mmstring,cigar[1],seqs[0],seqs[1]);
                        }
                    }
                } else {
                    pass2=pass1;
                    copyRow(prevRow,row,wc);
                }
            } else {
                pass2=pass1;
                copyRow(prevRow,row,wc);
            }
        }
        for(i=0;i<MAXROW;i++) free(prevRow[i]);
    } else if(singleEnd) {
        fprintf(stderr,"Entering single end mode\n");
        //while ((wc = lineFileChop(lf, row)) != 0) {
        while ((wc = myLineFileChop(row)) != 0) {
	    nReads++;
            int pass1 = passesFilters(row,wc);
            if(pass1==1){
                goodReads++;
                long coord = (long)atoi(row[3]);
                int chr_offset=0;
                if(row[2][0]=='c'){chr_offset=3;}
		int bitFlag = atoi(row[1]);
		char strnd = (bitFlag & 16) == 0 ? '+' : '-';
                //printf("%s\t%ld\t%ld\t%s\t%s\t%c\t%s\n",row[2]+chr_offset,coord,coord+getLenFromCigar(row[5])-1,row[5],row[9],strnd,row[0]);
                printf("%s\t%ld\t%ld\t%s\t%s\n",row[2]+chr_offset,coord,coord+getLenFromCigar(row[5])-1,row[5],row[9]);
            }
        }
    }
    //lineFileClose(&lf);
    fprintf(stderr,"\n%ld reads found\n",nReads);
    //fprintf(stderr,"%ld reads missing tags\n",missingTags);
    fprintf(stderr,"%ld reads removed\n",nReads - goodReads);
    fprintf(stderr,"%ld reads remaining\n\n",goodReads);
}

int main(int argc, char *argv[]) {
    int i;
    buffer = (char*)calloc(MAXFIELD*100+1,sizeof(char));
    op = (char*)calloc(100,sizeof(char));
    num = (char*)calloc(100,sizeof(char));
    tag = (char *) calloc(MAXSTR,sizeof(char));
    val = (char *) calloc(MAXSTR,sizeof(char));
    strBuf = (char *) calloc(MAXSTR,sizeof(char));
    if(op==NULL || num==NULL || tag==NULL || val==NULL || strBuf==NULL){fprintf(stderr, "memory allocation failure!\n"); return -1;}
    /* Process command line. */
    //optionInit(&argc, argv, options);
    if (argc == 1){ usage(); return -1;}
    //readSize = optionInt("readSize",75);
    //minQual      = optionInt("minQual",0);
    //minInsert    = optionInt("minInsert",-1);
    //maxInsert    = optionInt("maxInsert",270000000);
    //maxMismatch  = optionInt("maxMismatch",1000);
    //maxGapOpen   = optionInt("maxGapOpen",1000);
    //maxBestHit   = optionInt("maxBestHit",1000);
    //oFormat    = optionVal("oFormat","sam");
    //hashSize   = optionInt("hashSize",25);
    //skipMissing  = optionVal("skipMissing","F");
    //singleEnd    = optionExists("singleEnd");
    //bowtie       = optionExists("bowtie");
    minQual      = 0;
    minInsert    = -1;
    maxInsert    = 270000000;
    maxMismatch  = 1000;
    maxGapOpen   = 1000;
    maxBestHit   = 1000;
    singleEnd    = 0;
    bowtie       = 1;
    char* skipMissing0; skipMissing0  = calloc(10, sizeof(char));
    skipMissing0[0] = 'F';
    skipMissing0[1] = '\0';
    skipMissing = skipMissing0;

    for(i=0; i<argc-1; i++){if(strcmp(argv[i],"-minQual"    )==0){minQual    =atoi(argv[i+1]); break;}}
    for(i=0; i<argc-1; i++){if(strcmp(argv[i],"-minInsert"  )==0){minInsert  =atoi(argv[i+1]); break;}}
    for(i=0; i<argc-1; i++){if(strcmp(argv[i],"-maxInsert"  )==0){maxInsert  =atoi(argv[i+1]); break;}}
    for(i=0; i<argc-1; i++){if(strcmp(argv[i],"-maxMismatch")==0){maxMismatch=atoi(argv[i+1]); break;}}
    for(i=0; i<argc-1; i++){if(strcmp(argv[i],"-maxGapOpen" )==0){maxGapOpen =atoi(argv[i+1]); break;}}
    for(i=0; i<argc-1; i++){if(strcmp(argv[i],"-maxBestHit" )==0){maxBestHit =atoi(argv[i+1]); break;}}
    for(i=0; i<argc-1; i++){if(strcmp(argv[i],"-skipMissing")==0){skipMissing=     argv[i+1] ; break;}}
    for(i=0; i<argc;   i++){if(strcmp(argv[i],"-singleEnd"  )==0){singleEnd  =1; break;}}
    for(i=0; i<argc;   i++){if(strcmp(argv[i],"-bowtie"     )==0){bowtie     =1; break;}}
     //sortedByName = optionExists("sortedByName");
    //fragment   = optionExists("fragment");
    checkClArgs();
    //parseOutputFormat(oFormat);
    
    
    fprintf(stderr, "QC Filters: minQual=%d minInsert=%d maxInsert=%d maxMismatch=%d maxGapOpen=%d maxBestHit=%d oFormat%s skipMissing=%s\n", minQual,minInsert, maxInsert,maxMismatch,maxGapOpen,maxBestHit,oFormat,skipMissing);
        
        
        
        
    
    qcFilterBam(argv[1]);

    fprintf(stderr, "done.\n\n");
    
    free(op);
    free(num);
    //free(tag);
    //free(val);
    //free(strBuf);
    return 0;
}












