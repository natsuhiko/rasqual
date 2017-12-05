#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <zlib.h>

int parseCigarPoint(long start, char* c1, char* seq, long* at, long K, char* A0, char* A1, long* Y, char* num);
void parseCigarPoint2(long start, char* c1, char* seq, long* at, long K, char* A0, double* Y, char* num, int flag);


char dim(gzFile f, char** pbuf, long sizeOfBuf, long* pN, long* pP, long skip){// N x P matrix
        char* buf;
        buf = *pbuf;
        long N=0;
        long P=0;
        int i;
        char sep;
        //char* buf;
        int ret;
        //buf = (char*)calloc(1000,sizeof(char));
        while((ret = gzread(f, buf, sizeOfBuf)) != 0 && ret!=-1){
//              fprintf(stderr, "%s   %ld %d\n", buf, sizeof(buf), ret);
                for(i=0; i<ret/sizeof(char); i++){
                        if(N==0+skip && (buf[i]=='\n'||buf[i]==' '||buf[i]=='\t'||buf[i]==',')){if(buf[i]!='\n'){sep=buf[i];} P++;}
                        if(buf[i]=='\n'){N++;}
                }
        }
        gzseek(f, 0L, SEEK_SET);
        (*pN)=N;
        (*pP)=P;
        return sep;
}


void print(long* Y, long n){
	long i;
	for(i=0; i<n; i++){fprintf(stderr,"%ld ",Y[i]);}
	fprintf(stderr,"\n");
}

char* op;

int main(int argc, char** argv){// 1000G.marker.gz -> rs pos A1 A2
	op = (char*)calloc(100, sizeof(char));
	if(argc==1){
		fprintf(stderr,"Usage: tabix as_input_file CHR | %s bgzipped_vcf_file [-printFrag] [-overlap] > as_output\n\n", argv[0]);
		fprintf(stderr,"AS input file format consists of 5 columns separated by TABs:\n");
		fprintf(stderr,"  - chromosome\n  - start position (1-based)\n  - end position (1-based (only count length for [MDN=PX] CIGAR operators))\n  - CIGAR\n  - read sequence\n");
		return -1;
	}
	
	int i;
	char* buf; buf = (char*)calloc(1000,sizeof(char));
	long printFrag=0;
	long overlap=0;
	//long ascount;
	for(i=0; i<argc; i++){if(strcmp(argv[i],"-printFrag")==0){printFrag=1;}}
	for(i=0; i<argc; i++){if(strcmp(argv[i],"-overlap")==0){overlap=1;}}
	//for(i=0; i<argc; i++){if(strcmp(argv[i],"-ascount")==0){ascount=1;}}
	
	
	long N=0;
	gzFile fmarker;
	
	if((fmarker = gzopen(argv[1], "rb6f"))==NULL){
		fprintf(stderr, "Error: gzread failed to read \"%s\"!\n", argv[1]);
		return -1;
	}

	long P=0;
	//N = getLine(fmarker);
	dim(fmarker, &buf, 1000*sizeof(char), &N, &P, 0);
fprintf(stderr,"N of markers: %ld %ld\n", N, P);
	
	char** rs;
	long* pos;
	char* cpos;
	char* A1;
	char* A2;
	long* Y;
	double* ds;
	rs = (char**)calloc(N,sizeof(char*));  for(i=0;i<N;i++){ rs[i]=(char*)calloc(100,sizeof(char)); }
	cpos = (char*)calloc(100,sizeof(char));
	pos = (long*)calloc(N,sizeof(long));
	ds = (double*)calloc(N,sizeof(double));
	A1 = (char*)calloc(N,sizeof(char));
	A2 = (char*)calloc(N,sizeof(char));
	Y = (long*)calloc(N*3,sizeof(long));
	
	for(i=0; i<N; i++){ds[i]=0.0;}
	for(i=0; i<N*3; i++){Y[i]=0;}
	int ret;
	//char* buf;
	int l=0,cl=0,ncol=0;
	//buf = (char*)calloc(1000,sizeof(char));
	while((ret = gzread(fmarker, buf, sizeof(buf))) != 0 && ret!=-1){
		for(i=0;i<ret/sizeof(char);i++){
			if(buf[i]=='\n'){
				l++;
				ncol=cl=0;
			}else if(buf[i]==' ' || buf[i]=='\t' || buf[i]=='\n'){
				if(ncol==2){rs[l][cl]='\0';}else if(ncol==1){cpos[cl]='\0'; pos[l]=atoi(cpos);}
				if(buf[i]=='\n'){l++; ncol=cl=0;}else{cl=0; ncol++;}
			}else{
				if(ncol==2){
					rs[l][cl++]=buf[i];
				}else if(ncol==1){
					cpos[cl++]=buf[i];
				}else if(ncol==3){
					if(cl==0){A1[l]=buf[i];}else{A1[l]='0';}
					cl++;
				}else if(ncol==4){
					if(cl==0){
						A2[l]=buf[i];
						if(A2[l]=='-'){A2[l]='1';A1[l]='0';}
					}else{// allele length > 1
						A2[l]='1';A1[l]='0';
					}
					cl++;
				}
			}
		}
	}
	fprintf(stderr, "marker position: [%ld..%ld]\n",pos[0],pos[N-1]);
	//fprintf(stderr, "marker position: [%c %c %c %c..%c %c %c %c]\n",A1[0],A1[1],A1[2],A1[3],A1[N-4],A1[N-3],A1[N-2],A1[N-1]);
	//fprintf(stderr, "marker position: [%c %c %c %c..%c %c %c %c]\n",A2[0],A2[1],A2[2],A2[3],A2[N-4],A2[N-3],A2[N-2],A2[N-1]);
	
	char* chr;
	char* c1;
	char* seq;
	char* num;
	chr = (char*)calloc(100, sizeof(char));
	seq = (char*)calloc(1000, sizeof(char));
	c1 =  (char*)calloc(1000, sizeof(char));
	num = (char*)calloc(100, sizeof(char));
	if(num==NULL || c1==NULL || seq==NULL || chr==NULL){fprintf(stderr, "memory allocation error!\n"); return 0;}
	long left, right;
	long* at;
	at = (long*)calloc(100, sizeof(long));
	long j,j0=0;
	long l0=0, l1=0;
	while(scanf("%s\t%ld\t%ld\t%s\t%s", chr, &left, &right, c1, seq)!=EOF){
		l=0;
		int flag=0;
//fprintf(stderr, "%s\t%ld\t%ld\t%s\t%s\n", chr, left, right, c1, seq);
		for(j=j0; j<N; j++){
			if(pos[j]>right){break;}
			if(pos[j]>=left){l++;}else{j0=j+1;}
		}
		if(l>0){
			// flag: number of snps existing in the fragment
			flag=parseCigarPoint(left, c1, seq, pos+j0, l, A1+j0, A2+j0, Y+3*j0, num);
			    parseCigarPoint2(left, c1, seq, pos+j0, l, A1+j0, ds+j0, num, flag);
		}
		if(printFrag==1){// tabix out put
			if( (overlap==1&&flag>0) || (overlap==0&&flag==0) ){
				printf("%s\t%ld\t%ld\t%s\t%s\n", chr, left, right, c1, seq);
			}
		}
		if(flag==0){// no snp
			l1++;
		}
		//	//printf("%s\t%ld\t%ld\t%s\t%s\n", chr, left, right, c1, seq);
		//}else{// more than one snp | flag captures the number of snps
		//	printf("%d\n", flag);
		//}
		l0++;
	}
	if(printFrag==0){
		for(j=0; j<N; j++){
			printf("%s\t%ld\t%c\t%c\t%ld\t%ld\t%ld\t%lf\n", rs[j], pos[j], A1[j], A2[j], Y[j*3], Y[j*3+1], Y[j*3+2], ds[j]);
		}
	}
	//FILE* fcount;
	//fcount=fopen("count.bin","wb");
	//fwrite(Y,sizeof(long),N*3,fcount);
	//fclose(fcount);
	fprintf(stderr, "%ld out of %ld fragments overlapped with SNPs\n", l0-l1, l0);
	return 0;
}


int parseCigarPoint(long start, char* c1, char* seq, long* at, long K, char* A0, char* A1, long* Y, char* num){
        //char* op; // operation
        //op = (char*)calloc(100, sizeof(char));
        long len;// length in cigar
        long start0=start;
        long nseq=0;
        char* pc1;
        long nc1=0;
	long k;
	long* count;
        //pc1 = c1;
        long end;
	int flag=0;
	int flagin=0;
	long nchar=0;
	//fprintf(stderr, "%ld\n", at[0]);
        while((sscanf(c1+nc1, "%ld%[MIDNSH=PX]", &len, op))>0){
//fprintf(stderr, "<%ld %s %ld>\n", len, op, nseq);
                nchar = (long)sprintf(num, "%ld", len);
                nc1 += nchar+1;
                //nc1 += strlen(num)+1;
                //pc1 += strlen(num)+1;
		flagin=0;
                if(op[0]=='M' || op[0]=='=' || op[0]=='X'){
                        end = start+len-1;
			long nseq0=nseq;
			for(k=0; k<K; k++){
				count = Y+k*3;
				nseq=nseq0;
                        	if(start <= at[k] && at[k] <= end && A0[k]!='0'){
					flagin++;
                                	nseq += (at[k] - start);
                                	if(seq[nseq]==A0[k]){
                                        	count[0]++;
                                	}else if(seq[nseq]==A1[k]){
                                        	count[1]++;
                                	}else{
                                        	count[2]++;
                                	}
                                	//break;
                        	}
			}
                        start += len;
                        nseq = nseq0+len;
                }else if(op[0]=='D' || op[0]=='N' || op[0]=='P'){
                        start += len;
                }else{// S, H and I
//fprintf(stderr, "<%ld ", nseq);
                        nseq += len;
//fprintf(stderr, "%ld>\n", nseq);
                }
		if(flagin>0){flag++;}
        }

        return flag;
}

// counting expression dosage
void parseCigarPoint2(long start, char* c1, char* seq, long* at, long K, char* A0, double* Y, char* num, int flag){
        long len;// length in cigar
        long start0=start;
        long nseq=0;
        char* pc1;
        long nc1=0;
        long k;
        long end;
        long nchar=0;
        while((sscanf(c1+nc1, "%ld%[MIDNSH=PX]", &len, op))>0){
                nchar = (long)sprintf(num, "%ld", len);
                nc1 += nchar+1;
                if(op[0]=='M' || op[0]=='=' || op[0]=='X'){
                        end = start+len-1;
                        long nseq0=nseq;
                        for(k=0; k<K; k++){
                                nseq=nseq0;
                                if(start <= at[k] && at[k] <= end && A0[k]!='0'){
                                        Y[k]+= 1.0/(double)flag;
                                }
                        }
                        start += len;
                        nseq = nseq0+len;
                }else if(op[0]=='D' || op[0]=='N' || op[0]=='P'){
                        start += len;
                }else{// S, H and I
                        nseq += len;
                }
        }
}
