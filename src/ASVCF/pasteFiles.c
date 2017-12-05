#include <stdio.h>
#include <stdlib.h>
#include <zlib.h>
#include <string.h>
#include <math.h>

#define MAXF 1000
#define BUFSIZE 1000

long min(long a, long b){if(a>b){return b;}else{return a;}}
long max(long a, long b){if(a>b){return a;}else{return b;}}

int readLine(gzFile* f, char** cells, char* buf, int i0, int* ret0){
	int colnum=0;
	int breakFlag=0;
	int l=0, i;
	int ret = (*ret0);
	//if(ret<=0){return 0;}
	for(i=i0;i<ret/sizeof(char);i++){
		if(buf[i]=='\n'){
			cells[colnum][l]='\0';
			return i+1;
		}else if(buf[i]==' ' || buf[i]=='\t'){
			cells[colnum++][l]='\0';
			l=0;
		}else{
			if(l<MAXF-1){cells[colnum][l++]= buf[i]; }else{fprintf(stderr,"%s truncated!\n", cells[colnum]); };
		}
	}
	while((ret = gzread(*f, buf, BUFSIZE)) != 0 && ret!=-1){
		*ret0 = ret;
		for(i=0;i<ret/sizeof(char);i++){
			if(buf[i]=='\n'){
				cells[colnum][l]='\0';
				breakFlag=1;break;
			}else if(buf[i]=='\t'){
				cells[colnum++][l]='\0';
				l=0;
			}else{
				if(l<MAXF-1){ cells[colnum][l++]= buf[i]; }else{fprintf(stderr, "%s truncated!\n", cells[colnum]); };
			}
		}
		if(breakFlag>0){break;}
	}
	
	return i+1;	
} 

int main(int argc, char** argv){
	int i,j;
	if(argc==1){
		fprintf(stderr, "Usage: %s <file_1.gz, ..., file_K.gz> > output\n", argv[0]);
		return -1;
	}
	
	int K = argc-1;
	gzFile* fs;
	fs = (gzFile*)calloc(K, sizeof(gzFile));
	for(i=0; i<K; i++){
		if((fs[i] = gzopen(argv[i+1], "rb6f"))==NULL){
			fprintf(stderr, "Error: gzread failed to read \"%s\" (%dth file)!\n", argv[i+1], i+1);
			return -1;
		}
	}
	fprintf(stderr,"%d files are read\n", K);
	
	int ret;
	char* buf;
	buf=(char*)calloc(BUFSIZE,sizeof(char));
	int breakFlag=0;
	int P=0;
	int L=0;
	
	// counting ncol
	int header=0;
	int M=0;
	while((ret = gzread(fs[0], buf, BUFSIZE*sizeof(char))) != 0 && ret!=-1){
		for(i=0; i<ret/sizeof(char); i++){
			if(buf[i]=='#'){header=1;}
			if(buf[i]=='\n'){
				if(header>0){
					M++;
				}else{
					if(L==0)P++;
					L++;
				}
				header=0;
			}else if(buf[i]=='\t'){
				if(header==0){if(L==0)P++;}
			}
		}
	}
	gzseek(fs[0], 0L, SEEK_SET);
	fprintf(stderr, "%d Header lines and %d Data points x %d cells in %s\n", M, L, P, argv[1]);
	
	// print all headers
	if(M>0){
		int m=0;
		while((ret = gzread(fs[0], buf, sizeof(char))) != 0 && ret!=-1){
			printf("%c", buf[0]);
			if(buf[0]=='\n'){
				m++;
				if(m==M-1){
					printf("##FORMAT=<ID=AS,Number=.,Type=String,Description=\"Allele-specific counts from sequenced reads\">\n");
				}else if(m==M){
					break;
				}
			}
		}
	}else{
		printf("##FORMAT=<ID=AS,Number=.,Type=String,Description=\"Allele-specific counts from sequenced reads\">\n");	
	}
	// K times P string matrix
	char*** lines;
	lines=(char***)calloc(K,sizeof(char**));
	for(i=0; i<K; i++){
		lines[i]=(char**)calloc(P, sizeof(char*));
		for(j=0;j<P;j++){
			lines[i][j] = (char*)calloc(MAXF,sizeof(char));
		}
	}
	int p,k;
	int* is0;  is0 =(int*)calloc(K, sizeof(int));
	int* ret0; ret0=(int*)calloc(K, sizeof(int));
	for(i=0;i<K;i++){ret0[k]=is0[k]=0;}
	char** bufs;
	bufs = (char**)calloc(K, sizeof(char*));
	for(k=0; k<K; k++){bufs[k] = (char*)calloc(BUFSIZE, sizeof(char));}
	for(i=0; i<L; i++){
		is0[0] = readLine(fs,   lines[0],   bufs[0], is0[0], ret0);
		is0[1] = readLine(fs+1, lines[1]+9, bufs[1], is0[1], ret0+1);
		for(p=0; p<P; p++){
			for(k=0; k<K; k++){
				if(k<K-1){// body
					if(p<8){
						printf("%s\t",lines[k][p]);
					}else if(p==8){
						printf("%s:AS\t",lines[k][p]);
					}else{
						printf("%s:",lines[k][p]);
					}
				}else{// as
					if(p>=9){
						if(p<P-1){
							printf("%s\t",lines[k][p]);
						}else{
							printf("%s\n",lines[k][p]);
						}
					}
				}
			}
		}
	}
	return 0;
}


