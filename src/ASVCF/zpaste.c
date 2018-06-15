#include <stdio.h>
#include <stdlib.h>
#include <zlib.h>
#include <string.h>
#include <math.h>

#define MAXF 1040329
//#define MAXF 100
#define BUFSIZE 1000000

long min(long a, long b){if(a>b){return b;}else{return a;}}
long max(long a, long b){if(a>b){return a;}else{return b;}}

int readLine(gzFile* f, char* cells, char* buf, int i0, int* ret0){
	int breakFlag=0;
	int l=0, i;
	int ret = (*ret0);
	//if(ret<=0){return 0;}
	for(i=i0;i<ret/sizeof(char);i++){
		if(buf[i]=='\n'){
			cells[l]='\0';
			return i+1;
		}else{
			if(l<MAXF-1){cells[l++]= buf[i];}else{fprintf(stderr,"truncation occurs! %d\n", l);}
		}
	}
	while((ret = gzread(*f, buf, BUFSIZE)) != 0 && ret!=-1){
		*ret0 = ret;
		for(i=0;i<ret/sizeof(char);i++){
			if(buf[i]=='\n'){
				cells[l]='\0';
				breakFlag=1;break;
			}else{
				if(l<MAXF-1){cells[l++]= buf[i];}else{fprintf(stderr,"truncation occurs! %d\n", l);}
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
	
	while((ret = gzread(fs[0], buf, BUFSIZE*sizeof(char))) != 0 && ret!=-1){
		for(i=0;i<ret/sizeof(char);i++){
			if(buf[i]=='\n'){
				if(L==0)P++;
				L++;
			}else if(buf[i]==' ' || buf[i]=='\t'){
				if(L==0)P++;
			}
		}
	}
	gzseek(fs[0], 0L, SEEK_SET);
	fprintf(stderr, "%d times %d cells in %s\n", L, P, argv[1]);
	
	// K times P string matrix
	char** lines;
	lines=(char**)calloc(K,sizeof(char*));
	for(i=0; i<K; i++){
		lines[i]=(char*)calloc(MAXF, sizeof(char));
	}
	int p,k;
	int* is0; is0 = (int*)calloc(K, sizeof(int));
	int* ret0; ret0=(int*)calloc(K, sizeof(int));
	for(i=0;i<K;i++){ret0[k]=is0[k]=0;}
	char** bufs;
	bufs = (char**)calloc(K, sizeof(char*));
	for(k=0; k<K; k++){bufs[k] = (char*)calloc(BUFSIZE, sizeof(char));}
	for(i=0; i<L; i++){
		for(k=0; k<K; k++){is0[k] = readLine(fs+k, lines[k], bufs[k], is0[k], ret0+k);}
		for(k=0; k<K; k++){
			if(k<K-1){
				printf("%s\t",lines[k]);
			}else{
				printf("%s\n",lines[k]);
			}
		}
	}
	return 0;
}


