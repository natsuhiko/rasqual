#include <math.h>
#include <stdio.h>
#include <time.h>
#include <stdlib.h>
#include "sort.h"
#include "util.h"
#include "parseVCF.h"
//#include "sort.h"
#include "nbem.h"
#include "nbglm.h"




void printDouble(double* v, long n){
	long i;
		for(i=0; i<n; i++){fprintf(stderr, "%lf ", v[i]);}
		fprintf(stderr, "\n");
}
int isExon(long pos, long* starts, long* ends, long nexon){
	int i;
	if(nexon==0){return 1;}
	for(i=0; i<nexon; i++){
		if(starts[i]<=pos && pos<=ends[i]){
			return 1;
		}
	}
	return 0;
}

long countFields(char* posStr){
	long i;
	long n=0;
	if(strcmp(posStr,"NA")==0){return 0;}
	for(i=0; i<strlen(posStr); i++){
		if(posStr[i]==','){
			n++;
		}
	}
	return n+1;
}
long splitCSV(char* posStr, long n, long* pos){
	long i;
	long offset=0;
	char posChar[100];
	for(i=0; i<n-1; i++){
		sscanf(posStr+offset, "%ld,", pos+i);
		sprintf(posChar, "%ld", pos[i]);
		offset += strlen(posChar)+1;
	}
	sscanf(posStr+offset, "%ld", pos+i);
	sprintf(posChar, "%ld", pos[i]);
	return 0;
}

int isTestReg(char* chr, char* chr0, long TSS, long TSSPROX, long pos){// chr: chrom in which test SNP exists  chr0: feature chrom
	if(transQTL>0){
		if((strcmp(chr, chr0)!=0) || (strcmp(chr, chr0)==0 && (pos <= TSS-TSSPROX || pos >= TSS+TSSPROX))){
			return 1;
		}
	}else{
		if(pos >= TSS-TSSPROX && pos <= TSS+TSSPROX){
			return 1;
		}
	}
	return 0;
}
int isSameChr(char* chr, char* chr0){
	if(transQTL>0){
		if(strcmp(chr, chr0)==0){
			return 1;
		}else{
			return 0;
		}
	}
	return 1;
}

void usage();
/*void usage(){
	FILE *f;
	int flag;
	char ch;
	f = fopen("usage.txt","r");
	while (( ch = getc(f)) != EOF ) {
		if(ch=='#'){flag=1;}
		if(flag==0){
			putchar(ch);
		}else if(ch=='\n'){
			flag=0;
		}
	}   
}*/

int main(int argc, char** argv){
		
	long i, m, l, j;
	
	// verbose level
	verbose=0;
	verbose2=0;
	verbose3=0;
	for(i=0; i<argc; i++){if(strcmp(argv[i],"-v")==0){verbose2=100;verbose3=1;}}//atoi(argv[i+1]);}}
	for(i=0; i<argc; i++){if(strcmp(argv[i],"-V")==0){verbose2=100;verbose3=2;}}//atoi(argv[i+1]);}}
	for(i=0; i<argc; i++){if(strcmp(argv[i],"-Vv")==0){verbose2=100;verbose3=3;}}//atoi(argv[i+1]);}}
	for(i=0; i<argc; i++){if(strcmp(argv[i],"-VV")==0){verbose2=100;verbose3=4;}}//atoi(argv[i+1]);}}
	
	
	// initilization for loading VCF files
	if(init()==0){fprintf(stderr, "error on parseVCF\n"); return -1;};
	
	
	
	if(verbose2>10){
		fprintf(stderr, "\n");
		fprintf(stderr, "#################################\n");
		fprintf(stderr, "#		 RASQUAL v 1.0		 #\n");
		fprintf(stderr, "#################################\n");
		fprintf(stderr, "\n");
	}
	
	// usage
	if(argc==1){usage(); return -1;}//else{for(i=0;i<argc;i++){fprintf(stderr, "%s.", argv[i]);}}
	
    gzFile postVCF=NULL;
    for(i=0; i<argc-1; i++){if(strcmp(argv[i],"--vcf-file")==0){postVCF=gzopen(argv[i+1], "ab6f");break;}}
    FILE* fysim=NULL;
    for(i=0; i<argc-1; i++){if(strcmp(argv[i],"--bin-file")==0){fysim=fopen(argv[i+1], "ab");break;}}
    char* csnp=NULL;
#ifdef SIM
	double simp[5];
    csnp = (char*)calloc(100, sizeof(char));
	for(i=0; i<argc-1; i++){if(strcmp(argv[i],"--simulate")==0){sscanf(argv[i+1], "%[^,],%lf,%lf,%lf,%lf,%lf", csnp, simp, simp+1, simp+2, simp+3, simp+4);break;}}
    if(verbose3>0){printf("%s %lf %lf %lf %lf %lf\n", csnp, simp[0], simp[1], simp[2], simp[3], simp[4]);}
#endif
    
	// init hyper parameter setting
	ab=10.001;
	phiab=10.001;
	ad = 1.01;
	bd = 1.99;
	kappa = 1.01*2.0; 
	omega = 0.1*2.0; 
	sigma2 = 10000.0;
    beta0 = 100000.0;
    
	
	for(i=0; i<argc-1; i++){if(strcmp(argv[i],"-SIGMA")==0 && atof(argv[i+1])>0.){sigma2=(double)atof(argv[i+1]);break;}}	
	for(i=0; i<argc-1; i++){if(strcmp(argv[i],"-BETA")==0  && atof(argv[i+1])>0.){beta0 =(double)atof(argv[i+1]);break;}}
	for(i=0; i<argc-1; i++){if(strcmp(argv[i],"-ABPHI")==0 && atof(argv[i+1])>0.){phiab=(double)atof(argv[i+1]);break;}}	
	for(i=0; i<argc-1; i++){if(strcmp(argv[i],"-ABPI")==0  && atof(argv[i+1])>0.){ab   =(double)atof(argv[i+1]);break;}}	
	for(i=0; i<argc-1; i++){if(strcmp(argv[i],"-KAPPA")==0 && atof(argv[i+1])>0.){kappa=(double)atof(argv[i+1]);break;}}	
	for(i=0; i<argc-1; i++){if(strcmp(argv[i],"-OMEGA")==0 && atof(argv[i+1])>0.){omega=(double)atof(argv[i+1]);break;}}	
	for(i=0; i<argc-1; i++){if(strcmp(argv[i],"-ADELTA")==0 && atof(argv[i+1])>0.){ad=(double)atof(argv[i+1]);break;}}	
	for(i=0; i<argc-1; i++){if(strcmp(argv[i],"-BDELTA")==0 && atof(argv[i+1])>0.){bd=(double)atof(argv[i+1]);break;}}	
	
	if(verbose3>0){
		fprintf(stderr, "Hyper-parameters :\n");
		fprintf(stderr, " ab_phi: %lf\n", phiab);
		fprintf(stderr, " ab_pi : %lf\n", ab);
		fprintf(stderr, " kappa : %lf\n", kappa);
		fprintf(stderr, " omega : %lf\n\n", omega);
        fprintf(stderr, " beta0 : %lf\n\n", beta0);
        fprintf(stderr, " sigma2: %lf\n\n", sigma2);
	}
	
	
	// fix parameters
	fixParam=(int*)calloc(10,sizeof(int)); for(i=0; i<10; i++){fixParam[i]=0;}
	for(i=0; i<argc; i++){if(strcmp(argv[i],"--fix-delta")==0){fixParam[1]=1; break;}}
	for(i=0; i<argc; i++){if(strcmp(argv[i],"--fix-phi")==0){fixParam[2]=1; break;}}
	for(i=0; i<argc; i++){if(strcmp(argv[i],"--fix-theta")==0){fixParam[4]=1; break;}}
	
	//
	// program arguments
	//
	// het type (0: ours   1: sun   2: mcv)
	hetType=0;
	for(i=0; i<argc-1; i++){if(strcmp(argv[i],"-hetType")==0){hetType=atoi(argv[i+1]);}}
	// force to fit model
	int forced=0;
	for(i=0; i<argc; i++){if(strcmp(argv[i],"--force")==0){forced=1;}}//atoi(argv[i+1]);}}
	// print result
	int printAll=1;
	for(i=0; i<argc; i++){if(strcmp(argv[i],"--lead-snp")==0 || strcmp(argv[i],"-t")==0){printAll=0; break;}}
    int printRandom=0;
	for(i=0; i<argc; i++){if(strcmp(argv[i],"--random-snp")==0){printAll=0; printRandom=1; break;}}
    for(i=0; i<argc-1; i++){if(strcmp(argv[i],"--rsnp")==0){printAll=0; csnp=argv[i+1]; break;}}
	testImprinting=0;
	for(i=0; i<argc; i++){if(strcmp(argv[i],"--imprint")==0  || strcmp(argv[i],"-i")==0){testImprinting=1; break;}}
	// print VCF
	int printVCF=0;
	for(i=0; i<argc-1; i++){if(strcmp(argv[i],"--posterior-genotype")==0  || strcmp(argv[i],"-g")==0){printVCF=1; postVCF=gzopen(argv[i+1], "ab6f");break;}}
	// for likelihood ratio ties
	srand((unsigned)(time(NULL)+getpid()));
	int rand0=rand(), rand1=0;
	randomize=0;
	for(i=0; i<argc; i++){if(strcmp(argv[i],"-r")==0 || strcmp(argv[i],"--random-permutation")==0){randomize=1;}}
	// Usse ASE or not
	ASE=1;
	for(i=0; i<argc; i++){if(strcmp(argv[i],"--population-only")==0){ASE=0;}}
	for(i=0; i<argc; i++){if(strcmp(argv[i],"--as-only")==0){ASE=2;}}
	// fit Null model only
	Null=0;
	for(i=0; i<argc; i++){if(strcmp(argv[i],"--null")==0){Null=1;}}
	// trans QTL mapping
	transQTL=0;
	char* chr0;
	for(i=0; i<argc-1; i++){if(strcmp(argv[i],"-trans")==0){transQTL=1; chr0=argv[i+1]; if(verbose3>0){fprintf(stderr, "Trans QTL Mapping\nChromosome=%s\n\n", chr0);} break;}}
	// allelic prob from RSQ
	allelicProbEstByErrorRate=0;
	for(i=0; i<argc; i++){if(strcmp(argv[i],"-z")==0 || strcmp(argv[i],"--convert-imputation-score")==0){allelicProbEstByErrorRate=1; break;}}
	// GL relaxation
	double relaxGL=0.001;
	for(i=0; i<argc; i++){if(strcmp(argv[i],"--nominal-allelic-probability")==0){relaxGL=0.0; break;}}
	// No prior genotype
	int noPriorGenotype=0;
	for(i=0; i<argc; i++){if(strcmp(argv[i],"--population-allele-frequency")==0){noPriorGenotype=1; break;}}
	// No genotype likelihood
	noGL=0;
	for(i=0; i<argc; i++){if(strcmp(argv[i],"--fix-genotype")==0){noGL=1; break;}}
	// posterior update
	noPosteriorUpdate=0;
	for(i=0; i<argc; i++){if(strcmp(argv[i],"--no-posterior-update")==0){noPosteriorUpdate=1; break;}}
	if(verbose2>10){
		if(noPosteriorUpdate>0)fprintf(stderr, "No posterior probability update\n");
		if(allelicProbEstByErrorRate>0)fprintf(stderr, "Allelic probability has been estimated by error rate\n\n");
	}
	
	
	
	
	
	// data loading
	long N=0;
	for(i=0; i<argc-1; i++){if(strcmp(argv[i],"-n")==0 || strcmp(argv[i],"--sample-size")==0){N=atoi(argv[i+1]); break;}}
	if(N==0){fprintf(stderr, "Sample size is inappropriate\n"); return -1;}
	
	// Sample weights
	double* w;
	w=(double*)calloc(N, sizeof(double));
	for(i=0; i<N; i++){w[i]=1.0; }

	//if(verbose>10){printDouble(w,N);}
	
	// initial #snps
	long M=0, L=0, K=0;
	for(i=0; i<argc-1; i++){if(strcmp(argv[i],"-m")==0 || strcmp(argv[i],"--number-of-fsnps")==0){M=atoi(argv[i+1]); break;}} // # of exon SNPs
	for(i=0; i<argc-1; i++){if(strcmp(argv[i],"-l")==0 || strcmp(argv[i],"--number-of-testing-snps")==0){L=atoi(argv[i+1]); break;}} // # of Test SNPs
	for(i=0; i<argc-1; i++){if(strcmp(argv[i],"-j")==0 || strcmp(argv[i],"--feature-id")==0){K=atoi(argv[i+1]); break;}} // Region ID {1..50K}
	char* gid=NULL;
	for(i=0; i<argc-1; i++){if(strcmp(argv[i],"-f")==0 || strcmp(argv[i],"--feature-name")==0){gid=argv[i+1];}}
	if(gid==NULL){gid=(char*)calloc(100,sizeof(char)); sprintf(gid, "%ld", K);}
	
	//if(M==0 || L==0 || K==0){fprintf(stderr, "%s Numbers are inappropriate M=%ld L=%ld K=%ld!", gid, M, L, K); if(forced==0){fprintf(stderr, "\n"); return -1;}}
	
	if(verbose2>10){
		fprintf(stderr, "Sample Size		  : %ld\n", N);
		fprintf(stderr, "Init No. fSNPs	   : %ld\n", M);
		fprintf(stderr, "Init No. tested SNPs : %ld\n", L);
		fprintf(stderr, "Kth feature		  : %ld\n\n", K);
	}
	
	FILE* fy=NULL; // total fragment counts
	FILE* fk=NULL; // offset for Negative Binomial
	FILE* fx=NULL; // covariates
	
	for(i=0; i<argc-1; i++){if(strcmp(argv[i],"-y")==0 || strcmp(argv[i],"--feature-counts")==0){fy=fopen(argv[i+1],"rb"); break;}}
	for(i=0; i<argc-1; i++){if(strcmp(argv[i],"-k")==0 || strcmp(argv[i],"--sample-offsets")==0){fk=fopen(argv[i+1],"rb"); break;}}
	for(i=0; i<argc-1; i++){if(strcmp(argv[i],"-x")==0 || strcmp(argv[i],"--covariates")==0){fx=fopen(argv[i+1],"rb"); break;}}
	
	if(fy==NULL || fk==NULL){fprintf(stderr, "input files are not specified!\n"); return -1;}
	int freadres;
	long P=1;
	double* X;
	double tmp;
	if(fx==NULL){
		X = (double*)calloc(N, sizeof(double));
	}else{
		P = ncol(fx, N)+1;
		X = (double*)calloc(N*P, sizeof(double));
		freadres = fread(X+N, sizeof(double), N*(P-1), fx);
		for(i=1; i<P; i++){
			double m = mean(X+N*i, N);
			//double ss = sd(X+N*i, N);
			for(j=0; j<N; j++){
				X[i*N+j] -= m;
				//X[i*N+j] /= ss;
			}
		}
	}
	for(i=0; i<N; i++){ X[i] = 1.0; }
	for(i=0; i<argc-1; i++){if(strcmp(argv[i],"-p")==0 || strcmp(argv[i],"--number-of-covariates")==0){P=atoi(argv[i+1])+1;}} // # of covs after change
	
	double* y;
	y = (double*)calloc(N, sizeof(double));
	fseek(fy, N*(K-1)*sizeof(double), SEEK_SET);
	freadres = fread(y, sizeof(double), N, fy);
	
	double* ki;
	ki= (double*)calloc(N+randomize*M*N, sizeof(double));
	fseek(fk, N*(K-1)*sizeof(double), SEEK_SET);
	freadres = fread(ki, sizeof(double), N, fk);
	
	double totki = sum(ki, N)/(double)N;
	for(i=0; i<N; i++){ki[i] /= totki;}
	double tot = iwsum(y, ki, N);
	
	
	
	// feature region(s)
	char* sa=NULL;
	char* sb=NULL;
	for(i=0; i<argc-1; i++){if(strcmp(argv[i],"-s")==0 || strcmp(argv[i],"--feature-starts")==0){sa=argv[i+1];}}
	for(i=0; i<argc-1; i++){if(strcmp(argv[i],"-e")==0 || strcmp(argv[i],"--feature-ends")==0){sb=argv[i+1];}}
	int nexon=0;
	long* starts;
	long* ends;
	double cdnLen=0.0;
	if(sa!=NULL){
		nexon=countFields(sa);
	}else{fprintf(stderr, "no feature region specified\n"); return -1;}
	
	starts = (long*)calloc(nexon, sizeof(long));
	ends   = (long*)calloc(nexon, sizeof(long));
	splitCSV(sa, nexon, starts);
	splitCSV(sb, nexon, ends);
	for(i=0;i<nexon;i++){ cdnLen+=(double)(ends[i]-starts[i]); }
	if(verbose3>0)printLong(starts, nexon);
	if(verbose3>0)printLong(ends, nexon);
	
	// TSS
	long TSS = (starts[0]+ends[nexon-1])/2;
	for(i=0; i<argc-1; i++){if(strcmp(argv[i],"-c")==0 || strcmp(argv[i],"--cis-midpoint")==0){TSS=atol(argv[i+1]);}}
	long TSSPROX=270000000*2;
	for(i=0; i<argc-1; i++){if(strcmp(argv[i],"-w")==0 || strcmp(argv[i],"--cis-window-size")==0){TSSPROX=atol(argv[i+1])/2;}}
	
		
	for(i=0; i<argc-1; i++){if(strcmp(argv[i],"--effective-feature-length")==0){cdnLen = (double)atol(argv[i+1]);}}
	if(verbose2>10){fprintf(stderr, "Effective feature length			: %.0lf\n\n", cdnLen);}
	
	// SNP filter
	double MAF=0.05;
	double HWE=Qchisq(1.0e-8, 1.0);
	double RSQ=0.7;
	for(i=0; i<argc-1; i++){if(strcmp(argv[i],"-a")==0 || strcmp(argv[i],"--minor-allele-frequency")==0){MAF=(double)atof(argv[i+1]); break;}}
	for(i=0; i<argc-1; i++){if(strcmp(argv[i],"-h")==0 || strcmp(argv[i],"--hardy-weinberg-pvalue")==0){HWE=Qchisq((double)atof(argv[i+1]), 1.0); break;}}
	for(i=0; i<argc-1; i++){if(strcmp(argv[i],"-q")==0 || strcmp(argv[i],"--imputation-quality")==0){RSQ=(double)atof(argv[i+1]); break;}}
	double coverageDepth=0.05;
	for(i=0; i<argc-1; i++){if(strcmp(argv[i],"-d")==0 || strcmp(argv[i],"--min-coverage-depth")==0){coverageDepth=(double)atof(argv[i+1]);break;}}
	
	
	// loading VCF files
	VCF_info vinfo;
	int* dip;
	double* gen;
	double* gl;
	double* ap;
	long* ase;
	long pos;
	char* rs;
	char* al;
	dip=(int*)calloc(N*2, sizeof(int));
	gen=(double*)calloc(N, sizeof(double));
	gl =(double*)calloc(N*3, sizeof(double));
	ap =(double*)calloc(N*2, sizeof(double));
	ase=(long*)calloc(N*2, sizeof(long));
	al = (char*)calloc(2, sizeof(char));
	rs=(char*)calloc(1000, sizeof(char));
	int genType=1;
	
	double* Z=NULL;
	Z = (double*)calloc(N*(L+M+1)*2, sizeof(double));
	double* Y=NULL;
	Y = (double*)calloc(N*M*2, sizeof(double));
	double* km=NULL;
	km = (double*)calloc((M+1), sizeof(double));
	double* exon;
	exon = (double*)calloc(L, sizeof(double));
	
	long* poss;	     poss = (long*)calloc(L+2, sizeof(long));
	double* afs;     afs= (double*)calloc(L+2, sizeof(double));
	double* hwes;    hwes= (double*)calloc(L+2, sizeof(double));
	double* ias;     ias= (double*)calloc(L+2, sizeof(double));
	double* rsq;     rsq= (double*)calloc(L+2, sizeof(double));
	double* crs;     crs= (double*)calloc(L+2, sizeof(double));
	char** rss=NULL; rss = (char**)calloc(L+2, sizeof(char*));
	char** als;	     als = (char**)calloc(L+2, sizeof(char*));
	char** chrs;	 chrs =(char**)calloc(L+2, sizeof(char*));
	for(l=0; l<L; l++){//fprintf(stderr, "%ld ", l);
		als[l] = (char*)calloc(2, sizeof(char)); 
		rss[l] = NULL;  
		rss[l] = (char*)calloc(100, sizeof(char));
		if(rss[l]==NULL){return -1;}
		chrs[l] = (char*)calloc(100, sizeof(char)); 
	}	
	if(Z==NULL || rss==NULL || Y==NULL){fprintf(stderr, "memory allocation error!");return -1;}
	if(verbose2>10)fprintf(stderr, "\nMemory has been allocated...\n");
	
	m=l=0;
	numOfLoci=0;
	double asgaf[2];
	double asaf[2];
	double Y1=0.0;
	double Y0 = mean(y, N);
	long l0=0;
	double ep;
    maxCsnp=-1;
	while(parseLine(chrs[l], &pos, rss[l], als[l], &vinfo, dip, gen, gl, ap, ase, N, genType)>0){
        
        if(csnp!=NULL){ if(strcmp(rss[l],csnp)==0){ maxCsnp=l;} }  

        
        //simulation end
		
		afs[l] = getAF(gen,w,N);
        if(allelicProbEstByErrorRate>0 || noGL>0){
            hwes[l] = getHWE(dip,w,N);
        }else{
            hwes[l] = getHWEfromAP(ap, dip, w, N);
		}
        ias[l] = getIAfromAP(ap, gl, w, N); 
        rsq[l] = vinfo.RSQ;
		if(rsq[l]<0.0){rsq[l]=ias[l];}
		crs[l] = getCR(gen,w,N); //	 fprintf(stderr,"%lf ",crs[l]);
		poss[l]=pos;
		
        
		ep = afs[l]*(1-afs[l])*(1.0-rsq[l]);  //if(ep<0.01){ep=0.01;}
		if(noGL>0){ep=0.0;}
		
		if(allelicProbEstByErrorRate>0 || noGL>0){
			for(i=0; i<N; i++){
				Z[l*2*N+i*2]   = dip[i*2]==0  ?ep:(1.0-ep);
				Z[l*2*N+i*2+1] = dip[i*2+1]==0?ep:(1.0-ep);
			}
            
		}else{
			cblas_dcopy(2*N, ap, 1, Z+l*2*N, 1);
			//for(i=0; i<N*2; i++){if(Z[l*2*N+i]>1-ep){Z[l*2*N+i]=1-ep;}else if(Z[l*2*N+i]<ep){Z[l*2*N+i]=ep;}} // ep adjust
			//fprintf(stderr, "%s ", rss[l]); print(ap,2*N);
		}
		if(relaxGL>0.0){for(i=0; i<2*N; i++){
			if(Z[l*2*N+i]>1.0-relaxGL){
				Z[l*2*N+i] = 1.0 - relaxGL;
			}else if(Z[l*2*N+i]<relaxGL){
				Z[l*2*N+i] = relaxGL;
			}
		}}
		if(noPriorGenotype>0){
			for(i=0; i<N*2; i++){Z[l*2*N+i] = afs[l];}
		}
        
        // testing SNPs
        //if(l==maxCsnp){if(verbose3>0){fprintf(stderr, "%ld %s AF=%lf RSQ=%lf HWE=%lf %s %s %ld %ld %ld\n", l, rss[l], afs[l], rsq[l], hwes[l], chrs[l], chr0, TSS, TSSPROX, pos);}}
        if(afs[l]>MAF && afs[l]<1.0-MAF && rsq[l]>RSQ && (hwes[l]<HWE||noPriorGenotype==1) && isTestReg(chrs[l], chr0, TSS, TSSPROX, pos)){// tested
            exon[l]=1.0; numOfLoci++;
		}else{exon[l]=-1.0;}// not tested
	
		//if(strcmp("rs9945088",rss[l])==0){fprintf(stderr, "%s ", rss[l]); print(ap,2*N); print(Z+l*2*N,2*N);}
		//if(strcmp("rs2296475",rss[l])==0){fprintf(stderr, "%s ", rss[l]); print(ap,2*N);}
		if(ASE>0 && isExon(pos, starts, ends, nexon) && afs[l]>0.005 && afs[l]<0.995 && rsq[l]>0. && isSameChr(chrs[l],chr0)>0){
			km[m] = asgaf[0] = asgaf[1] = asaf[0] = asaf[1] = 0.0;
			for(i=0; i<N; i++){
				Y[m*2*N+i*2]   = (double)ase[i*2];
				Y[m*2*N+i*2+1] = (double)ase[i*2+1];
				asgaf[0] += gen[i]*(Y[m*2*N+i*2]+Y[m*2*N+i*2+1])/ki[i];
				asgaf[1] += 2.0*(Y[m*2*N+i*2]+Y[m*2*N+i*2+1])/ki[i];
				asaf[0] += Y[m*2*N+i*2];
				asaf[1] += Y[m*2*N+i*2+1];
				km[m] += (Y[m*2*N+i*2]+Y[m*2*N+i*2+1])/ki[i]/((tot*150)/cdnLen);
				//Y1 += (Y[m*2*N+i*2]+Y[m*2*N+i*2+1])/(double)N;
			}
			//if((asaf[0]+asaf[1])>0.0){asaf[1]=asaf[1]/(asaf[0]+asaf[1]);}else{asaf[1]=0.0;};
			asgaf[0]/=asgaf[1];
			asaf[0]/=(asaf[0]+asaf[1]);
			//km[m] /= tot;
			if(verbose3>0){fprintf(stderr, "%s %ld km=%lf af=%lf hwe=%lf asgaf=%lf  asaf=%lf\n", rss[l], pos, km[m], afs[l], hwes[l], asgaf[0], asaf[0]);}
			if(km[m]>coverageDepth && asgaf[0]>0.005 && asgaf[0]<0.995 && asaf[0]>0.0 && asaf[0]<1.0){
                for(i=0; i<N; i++){Y1 += (Y[m*2*N+i*2]+Y[m*2*N+i*2+1])/(double)N;}
				if(verbose3>0){fprintf(stderr, "RSQ=%lf EP=%lf\n", rsq[l], ep);}
				//if(m==0){ hoge++; }else{ if(poss[l]-poss[l0]<150){hoge += ((double)(poss[l]-poss[l0]))/150.0; }else{hoge++;} }
				km[m]=1.0;
				m++; // num of eSNPs
				exon[l]*=2.0;
				l0 = l;
			}
		}
		if(exon[l]!=(-1)){l++;}
	}
    //if(verbose3>0)fprintf(stderr, "%lf %ld %s\n", exon[maxCsnp], maxCsnp, csnp);
	asNonasRatio0=0.0;
    //if(m>0){ if(Y1>Y0){asNonasRatio0=0.99/(double)m;}else{asNonasRatio0 = Y1/Y0/(double)(m);} }else{ asNonasRatio0 = 0.0; }
#ifdef SIM
    int nrowZ = l;
    //printf("%s\t%lf\t%lf\t%lf\t%ld\n",gid,Y1,Y0,Y1/Y0,m);
    //return 0;
	if(m>0){ if(Y1>Y0){asNonasRatio0=0.99/(double)m;}else{asNonasRatio0 = Y1/Y0/(double)(m);} }else{ asNonasRatio0 = 0.0; }
    for(i=0; i<m; i++){km[i]=asNonasRatio0;}
    km[m] = 1.0-asNonasRatio0*((double)m);
    //print(km,m+1);
#endif
	double effrlen = 150.0;
	//if(m>0){ if(hoge*effrlen>cdnLen){asNonasRatio0=0.99;}else{asNonasRatio0=hoge*effrlen/cdnLen;} }else{ asNonasRatio0=0.0; }
	//if(m>0){ if(((double)m)*75.0>(cdnLen-74.0)){asNonasRatio=0.99/(double)(m);}else{asNonasRatio=75.0/(cdnLen-74.0);} }else{ asNonasRatio=0.0; }
	//for(l=0;l<m;l++){km[l]=asNonasRatio0;}asNonasRatio0=0.0;
	//for(l=0;l<m;l++){if(m<6.0){km[l]=m;}else{km[l]=5.0;}}
	//for(l=0;l<m;l++){ if(Y1>Y0){km[l]=0.99/(double)(m);}else{ km[l]=Y1/Y0/(double)(m);} }
	//print(km, m);
	//for(l=0;l<m;l++){for(i=0; i<N; i++){Y[i*2] += Y[m*2*N+i*2]; Y[i*2+1] += Y[m*2*N+i*2+1];   }}
	//m = 1;
	if(verbose2>10){
		fprintf(stderr, "\nData has been loaded...\n\n\n");
		fprintf(stderr, "No. fSNPs	   : %ld\n", m);
		fprintf(stderr, "No. testing SNPs: %ld\n", numOfLoci);
		fprintf(stderr, "No ase/ase ratio: %lf\n", asNonasRatio0);
		//fprintf(stderr, "Effective no of eSNPs: %lf\n", hoge);
	}//fprintf(stderr, "No. eSNPs : %d\n\n", m);
	
	
	// average number of joint tests
	NOfSigLoci=(double)numOfLoci;
	for(i=0; i<argc-1; i++){if(strcmp(argv[i],"--number-of-significant-loci")==0){NOfSigLoci=(double)atof(argv[i+1]);break;}} 
	//fprintf(stderr, "num of sig=%lf\n", NOfSigLoci); 
	
	
	
	
	
	
	// ASEQTL
	double* lkhdDiff;
	double* randlkhdDiff;
	double* ppi;
	double* pdelta;
	double* pphi;
	double* pbeta;
	double* ptheta;
	double* pasr;
	double* ptval;
	double* pkld;
	int* pitr;
	int* pbound;
	lkhdDiff = (double*)calloc(L+2, sizeof(double)); 
	randlkhdDiff = (double*)calloc(L+2, sizeof(double)); for(i=0; i<L+1; i++){randlkhdDiff[i]=rand();}
	ppi	  = (double*)calloc(L+2, sizeof(double));
	pdelta   = (double*)calloc(L+2, sizeof(double));
	pphi	 = (double*)calloc(L+2, sizeof(double));
	pbeta   = (double*)calloc(L+2, sizeof(double));
	ptheta   = (double*)calloc(L+2, sizeof(double));
	pasr   = (double*)calloc(L+2, sizeof(double));
	ptval   = (double*)calloc(L+2, sizeof(double));
	pkld   = (double*)calloc(2*(L+2)+1, sizeof(double));
	pitr	 = (int*)calloc(L+2, sizeof(int));
	pbound   = (int*)calloc(L+2, sizeof(int));
	//fprintf(stderr,"%s ", gid);
	ones = (double*)calloc(m, sizeof(double)); for(i=0; i<m; i++){ ones[i]=1.0; }
	//if(randomize>0){randomise4(y, Y, ki, X, N, m, P);}
    
    
    
    
    // 
    // Simulate fragment counts
    //
#ifdef SIM
    if(simp[0]<0.0){
        gzclose(postVCF);
        fwrite(y, sizeof(double), N, fysim);
        fclose(fysim);
        return 0;
    }
    init_gsl_rand();
    double *yx, *zx, *zc, *zct, *zxt;
    yx = calloc(N*2, sizeof(double));
    zx = Z;
    zc = Z+N*maxCsnp*2;
    zct = Z+N*2*(L+M);
    zxt = Z+N*2*L;
    double H[10];
    double Ktmp[8];
    double K2[20];
    int* Yk;    Yk = (int*)calloc((m+1)*N, sizeof(int));
    int yxtmp;
    double thij, thijA, thijB;
    getK(simp[0], simp[1], simp[2], Ktmp, K2);
    int ll=0;
    for(i=0; i<N; i++){
        rG(zc+i*2, zct+i*2);
        thij = ( (1.0-zct[i*2])*(1.0-zct[i*2+1])*2.0*(1.0-simp[0]) + zct[i*2]*(1.0-zct[i*2+1]) + (1.0-zct[i*2])*zct[i*2+1] + zct[i*2]*zct[i*2+1]*2.0*simp[0] )*ki[i];//*(1.0-asNonasRatio0*((double)m));
        y[i] = (double)rbinom(exp(simp[4])*thij, simp[3]*thij);
        // dirichlet multinom
        //for(j=0; j<m+1; j++){km[(m+1)*2+j]=km[j]*simp[3]*thij*10000.0;}; rdirmultinom(Yk+i*(m+1), m+1, (int)y[i], km+m+1, km+(m+1)*2);
        // multinom
        rmultinom(Yk+i*(m+1), m+1, (int)y[i], km);
        //printInt(Yk+i*(m+1), m+1); for(j=0; j<m; j++){fprintf(stderr,"%.0lf ",Y[N*j*2+i*2]+Y[N*j*2+i*2+1]);}; fprintf(stderr,"\n");
    }
    for(l=0;l<nrowZ;l++){
            zx = Z+l*N*2;
            
            if(fabs(exon[l])>1.5){// feature snp
                for(i=0; i<N; i++){
                    if(strcmp(rss[l],csnp)==0){// causal snp
                        getPrior1(H, zct+i*2, 1.0);
                    }else{
                        rG(zx+i*2, zxt+i*2);
                        getPrior( H, zct+i*2, zxt+i*2, 1.0);
                    }
                    thijA = thijB = 0.0;
                    for(j=0; j<10; j++){
                        thijA += H[j]*K2[j*2];
                        thijB += H[j]*K2[j*2+1];
                    }
                    //yx[i*2]   = rbinom(exp(simp[4])*ki[i]*asNonasRatio0*thijA, simp[3]*ki[i]*asNonasRatio0*thijA); //Y[N*2*ll+i*2];
                    //yx[i*2+1] = rbinom(exp(simp[4])*ki[i]*asNonasRatio0*thijB, simp[3]*ki[i]*asNonasRatio0*thijB); //Y[N*2*ll+i*2+1];
                    yxtmp = Yk[i*(m+1)+ll]; //rbinom(exp(simp[4])*ki[i]*asNonasRatio0*(thijA+thijB), simp[3]*ki[i]*asNonasRatio0*(thijA+thijB)); 
                    yx[i*2]   = (double)rbetabinom(simp[3]*ki[i]*thijA, simp[3]*ki[i]*thijB, yxtmp); 
                    yx[i*2+1] = (double)(yxtmp - (int)yx[i*2]);
                    //y[i] += yx[i*2]+yx[i*2+1];
                }
                ll++;
            }else{// other snp
                for(i=0; i<N; i++){ yx[i*2] = yx[i*2+1] = 0.0; }
            }
            gzprintf(postVCF, "%s\t%ld\t%s\t%c\t%c\t.\t.\tRSQ=%lf;GID=%s\tGT:GP:AS\t", chrs[l], poss[l], rss[l], als[l][0], als[l][1], rsq[l], gid);
            for(i=0; i<N; i++){
                gzprintf(postVCF, "%d|%d:%lf,%lf,%lf:%.0lf,%.0lf", 
                         zx[i*2]>0.5?1:0, zx[i*2+1]>0.5?1:0 , 
                         (1.0-zx[i*2])*(1.0-zx[i*2+1]), zx[i*2]*(1.0-zx[i*2+1])+(1.0-zx[i*2])*zx[i*2+1], zx[i*2]*zx[i*2+1], 
                         //zx[i*2], zx[i*2+1], 
                         yx[i*2], yx[i*2+1]); 
                if(i<N-1){gzprintf(postVCF, "\t");}else{gzprintf(postVCF, "\n");}
            }
    }
    gzclose(postVCF);
    fwrite(y, sizeof(double), N, fysim);
    fclose(fysim);
    return 0;
#endif
    //
    //  Simulation end
    //
    
    if((m+1)*numOfLoci>3e4 && forced==0){// stop when N of fSNPs x rSNPs too large
        fprintf(stderr, "Estimated computational time is too long for %s (Sample Size=%ld, NofrSNPs=%ld, NoffSNPs=%ld)...aborted.\n", gid, N, numOfLoci, m);
        printf("%s\tSKIPPED\t%s\t-1\tN\tN\t-1.0\t-1.0\t-1.0\t-1.0\t0.0\t-1.0\t-1.0\t-1.0\t-1.0\t-1.0\t%ld\t%ld\t-1\t-1\t-1\t0.0\t0\t-1.0\t-1.0\n", 
               gid, chrs[0], 
               m, numOfLoci);
        return 0;
    }
    
    if(tot<=0.0){ printf("%s\trs\t%s\t0\tX\tX\t-1.0\t-1.0\t-1.0\t-1.0\t-1.0\t-1.0\t-1.0\t-1.0\t-1.0\t0\t%ld\t%ld\t-1\t-1\t-1\t-1\t-1\t0\t0.0\n", gid, "0", m,numOfLoci); return 0;}
	if(numOfLoci==0){printf("%s\trs\t%s\t0\tX\tX\t-1.0\t-1.0\t-1.0\t-1.0\t-1.0\t-1.0\t-1.0\t-1.0\t-1.0\t0\t%ld\t%ld\t-1\t-1\t-1\t-1\t-1\t0\t0.0\n", gid, "0", m,numOfLoci); return 0;}
	
    if(printRandom>0){
        init_gsl_rand();
        maxCsnp = (long)floor(runif()*((double)numOfLoci));
    }
    
	//if(m>0 || forced>0){
	ASEQTLALL(y, Y, Z, X, P, ki, w, km, exon, L, N, m, rss, lkhdDiff, ppi, pdelta, pphi, pbeta, ptheta, pasr, pitr, pbound, ptval, pkld);
	if(verbose3>0){fprintf(stderr, "ASEQTLALL finished.\n");}
	//}else{fprintf(stderr,"No marker tested.\n"); return -1;}
	
	
	//-------------- VCF start --------------
	if(printVCF>0){
		double *yx, *zx, *zc;
		double *Zx;
		Zx = Z+N*L*2;
		zc = Zx+N*m*2;
		int ll=0;
		for(l=0;l<L;l++){
			if(fabs(exon[l])>1.5){// feature snp
				if(l==maxCsnp){
					gzprintf(postVCF, "%s\t%ld\t%s\t%c\t%c\t.\tLEAD_QTL_SNP\t%s\tGT:GP:AP:AS\t", chrs[l], poss[l], rss[l], als[l][0], als[l][1], gid);
				}else{
					gzprintf(postVCF, "%s\t%ld\t%s\t%c\t%c\t.\t.\t%s\tGT:GP:AP:AS\t", chrs[l], poss[l], rss[l], als[l][0], als[l][1], gid);
				}
				zx = Zx+N*2*ll;
				yx = Y+N*2*ll;
				for(i=0; i<N; i++){
					gzprintf(postVCF, "%d|%d:%lf,%lf,%lf:%lf,%lf:%.0lf,%.0lf", 
							 zx[i*2]>0.5?1:0, zx[i*2+1]>0.5?1:0 , 
							 (1.0-zx[i*2])*(1.0-zx[i*2+1]), zx[i*2]*(1.0-zx[i*2+1])+(1.0-zx[i*2])*zx[i*2+1], zx[i*2]*zx[i*2+1], 
							 zx[i*2], zx[i*2+1], yx[i*2], yx[i*2+1]);
					if(i<N-1){gzprintf(postVCF, "\t");}else{gzprintf(postVCF, "\n");}
				}
				ll++;
			}else if(l==maxCsnp){
				gzprintf(postVCF, "%s\t%ld\t%s\t%c\t%c\t.\tLEAD_QTL_SNP\t%s\tGT:GP:AP:AS\t", chrs[l], poss[l], rss[l], als[l][0], als[l][1], gid);
				for(i=0; i<N; i++){
					gzprintf(postVCF, "%d|%d:%lf,%lf,%lf:%lf,%lf:%.0lf,%.0lf", 
							 zc[i*2]>0.5?1:0, zc[i*2+1]>0.5?1:0 , 
							 (1.0-zc[i*2])*(1.0-zc[i*2+1]), zc[i*2]*(1.0-zc[i*2+1])+(1.0-zc[i*2])*zc[i*2+1], zc[i*2]*zc[i*2+1], 
							 zc[i*2], zc[i*2+1], 0.0, 0.0); 
					if(i<N-1){gzprintf(postVCF, "\t");}else{gzprintf(postVCF, "\n");}
				}
			}
		}
		gzclose(postVCF);
	}
	
	//--------------- VCF end ---------------
	
	
	
	double maxld=0.0;
	int maxl  = -1;
	int maxlr = -1;
	int numOfTies=0;
	//fprintf(stderr,"TSS=%ldn",TSS);
	afs[L] = afs[L+1] = 0.5; rsq[L] = rsq[L+1]=1.0; hwes[L]=hwes[L+1]=0.0;
	if(Null==0){for(l=0; l<L; l++){
	//if(exon[l]>=0.0)printf("%s\t%s\t%s\t%ld\t%c\t%c\t%lf\t%lf\t%lf\t%lf\t%3.22lf\t%lf\t%ld\t%ld\n", gid, rss[l], chr, poss[l], als[l][0], als[l][1], afs[l], hwes[l], ias[l], rsq[l], lkhdDiff[l], ppi[l], m, l);
		if(afs[l]>MAF && afs[l]<1.0-MAF && rsq[l]>RSQ && (hwes[l]<HWE||noPriorGenotype==1) && isTestReg(chrs[l], chr0, TSS, TSSPROX, poss[l]) ){// && (pbound[l+1]+pbound[l+1+L])<1 && ptheta[l+1]<=ptheta[l+1+L]){
			//printf("%ld %lf\n", poss[l],lkhdDiff[l]);
			if(round(maxld*100000)<round(lkhdDiff[l]*100000)){
				numOfTies=0;
				rand0=rand();
				maxld=lkhdDiff[l];
				maxl=l;
				maxlr=l;
			}else if(round(maxld*100000)==round(lkhdDiff[l]*100000)){
				numOfTies++;
				if(abs(poss[l]-TSS)<abs(poss[maxl]-TSS)){
					//maxld=lkhdDiff[l];
					maxl=l;	
				}
				rand1=rand();
				//if(randlkhdDiff[maxlr]<randlkhdDiff[l]){
				if(rand0<rand1){
					maxlr=l;
					rand0=rand1;
				}
			}
		}else{lkhdDiff[l]=0.0;}
		//printf("%s\t%ld\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n", rss[l], poss[l], afs[l], hwes[l], ias[l], rsq[l], lkhdDiff[l], ppi[l]);
	}
	l=maxl; }else{l=L;}//null
	double minBH = getMinBH(lkhdDiff, L);
	if(printAll>0){
		l=0;
		for(l0=0; l0<L; l0++){
			if(abs(exon[l0])>0.5){printf("%s\t%s\t%s\t%ld\t%c\t%c\t%lf\t%lf\t%lf\t%lf\t%.10lf\t%lf\t%lf\t%lf\t%lf\t%ld\t%ld\t%ld\t%d\t%d\t%ld\t%lf\t%d\t%lf\t%lf\n", gid, rss[l], chrs[l], poss[l], als[l][0], als[l][1], afs[l], hwes[l], ias[l], rsq[l], lkhdDiff[l],
					ppi[l], pdelta[l], pphi[l], ptheta[l], l, m, numOfLoci, pitr[l], pitr[0], maxlr<0?-1:poss[maxlr], lkhdDiff[L+1], pbound[l]+pbound[L], pkld[l], pkld[l+L+2]);
				l++;
			}
		}
    }else{
		if(l<0){// no max lkhd
			printf("%s\trs\t%s\t0\tN\tN\t-1.0\t-1.0\t-1.0\t-1.0\t-1.0\t-1.0\t-1.0\t-1.0\t-1.0\t%ld\t%ld\t%ld\t-1\t-1\t-1\t-1\t-1\t0.0\n", gid, chrs[0], l,m,numOfLoci);
		}else if(l<L){// normal
			printf("%s\t%s\t%s\t%ld\t%c\t%c\t%lf\t%lf\t%lf\t%lf\t%.10lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%ld\t%ld\t%d\t%d\t%ld\t%lf\t%d\t%lf\t%lf\n", 
				   gid, rss[l], chrs[l], poss[l], als[l][0], als[l][1], afs[l], hwes[l], ias[l], 
                   rsq[l], 
				   lkhdDiff[l], 
//minBH, 
                   ppi[l], pdelta[l], pphi[l], ptheta[l], pbeta[l], m, numOfLoci, pitr[l], pitr[0], 
				   maxlr<0?-1:poss[maxlr], lkhdDiff[L+1], pbound[l]+pbound[L], pkld[l], pkld[l+L+2]);
		}else if(l==L){// null lkhd
			printf("%s\tAI\t%s\t-1\tN\tN\t-1.0\t-1.0\t-1.0\t%lf\t%.10lf\t%lf\t%lf\t%lf\t%lf\t%ld\t%ld\t%ld\t%d\t%d\t%ld\t%lf\t%d\t%lf\t%lf\t%lf\n", 
                   gid, chrs[l], asNonasRatio0, 
                   lkhdDiff[L],  ppi[L], pdelta[L], pphi[L], ptheta[L], l, m, numOfLoci, pitr[L], pitr[L], 
                   maxlr<0?-1:poss[maxlr], lkhdDiff[L+1], pbound[L], pkld[L], pkld[L], pkld[2*L+2]);
		}
	}
	if(testImprinting>0){
		printf("%s\tIMP\t%s\t-1\tN\tN\t-1.0\t-1.0\t-1.0\t-1.0\t%.10lf\t%lf\t%lf\t%lf\t%lf\t%ld\t%ld\t%ld\t%d\t%d\t%ld\t%lf\t%d\t%lf\t%lf\n", 
               gid, chrs[0], 
               lkhdDiff[L+1], ppi[L+1], pdelta[L+1], pphi[L+1], ptheta[L+1], l, m, numOfLoci, pitr[L+1], pitr[L+1], 
               maxlr<0?-1:poss[maxlr], lkhdDiff[L+1], pbound[L+1]+pbound[L], pkld[L+1], pkld[2*L+3]);
	}
	return 0;
}









