#include <math.h>
#include <float.h>
#include <stdio.h>
#include <time.h>
#include <stdlib.h>
#include "sort.h"
#include "util.h"
#include "parseVCF.h"
//#include "sort.h"
#include "nbem.h"
#include "nbglm.h"

int parseScond(char* scond, char*** cs0, double** ce0){
    int offs=0, offs1=0;
    int ncs=0, i;
    if(scond==NULL){return ncs;}
    ncs++;
    for(i=0; i<strlen(scond); i++){
        if(scond[i]==','){ncs++;}
    }
    char **cs = (char**)calloc(ncs, sizeof(char*));
    double *ce = (double*)calloc(ncs, sizeof(double));
    for(i=0; i<ncs; i++){ cs[i] = (char*)calloc(1000, sizeof(double)); }
    for(i=0; i<ncs; i++){
        offs += offs1;
        sscanf(scond+offs, "%[^,:]:%lf%n", cs[i], ce+i, &offs1);
        offs1++;
    }
    (*cs0)=cs; (*ce0)=ce;
    return ncs;
}


double trunc(double x){
	return x<1.0e-5 ? 1.0e-5 : x;
}

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
        fprintf(stderr, "#          RASQUAL v 1.0        #\n");
        fprintf(stderr, "#################################\n");
        fprintf(stderr, "\n");
    }

    // usage
    if(argc==1){usage(); return -1;}//else{for(i=0;i<argc;i++){fprintf(stderr, "%s.", argv[i]);}}

    gzFile postVCF=NULL;
    for(i=0; i<argc-1; i++){if(strcmp(argv[i],"--vcf-file")==0){postVCF=gzopen(argv[i+1], "ab6f");break;}}
    char* csnp=NULL;
#ifdef SIM
    FILE* fysim=NULL;
    for(i=0; i<argc-1; i++){if(strcmp(argv[i],"--bin-file")==0){fysim=fopen(argv[i+1], "ab");break;}}
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
    // N of threads
    nthreads=1;
    for(i=0; i<argc-1; i++){if(strcmp(argv[i],"--n-threads")==0){nthreads=atoi(argv[i+1]); break;}}
    if(nthreads>1)pthread_mutex_init( &mutex, NULL );
    // het type (0: ours   1: sun   2: cht)
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
    //srand((unsigned)193849);
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
    char* chr0=NULL;
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
    double* w=NULL;	w=(double*)calloc(N, sizeof(double)); for(i=0; i<N; i++){w[i]=1.0;}

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
        fprintf(stderr, "Sample Size		   : %ld\n", N);
        fprintf(stderr, "Init No. fSNPs	       : %ld\n", M);
        fprintf(stderr, "Init No. testint SNPs : %ld\n", L);
        fprintf(stderr, "Kth feature		   : %ld\n\n", K);
    }

    FILE* fy=NULL; // total fragment counts
    FILE* fk=NULL; // offset for Negative Binomial
    FILE* fx=NULL; // covariates
    char* scond=NULL;
    for(i=0; i<argc-1; i++){if(strcmp(argv[i],"-y")==0 || strcmp(argv[i],"--feature-counts")==0){fy=fopen(argv[i+1],"rb"); break;}}
    for(i=0; i<argc-1; i++){if(strcmp(argv[i],"-k")==0 || strcmp(argv[i],"--sample-offsets")==0){fk=fopen(argv[i+1],"rb"); break;}}
    for(i=0; i<argc-1; i++){if(strcmp(argv[i],"-k2")==0 || strcmp(argv[i],"--conditional-snps")==0){ scond=argv[i+1]; break;}}
    for(i=0; i<argc-1; i++){if(strcmp(argv[i],"-x")==0 || strcmp(argv[i],"--covariates")==0){if((fx=fopen(argv[i+1],"rb"))==NULL){fprintf(stderr, "Covariate file does not exists.\n"); return 1;}; break;}}

    if(fy==NULL || fk==NULL){fprintf(stderr, "input files are not specified!\n"); return -1;}
    int freadres;
    long P=1;
    double* X;
    //double tmp;
    
    if(fx==NULL){
        X = (double*)calloc(N*P, sizeof(double));
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
	// fragment count
	double* y=NULL;   y = (double*)calloc(N, sizeof(double));
	fseek(fy, N*(K-1)*sizeof(double), SEEK_SET);
	freadres = fread(y, sizeof(double), N, fy);
	
    // allele specific offsets (k2[i*2]+k2[i*2+1])/ 2 -> sample offset
    double* ki2=NULL;  ki2 = (double*)calloc(N*2+randomize*M*N*2, sizeof(double));
    for(i=0; i<2*N; i++){ki2[i] = 1.0;}
	char** condSnp=NULL;
	double* condEff=NULL;
    int nCondSnp = parseScond(scond, &condSnp, &condEff);
        
        // sample offsets
	double* ki=NULL;  ki= (double*)calloc(N+randomize*M*N+N, sizeof(double));
	fseek(fk, N*(K-1)*sizeof(double), SEEK_SET);
	freadres = fread(ki, sizeof(double), N, fk);
	if(verbose3>10){fprintf(stderr, "%d\n", freadres);}
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
	
	starts = (long*)calloc(nexon+1, sizeof(long));
	ends   = (long*)calloc(nexon+1, sizeof(long));
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
	
    double fMAF=0.0;
	double fHWE=DBL_MAX;
	double fRSQ=0.0;
	double coverageDepth=0.05;
	for(i=0; i<argc-1; i++){if(strcmp(argv[i],"--minor-allele-frequency-fsnp")==0){fMAF=(double)atof(argv[i+1]); break;}}
	for(i=0; i<argc-1; i++){if(strcmp(argv[i],"--hardy-weinberg-pvalue-fsnp")==0){ fHWE=Qchisq((double)atof(argv[i+1]), 1.0); break;}}
	for(i=0; i<argc-1; i++){if(strcmp(argv[i],"--imputation-quality-fsnp")==0){    fRSQ=(double)atof(argv[i+1]); break;}}
	for(i=0; i<argc-1; i++){if(strcmp(argv[i],"-d")==0 || strcmp(argv[i],"--min-coverage-depth")==0){coverageDepth=(double)atof(argv[i+1]);break;}}
	
	//
	// loading VCF files
    //
	VCF_info vinfo;
	int* dip=NULL;     dip=(int*)calloc(N*2, sizeof(int));
	double* gen=NULL;  gen=(double*)calloc(N, sizeof(double));
	double* gl=NULL;   gl =(double*)calloc(N*3, sizeof(double));
	double* ap=NULL;   ap =(double*)calloc(N*2, sizeof(double));
	long* ase=NULL;    ase=(long*)calloc(N*2, sizeof(long));
	long pos=0;     
	//char* rs=NULL;     rs=(char*)calloc(1000, sizeof(char));
	//char* al=NULL;     al = (char*)calloc(2, sizeof(char));
	int genType=1;
	
	double* Z=NULL;
	Z = (double*)calloc(N*L*2+N*(M+1+2)*2+N*(M+2)*2, sizeof(double));
	double* Y=NULL;
	Y = (double*)calloc(N*(M+2)*2, sizeof(double));
	double* km=NULL;
	km = (double*)calloc((M+1+2), sizeof(double));
    if(Z==NULL || Y==NULL || km==NULL){fprintf(stderr, "memory allocation error (Z, Y, km)...aborted.\n");return 1;}
	
    double* exon=NULL;    exon = (double*)calloc(L+2, sizeof(double));
	long* poss=NULL;	  poss = (long*)calloc(L+2, sizeof(long));
	double* afs=NULL;     afs= (double*)calloc(L+2, sizeof(double));
	double* hwes=NULL;    hwes= (double*)calloc(L+2, sizeof(double));
	double* ias=NULL;     ias= (double*)calloc(L+2, sizeof(double));
	double* rsq=NULL;     rsq= (double*)calloc(L+2, sizeof(double));
	double* crs=NULL;     crs= (double*)calloc(L+2, sizeof(double));
	char** rss=NULL;      rss = (char**)calloc(L+2, sizeof(char*));
	char** als=NULL;	  als = (char**)calloc(L+2, sizeof(char*));
	char** chrs=NULL;	  chrs =(char**)calloc(L+2, sizeof(char*));
    if(exon==NULL || poss==NULL || afs==NULL || hwes==NULL || ias==NULL || rsq==NULL || crs==NULL || rss==NULL || als==NULL || chrs==NULL){fprintf(stderr, "memory allocation error (Z, Y, km)...aborted.\n");return 1;}
	for(l=0; l<L+2; l++){//fprintf(stderr, "%ld ", l);
		als[l] = (char*)calloc(2, sizeof(char)); 
		rss[l] = NULL;  
		rss[l] = (char*)calloc(100, sizeof(char));
		if(rss[l]==NULL){return -1;}
		chrs[l] = (char*)calloc(100, sizeof(char)); 
	}
    
	if(verbose2>10)fprintf(stderr, "\nMemory has been allocated...\n");
	
	m=l=0;
	numOfLoci=0;
	double asgaf[2];
	double asaf[2];
	double Y1=0.0;
	//double Y0 = mean(y, N);
	long l0=0;
	double ep;
    maxCsnp=-1;
    int nco=0;
    //FILE* vcffile;
    //vcffile=fopen("vcffile.txt","r");
    while(parseLine(chrs[l], &pos, rss[l], als[l], &vinfo, dip, gen, gl, ap, ase, N, genType)>0){
        if(csnp!=NULL){ if(strcmp(rss[l],csnp)==0){ maxCsnp=l;} }

        afs[l] = getAF(gen,w,N);
        if(allelicProbEstByErrorRate>0 || noGL>0){
            hwes[l] = getHWE(dip,w,N);
        }else{
            hwes[l] = getHWEfromAP(ap, dip, w, N);
        }
        ias[l] = getIAfromAP(ap, gl, w, N); 
        rsq[l] = vinfo.RSQ;
        // jkl
        if(verbose3>1){fprintf(stderr, "%lf %lf %lf %lf\n", afs[l], hwes[l], ias[l], rsq[l]);}
        if(rsq[l]<0.0){rsq[l]=ias[l];}
        crs[l] = getCR(gen,w,N); 
        poss[l]=pos;

        ep = afs[l]*(1-afs[l])*(1.0-rsq[l]);
        if(noGL>0){ep=0.0;}
            if(allelicProbEstByErrorRate>0 || noGL>0){
            for(i=0; i<N; i++){
                Z[l*2*N+i*2]   = dip[i*2]==0  ?ep:(1.0-ep);
                Z[l*2*N+i*2+1] = dip[i*2+1]==0?ep:(1.0-ep);
            }
        }else{
            cblas_dcopy(2*N, ap, 1, Z+l*2*N, 1);
        }
        if(relaxGL>0.0){
            for(i=0; i<2*N; i++){
                if(Z[l*2*N+i]>1.0-relaxGL){
                    Z[l*2*N+i] = 1.0 - relaxGL;
                }else if(Z[l*2*N+i]<relaxGL){
                    Z[l*2*N+i] = relaxGL;
                }
            }
        }
        if(noPriorGenotype>0){
            for(i=0; i<N*2; i++){Z[l*2*N+i] = afs[l];}
        }
        
        // testing SNPs
        if(afs[l]>MAF && afs[l]<1.0-MAF && rsq[l]>RSQ && (hwes[l]<HWE||noPriorGenotype==1) && isTestReg(chrs[l], chr0, TSS, TSSPROX, pos)){// tested
            exon[l]=1.0; numOfLoci++;
        }else{
            exon[l]=-1.0;
        }// not tested
        
        // conditional SNPs
        for(nco=0; nco<nCondSnp; nco++){
            if(strcmp(rss[l], condSnp[nco])==0){
                //printM(Z+l*2*N, 2, N);
                for(i=0; i<N*2; i++){
                    ki2[i] *= 2.0*( Z[l*2*N+i]*condEff[nco] + (1.0-Z[l*2*N+i])*(1.0-condEff[nco]) );
                }
            }
        }
        
        // feature SNPs
        if(ASE>0 && isExon(pos, starts, ends, nexon) && afs[l]>fMAF && afs[l]<(1.0-fMAF) && (hwes[l]<fHWE||noPriorGenotype==1) && rsq[l]>fRSQ && isSameChr(chrs[l],chr0)>0){
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
                km[m]=1.0;
                m++; // num of eSNPs
                if(m>M){fprintf(stderr, "\nThe number of fSNPs is greater than that specified by \'-m\' option.\nAborted...\n"); return -1;}
                exon[l]*=2.0;
                l0 = l;
            }
        }
        if(exon[l]!=(-1)){
            l++; 
            if(l>L){fprintf(stderr, "\nThe number of lines from STDIN is greater than the number of loci specified by \'-l\' option.\nAborted...\n"); return -1;}
        }
        // 2: fSNP & rSNP; 1: rSNP; -2: fSNP
    }
    //printM(ki, 1, N);
	for(i=0; i<N; i++){ 
        double ki2tot = (ki2[i*2]+ki2[i*2+1])/2.0;
        ki[i] *= ki2tot; 
        ki2[i*2+0] = ki2[i*2+0]/ki2tot;
        ki2[i*2+1] = ki2[i*2+1]/ki2tot;
        //fprintf(stderr, "%lf %lf %lf\n", ki2[i*2]+ki2[i*2+1], ki2[i*2], ki2[i*2+1]); 
    }
    //printM(ki2, 2, 10);
    //return 0;
    //for(i=0; i<20; i++){ki2[i]/=2.0;}
    //printM(ki, 1, N);
    //for(i=0; i<N; i++){ X[N*(P-1)+i] = log(ki2[i*2]+ki2[i*2+1]); fprintf(stderr, "%lf %lf %lf\n", ki2[i*2]+ki2[i*2+1], ki2[i*2], ki2[i*2+1]); }
	asNonasRatio0=0.0;
    
    //totki = sum(ki, N)/(double)N;
	//for(i=0; i<N; i++){ki[i] /= totki;}
    
    
	if(verbose2>10){
		fprintf(stderr, "\nData has been loaded...\n\n\n");
		fprintf(stderr, "No. fSNPs	   : %ld\n", m);
		fprintf(stderr, "No. all SNPs      : %ld %ld\n", l, L);
		fprintf(stderr, "No. testing SNPs: %ld\n", numOfLoci);
		fprintf(stderr, "No ase/ase ratio: %lf\n", asNonasRatio0);
	}
	
	// average number of joint tests
	NOfSigLoci=(double)numOfLoci;
	for(i=0; i<argc-1; i++){if(strcmp(argv[i],"--number-of-significant-loci")==0){NOfSigLoci=(double)atof(argv[i+1]);break;}} 
	
	if(((m+1)*numOfLoci>3e4 && forced==0) || tot<=0.0 || numOfLoci==0){// stop when N of fSNPs x rSNPs too large
        if((m+1)*numOfLoci>3e4 && forced==0)fprintf(stderr, "Estimated computational time is too long for %s (Sample Size=%ld, NofrSNPs=%ld, NoffSNPs=%ld)...aborted.\n", gid, N, numOfLoci, m);
        printf("%s\tSKIPPED\t%s\t-1\tN\tN\t-1.0\t-1.0\t-1.0\t0.0\t0.0\t-1.0\t-1.0\t-1.0\t-1.0\t-1.0\t%ld\t%ld\t-1\t-1\t-1\t0.0\t0\t-1.0\t-1.0\n", 
               gid, chrs[0], m, numOfLoci);
	if(printVCF==1){gzclose(postVCF);}
        return 0;
    }
    if(printRandom>0){
        init_gsl_rand();
        maxCsnp = (long)floor(runif()*((double)numOfLoci));
    }
    //
    // loading VCF end
    //
    
    
    
    
    
    
    
	// ASEQTL
	double* lkhdDiff=NULL;        lkhdDiff     = (double*)calloc(L+4, sizeof(double)); 
	double* randlkhdDiff=NULL;    randlkhdDiff = (double*)calloc(L+4, sizeof(double)); for(i=0; i<L+1; i++){randlkhdDiff[i]=rand();}
	double* ppi=NULL;             ppi	       = (double*)calloc(L+2, sizeof(double));
	double* pdelta=NULL;          pdelta       = (double*)calloc(L+2, sizeof(double));
	double* pphi=NULL;            pphi	       = (double*)calloc(L+2, sizeof(double));
	double* pbeta=NULL;           pbeta        = (double*)calloc(L+2, sizeof(double));
	double* ptheta=NULL;          ptheta       = (double*)calloc(L+2, sizeof(double));
	double* pasr=NULL;            pasr         = (double*)calloc(L+2, sizeof(double));
	double* ptval=NULL;           ptval        = (double*)calloc(L+2, sizeof(double));
	double* pkld=NULL;            pkld         = (double*)calloc(2*(L+2)+1, sizeof(double));
	int* pitr=NULL;               pitr	       = (int*)calloc(L+2, sizeof(int));
	int* pbound=NULL;             pbound       = (int*)calloc(L+2, sizeof(int));
    int* tested=NULL;             tested       = (int*)calloc(L+2, sizeof(int));
    if(lkhdDiff==NULL || ppi==NULL || pdelta==NULL || pphi==NULL || pbeta==NULL || ptheta==NULL || pasr==NULL || ptval==NULL || pkld==NULL || pitr==NULL || pbound==NULL){fprintf(stderr, "memory allocation error (Z, Y, km)...aborted.\n");return 1;}
    
	ones = (double*)calloc(m, sizeof(double)); for(i=0; i<m; i++){ ones[i]=1.0; }
    
    
    // moment est (NULL)
    nbbem_init0(y, ki, N, pbeta+L, ptheta+L);
    
    // glm est with covariates (NULL) -> ki updated
    double* betaGlm; betaGlm=(double*)calloc(P, sizeof(double)); clear1(betaGlm, P);
    double* work; work = (double*)calloc((2*N+3*P)*(P+2), sizeof(double));
    double* dki;  dki = ki+N+randomize*m*N; //(double*)calloc(N, sizeof(double));
    betaGlm[0] = pbeta[L];
    double lkhdNullglm = nbglm(y, X, ki, dki, w, N, P, betaGlm, ptheta+L, 0, work);
    if(verbose3>0){
        fprintf(stderr, "Glm. results:\n");
        fprintf(stderr, "Lkhd GLM=%lf Theta=%lf\n", lkhdNullglm, ptheta[L]);
        printM(betaGlm, 1, P);
        printM(y,1,84);	
        printM(ki,1,84);
        printM(dki,1,84);
        fprintf(stderr, "\n");
    }
    pbeta[L] = betaGlm[0];
    double mdki = mean(dki, N);
    for(i=0; i<N; i++){dki[i] /= mdki;}
     
    //randomization
    if(randomize>0){ randomPerm(y, Y, Z, ki, ki2, exon, L, N, m); }
    
	// Null
	ASEQTLALL(y, Y, Z, X, P, ki, ki2, w, km, exon, L, N, m, rss, lkhdDiff, ppi, pdelta, pphi, pbeta, ptheta, pasr, pitr, pbound, ptval, pkld);
	if(verbose3>0){fprintf(stderr, "Fitting null hypothesis finished.\n");}
    
    // Alternative    
    if(nthreads>1 && maxCsnp<0){
        RASQUAL_INP* ri; ri=(RASQUAL_INP*)calloc(nthreads, sizeof(RASQUAL_INP));
        pthread_t* pid; pid=(pthread_t*)calloc(nthreads, sizeof(pthread_t));
        if(verbose3>0){fprintf(stderr, "\nN of threads: %ld\n\n", nthreads);}
        for(i=0; i<nthreads; i++){
            ri[i].y = y;
            ri[i].Y = Y;
            ri[i].Z = Z;
            ri[i].X = X;
            ri[i].P = P;
            ri[i].ki = ki;
            ri[i].ki2 = ki2;
            ri[i].w = w;
            ri[i].km = km;
            ri[i].exon = exon;
            ri[i].L = L;
            ri[i].N = N;
            ri[i].Lx = m;
            ri[i].rss = rss;
            ri[i].lkhdDiff = lkhdDiff;
            ri[i].ppi = ppi;
            ri[i].pdelta=pdelta;
            ri[i].pphi = pphi;
            ri[i].pbeta = pbeta;
            ri[i].ptheta = ptheta;
            ri[i].pasr = pasr;
            ri[i].pitr = pitr;
            ri[i].pbound = pbound;
            ri[i].ptval = ptval;
            ri[i].pkld = pkld;
            //ri[i].rsnp_start = maxCsnp<0 ? 0 : maxCsnp;
            //ri[i].rsnp_end = maxCsnp<0 ? L : maxCsnp+1;
            ri[i].rsnp_start = 0; //(L/nthreads)*i;
            ri[i].rsnp_end   = L; //(L/nthreads)*(i+1);
            //ri[i].rsnp_end = ri[i].rsnp_end>L ? L : ri[i].rsnp_end;
            ri[i].tested = tested;
            ri[i].tid = i+1;
            int pthflag;
            if( (pthflag = pthread_create(pid+i, NULL, (void*)ASEQTLALL_MP, (void*)(ri+i))) != 0){fprintf(stderr, "Thread not created...aborted.\n");return 1;};
        }
        for(i=0; i<nthreads; i++){
            int pthflag;
            if( (pthflag = pthread_join(pid[i], NULL)) !=0 ){fprintf(stderr, "Thread not joined...aborted.\n"); return 1;};
        }
    }else{
        //fprintf(stderr, "main=%lf\n", ptheta[2]);
        ASEQTLALL_ALT(y, Y, Z, X, P, ki, ki2, w, km, exon, L, N, m, rss, lkhdDiff, ppi, pdelta, pphi, pbeta, ptheta, pasr, pitr, pbound, ptval, pkld, maxCsnp<0 ? 0 : maxCsnp, maxCsnp<0 ? L : maxCsnp+1, tested, 1);
    }
    
	
	if(verbose3>0){fprintf(stderr, "ASEQTLALL finished.\n");}
    
	
	
	
	
	
	double maxld=0.0;
	int maxl  = -1;
	int maxlr = -1;
	int numOfTies=0;
	//fprintf(stderr,"TSS=%ldn",TSS);
	afs[L] = afs[L+1] = 0.5; rsq[L] = rsq[L+1]=1.0; hwes[L]=hwes[L+1]=0.0;
	if(Null==0){for(l=0; l<L; l++){
        //fprintf(stderr, "main %ld\n", l);
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
	//double minBH = getMinBH(lkhdDiff, L);
    double* qval;
    qval = (double*)calloc(L+2, sizeof(double));
    double minBH = getBH(lkhdDiff, qval, L);
    //fprintf(stderr, "main minBH=%lf\n", minBH);
	if(printAll>0){
		l=0;
		for(l0=0; l0<L; l0++){
			if(abs(exon[l0])>0.5){
                //fprintf(stderr, "main %ld\n", l0);
                if(exon[l0]>0.5){
                    printf("%s\t%s\t%s\t%ld\t%c\t%c\t%lf\t%lf\t%lf\t%.10lf\t%.10lf\t%lf\t%lf\t%lf\t%lf\t%ld\t%ld\t%ld\t%d\t%d\t%ld\t%lf\t%d\t%lf\t%lf\n", gid, rss[l], chrs[l], poss[l], als[l][0], als[l][1], afs[l], hwes[l], //ias[l], 
                           rsq[l], 
                           qval[l],
                           lkhdDiff[l],
                           ppi[l], pdelta[l], pphi[l], ptheta[l], l, m, numOfLoci, pitr[l], pitr[0], maxlr<0?-1:poss[maxlr], lkhdDiff[L], pbound[l]+pbound[L], pkld[l*2], pkld[l*2+1]);
                }
				l++;
			}
		}
    }else{
		if(l<0){// no max lkhd
			printf("%s\trs\t%s\t0\tN\tN\t-1.0\t-1.0\t-1.0\t0.0\t-1.0\t-1.0\t-1.0\t-1.0\t-1.0\t%ld\t%ld\t%ld\t-1\t-1\t-1\t-1\t-1\t0.0\n", gid, chrs[0], l,m,numOfLoci);
		}else if(l<L){// normal
			printf("%s\t%s\t%s\t%ld\t%c\t%c\t%lf\t%lf\t%lf\t%.10lf\t%.10lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%ld\t%ld\t%d\t%d\t%ld\t%lf\t%d\t%lf\t%lf\n", 
				   gid, rss[l], chrs[l], poss[l], als[l][0], als[l][1], afs[l], hwes[l], //ias[l], 
                   rsq[l], 
                   minBH,
				   lkhdDiff[l], 
                   ppi[l], pdelta[l], pphi[l], ptheta[l], pbeta[l], m, numOfLoci, pitr[l], pitr[0], 
				   maxlr<0?-1:poss[maxlr], lkhdDiff[L], pbound[l]+pbound[L], pkld[l*2], pkld[l*2+1]);
		}else if(l==L){// null lkhd
			printf("%s\tAI\t%s\t-1\tN\tN\t-1.0\t-1.0\t-1.0\t0.0\t%.10lf\t%lf\t%lf\t%lf\t%lf\t%ld\t%ld\t%ld\t%d\t%d\t%ld\t%lf\t%d\t%lf\t%lf\n", 
                   gid, chrs[l],  
                   0.0,  ppi[L], pdelta[L], pphi[L], ptheta[L], l, m, numOfLoci, pitr[L], pitr[L], 
                   maxlr<0?-1:poss[maxlr], lkhdDiff[L], pbound[L], pkld[L*2], pkld[L*2+1]);
		}
	}
	if(testImprinting>0){
		printf("%s\tIMP\t%s\t-1\tN\tN\t-1.0\t-1.0\t-1.0\t0.0\t%.10lf\t%lf\t%lf\t%lf\t%lf\t%ld\t%ld\t%ld\t%d\t%d\t%ld\t%lf\t%d\t%lf\t%lf\n", 
               gid, chrs[0], 
               lkhdDiff[L+1], ppi[L+1], pdelta[L+1], pphi[L+1], ptheta[L+1], l, m, numOfLoci, pitr[L+1], pitr[L+1], 
               maxlr<0?-1:poss[maxlr], lkhdDiff[L], pbound[L+1]+pbound[L], pkld[(L+1)*2], pkld[(L+1)*2+1]);
	}
    if(nthreads>1)pthread_mutex_destroy( &mutex );
    
    
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
					gzprintf(postVCF, "%s\t%ld\t%s\t%c\t%c\t.\tLEAD_QTL_SNP\t%s;RSQ=1.0\tGT:GP:AP:AS\t", chrs[l], poss[l], rss[l], als[l][0], als[l][1], gid);
				}else{
					gzprintf(postVCF, "%s\t%ld\t%s\t%c\t%c\t.\t.\t%s;RSQ=1.0\tGT:GP:AP:AS\t", chrs[l], poss[l], rss[l], als[l][0], als[l][1], gid);
				}
				zx = Zx+N*2*ll;
				yx = Y+N*2*ll;
				for(i=0; i<N; i++){
					gzprintf(postVCF, "%d|%d:%lf,%lf,%lf:%lf,%lf:%.0lf,%.0lf", 
							 zx[i*2]>0.5?1:0, zx[i*2+1]>0.5?1:0 , 
							 (1.0-zx[i*2])*(1.0-zx[i*2+1]), zx[i*2]*(1.0-zx[i*2+1])+(1.0-zx[i*2])*zx[i*2+1], zx[i*2]*zx[i*2+1], 
							 log10(trunc(zx[i*2])), log10(trunc(zx[i*2+1])), yx[i*2], yx[i*2+1]);
					if(i<N-1){gzprintf(postVCF, "\t");}else{gzprintf(postVCF, "\n");}
				}
				ll++;
			}else if(l==maxCsnp){
				gzprintf(postVCF, "%s\t%ld\t%s\t%c\t%c\t.\tLEAD_QTL_SNP\t%s;RSQ=1.0\tGT:GP:AP:AS\t", chrs[l], poss[l], rss[l], als[l][0], als[l][1], gid);
				for(i=0; i<N; i++){
					gzprintf(postVCF, "%d|%d:%lf,%lf,%lf:%lf,%lf:%.0lf,%.0lf", 
							 zc[i*2]>0.5?1:0, zc[i*2+1]>0.5?1:0 , 
							 (1.0-zc[i*2])*(1.0-zc[i*2+1]), zc[i*2]*(1.0-zc[i*2+1])+(1.0-zc[i*2])*zc[i*2+1], zc[i*2]*zc[i*2+1], 
							 log10(trunc(zc[i*2])), log10(trunc(zc[i*2+1])), 0.0, 0.0); 
					if(i<N-1){gzprintf(postVCF, "\t");}else{gzprintf(postVCF, "\n");}
				}
			}
		}
		gzclose(postVCF);
	}
	//--------------- VCF end ---------------
    
    
	return 0;
}









