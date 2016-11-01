#include "util.h"



void init_gsl_rand(){
	gsl_rng_env_setup();
	rngT = gsl_rng_default;
	rng = gsl_rng_alloc (rngT);
	gsl_rng_set (rng, (unsigned int)time( NULL )+getpid());
}

void rGrand(double* z, double* zt){
	zt[0] = (z[0]>0.5 ? 1.0 : 0.0);
	zt[1] = (z[1]>0.5 ? 1.0 : 0.0);
}
void rG(double* z, double* zt){
	int a = gsl_ran_binomial (rng, z[0], 1);
	zt[0] = (a==1 ? 1.0 : 0.0);
	a = gsl_ran_binomial (rng, z[1], 1);
	zt[1] = (a==1 ? 1.0 : 0.0);
}

double runif(){
	return gsl_ran_flat (rng, 0.0, 1.0);
}

void rmultinom(int* ys, int m, int Y, double* probs){
	gsl_ran_multinomial(rng, m, Y, probs, ys);
}


void rdirmultinom(int* ys, int m, int Y, double* probs, double* alpha){
	gsl_ran_dirichlet(rng, m, alpha, probs);
	gsl_ran_multinomial(rng, m, Y, probs, ys);
}

int rbetabinom(double thA, double thB, int n){
	int k = gsl_ran_binomial (rng, gsl_ran_beta (rng, thA, thB), n);
	return k;
}
int rbinom(double mu, double th){
	int k = gsl_ran_poisson (rng, mu*gsl_ran_gamma (rng, th, 1.0/th));
	return k;
}

// len(IPIV) N+1
// len(WORK) N*N
void solveLapack(double* A, double* b, int N, integer* ipiv, double* work)
{
    // solving A^T x = b
    integer lwork = 2*N*N;
    integer info;
    char uplo = 'L';
    integer n = N;
    integer lda = 10;
    integer ldb = 10;
    integer nrhs = 1;
    dsytrf_(&uplo, &n, A, &lda, (integer*)ipiv, work, &lwork, &info);
    //fprintf(stderr, "%ld ", info);
    dsytrs_(&uplo, &n, &nrhs, A, &lda, (integer*)ipiv, b, &ldb, &info);
    //fprintf(stderr, "%ld \n", info);
    //dgetri_(&N,A,&N,IPIV,WORK,&LWORK,&INFO);
}

void inverseLapack(double* A, int N, integer* ipiv, double* work){
    integer lwork = N*N;
    integer info;
    char uplo = 'L';
    integer n = N;
    integer lda = 10;
    //integer nrhs = 1;
    dsytrf_(&uplo, &n, A, &lda, (integer*)ipiv, work, &lwork, &info);
    dsytri_(&uplo, &n, A, &lda, (integer*)ipiv, work, &info);
}

void getStepSize5(double* H, double* g, double* step, integer* ipiv){
	cblas_dcopy(5, g, 1, step, 1);
	solveLapack(H, step, 5, ipiv, H+50);
}

void getStepSize6(double* H, double* g, double* step, integer* ipiv){
    cblas_dcopy(6, g, 1, step, 1);
    solveLapack(H, step, 6, ipiv, H+36);
}

double pchisq(double q, double k){// lower = F
	return  gsl_sf_gamma_inc_Q(k/2.0, q/2.0); 
}
double Qchisq(double p0, double k){
	double a=0.0;
	double b=10000.0;
	int i;
	if(pchisq(b,k)>p0){return(b);}
	for(i=0; i<100; i++){
		if(pchisq((a+b)/2.0,k)>p0){
			a=(a+b)/2.0;
		}else{
			b=(a+b)/2.0;
		}
		if(fabs(a-b)<1e-8){break;}
	}
	return (a+b)/2.0;
}

double inverse(double* a, double* inv){
    //
    //  a[0] a[1] a[2]
    //       a[3] a[4]
    //            a[5]
    //
    double det = a[0]*a[3]*a[5]+a[1]*a[4]*a[2]+a[2]*a[1]*a[4]-a[0]*a[4]*a[4]-a[2]*a[3]*a[2]-a[1]*a[1]*a[5];
    inv[0] = (a[3]*a[5] - a[4]*a[4])/det;
    inv[1] = (a[2]*a[4] - a[1]*a[5])/det;
    inv[2] = (a[1]*a[4] - a[2]*a[3])/det;
    inv[3] = (a[0]*a[5] - a[2]*a[2])/det;
    inv[4] = (a[2]*a[1] - a[0]*a[4])/det;
    inv[5] = (a[0]*a[3] - a[1]*a[1])/det;
    return det;
}

double getStepSize(double* H, double* g, int fixPi){
    //
    //  H[0] H[1] H[2]
    //       H[3] H[4]
    //            H[5]
    //
	double det;
	if(fixPi>0){
		det = H[3]*H[5]-H[4]*H[4];
		double s2 = ( H[5]*g[1]-H[4]*g[2])/det;
		double s3 = (-H[4]*g[1]+H[3]*g[2])/det;
		g[0]=0.0; g[1]=s2; g[2]=s3;
	}else{
		double* inv;
		inv = H+6;//(double*)calloc(6, sizeof(double));
		det = inverse(H, inv);
		double s1 = inv[0]*g[0] + inv[1]*g[1] + inv[2]*g[2];
		double s2 = inv[1]*g[0] + inv[3]*g[1] + inv[4]*g[2];
		double s3 = inv[2]*g[0] + inv[4]*g[1] + inv[5]*g[2];
		g[0]=s1; g[1]=s2; g[2]=s3;
		//free(inv);
	}
	return det;
}

double pf1(double fval, int df1, int df2){
	double a,b,x;
	a = ((double)df1)/2.0;
	b = ((double)df2)/2.0;
	x = ((double)df1)*fval/(((double)df1)*fval+(double)df2);
	return 1.0-gsl_sf_beta_inc (a, b, x);
}


void pf(double* fvals, int N, int df1, int df2){
	int i;
	for(i=0; i<N; i++){fvals[i]=pf1(fvals[i], df1, df2);}
}

double sign1(double g){
	if(g>=0.0){
		return 1.0;
	}else{
		return -1.0;
	}
}

void dcopy(long n, double* x, double* y){
    long i;
    for(i=0; i<n; i++){
        y[i]=x[i];
    }
}

long ncol(FILE* fp, long L){
	fseek(fp, 0L, SEEK_END);
	long sz = ftell(fp);
	fseek(fp, 0L, SEEK_SET);
    if(sz % (8*L) == 0){
        return sz/8/L;
    }else{
        return -1;
    }
}

long zncol(gzFile fp, long L){
    long sz=0;
    double tmp[1];
    //while(gzseek(fp, 16777216, SEEK_CUR)>0){ if(sz<0){return -1;}};
    while(gzread(fp, tmp, sizeof(double))>0){ if(sz<0){return -1;}};
	sz = gztell(fp);
    //printf("sz=%ld\n",sz/8);
	gzseek(fp, 0L, SEEK_SET);
    if(sz % (8*L) == 0){
        return sz/8/L;
    }else{
        return -1;
    }
}

long Fread(double* x, long sizeOf, long N, gzFile f){
    long pos=0;
    long l=0;
    //long i;
    l=gzread(f, x, sizeOf*N);
    pos=l/sizeOf;
    //for(i=0; i<N; i++){
    //    l=gzread(f, x+pos, sizeOf);
    //    if(l<0){return -1;}
    //    pos+=l/sizeOf;
    //}
    return pos;
}


/*
 void isna(double* v, long n, char tex[]){
 long i;
 for(i=0; i<n; i++){if(isnan(v[i])>0){fprintf(tex);printf("\n");break;}}	
 }
 */

long imax(double* x, long n){
    long i;
    long im=-1;
    double th = -(1.0e31);
    for(i=0; i<n; i++){
        if(x[i]>th && isnan(x[i])==0 && isinf(x[i])==0){
            im = i;
            th = x[i];
        }
    }
    return im;
}
/*double max(double* x, long n){
 long i;
 long im=-1;
 double th = -(1.0e31);
 for(i=0; i<n; i++){
 if(x[i]>th && isnan(x[i])==0 && isinf(x[i])==0){
 im = i;
 th = x[i];
 }
 }
 return th;
 }*/

double amax(double* x, long n){
    long i;
    double th = 0.0;
    for(i=0; i<n; i++){
        if(th<fabs(x[i])){
            th=fabs(x[i]);
        }
    }
    return th;
}

double doub(double x){
    return x*x;
}
/*
 double min(double a, double b){
 if(a>b){return b;}else{return a;}
 }
 */

long minL(long i, long j){
    if(i>j){return j;}else{return i;}
}


void scale(double* x, double* w, long L){
    long i;
    double m=0.0;
    double s=0.0;
    double n=0.0;
    for(i=0; i<L; i++){
        if(isnan(x[i])==0){
            m += w[i]*x[i];
            s += w[i]*x[i]*x[i];
            n += w[i];
        }
    }
    m /= n;
    s /= n;
    s -= m*m;
    s = sqrt(s);
    for(i=0; i<L; i++){
        if(isnan(x[i])==0){
            x[i] = (x[i]-m)/s;
        }else{
            x[i] = 0.0;
        }
    }
}

void progress(long j, long L){
    //long i;
    //double th=(double)L/100.0;
    //if(floor((double)(i+1)/th) > floor((double)i/th)){fprintf(stderr,"=");}
}

double lgamma(double x){
	return gsl_sf_lngamma(x);
}

double lbeta(double a, double b){
	return gsl_sf_lnbeta(a, b);
}

//double lgamma_a(double z){
//	return gsl_sf_lngamma(z);
//return z*log(z) - z -log(z/2.0/M_PI) + 1.0/12.0/z - 1.0/360.0/exp(3.0*log(z)) + 1.0/1260.0/exp(5.0*log(z));
//}

double digamma(double x){
	return gsl_sf_psi(x);
}

double trigamma(double x){
	return gsl_sf_psi_1(x);
}

double psi_2(double x){
	return gsl_sf_psi_n(2,x);
}

double ltrace(double* R, long n, long ldR){
    double res = 0.0;
    long i;
    for(i=0; i<n; i++){
        if(R[i+i*ldR]>0.0){res += log(R[i+i*ldR]);}
    }
    return res;
}


double lsum(double* x, long n){
    double res = 0.0;
    long i;
    for(i=0; i<n; i++){
        if(x[i]>0.0){res += log(x[i]);}
    }
    return res;
}

void WriteBin(double* x, long n, char* fname){
    FILE* fout;
    fout = fopen(fname, "wb");
    fwrite(x, sizeof(double), n, fout);
    fclose(fout);
}



void printInt(int* v, long n){
    long i;
    for(i=0; i<n; i++){fprintf(stderr,"%d ", v[i]);}
    fprintf(stderr, "\n");
}


void printLong(long* v, long n){
    long i;
    for(i=0; i<n; i++){fprintf(stderr,"%ld ", v[i]);}
    fprintf(stderr, "\n");
}

void print2(double* v, long n, long offs){
	long i;
	for(i=0; i<n; i++){fprintf(stderr, "%lf,", v[i*offs]);}
	fprintf(stderr, "\n");
}	

void print(double* v, long n){
	long i;
	for(i=0; i<n; i++){fprintf(stderr, "%lf,", v[i]);}
	fprintf(stderr, "\n");
}

void printSep(double* v, long n, char sep){
    long i;
    for(i=0; i<n; i++){fprintf(stderr, "%ld %lf%c", i, v[i], sep);}
    fprintf(stderr, "\n");
}

void lprint(double* v, long n){
    long i;
    for(i=0; i<n; i++){fprintf(stderr, "%lf,", log(v[i]));}
    fprintf(stderr, "\n");
}

void printL(double* m, long n, long p){
    long i,j,l=0;
    for(i=0; i<n; i++){
        for(j=0; j<p; j++){
            if(j>=i){
                fprintf(stderr, "%lf ", m[l++]);
            }else{
                fprintf(stderr, "*********** ");
            }
        }
        fprintf(stderr, "\n");
    }
}


void printM(double* m, long n, long p){
    long i,j;
    for(i=0; i<n; i++){
        for(j=0; j<p; j++){
            fprintf(stderr, "%lf ", m[i+j*n]);
        }
        fprintf(stderr, "\n");
    }
}

double ifelse(long a, long b, double* Beta){
	return ((a) > (b) ? (Beta[a]) : (0.0));
	//if(a>b){return 1.0;}else{return 0.0;}
	//return -1.0;
}

double nrm2(double* x, long n, long offset, long inc){
	return cblas_dnrm2(n, x+offset, inc);
}

double mean(double* x, long n){
    long i;
    double res=0.;
    double nef=0.;
    for(i=0; i<n; i++){
        if(isnan(x[i])==0){
            res += x[i];
            nef++;
        }
    }
    return res/nef;
}


double lmean(double* x, long n){
    long i;
    double res=0.;
    double nef=0.;
    for(i=0; i<n; i++){
        if(isnan(x[i])==0 && x[i]>0.0){
            res += log(x[i]);
            nef++;
        }
    }
    return res/nef;
}

double sd(double *x, long n){
    long i;
    double res=0.;
    double nef=0.;
    double m = mean(x, n);
    for(i=0; i<n; i++){
        if(isnan(x[i])==0){
            res += pow(x[i],2.0);
            nef++;
        }
    }
    return sqrt(res/nef-m*m);
}


double lsd(double *x, long n){
    long i;
    double res=0.;
    double nef=0.;
    double m = lmean(x, n);
    for(i=0; i<n; i++){
        if(isnan(x[i])==0 && x[i]>0.0){
            res += pow(log(x[i]),2.0);
            nef++;
        }
    }
    return sqrt(res/nef-m*m);
}

double wmean(double* x, double* w, long n){
    long i;
    double res=0.;
    double nef=0.;
    for(i=0; i<n; i++){
        if(isnan(x[i])==0){
            res += w[i]*x[i];
            nef += w[i];
        }
    }
    return res/nef;
}
double wvar(double* x, double* w, long n){
    long i;
    double res=0.;
    double nef=0.;
    double mm = wmean(x, w, n);
    for(i=0; i<n; i++){
        if(isnan(x[i])==0){
            res += w[i]*pow(x[i]-mm, 2.0);
            nef += w[i];
        }
    }
    return res/nef;
}

double asum(double *x, long n){
	return asum0(x, n, 0, 1);
}

double esum(double *x, long n, double y){
    long i;
    double res=0.0;
    for(i=0; i<n; i++){res += exp(x[i])*y;}
    return res;
}

double wsum(double* x, double* k, long n){
    long i;
    double res=0.0;
    for(i=0; i<n; i++){res += x[i]*k[i];}
    return res;
}

void wshift(double* x, double* k, long n){
    double mm = wsum(x, k, n);
    long i;
    for(i=0; i<n; i++){ x[i] -= mm; }
}

double iwsum(double* x, double* k, long n){
    long i;
    double res=0.0;
    for(i=0; i<n; i++){res += x[i]/k[i];}
    return res;
}

double iwlsum(long* x, double* k, long n){
    long i;
    double res=0.0;
    for(i=0; i<n; i++){res += ((double)(x[i*2]+x[i*2+1]))/k[i];}
    return res;
}

double sum(double *x, long n){
    long i;
    double res=0.0;
    for(i=0; i<n; i++){res += x[i];}
    return res;
}

long sumLong(long* x, long n){
    long i;
    long res=0;
    for(i=0; i<n; i++){res += x[i];}
    return res;
}

double asum0(double* x, long n, long offset, long inc){
    long i;
    double res=0.0;
    for(i=0; i<n; i++){res += fabs(x[i*inc+offset]);}
    return res;
	//return cblas_dasum(n, x+offset, inc);
}

double ddot(double* x, double* y, long n, long offset, long inc){
	return cblas_ddot(n, x+offset, inc, y+offset, inc);
}

void clear(double* v, long N){
    long i;
    for(i=0; i<N; i++){v[i]=0.0;}
    //cblas_dscal(N, 0.0, v, 1);
}

void clear1(double* v, long N){
	long i;
	for(i=0; i<N; i++){v[i]=0.0;}
	//cblas_dscal(N, 0.0, v, 1);
}

void clearLong(long* v, long N){
	long i;
	for(i=0; i<N; i++){v[i] = 0;}
}

void clearInt(int* v, int N){
	int i;
    for(i=0; i<N; i++){v[i] = 0;}
}


void apply1(double* x, long n, long p, double* y){
	long i,j;
	for(i=0; i<n; i++){
		for(j=0; j<p; j++){
			y[i] += x[i*n+j];
		}
	}
}


void apply2(double* x, long n, long p, double* y){
	long i,j;
	for(i=0; i<n; i++){
		for(j=0; j<p; j++){
			y[j] += x[i*n+j];
		}
	}
}

double logit(double x){
    if(x<=1.0&&x>=0.0){
        return log(x/(1-x));
    }else if(x<=0.0){
        return -30.0;
    }else if(x>=1.0){
        return 30.0;
    }
    return 0.0;
}

double expit(double x){
    if(x>0){
    	return 1.0/(1.0+exp(-x));
    }else{
        return exp(x)/(1.0+exp(x));
    }
}

void BackSolveOld(double* R, double* y, double* x, long p){
    // To solve R x = y
	long i,j;
	x[p-1] = y[p-1]/R[p*p-1];
	for(i=p-2; i>=0; i--){
		double tm = 0.0;
		for(j=p-1; j>i; j--){
			tm += R[p*j+i]*x[j];
		}
		x[i] = (y[i]-tm)/R[i*p+i];
	}
}

void QRold(double* X, long n, long p, double* R){
    // QR decomp using Gram-Schmidt orthogonalization
    // X is replaced by Q
    
    long i, j;
    
    clear1(R, p*p);
    R[0] = cblas_dnrm2(n, X, 1);
    cblas_dscal(n, 1.0/R[0], X, 1);
    
    for(i=1; i<p; i++){
        for(j=0; j<i; j++){
            R[i*p+j] = cblas_ddot(n, X+n*j, 1, X+n*i, 1);
            cblas_daxpy(n, -R[i*p+j], X+n*j, 1, X+n*i, 1);
        }
        R[i*p+i] = cblas_dnrm2(n, X+n*i, 1);
        cblas_dscal(n, 1.0/R[i*p+i], X+n*i, 1);
    }
}


void sum2one(double* x, int n){
	double tot=0.0;
	int i;
	for(i=0;i<n;i++){
		tot+=x[i];
	}
	for(i=0;i<n;i++){
		x[i]/=tot;
	}
}






int rD(double* zc, double* zx){
	int zc1 = gsl_ran_binomial(rng, zc[0], 1);
	int zc2 = gsl_ran_binomial(rng, zc[1], 1);
	int zx1 = gsl_ran_binomial(rng, zx[0], 1);
	int zx2 = gsl_ran_binomial(rng, zx[1], 1);
	if((zc1+zc2)==0){
		if((zx1+zx2)==0){
			return 0;
		}else if((zx1+zx2)==1){
			return 1;
		}else{
			return 2;
		}
	}else if((zc1+zc2)==1){
		if((zx1+zx2)==0){
			return 3;
		}else if((zx1+zx2)==1){
			if(zc1==zx2){
				return 4;
			}else{
				return 5;
			}
		}else{
			return 6;
		}
	}else{
		if((zx1+zx2)==0){
			return 7;
		}else if((zx1+zx2)==1){
			return 8;
		}else{
			return 9;
		}
	}
}


