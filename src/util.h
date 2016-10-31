#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
//#include <cblas.h>
#include <blaswrap.h>
#include <gsl_cblas.h>
#include <gsl_rng.h>
#include <gsl_randist.h>
#include <f2c.h>
#include <clapack.h>
#include <gsl_sf_gamma.h>
#include <gsl_sf_psi.h>
#include <zlib.h>
#include <time.h>

//#include <R.h>


const gsl_rng_type * rngT;
gsl_rng * rng;
void init_gsl_rand();
double runif();
void rdirmultinom(int* ys, int m, int Y, double* probs, double* alpha);
void rmultinom(int* ys, int m, int Y, double* probs);
int rbinom(double mu, double th);
int rbetabinom(double thA, double thB, int n);
void rG(double* z, double* zt);
int rD(double* zc, double* zx);


double pchisq(double q, double k);
double Qchisq(double p0, double k);


double inverse(double* a, double* inv);
double getStepSize(double* H, double* g, int fixPi);
void inverseLapack(double* A, int N, integer* ipiv, double* work);
void solveLapack(double* A, double* b, int N, integer* ipiv, double* work);
void getStepSize5(double* H, double* g, double* step, integer* ipiv);
void getStepSize6(double* H, double* g, double* step, integer* ipiv);

double pf1(double fval, int df1, int df2);
void pf(double* fvals, int N, int df1, int df2);

double sign1(double g);
void dcopy(long n, double* x, double* y);

long zncol(gzFile fp, long L);
long Fread(double* x, long sizeOf, long N, gzFile f);

long imax(double* x, long n);
//double max(double* x, long n);
double doub(double x);

double esum(double *x, long n, double y);

//double lgamma_a(double z);
double lgamma(double x);
double lbeta(double a, double b);
double digamma(double x);
double trigamma(double x);
double psi_2(double x);

long ncol(FILE* fp, long L);
void scale(double* x, double* v, long L);

void progress(long j, long L);
void isna(double* v, long n, char tex[]);
void print(double* v, long n);
void printSep(double* v, long n, char sep);
void lprint(double* v, long n);
void print2(double* v, long n, long offs);
void printLong(long* v, long n);
void printInt(int* v, long n);
void printM(double* m, long n, long p);
void printL(double* m, long n, long p);
void clear(double* v, long N);
void clear1(double* v, long N);
void clearInt(int* v, int N);
void clearLong(long* v, long N);
void apply1(double* x, long n, long p, double* y);
void apply2(double* x, long n, long p, double* y);
void BackSolveOld(double* R, double* y, double* x, long p);
void QRold(double* X, long n, long p, double* R);
double ifelse(long a, long b, double* Beta);
double nrm2(double* x, long n, long offset, long inc);
double asum0(double* x, long n, long offset, long inc);
double asum(double* x, long n);
double sum(double *x, long n);
void wshift(double *x, double *k, long n);
double wsum(double *x, double *k, long n);
double iwsum(double *x, double *k, long n);
double iwlsum(long *x, double *k, long n);
long sumLong(long *x, long n);
double ddot(double* x, double* y, long n, long offset, long inc);

double ltrace(double* R, long n, long ldR);
double lsum(double* x, long n);

void WriteBin(double* x, long n, char* fname);

double logit(double x);
double expit(double x);


long minL(long i, long j);
//double min(double a, double b);

double mean(double* x, long n);
double lmean(double* x, long n);
double lsd(double* x, long n);
double sd(double* x, long n);
double wmean(double* x, double* w, long n);
double wvar(double* x, double* w, long n);

void sum2one(double* x, int n);
