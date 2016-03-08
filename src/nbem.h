#include <pthread.h>

#define GZ_WMODE "wb6f"
#define GZ_RMODE "rb6f"

#define MAXITR 50
#define MAXITR2 30
#define MAXITR3 30
#define PRESTEPS 1

pthread_mutex_t mutex;
typedef struct{
	double* y;
	double* Y;
	double* Z;
	double* X;
	double* ki;
	double* ki2;
	double* w;
	double* km;
	double* exon;
	long P;
	long L;
	long N;
	long Lx;
	char** rss;
	double* lkhdDiff;
	double* ppi;
	double* pdelta; 
	double* pphi;
	double* pbeta;
	double* ptheta; 
	double* pasr; 
	int* pitr;
	int* pbound; 
	double* ptval;
	double* pkld;
	long rsnp_start;
	long rsnp_end;
	int* tested;
	int tid;
} RASQUAL_INP;


long maxCsnp;
long nthreads;

int verbose;
int verbose2;
int verbose3;

int ASE;
int transQTL;
int noPosteriorUpdate;
int noGL;
int allelicProbEstByErrorRate;
long Null;
int randomize;
double NOfSigLoci;
int* fixParam;
int hetType;
int testImprinting;

double theta2;

long numOfLoci;

// hyper parameters
double ab;//=1.501;
double phiab;//=1.501;
double ad;// = 1.1;
double bd;// = 57.0;
double kappa;// = .0700000613259765*2.;
double omega;// = 0.002000092291342*2.;
double sigma2;

double* ones;

double asNonasRatio0;

double beta0;
double theta0;

void* ASEQTLALL_MP(void *args);
long ASEQTLALL(double* y, double* Y, double* Z, double* X, long P, double* ki, double* ki2, double* w, double* km, double* exon, long L, long N, long M, char** rss, double* lkhdDiff, double* ppi, double* pdelta, double* pphi, double* pbeta, double* ptheta, double* pasr, int* pitr, int* pbound, double* ptval, double* pkld);
long ASEQTLALL_ALT(double* y, double* Y, double* Z, double* X, long P, double* ki, double* ki2, double* w, double* km, double* exon, long L, long N, long Lx, char** rss, double* lkhdDiff, double* ppi, double* pdelta, double* pphi, double* pbeta, double* ptheta, double* pasr, int* pitr, int* pbound, double* ptval, double* pkld, long csnp_start, long csnp_end, int* tested, int tid);

void getInformation(double* y, double* Y, double* h0, double* H0, double* h1, double* H1, double* H2, double* ki, double* dki, double* ki2, double* K0, double* K2, double* K, double* km, long Lx, long N0, long J0, long J, double beta, double th, double pi, double delta, double phi, double asr, double* work, double* hess, integer* ipiv, double* a, double* A);

double ASEQTL(double* y, double* Y, double* h, double* h0, double* ki, double* dki, double* ki2, double* km, long L, long N, long M, long csnp, double* work, double* ppi, double* pdelta, double* pphi, double* pbeta, double* ptheta, double* pasr, int* pitr, int* pbound, double* tval, double* pkld);

double getLkhdAll(double* y, double* Y, double* h0, double* h, double* H0, double* H, double* ki, double* dki, double* ki2, double* kj0, double* kj, double* ktmp, double* km, long L, long N0, long J0, long J, double beta, double* th, double* pi, double delta, double phi, double asr, double* cls, double* hYg, double* a, double* A, double* H2, double* grad, double* hess);

double getQval(double* y, double* Y, double* h1, double* H1, double* h0, double* H0, double* H2, double* ki, double* dki, double* ki2, double* kj0, double* kj, double* ktmp, double* km, long L, long N0, long J0, long J, double* th, double beta, double* pi, double delta, double phi, double asr);

long nbbem_init0(double* y, double* ki, long N, double* pbeta, double* ptheta);

long getPrior(double* h0, double* z1, double* z2, double wi);
long getPrior1(double* h0, double* z, double wi);
long getPriorNull(double* h0, double* z, double wi);
long getPriorIP(double* h0, double* z2, double wi);
void getK0s(double pi, double* K0);
void getK2s(double pi, double delta, double phi, double* ktmp, double* w, double* K2);
void getK(double pi, double delta, double phi, double* K, double* w, double* K2);
void getdKdPi(double pi, double delta, double phi, double* K, double* w, double* K2);
void getd2KdPi2(double pi, double delta, double phi, double* K, double* w, double* K2);
void getdKdDelta(double pi, double delta, double phi, double* K, double* w, double* K2);
void getd2KdDelta2(double pi, double delta, double phi, double* K, double* w, double* K2);
void getdKdPhi(double pi, double delta, double phi, double* K, double* w, double* K2);
void getd2KdPhi2(double pi, double delta, double phi, double* K, double* w, double* K2);

double getGradHess(double* y, double* Y, double* h0, double* h, double* H0, double* H, double* ki, double* dki, double* ki2, double* K0, double* K2, double* ktmp, double* km, long Lx, long N0, long J0, long J, double beta, double* th, double pi, double delta, double phi, double asr, double* cls, double* pYg, double* a, double* A, double* H2, double* grad, double* hess);
double getGrad(double* A, double* K2, double* K2p);
double getHess1(double* A, double* K2, double* K2p, double* K2pp);
double getHess2(double* A, double* K2, double* K2p, double* K2d, double* K2pd);
double getHess3(double* A, double* K2, double* K2p);


void getd2KdPidDelta(double pi, double delta, double phi, double* K, double* w, double* K2);
void getd2KdPidPhi(double pi, double delta, double phi, double* K, double* w, double* K2);
void getd2KdDeltadPhi(double pi, double delta, double phi, double* K, double* w, double* K2);

void expand(double* s1, double* s2, long N, double* K, double* w, double* K2);
double vsum(double* x, double* y, long N);
double vsum2(double* x, double* y, long N);
double vsum3(double* x, double* y, double* z, long J);
double vsum4(double* x, double* y1, double* y2, long J, long M);
double vsum5(double* x, double* y, long N);

double getCE(double* A, double* H1, double* h1, long j0);
double getCE2(double* K2, double* K2p, double* A, double* H1, double* h1, long j0);
double getCV(double* A, double* H1, double* h1, long j0);
double getCV2(double* K2, double* K2p, double* A, double* H1, double* h1, long j0);
double getCCV(double* K2, double* K2p, double* A, double* H1, double* h1, long j0);
double getCCV2(double* K2, double* K2p, double* K2d, double* A, double* H1, double* h1, long j0);

double getIA1(double* gl, long N);
double getIA2(double* gl, long N);

void gencall(double* H1, double* Z, double* Zx, long N, long L, double* exon);
void gencallAlt(double* H1, double* h1, double* Z, double* Zx, long N, long L, long Lx, double* exon, long csnp);

void H3toZ(double* h3, double* z);
void H3toZx(double* H1, double* zc, double* zx, double* work);
void H3toZx2(double* H1, double* zc, double* zx, double* work);
double getCov(double* h0, double* h1, long N0);
double getCov2(double* h0, double* h1, long N0);

double getAA(double* gl);
double getAR(double* gl);

void printMapGen(double* z, long N, long L);

void randomPerm(double* y, double* Y, double* Z, double* ki, double* ki2, double* exon, int L, int N, int m);
