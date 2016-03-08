#include "util.h"
#include "nbem.h"
#include "nbglm.h"
#include "sort.h"

#define upper 1.0
#define lower 0.0

int ig[10] = {0,0,0,1,1,1,1,2,2,2};
int jg[3]  = {0, 3, 7};
int kg[3]  = {3, 4, 3};

double d0(double delta){return (delta-lower)/(upper-lower);}




long asemethod = 0;

void setASEMethod(int* i){
    asemethod = (long)(*i);
}

double oasr = 2.0;

double lambda=0.01;

double dipc[20] = { 0.,0.,0.,0.,0.,0., 1.,0.,1.,0.,1.,0.,1.,0., 1.,1.,1.,1.,1.,1.};
double dipx[20] = { 0.,0.,0.,1.,1.,1., 0.,0.,0.,1.,1.,0.,1.,1., 0.,0.,0.,1.,1.,1.};


// y : N
// X : N x P
// Y : 2 x N x Lx
// Z : 2 x N x L
long ASEQTLALL(double* y, double* Y, double* Z, double* X, long P, double* ki, double* ki2, double* w, double* km, double* exon, long L, long N, long Lx, char** rss, double* lkhdDiff, double* ppi, double* pdelta, double* pphi, double* pbeta, double* ptheta, double* pasr, int* pitr, int* pbound, double* ptval, double* pkld){
    
    double* work=NULL;
    //double* z;
    double* z2;
    double* h0;
    double* H0;
    
    
    
    double tval;
    
    //double maxLR=-1.0e31;
    //maxCsnp=-1;
    
    long i, l;
    long top=0;
    long ll;
    h0=(double*)calloc(N*3*2,sizeof(double));
    H0=(double*)calloc(N*10*Lx*2,sizeof(double));
    long lwork = 78+N*Lx*10+3*N+300+2*(Lx+1)+53+Lx+3*Lx+N*10;// + 2*(N+4*P)*(P+4);//(N+P)*P + P*P + 2*(N+P) + 2*N + 4*P+;
    work=(double*)calloc(lwork,sizeof(double));
    if(work==NULL){fprintf(stderr, "memory allocation failur (work space is not allocated)\n"); return -1;}
    
    
    
    
    
    ///double* Zx;// posterior allelic prob
    //Zx = Z+L*N*2;
    //double* H1null;
    double* H1;
    double* h1;
    //H1null = (double*)calloc(N*Lx*10, sizeof(double));
    h1 = h0 + N*3;
    H1 = H0 + N*10*Lx;
    
    ppi[L]=pphi[L]=0.5;
    pdelta[L]=(ad-1.0)/(ad+bd-2.0);
    
    //double lkhdNull;
    double lkhdR;
    int itr_rand;
    double pitmp, deltatmp, phitmp, betatmp, thetatmp;
    double kldtmp[2];
    int itrtmp, boundtmp;
    
    double deltaR[5] = {(ad-1.0)/(ad+bd-2.0), (ad-1.0)/(ad+bd-2.0), 0.5, 0.5, 0.5};
    double phiR[5]   = {0.1, 0.9, 0.1, 0.9, 0.5};
    
    
    //lkhdNull = -1.0/0.0;
    
    double nasRat = 0.0;
    
    
    double* dki;  dki = ki+N+randomize*Lx*N;
    double* zr;
    double* Zr=Z+N*(L+Lx+1)*2;
    
    
    if(fixParam[4]==1){theta0=ptheta[L]=(kappa-2.0)/omega;}
    //###################
    //#  Null without AS
    //###################
    for(i=0;i<N;i++){h0[i*3+1] = h1[i*3+1] = 0.5*w[i];   h0[i*3+2] = h1[i*3+2] = h0[i*3+0] = h1[i*3+0] = 0.25*w[i];}
    fixParam[0]=1;
    if(verbose3>0)fprintf(stderr, "Grand Null w/o AS ");
    
    double lkhdNull0 = lkhdDiff[L] = lkhdDiff[L+2] = ASEQTL(y, Y, h0, NULL, ki, dki, ki2, km, L, N, 0, -1, work, ppi+L, pdelta+L, pphi+L, pbeta+L, ptheta+L, &nasRat, pitr+L, pbound+L, &tval, pkld+2*L);
    for(l=0;l<L;l++){ppi[l]=0.5; pdelta[l]=(ad-1.0)/(ad+bd-2.0); pphi[l]=0.5; pbeta[l]=pbeta[L]; ptheta[l]=ptheta[L];}
    
    
    //###################
    //#  Null with AS
    //###################
    double lkhdNull0as = -1.0e20;
    if(Lx>0){
        ll=0;
        
        for(l=0;l<L;l++){
            z2=Z+N*2*l;
            if(fabs(exon[l])>1.5){
                if(randomize>0){
                    zr = Zr+N*2*ll;
                    for(i=0;i<N;i++){ getPriorNull(H0+10*i+N*10*ll, zr+i*2, w[i]); }
                    for(i=0;i<N;i++){ getPriorNull(H1+10*i+N*10*ll, zr+i*2, w[i]);}
                }else{
                    for(i=0;i<N;i++){ getPriorNull(H0+10*i+N*10*ll, z2+i*2, w[i]); }
                    for(i=0;i<N;i++){ getPriorNull(H1+10*i+N*10*ll, z2+i*2, w[i]);}
                }
                ll++;
            }
        }
        ppi[L]    = 0.5;
        pphi[L]   = 0.5;
        pdelta[L] = (ad-1.0)/(ad+bd-2.0);
        pasr[L]   = asNonasRatio0;
        if(verbose3>0)fprintf(stderr, "Grand Null with AS ");
        lkhdNull0as  = lkhdDiff[L] = lkhdDiff[L+3] = ASEQTL(y, Y, h0, H0, ki, dki, ki2, km, L, N, Lx, -1, work, ppi+L, pdelta+L, pphi+L, pbeta+L, ptheta+L, pasr+L, pitr+L, pbound+L, &tval, pkld+2*L); 
        
        //cblas_dcopy(N*Lx*10, H1, 1, H1null, 1);
        //gencall(H1null, Z, Zx, N, L, exon);
        
        // multi starting
        for(itr_rand=0; itr_rand<5; itr_rand++){
            pitmp    = 0.5;
            deltatmp = fixParam[1]==1 ? pdelta[L] : deltaR[itr_rand];
            phitmp   = fixParam[2]==1 ? 0.5       : phiR[itr_rand];
            betatmp  = pbeta[L];
            thetatmp = ptheta[L];
            
            if(verbose3>0)fprintf(stderr, "Grand Null with AS ");
            lkhdR  = ASEQTL(y, Y, h0, H0, ki, dki, ki2, km, L, N, Lx, -1, work, &pitmp, &deltatmp, &phitmp, &betatmp, &thetatmp, pasr+L, &itrtmp, &boundtmp, &tval, kldtmp);
            if(lkhdR>lkhdNull0as && boundtmp==0){
                lkhdNull0as = lkhdDiff[L] = lkhdDiff[L+3] = lkhdR;
                ppi[L]    = 0.5;
                pdelta[L] = deltatmp;
                pphi[L]   = phitmp;
                pbeta[L]  = betatmp;
                ptheta[L] = thetatmp;
                
                pitr[L]   = itrtmp;
                pbound[L] = boundtmp;
                pkld[2*L]   = kldtmp[0];
                pkld[2*L+1] = kldtmp[0];
                
                //cblas_dcopy(N*Lx*10, H1, 1, H1null, 1);
                //gencall(H1null, Z, Zx, N, L, exon);
            }
        }
        // multi starting end
    }
    //lkhdDiff[L]= 0.0*2.0*(lkhdNull0as);
    if(verbose3>1){fprintf(stderr, "L0=%lf L0AS=%lf\n", lkhdNull0, lkhdNull0as);}
    
    // alt here
    
    //##############
    // Imprinting
    //##############
    if(Lx>0 && testImprinting){
        ll=0;
        for(i=0;i<N;i++){h0[i*3+1] = 1.0*w[i];   h0[i*3+2] = h0[i*3+0] = 0.0*w[i];}
        for(l=0;l<L;l++){
            z2=Z+N*2*l;
            if(fabs(exon[l])>1.5){
                if(randomize>0){
                    zr = Zr+N*2*ll;
                    for(i=0;i<N;i++){ getPriorIP(H0+10*i+N*10*ll, zr+i*2, w[i]); }// xSNP = cSNP
                }else{
                    for(i=0;i<N;i++){ getPriorIP(H0+10*i+N*10*ll, z2+i*2, w[i]); }// xSNP = cSNP
                }
                ll++;
            }
        }
        pasr[L+1]   = asNonasRatio0;
        
        ppi[L+1]    = 0.9;
        pdelta[L+1] = pdelta[L]; //(ad-1.0)/(ad+bd-2.0);
        pphi[L+1]   = pphi[L]; //0.5;
        pbeta[L+1]  = pbeta[L];
        ptheta[L+1] = ptheta[L];
        fixParam[0] = 0;
        if(verbose3>0)fprintf(stderr, "Imprinting         ");
        lkhdDiff[L+1]=2.0*ASEQTL(y, Y, h0, H0, ki, dki, ki2, km, L, N, Lx, -1, work, ppi+L+1, pdelta+L+1, pphi+L+1, pbeta+L+1, ptheta+L+1, pasr+L+1, pitr+L+1, pbound+L+1, &tval, pkld+2*(L+1))-2*lkhdNull0as;
    }
    
    
    
    //fprintf(stderr, "lkhd00000=%lf\n", lkhdDiff[L]);
    free(work);
    free(h0);
    //if(Lx>0)free(H1null);
    if(Lx>0)free(H0);
    return top;
}


void* ASEQTLALL_MP(void *args){
    RASQUAL_INP ri = *(RASQUAL_INP *) args;
    //fprintf(stderr, "ASEQTLALL_MP: init %ld %ld\n", ri->rsnp_start, ri->rsnp_end);
    ASEQTLALL_ALT(ri.y, ri.Y, ri.Z, ri.X, ri.P, ri.ki, ri.ki2, ri.w, ri.km, ri.exon, ri.L, ri.N, ri.Lx, ri.rss, ri.lkhdDiff, 
                  ri.ppi, ri.pdelta, ri.pphi, ri.pbeta, ri.ptheta, ri.pasr, ri.pitr, ri.pbound, ri.ptval, ri.pkld, ri.rsnp_start, ri.rsnp_end, ri.tested, ri.tid);
    //pthread_exit(NULL);
    return args;
}

long ASEQTLALL_ALT(double* y, double* Y, double* Z, double* X, long P, double* ki, double* ki2, double* w, double* km, double* exon, long L, long N, long Lx, char** rss, double* lkhdDiff, double* ppi, double* pdelta, double* pphi, double* pbeta, double* ptheta, double* pasr, int* pitr, int* pbound, double* ptval, double* pkld, long csnp_start, long csnp_end, int* tested, int tid){
    
    double* z;
    double* z2;
    
    //double maxLR;
    double tval;
    
    //if(verbose3>0)fprintf(stderr, "ASEQTLALL_ALT: rSNP [%ld, %ld]\n", csnp_start, csnp_end);
    
    double lkhdNull0=lkhdDiff[L+2], lkhdNull0as=lkhdDiff[L+3];
    
    int i, l, ll;
    
    
    double* work=NULL;
    long lwork = 78+N*Lx*10+3*N+300+2*(Lx+1)+53+Lx+3*Lx+N*10;// + 2*(N+4*P)*(P+4);//(N+P)*P + P*P + 2*(N+P) + 2*N + 4*P+;
    work=(double*)calloc(lwork, sizeof(double));
    
    double maxLR = -1.0e-20;
    double* Zx = Z+N*L*2;
    
    
    //double* Zx;// posterior allelic prob
    //Zx = Z+L*N*2; length: 2*N*(Lx+1)
    double* h0=NULL;
    double* H0=NULL;
    double* h1;
    double* H1;
    //h0=work + 78+N*Lx*10+3*N+300+2*(Lx+1)+53+Lx+3*Lx+N*10; //
    //H0=work + 78+N*Lx*10+3*N+300+2*(Lx+1)+53+Lx+3*Lx+N*10 + N*6; //
    h0=(double*)calloc(N*3*2,sizeof(double));
    H0=(double*)calloc(N*10*Lx*2,sizeof(double));
    h1 = h0 + N*3;
    H1 = H0 + N*10*Lx;
    
    if(work==NULL || h0==NULL || H0==NULL){fprintf(stderr, "memory allocation failur (work space is not allocated)\n"); return -1;} //else{fprintf(stderr, "ASEQTLALL_ALT: memory alloc\n");}
    
    double* dki;  dki = ki+N+randomize*Lx*N;
    double* zr;
    double* Zr=Z+N*(L+Lx+1)*2;
    
    
    
    if(verbose3>0){fprintf(stderr, "\n\nAlternative Start\n\n");}
    //#################
    //#   Alternative
    //#################
    fixParam[0]=0;
    long csnp;
    //csnp_start = (maxCsnp<0 ? 0 : maxCsnp);
    //csnp_end   = (maxCsnp<0 ? L : maxCsnp+1);
    //fprintf(stderr, "ASEQTL start=%ld end%ld\n", csnp_start, csnp_end);
    for(csnp=csnp_start; csnp<csnp_end; csnp++){
        if(exon[csnp]>0.5){
            pthread_mutex_lock( &mutex );
            if(tested[csnp]==0){
                tested[csnp] = tid;
            }
            pthread_mutex_unlock( &mutex );
            if(tested[csnp] == tid){
                if(verbose3>1)fprintf(stderr, "rSNP %ld is taken by thread %d\n", csnp, tid);
                clear1(work, lwork);
                z=Z+N*2*csnp;
                for(i=0;i<N;i++){
                    h0[i*3]   = h1[i*3]   = w[i] * (1.0-z[i*2])*(1.0-z[i*2+1]);
                    h0[i*3+2] = h1[i*3+2] = w[i] *       z[i*2]*z[i*2+1];
                    h0[i*3+1] = h1[i*3+1] = w[i] * (1.0-h0[i*3]-h0[i*3+2]);
                }
                pasr[csnp]=0.0;
                if(verbose3>0)fprintf(stderr, "%ld\t%s %lf ", csnp, rss[csnp], ptheta[csnp]);
                lkhdDiff[csnp]=2.0*ASEQTL(y, Y, h0, H0, ki, dki, ki2, km, L, N, 0, csnp, work, ppi+csnp, pdelta+csnp, pphi+csnp, pbeta+csnp, ptheta+csnp, 
                                          pasr+csnp, pitr+csnp, pbound+csnp, &tval, pkld+csnp*2)-2.0*lkhdNull0;

                // Test with AS
                if(Lx>0 && (lkhdDiff[csnp]>Qchisq(min(1.0, NOfSigLoci/(double)numOfLoci), 1.0) || ASE==2)){
                    if(verbose3>10){fprintf(stderr, "\n\nAlternative with AS\n\n");}
                    // init H0
                    ll=0;
                    
                    for(l=0;l<L;l++){
                        z2=Z+N*2*l;  // feature snp
                        if(fabs(exon[l])>1.5){
                            if(randomize>0){
                                zr = Zr+N*2*ll;
                                for(i=0;i<N;i++){ getPrior( H0+10*i+N*10*ll, z+i*2, zr+i*2,  w[i]); }// xSNP = cSNP
                                for(i=0;i<N;i++){ getPrior( H1+10*i+N*10*ll, z+i*2, zr+i*2,  w[i]); }// xSNP = cSNP
                            }else{
                                if(l==csnp){// p(G|G)=1
                                    for(i=0;i<N;i++){ getPrior1(H0+10*i+N*10*ll, z+i*2, w[i]); }// xSNP = cSNP
                                    for(i=0;i<N;i++){ getPrior1(H1+10*i+N*10*ll, z+i*2, w[i]); }// xSNP = cSNP
                                }else{
                                    for(i=0;i<N;i++){ getPrior( H0+10*i+N*10*ll, z+i*2, z2+i*2, w[i]); }
                                    for(i=0;i<N;i++){ getPrior( H1+10*i+N*10*ll, z+i*2, z2+i*2, w[i]); }
                                }
                            }
                            ll++;
                        }
                    }
                    //ppi[csnp+1]    = 0.49999;
                    pdelta[csnp] = pdelta[L];
                    pphi[csnp]   = pphi[L];
                    pbeta[csnp]  = pbeta[L];
                    ptheta[csnp] = ptheta[L];
                    pasr[csnp]   = asNonasRatio0;
                    if(verbose3>0)fprintf(stderr, "%ld\t%s ", csnp, rss[csnp]);
                    
                    lkhdDiff[csnp]=2.0*ASEQTL(y, Y, h0, H0, ki, dki, ki2, km, L, N, Lx, csnp, work, ppi+csnp, pdelta+csnp, pphi+csnp, pbeta+csnp, ptheta+csnp, pasr+csnp, pitr+csnp, pbound+csnp, &tval, pkld+csnp*2)-2.0*lkhdNull0as;
                    
                    
                }
                if(lkhdDiff[csnp]>maxLR && nthreads==1){
                    maxLR=lkhdDiff[csnp];
                    maxCsnp = csnp;
                    gencallAlt(H1, h1, Z, Zx, N, L, Lx, exon, csnp);
                }
                ptval[csnp] = tval;
            }
        }
    }
    
    
    
    free(work);
    free(h0);
    if(Lx>0)free(H0);
    
    return 0;
}


int isConv(double* grad, int n, double th){
    int res=0;
    int geta = 1;
    int i;
    for(i=n-1; i>=0; i--){
        res += fabs(grad[i])<th ? 0 : geta*(1-fixParam[i]);
        geta *= 10;
    }
    return res;
}





/*
 
 
 pi    0 5 10 15 20
 delta 1 6 11 16 21
 phi   2 7 12 17 22
 beta  3 8 13 18 23
 theta 4 9 14 19 24
 
 
 */

// h0: genotype likelihood at csnp
// H0: diplotype likelihood between csnp and xsnp
// h1: posterior probability given fragment count
// L : num of snps
// N : sample size
// M : num of xSNP
double ASEQTL(double* y, double* Y, double* h0, double* H0, double* ki, double* dki, double* ki2, double* km, long L, long N0, long Lx, long csnp, double* work, double* ppi, double* pdelta, double* pphi, double* pbeta, double* ptheta, double* pasr, int* pitr, int* pbound, double* tval, double* pkld){
    
    integer* ipiv=NULL;
    ipiv = (integer*)calloc(11, sizeof(integer));
    if(ipiv==NULL){fprintf(stderr, "ASEQTL: memory allocation error.\n");}
    
    
    
    long itr, itr_step, j;
    double* K;
    double* K2;
    double* K0;
    //double* K0p;
    //double* K0pp;
    double* K2p;
    double* K2d;
    double* K2h;
    double* K2pd;
    double* K2dh;
    double* K2ph;
    double* K2pp;
    double* K2dd;
    double* K2hh;
    
    int nofp = 5;
    
    
    double* a;
    double* A;
    double* h1;
    double* H1;
    double* H2;
    double* cls;
    double* pYg;
    
    double* grad;
    double* hess;
    double* step;
    
    
    double ssize=1.0;
    double beta1, phi1, asr1, delta1;
    double beta,  phi,  asr, delta;
    
    double pi[2];
    double pi1[2];
    double theta[2];
    double theta1[2];
    //double delta[2];
    //double delta1[2];
    if(transQTL==0){
        pi[0] = pi[1] = pi1[0] = pi1[1] = ppi[0];
    }else{
        pi[0] = pi1[0] = ppi[0]; pi[1] = pi1[1] = 0.5;
    }
    delta = pdelta[0];//ad/(ad+bd);
    phi   = pphi[0];//0.5;
    theta[0] = theta[1] = ptheta[0];
    beta  = pbeta[0];
    asr = asr1 = pasr[0];
    
    
    
    double qval, qvalNext;
    
    
    
    //integer ipiv[11];
    
    //asr=asr1=0.1;
    double lkhd;
    double lkhd_old;
    double lkhdDiff;
    
    h1 = h0+N0*3;
    H1 = H0+N0*Lx*10;
    
    K     =work;
    A     =work+8;
    //h1    =work+58;
    
    K2    =work+58+3*N0;
    K2p   =work+78+3*N0;
    K2d   =work+78+3*N0+20;
    K2h   =work+78+3*N0+40;
    K2pp  =work+78+3*N0+60;
    K2pd  =work+78+3*N0+80;
    K2ph  =work+78+3*N0+100;
    K2dd  =work+78+3*N0+120;
    K2dh  =work+78+3*N0+140;
    K2hh  =work+78+3*N0+160;
    
    
    grad = work+78+3*N0+180;
    step = work+78+3*N0+190;
    hess = work+78+3*N0+200; // length 50 = 25 * 2 for hess and inv hess
    //H1   = work+78+3*N0+300+2*(Lx+1);
    
    K0  = work+78+3*N0+300+2*(Lx+1)+N0*Lx*10; // K0p and K0pp should be allocated
    //K0p = work+78+3*N0+300+2*(Lx+1)+N0*Lx*10+3;
    //K0pp= work+78+3*N0+300+2*(Lx+1)+N0*Lx*10+6;
    a   = work+78+3*N0+300+2*(Lx+1)+N0*Lx*10+9;
    cls = work+78+3*N0+300+2*(Lx+1)+N0*Lx*10+53;
    pYg = work+78+3*N0+300+2*(Lx+1)+N0*Lx*10+53+Lx;
    
    
    H2  = work+78+3*N0+300+2*(Lx+1)+N0*Lx*10+53+Lx+3*Lx;// length N0 * 10
    //if(Lx>1){cblas_dgemv(CblasColMajor, CblasNoTrans, N0*10, Lx, 1.0, H0, N0*10, ones, 1, 0.0, H2, 1);}else if(Lx==1){cblas_dcopy(N0*10, H0, 1, H2, 1);}else if(Lx==0){clear1(H2,N0*10);}
    
    
    lkhd_old=lkhd=-1e16;
    
    if(verbose3>1){fprintf(stderr, "beta=%lf pi=%lf delta=%lf phi=%lf theta=%lf\n", beta, pi[0], delta, phi, theta[0]);}
    
    
    double xi=.0;
    
    int stuck=0;
    
    
    for(itr=0; itr<MAXITR3+PRESTEPS+200; itr++){
        
        //fprintf(stderr, "E-step");
        // E-step
        lkhd = getLkhdAll(y, Y, h0, h1, H0, H1, ki, dki, ki2, K0, K2, K, km, Lx, N0, 3, 10, beta, theta, pi, delta, phi, asr, cls, pYg, a, A, H2, grad, hess);
        
        // M-step
        
        //fprintf(stderr, "M-step");
        
        
        int kk;
        // parameter fixation
        //printM(hess,10,6);
        for(j=0; j<5; j++){
            if(fixParam[j]>0){ // null
                grad[j] = 0.0; 
                //fprintf(stderr,"\n\n\nFix\n");
                for(kk=0; kk<5; kk++){
                    //fprintf(stderr, "%lf ", hess[kk*10+j]);
                    hess[kk*10+j] = 0.0;
                }
                for(kk=0; kk<10; kk++){
                    //fprintf(stderr, "%lf ", hess[j*10+kk]);
                    hess[j*10+kk] = 0.0;
                }
                //fprintf(stderr,"\nFix End\n\n");
            }
        }
        //bar
        //if(verbose3>1){fprintf(stderr, "\n");printM(hess,10,6);}
        //if(verbose3>1){fprintf(stderr, "grad> ");print(grad,5);}
        // prior
        
        grad[0] += (ab-1.0)*(1.0-2.0*pi[0]);
        hess[0] += -2.0*(ab-1.0)*pi[0]*(1.0-pi[0]);
        
        grad[1] += (ad-1.0)*(1.0-delta)-(bd-1.0)*delta;
        hess[11] += -((ad-1.0)+(bd-1.0))*delta*(1.0-delta);
        
        grad[2] += (phiab*1.-1.0)*(1.0-2.0*phi);        
        hess[22]+= -2.0*(phiab*1.-1.0)*phi*(1.0-phi);
        
        grad[3] += -theta[0]*(beta-beta0)/sigma2;
        grad[4] += 0.5-0.5*theta[0]*pow(beta-beta0,2.0)/sigma2     + (kappa/2.0-1.0) - omega/2.0*theta[0];
        
        hess[33]+= -theta[0]/sigma2;
        if(fixParam[4]==0){
            hess[34]+= -theta[0]*(beta-beta0)/sigma2;
        }
        hess[44]+= -theta[0]*0.5*pow(beta-beta0,2.0)/sigma2        - omega/2.0*theta[0];;
#ifdef DIFFTHETA
        grad[5] += (kappa/2.0-1.0) - omega/2.0*theta[1];
        hess[55] += - omega/2.0*theta[1];
#else
        grad[5]=hess[5]=hess[15]=hess[25]=hess[35]=hess[45]=hess[55]=0.0;
#endif
        
        //if(Lx>0){
        //    grad[5] += (ad-1.0)*(1.0-asr)-(bd-1.0)*asr;
        //    hess[55]+= -((ad-1.0)+(bd-1.0))*asr*(1.0-asr);
        //}
        //if(verbose3>1){fprintf(stderr, "hess>\n");printM(hess, 10,5);fprintf(stderr, "\n");}
        if(verbose3>1){fprintf(stderr, "grad2>");print(grad,5);}
        
        
        
        
        lkhdDiff=fabs(lkhd-lkhd_old);
        //if(fabs(lkhdDiff)<1e-4 && fabs(grad[0])<1e-5 && fabs(grad[1])<1e-5 && fabs(grad[2])<1e-5 && fabs(grad[3])<1e-5 && fabs(grad[4])<1e-5 && fabs(grad[5])<1e-5){break;}else{lkhd_old=lkhd;}
        
        
        
        if(fabs(lkhdDiff)<1e-4 && isConv(grad, nofp, 1.0e-5)==0){break;}else{lkhd_old=lkhd;}
        
        
        if(Lx==0){
            getStepSize5(hess, grad, step, ipiv);
        }else{
#ifdef DIFFTHETA
            getStepSize6(hess, grad, step, ipiv);
#else
            getStepSize5(hess, grad, step, ipiv);
#endif
        }
        //printM(hess, 5,5);fprintf(stderr, "\n");
        if(verbose3>1){fprintf(stderr, "step> ");print(step, 5);}
        
        stuck++;
        //if(stuck>1.5){step[1]=step[2]=step[3]=step[4]=0.0;}
        ssize=1.0;
        //qval = getLkhdAll(y, Y, h0, h1, H0, H1, ki, dki, ki2, K0, K2, K, km, Lx, N0, 3, 10, beta, theta, pi, delta, phi, cls, pYg, a, A);
        qval = getQval(y, Y, h1, H1, h0, H0, H2, ki, dki, ki2, K0, K2, K, km, Lx, N0, 3, 10, theta, beta, pi, delta, phi, asr);
        xi=0.1;
        for(itr_step=0; itr_step<50; itr_step++){
            pi1[0]    = fixParam[0]==0 ? expit(logit(pi[0])-ssize*step[0]) : 0.5;
            delta1    = fixParam[1]==0 ? expit(logit(delta)-ssize*step[1]) : (ad-1.0)/(ad+bd-2.0); 
            phi1      = fixParam[2]==0 ? expit(logit(phi)  -ssize*step[2]) : 0.5;
            beta1     = fixParam[3]==0 ? beta              -ssize*step[3]  : beta0;
            theta1[0] = fixParam[4]==0 ? exp(log(theta[0]) -ssize*step[4]) : (kappa-1.0)/(omega+0.0*pow(beta-beta0,2.0)/sigma2);
            //if(Lx>0){ asr1 = expit(logit(asr)  -ssize*step[5]); if(asr1   <= 1.0e-7 || asr1    >= 1.0-1.0e-7){asr1    = (asr1<0.5)    ? 1.0e-7 : 1.0-1.0e-7;}}else{asr1=0.0;};
            //fprintf(stderr,"asr=%lf\n", asr1);
            if(pi1[0] <= 1.0e-7 || pi1[0] >= 1.0-1.0e-7){pi1[0] = (pi1[0]<0.5) ? 1.0e-7 : 1.0-1.0e-7;}// ssize /= fabs(step[0]);}
            if(phi1   <= 1.0e-7 || phi1   >= 1.0-1.0e-7){phi1   = (phi1<0.5)   ? 1.0e-7 : 1.0-1.0e-7;}// ssize /= fabs(step[2]);}
            //if(phi1 <= 1.0e-7 || phi1   >= 1.0-1.0e-7){phi1   = 0.5;}
            if(delta1 <= 1.0e-7 || delta1 >= 1.0-1.0e-7){delta1 = (delta1<0.5) ? 1.0e-7 : 1.0-1.0e-7;}// ssize /= 100.0;}
            //if(delta1 <= 1.0e-7 || delta1 >= 1.0-1.0e-7){delta1 = 0.5;}
            if(beta1>100.0){   beta1 = 100.0;}  else if( beta1<(-100.)){beta1=-100.0;}// ssize /= fabs(step[3]);}
            if(theta1[0]>10000.0){theta1[0] = 10000.0;}else if( theta1[0]<1e-5 ){theta1[0]=1e-5;}// ssize /= fabs(step[4]);}
            
#ifdef DIFFTHETA
            theta1[1] = exp(log(theta[1])    -ssize*step[5]);
            if(theta1[1]>10000.0){theta1[1] = 10000.0;}else if( theta1[1]<1e-5 ){theta1[1]=1e-5;}// ssize /= fabs(step[4]);}
#else
            theta1[1] = theta1[0];
#endif
            if(transQTL==0){pi1[1] = pi1[0];}
            //qvalNext = getLkhdAll(y, Y, h0, h1, H0, H1, ki, dki, K0, K2, K, km, Lx, N0, 3, 10, beta1, theta1, pi1, delta1, phi1, cls, pYg, a, A);
            qvalNext = getQval(y, Y, h1, H1, h0, H0, H2, ki, dki, ki2, K0, K2, K, km, Lx, N0, 3, 10, theta1, beta1, pi1, delta1, phi1, asr1);
            if(verbose3>1)fprintf(stderr, "%ld: %lf %lf ssize=%lf %lf %lf %lf %lf %lf %lf\n", itr_step, qval, qvalNext, ssize, pi1[0], delta1, phi1, beta1, theta1[0],theta1[1]); 
            
            if(qvalNext + xi > qval){
                //if(verbose3>1)fprintf(stderr, "log qval diff = %lf\n", log(qvalNext-qval));
                //if(log(qvalNext + xi - qval)<(-20.0)){stuck2++;}else{stuck2=0;}
                stuck=0;
                qval  = qvalNext;
                pi[0] = pi1[0];
                pi[1] = pi1[1];
                phi   = phi1;
                delta = delta1;
                beta  = beta1;
                theta[0] = theta1[0];
                theta[1] = theta1[1];
                asr =asr1;
                break;
            }else if(ssize>0.0){ssize *= -1.0;}else{ssize/=-2.0;}//if(ssize>0.0){ssize*=(-1.0);}else{ssize/=-2.0;}}
        }
        
        
        //if(stuck2>3){xi=0.0;}else{xi=0.0;}
        
        if(stuck>0){
            if(verbose3>1)fprintf(stderr, "stuck=%d\n", stuck);
            //itr--;
            xi+=.1;
            //if(stuck>5){xi2=10.0;}else if(stuck>10){xi2=1000.0;}
            //xi*=2.0;
        }else{
            //fprintf(stderr, "Q=%lf %lf %lf\n", qvalNext, grad[0], pi);
            xi = 0.0;
            //xi2= 0.0;
        }
        
        
        
        if(verbose3>1)fprintf(stderr, "ASEQTL [%ld] pi=%lf delta=%lf phi=%lf beta=%lf theta=%lf theta2=%lf asr=%lf lkhd=%lf ldiff=%lf\n", itr, pi[0], delta, phi, beta, theta[0], theta[1], asr, lkhd,lkhdDiff);
    }
    
    //getInformation(y,Y,h0,H0,h1,H1,H2,ki,dki,ki2,K0,K2,K,km,Lx,N0,3,10,beta,theta,pi,delta,phi,K2pp,hess,ipiv,a,A);
    //inverseLapack(hess, 5, ipiv, hess+25);
    
    //updateHess2( y,  Y,  h0,  H0,  h1,  H1,  H2,  ki,  dki,  ki2, K0,  K2,  K,  km, Lx, N0, 3, 10, beta, theta, pi, delta, phi,  work,  hess, ipiv,  a,  A);
    
    //updateHess(y, Y, h1, H1, ki, dki, ki2, K0, dpdp, K2, K2p, K2d, K2h, K, km, N0, N0*Lx, 3, 10, 2, beta[0], theta[0], pi, delta, phi, d2pdp2, hess, ipiv);
    
    
    double kld  = getCov(h0, h1, N0);
    double kld2 = getCov2(H0, H1, N0*Lx);
    
    // K-L divergence
    //for(j=0; j<N0*3;     j++){if(h1[j]>0.0&&h0[j]>0.0)kld += h1[j]*log(h1[j]/h0[j]);}
    //for(j=0; j<N0*10*Lx; j++){if(H1[j]>0.0&&H0[j]>0.0)kld += H1[j]*log(H1[j]/H0[j]);}
    
    // !!!!!!!!!!!!!
    (*ppi)   = pi[0];
    (*pphi)  = phi;
    (*pdelta)= delta;
    (*ptheta)=theta[0];
    (*pbeta) =beta;
    (*pitr)  = itr;
    (*pasr) = asr;
    (*pbound)= (lkhdDiff<1e-4?0:1000000) + isConv(grad, nofp, 1.0e-5);
    
    
    (*tval) = -pow(logit(pi[0]),2.0)/hess[0];
    pkld[0] = kld2;//theta[1];
    //pkld[L+2] = getIA1(h1,N0);//theta[1];
    //pkld[L+2] = getIA2(H1,N0*Lx);//theta[1];
    pkld[1] = kld;
    theta2=theta[1];
    
    
    free(ipiv);
    
    
    if(verbose3>0)fprintf(stderr, "[%ld] pi=%lf delta=%lf phi=%lf beta=%lf theta=%lf theta2=%lf asr=%lf lkhd=%lf ldiff=%lf gpi=%lf gd=%lf gh=%lf gb=%lf gt=%lf gR=%lf tval=%lf kld=%lf\n", itr,pi[0], delta, phi, beta, theta[0], theta[1], asr, lkhd, lkhdDiff, grad[0], grad[1], grad[2], grad[3], grad[4], grad[5], *tval, kld);
    return lkhd;
    
}











double getLkhdAll(double* y, double* Y, double* h0, double* h, double* H0, double* H, double* ki, double* dki, double* ki2, double* K0, double* K2, double* ktmp, double* km, long Lx, long N0, long J0, long J, double beta, double* th, double* pi, double delta, double phi, double asr, double* cls, double* pYg, double* a, double* A, double* H2, double* grad, double* hess){
    double lkhd = 0.0;
    long i0,j0,m,i,j,l,i1;
    double thij, lmij, lcij, lci, cl0, Lci, kij;
    double mu = exp(beta);
    double denom;
    
    double dLx=(double)Lx;
    if(Lx==0){dLx=1.0;}
    
    double* K0p = K0+3;
    double* K0pp = K0+6;
    
    clear1(h, N0*3);
    clear1(H, N0*Lx*10);
    
    // as counts    
    
    getK0s(pi[0], K0);
    
    
    
    /*fprintf(stderr, "J0=%ld\n", J0);
     fprintf(stderr, "L =%ld\n", L);
     fprintf(stderr, "M =%ld\n", M);
     fprintf(stderr, "J =%ld\n", J);*/
    
    for(i0=0; i0<N0; i0++){
        getK2s(pi[1], delta, phi, ktmp, ki2+i0*2, K2);
        lci=0.0;
        for(j0=0; j0<J0; j0++){
            if(h0[i0*J0+j0]>0.0){
                kij = K0[j0];
#ifdef NOKJ
                thij = th[0]*ki[i0];
#else
                thij = th[0]*ki[i0]*kij;
#endif
                lmij = mu*ki[i0]*dki[i0]*kij;
                if(ASE==2){
                    h[i0*J0+j0] = log(h0[i0*J0+j0]);
                }else{
                    h[i0*J0+j0] = log(h0[i0*J0+j0]) + lgamma(thij+y[i0]) - lgamma(y[i0]+1.0) - lgamma(thij) + y[i0]*log(lmij) + (thij)*log(thij) - (y[i0]+thij)*log(lmij+thij);
                }
                lci += h[i0*J0+j0];
            }
        }
        if(lci<-10.0){lci=-10.0;}
        cl0 = lci/(double)J0;
        for(j0=0; j0<J0; j0++){
            if(h0[i0*J0+j0]>0.0){
                h[i0*J0+j0] -= cl0;
            }
        }
        // fuga
        // for each locus
        clear1(pYg, Lx*3);
        for(l=0;l<Lx;l++){
            i=i0+l*N0;
            i1=i0+(l+1)*N0*randomize;
            lci=0;
            for(j=0; j<J; j++){
                if(H0[i*J+j]>0.0 && h0[i0*J0+ig[j]]>0.0){
                    
                    lcij = log(H0[i*J+j]) - log(h0[i0*J0+ig[j]])
                                - log(Y[i*2]+Y[i*2+1]+1.0) - lbeta(Y[i*2]+1.0, Y[i*2+1]+1.0)
                                + lbeta(Y[i*2]+th[1]*ki[i1]*K2[j*2]*km[l], Y[i*2+1]+th[1]*ki[i1]*K2[j*2+1]*km[l]) 
                                - lbeta(       th[1]*ki[i1]*K2[j*2]*km[l],          th[1]*ki[i1]*K2[j*2+1]*km[l]);
                    H[i*J+j] = lcij;
                    lci += lcij;
                }
            }
            if(lci<-10.0){lci=-10.0;}
            cls[l]=lci/(double)J;
            
            for(j=0; j<J; j++){
                if(H0[i*J+j]>0.0){
                    H[i*J+j] -= cls[l];
                    pYg[l*3+ig[j]] += exp(H[i*J+j]);
                    
                }
            }
        }
        lkhd += cl0;
        for(l=0;l<Lx;l++){lkhd += cls[l];}
        for(l=0;l<Lx;l++){
            for(j0=0; j0<J0; j0++){ if(pYg[l*J0+j0]>0.0){
                h[i0*J0+j0] += log(pYg[l*J0+j0]);
            }}
        }
        Lci=0.0;
        for(j0=0; j0<J0; j0++){if(h0[i0*J0+j0]>0.0){
            Lci += exp(h[i0*J0+j0]);
        }}
        if(Lci>0){ lkhd += log(Lci); }//else{for(j0=0; j0<J0; j0++){fprintf(stderr, "%ld h0i=%lf hi=%lf\n", i0, h0[i0*J0+j0], h[i0*J0+j0]);}}
        
        // posterior calc
        // fSNPs
        for(l=0;l<Lx;l++){
            i=i0+l*N0;
            for(j=0; j<J; j++){if(H0[i*J+j]>0.0 && pYg[l*3+ig[j]]>0.0){
                H[i*J+j] += -log(pYg[l*3+ig[j]]) + h[i0*3+ig[j]];
            }}
            double Lcli=0.0;
            for(j=0; j<J; j++){if(H0[i*J+j]>0.0){
                Lcli += exp(H[i*J+j]);
            }}
            for(j=0; j<J; j++){if(H0[i*J+j]>0.0 && Lcli>0.0){
                if(noPosteriorUpdate>0){
                	H[i*J+j] = H0[i*J+j]; // no posterior update
                }else{
                    H[i*J+j] = ( exp(H[i*J+j]) )/(Lcli);
                }
            }else{H[i*J+j]=0.0;}}
        }
        // rSNP
        for(j0=0; j0<J0; j0++){
            if(h0[i0*J0+j0]>0.0 && Lci >0.0){
                if(noPosteriorUpdate>0){
                	h[i0*J0+j0] = h0[i0*J0+j0]; // no posteror update
                }else{
                    h[i0*J0+j0] = exp(h[i0*J0+j0])/Lci;
            	}
            }else{
                h[i0*J0+j0]=0.0;
            }
        }
    }
    // prior
    lkhd += kappa/2.0*log(omega/2.0) + (kappa/2.0-1.0)*log(th[0]) - omega*th[0]/2.0 - lgamma(kappa/2.0);
    lkhd += 0.5*log(th[0]) - 0.5*log(2.0*M_PI) - 0.5*log(sigma2) - th[0]*pow(beta-beta0,2.0)/sigma2/2.0;
    
#ifdef DIFFTHETA
    lkhd += kappa/2.0*log(omega/2.0) + (kappa/2.0-1.0)*log(th[1]) - omega*th[1]/2.0 - lgamma(kappa/2.0);
#endif
    
    lkhd += (ad-1.0)*log(delta) + (bd-1.0)*log(1.0-delta);
    lkhd += (ab-1.0)*(log(pi[0])+log(1.0-pi[0]));
    lkhd += (phiab-1.0)*(log(phi)+log(1.0-phi));
    
    //if(Lx>0)lkhd += (ad-1.0)*log(asr) + (bd-1.0)*log(1.0-asr);
    
    if(isnan(lkhd)>0){fprintf(stderr, "Lkhd nan! %lf %lf \n", pi[0], pi[1]);}
    
    // calculation for M-step
    clear1(grad, 20);
    clear1(hess, 100);
    //if(Lx>1){cblas_dgemv(CblasColMajor, CblasNoTrans, N0*10, Lx, 1.0, H, N0*10, ones, 1, 0.0, H2, 1);}else if(Lx==1){cblas_dcopy(N0*10, H, 1, H2, 1);}
    clear1(a, J0*5);
    mu=exp(beta);
    
    //double kmL = 0.0; if(L>0){kmL = sum(km, L);}
    double k_p, k_pp, k_d, k_h, k_pd, k_ph, k_dd, k_dh, k_hh;//, k_R, k_RR, k_pR, k_dR, k_hR;
    double abij, c2de, acd, bde;
    for(i=0; i<N0; i++){
        for(j=0; j<J0; j++){
            if(h0[i*J0+j]>0.0){
                
                kij = K0[j]; 
                k_p = K0p[j]/kij;
                k_pp= K0pp[j]/kij;
                k_d = k_dd = k_h = k_hh = k_pd= k_ph = k_dh = 0.0;
                
                lmij = mu*ki[i]*dki[i]*kij;
#ifdef NOKJ
                thij = th[0]*ki[i];
#else
                thij = th[0]*ki[i]*kij;
#endif
                //thij = th[0]*ki[i]*kij;
                denom = 1.0 + lmij/thij;
                a[j]      = h[i*J0+j] * (y[i]-lmij)/denom;
                a[j+J0]   = h[i*J0+j] * thij * (digamma( thij + y[i]) - digamma( thij) + log(thij) - log( lmij + thij) + (lmij - y[i])/(lmij + thij) );
                a[j+J0*2] = h[i*J0+j] * ( -lmij/denom - (y[i]-lmij)*(1.0+2.0*lmij/thij)/denom/denom );
                a[j+J0*3] = h[i*J0+j] * ( (y[i]-lmij)*lmij*thij/pow(thij+lmij,2.0) );
                a[j+J0*4] = h[i*J0+j] * pow(thij,2.0) * ( trigamma(thij + y[i]) - trigamma(thij) + 1.0/thij  - 1.0/(lmij + thij) - (lmij - y[i])/pow(lmij + thij, 2.0));
                
#ifdef NOKJ
                abij = (a[j]);
                c2de = (a[j+J0*2]);
                acd =  (a[j]      + a[j+J0*2]);
                bde =  (a[j+J0*3]);
#else
                abij = (a[j]      + a[j+J0]);
                c2de = (a[j+J0*2] +2.0*a[j+J0*3]+a[j+J0*4]);
                acd =  (a[j]      + a[j+J0*2] + a[j+J0*3]);
                bde =  (a[j+J0]   + a[j+J0*3] + a[j+J0*4]);
#endif
                if(ASE<2){
                    grad[0] += abij*k_p;
                    grad[1] += abij*k_d;
                    grad[2] += abij*k_h;
                    grad[3] += a[j];
                    grad[4] += a[j+J0];
                    //grad[5] += abij*k_R;
                    
                    //hoge fuga foo
                    
                    hess[0] += c2de * k_p*k_p + abij * k_pp;
                    hess[1] += c2de * k_p*k_d + abij * k_pd;
                    hess[2] += c2de * k_p*k_h + abij * k_ph;
                    hess[3] += acd  * k_p;
                    hess[4] += bde  * k_p;
                    //hess[5] += c2de * k_p*k_R + abij * k_pR;
                    
                    hess[11] += c2de * k_d*k_d + abij * k_dd;
                    hess[12] += c2de * k_d*k_h + abij * k_dh;
                    hess[13] += acd  * k_d;
                    hess[14] += bde  * k_d;
                    //hess[15] += c2de * k_d*k_R + abij * k_dR;
                    
                    hess[22] += c2de * k_h*k_h + abij * k_hh;
                    hess[23] += acd  * k_h;
                    hess[24] += bde  * k_h;
                    //hess[25] += c2de * k_h*k_R + abij * k_hR;
                    
                    hess[33] += a[j]+a[j+J0*2];
                    hess[34] += a[j+J0*3];
                    //hess[35]+= acd  * k_R;
                    
                    hess[44] += a[j+J0]+a[j+J0*4];
                    //hess[45]+= bde  * k_R;
                    
                    //hess[55]+= c2de * k_R*k_R + abij * k_RR;
                }
            }
        }
    }
    
    if(Lx>0){
        double dlAB, d2lAB, thijA, thijB;
        double* K2p   = K2+20;
        double* K2d   = K2+40;
        double* K2h   = K2+60;
        double* K2pp  = K2+80;
        double* K2pd  = K2+100;
        double* K2ph  = K2+120;
        double* K2dd  = K2+140;
        double* K2dh  = K2+160;
        double* K2hh  = K2+180;
        
        for(i0=0; i0<N0; i0++){
            clear1(A, 5*J);
            getK2s(pi[1],delta,phi,ktmp,ki2+i0*2,K2);
            
            for(l=0;l<Lx;l++){// for each xSNP l
                i=i0+l*N0;
                i1 = i%N0 + (i/N0+1)*N0*randomize;       //if(randomize>0){i1=i+N0;}else{i1=i%N0;}
                double AB=th[1]*ki[i1]*km[i/N0];
                double Yitot=Y[i*2]+Y[i*2+1];
                for(j=0; j<J; j++){// all possible 
                    thijA = AB*K2[j*2];
                    thijB = AB*K2[j*2+1];
                    
                    dlAB  = digamma(thijA+thijB)  - digamma(thijA+thijB + Yitot) ;
                    d2lAB = trigamma(thijA+thijB) - trigamma(thijA+thijB + Yitot) ;
                    //dlAB=d2lAB=0.0;
                    A[j]     += H[i*J+j] * thijA         * (digamma(thijA + Y[i*2])    - digamma(thijA)  + dlAB );
                    A[j+J]   += H[i*J+j] * thijB         * (digamma(thijB + Y[i*2+1])  - digamma(thijB)  + dlAB );
                    A[j+J*2] += H[i*J+j] * thijA * thijA * (trigamma(thijA + Y[i*2])   - trigamma(thijA) + d2lAB );
                    A[j+J*3] += H[i*J+j] * thijA * thijB * d2lAB ;
                    A[j+J*4] += H[i*J+j] * thijB * thijB * (trigamma(thijB + Y[i*2+1]) - trigamma(thijB) + d2lAB );
                    
                    //if(i1==0 && H[i*J+j]>0.5){ fprintf(stderr, "%d %d %d %lf %lf %lf %lf   %lf %lf\n", i%N0, i/N0, j, K2[j*2], K2[j*2+1], ki2[j*2], ki2[j*2+1], pi[0], pi[1]); }
                }
            }
            
            // theta

            grad[4] += sum(A, 10) + sum(A+10, 10);
            hess[44]+= sum(A+20, 10) + 2.0*sum(A+30, 10) + sum(A+40, 10) + sum(A, 10) + sum(A+10, 10);
            
            hess[4] += getHess3(A, K2, K2p);
            hess[14] += getHess3(A, K2, K2d);
            hess[24]+= getHess3(A, K2, K2h);
            
            // others
            if(transQTL==0)grad[0] += getGrad(A, K2, K2p);
            grad[1] += getGrad(A, K2, K2d);
            grad[2] += getGrad(A, K2, K2h);
            
            if(transQTL==0){
                hess[0] += getHess1(A, K2, K2p, K2pp);
                hess[1] += getHess2(A, K2, K2p, K2d, K2pd);
                hess[2] += getHess2(A, K2, K2p, K2h, K2ph);
            }
            
            hess[11] += getHess1(A, K2, K2d, K2dd);
            hess[12] += getHess2(A, K2, K2d, K2h, K2dh);
            
            hess[22]+= getHess1(A, K2, K2h, K2hh);
        }
    }
    
    return lkhd;
}






double getQval(double* y, double* Y, double* h1, double* H1, double* h0, double* H0, double* H2, double* ki, double* dki, double* ki2, double* K0, double* K2, double* ktmp, double* km, long Lx, long N0, long J0, long J, double* th, double beta, double* pi, double delta, double phi, double asr){
    double lkhd = 0.0;
    long i0,j0,m,i,j,l,i1,i2;
    double thij, lmij, kij;
    double mu = exp(beta);
    
    double dLx = (double)Lx;
    
    double AB, Aj, Bj;
    if(Lx==0){dLx=1.0;}
    
    K0[0] = 2.0*(1-pi[0]); K0[1] = 1.0; K0[2] = 2.0*pi[0];
    
    for(i0=0; i0<N0; i0++){
        for(j0=0; j0<J0; j0++){
            if(h0[i0*J0+j0]>0.0){
                kij = K0[j0];
                
#ifdef NOKJ
                thij = th[0]*ki[i0];
#else
                thij = th[0]*ki[i0]*kij;
#endif
                //thij = th[0]*ki[i0]*kij;
                lmij = mu*ki[i0]*dki[i0]*kij;
                if(ASE==2){
                    lkhd += h1[i0*J0+j0];
                }else{
                    lkhd += h1[i0*J0+j0] * ( lgamma(thij+y[i0]) - lgamma(y[i0]+1.0) - lgamma(thij) + y[i0]*log(lmij) + (thij)*log(thij) - (y[i0]+thij)*log(lmij+thij) );
                }
            }
        }
    }

    
    for(i0=0; i0<N0; i0++){
        getK(pi[1],delta,phi,ktmp,ki2+i0*2,K2);
        for(l=0;l<Lx;l++){// for each xSNP l
            i=i0+l*N0;
            i2=i*2;
            i1=i0 + (l+1)*N0*randomize;
            AB = th[1]*ki[i1]*km[l];
            for(j=0; j<J; j++){
                //if(H0[i*J+j]>0.0){
                //lkhd += H1[i*J+j] * ( - lgamma(th[1]*ki[i1]*km[l]*(K2[j*2]+K2[j*2+1])+Y[i2]+Y[i2+1]) + lgamma(Y[i2]+Y[i2+1]+1.0) + lgamma(th[1]*ki[i1]*km[l]*(K2[j*2]+K2[j*2+1])) );
                //for(m=0; m<2; m++){
                //    thij = th[1]*ki[i1]*K2[j*2+m]*km[l];
                //    lkhd += H1[i*J+j] * ( lgamma(thij+Y[i2+m]) - lgamma(Y[i2+m]+1.0) - lgamma(thij) );
                //}
                Aj = AB*K2[j*2];
                Bj = AB*K2[j*2+1];
                lkhd += H1[i*J+j] * ( -log(Y[i2]+Y[i2+1]+1.0) - lbeta(Y[i2]+1.0, Y[i2+1]+1.0) 
                                     + lbeta(Y[i2]+Aj, Y[i2+1]+Bj) 
                                     - lbeta(      Aj,         Bj) );
                //}
            }
        }
    }
    
    
    // prior
    lkhd += kappa/2.0*log(omega/2.0) + (kappa/2.0-1.0)*log(th[0]) - omega*th[0]/2.0 - lgamma(kappa/2.0);
    lkhd += 0.5*log(th[0]) - 0.5*log(2.0*M_PI) - 0.5*log(sigma2) - th[0]*pow(beta-beta0,2.0)/sigma2/2.0;
    
#ifdef DIFFTHETA
    lkhd += kappa/2.0*log(omega/2.0) + (kappa/2.0-1.0)*log(th[1]) - omega*th[1]/2.0 - lgamma(kappa/2.0);
#endif
    
    //lkhd += (ad-1.0)*log(delta-lower) + (bd-1.0)*log(upper-delta);
    lkhd += (ad-1.0)*log(delta)  + (bd-1.0)*log(1.0-delta);
    lkhd += (ab-1.0)*(log(pi[0]) + log(1.0-pi[0]));
    lkhd += (phiab-1.0)*(log(phi)+ log(1.0-phi));
    
    //if(Lx>0)lkhd += (ad-1.0)*log(asr) + (bd-1.0)*log(1.0-asr);
    
    if(isnan(lkhd)>0){fprintf(stderr, "Qval nan! %lf %lf\n", pi[0], pi[1]);}
    return lkhd;
}


















long nbbem_init0(double* y, double* ki, long N, double* pbeta, double* ptheta){
    long i;
    // init mean
    double mu = 0.0;
    double va = 0.0;
    double A=(double)(N);
    
    for(i=0; i<N; i++){// sample
        mu += y[i]/ki[i]/A;
    }
    pbeta[0] = log(mu);
    if(beta0>9999.0){beta0=pbeta[0];}
    for(i=0; i<N; i++){
        va += pow(y[i]/ki[i]-mu, 2.0)/A;
    }
    // init theta
    theta0 = mu*mu/(va-mu);
    if(theta0<0.0){theta0=34;}else if(theta0<0.001){theta0=0.001;}else if(theta0>1000.){theta0=1000.;}
    ptheta[0]=theta0;
    if(verbose3>0){
        fprintf(stderr, "\n\nva0=%lf mu0=%lf th0=%lf\n\n", va, mu, theta0);
    }
    
    return 0;
}





// z1 cSNP
// z2 xSNP

long getPrior(double* h0, double* z1, double* z2, double wi){
	double D;
	double u00, u01, u10, u11;
	double v00, v01, v10, v11;
	if(z1[0]<0.5){
		if(z2[0]<0.5){// 00
			D = (1.0-z1[0]) < (1.0-z2[0]) ? (1.0-z1[0])*z2[0] : (1.0-z2[0])*z1[0];
		}else{// 01
			D = (1.0-z1[0]) < z2[0] ? -(1.0-z1[0])*(1.0-z2[0]) : -z2[0]*z1[0];
		}
	}else{
		if(z2[0]<0.5){// 10
			D = z1[0] < (1.0-z2[0]) ? -z1[0]*z2[0] : -(1.0-z2[0])*(1.0-z1[0]);
		}else{// 11
			D = z1[0] < z2[0] ? z1[0]*(1.0-z2[0]) : z2[0]*(1.0-z1[0]);
		}
	}
    D=0.0;
	u00 = (1.0-z1[0])*(1.0-z2[0]) + D;
	u01 = (1.0-z1[0])*z2[0]       - D;
	u10 = z1[0]*(1.0-z2[0])       - D;
	u11 = z1[0]*z2[0]             + D;
	
	if(z1[1]<0.5){
		if(z2[1]<0.5){// 00
			D = (1.0-z1[1]) < (1.0-z2[1]) ? (1.0-z1[1])*z2[1] : (1.0-z2[1])*z1[1];
		}else{// 01
			D = (1.0-z1[1]) < z2[1] ? -(1.0-z1[1])*(1.0-z2[1]) : -z2[1]*z1[1];
		}
	}else{
		if(z2[1]<0.5){// 10
			D = z1[1] < (1.0-z2[1]) ? -z1[1]*z2[1] : -(1.0-z2[1])*(1.0-z1[1]);
		}else{// 11
			D = z1[1] < z2[1] ? z1[1]*(1.0-z2[1]) : z2[1]*(1.0-z1[1]);
		}
	}
    D=0.0;
	v00 = (1.0-z1[1])*(1.0-z2[1]) + D;
	v01 = (1.0-z1[1])*z2[1]       - D;
	v10 = z1[1]*(1.0-z2[1])       - D;
	v11 = z1[1]*z2[1]             + D;
	
	h0[0] = u00*v00;
	h0[1] = u00*v01 + v00*u01;
	h0[2] = u01*v01;
	h0[3] = u10*v00 + v10*u00;
	h0[4] = u10*v01 + v10*u01;
	h0[5] = u00*v11 + v00*u11;
	h0[6] = u11*v01 + v11*u01;
	h0[7] = u10*v10;
	h0[8] = u10*v11 + v10*u11;
	h0[9] = u11*v11;
#ifdef RECOMB
	h0[4]=h0[5]=(h0[4]+h0[5])/2.0;
#endif
    return 0;
}

long getPriorIndep(double* h0, double* z1, double* z2, double wi){
    h0[0] = (1.0-z1[0])*(1.0-z2[0])*(1.0-z1[1])*(1.0-z2[1])*wi;
    h0[1] = ((1.0-z1[0])*(1.0-z2[0])*(1.0-z1[1])*(z2[1])+(1.0-z1[0])*(z2[0])*(1.0-z1[1])*(1.0-z2[1]))*wi;
    h0[2] = (1.0-z1[0])*(z2[0])*(1.0-z1[1])*(z2[1])*wi;
    
    h0[3] = ((z1[0])*(1.0-z2[0])*(1.0-z1[1])*(1.0-z2[1])+(1.0-z1[0])*(1.0-z2[0])*(z1[1])*(1.0-z2[1]))*wi;
    h0[4] = ((z1[0])*(1.0-z2[0])*(1.0-z1[1])*(z2[1])+(1.0-z1[0])*(z2[0])*(z1[1])*(1.0-z2[1]))*wi;            // 10/01 + 01/10
    h0[5] = ((z1[0])*(z2[0])*(1.0-z1[1])*(1.0-z2[1])+(1.0-z1[0])*(1.0-z2[0])*(z1[1])*(z2[1]))*wi;            // 11/00 + 00/11
    h0[6] = ((z1[0])*(z2[0])*(1.0-z1[1])*(z2[1])+(1.0-z1[0])*(z2[0])*(z1[1])*(z2[1]))*wi;
    
    h0[7] = (z1[0])*(1.0-z2[0])*(z1[1])*(1.0-z2[1])*wi;
    h0[8] = ((z1[0])*(1.0-z2[0])*(z1[1])*(z2[1])+(z1[0])*(z2[0])*(z1[1])*(1.0-z2[1]))*wi;
    h0[9] = z1[0]*z2[0]*z1[1]*z2[1]*wi;
    return 1;
}
long getPrior1(double* h0, double* z, double wi){
    h0[0] = (1.0-z[0])*(1.0-z[1])*wi;
    h0[1] = 0.0;
    h0[2] = 0.0;
    
    h0[3] = 0.0;
    h0[4] = 0.0;
    h0[5] = ((1.0-z[0])*(z[1]) + (z[0])*(1.0-z[1]))*wi;
    h0[6] = 0.0;
    
    h0[7] = 0.0;
    h0[8] = 0.0;
    h0[9] = z[0]*z[1]*wi;
    return 1;
}
long getPriorNull(double* h0, double* z2, double wi){
	double z10, z11;
	z10 = z11 = 0.5;
    h0[0] = (1.0-z10)*(1.0-z2[0])*(1.0-z11)*(1.0-z2[1])*wi;
    h0[1] = ((1.0-z10)*(1.0-z2[0])*(1.0-z11)*(z2[1])+(1.0-z10)*(z2[0])*(1.0-z11)*(1.0-z2[1]))*wi;
    h0[2] = (1.0-z10)*(z2[0])*(1.0-z11)*(z2[1])*wi;
    
    h0[3] = ((z10)*(1.0-z2[0])*(1.0-z11)*(1.0-z2[1])+(1.0-z10)*(1.0-z2[0])*(z11)*(1.0-z2[1]))*wi;
    h0[4] = ((z10)*(1.0-z2[0])*(1.0-z11)*(z2[1])+(1.0-z10)*(z2[0])*(z11)*(1.0-z2[1]))*wi;
    h0[5] = ((z10)*(z2[0])*(1.0-z11)*(1.0-z2[1])+(1.0-z10)*(1.0-z2[0])*(z11)*(z2[1]))*wi;
    h0[6] = ((z10)*(z2[0])*(1.0-z11)*(z2[1])+(1.0-z10)*(z2[0])*(z11)*(z2[1]))*wi;
    
    h0[7] = (z10)*(1.0-z2[0])*(z11)*(1.0-z2[1])*wi;
    h0[8] = ((z10)*(1.0-z2[0])*(z11)*(z2[1])+(z10)*(z2[0])*(z11)*(1.0-z2[1]))*wi;
    h0[9] = z10*z2[0]*z11*z2[1]*wi;
    return 1;
}
long getPriorIP(double* h0, double* z2, double wi){
    h0[0] = 0.0;
    h0[1] = 0.0;
    h0[2] = 0.0;
    
    h0[3] = (1.0-z2[0])*(1.0-z2[1])*wi;
    h0[4] = 0.5*((z2[0])*(1.0-z2[1])+(1.0-z2[0])*(z2[1]))*wi;
    h0[5] = 0.5*((z2[0])*(1.0-z2[1])+(1.0-z2[0])*(z2[1]))*wi;
    h0[6] = z2[0]*z2[1]*wi;
    
    h0[7] = 0.0;
    h0[8] = 0.0;
    h0[9] = 0.0;
    return 1;
}








/*
 long getCondPrior(double* h0, double* z1, double* z2){
 h0[0] = h0[3] = h0[7] = (1.0-z2[0])*(1.0-z2[1]);
 h0[1] = h0[8] = (1.0-z2[0])*(z2[1])+(z2[0])*(1.0-z2[1]);
 h0[2] = h0[6] = h0[9] = (z2[0])*(z2[1]);
 
 if(((1.0-z1[0])*(z1[1]) + (z1[0])*(1.0-z1[1]))>0.0){
 h0[4] = ((z1[0])*(1.0-z2[0])*(1.0-z1[1])*(z2[1])+(1.0-z1[0])*(z2[0])*(z1[1])*(1.0-z2[1]))/((1.0-z1[0])*(z1[1]) + (z1[0])*(1.0-z1[1]));
 h0[5] = ((z1[0])*(z2[0])*(1.0-z1[1])*(1.0-z2[1])+(1.0-z1[0])*(1.0-z2[0])*(z1[1])*(z2[1]))/((1.0-z1[0])*(z1[1]) + (z1[0])*(1.0-z1[1]));
 }else{
 h0[4] = 0.0;
 h0[5] = 0.0;
 }
 
 return 1;
 }
 long getPriorHet(double* h0, double* z1, double* z2, double wi){
 h0[0] = 0.0;
 h0[1] = 0.0;
 h0[2] = 0.0;
 
 h0[3] = 0.0;
 h0[4] = ((z1[0])*(1.0-z2[0])*(1.0-z1[1])*(z2[1])+(1.0-z1[0])*(z2[0])*(z1[1])*(1.0-z2[1]))*wi;
 h0[5] = ((z1[0])*(z2[0])*(1.0-z1[1])*(1.0-z2[1])+(1.0-z1[0])*(1.0-z2[0])*(z1[1])*(z2[1]))*wi;
 h0[6] = 0.0;
 
 h0[7] = 0.0;
 h0[8] = 0.0;
 h0[9] = 0.0;
 return 1;
 }
 long getPrior1Het(double* h0, double* z, double wi){
 h0[0] = 0.0;
 h0[1] = 0.0;
 h0[2] = 0.0;
 
 h0[3] = 0.0;
 h0[4] = 0.0;
 h0[5] = ((1.0-z[0])*(z[1]) + (z[0])*(1.0-z[1]))*wi;;
 h0[6] = 0.0;
 
 h0[7] = 0.0;
 h0[8] = 0.0;
 h0[9] = 0.0;
 return 1;
 }
 long getCondPrior1(double* h0, double* z){
 h0[0] = h0[5] = h0[9] = 1.0;
 h0[1] = h0[2] = h0[3] = h0[4] = h0[6] = h0[7] = h0[8] = 0.0;
 return 1;
 }
 */






double getGrad(double* A, double* K2, double* K2p){
    long j;
    double res=0.0;
    for(j=0; j<10; j++){
        res += A[j]*K2p[j*2]/K2[j*2] + A[j+10]*K2p[j*2+1]/K2[j*2+1];
    }
    return res;
}
double getHess1(double* A, double* K2, double* K2p, double* K2pp){
    long j;
    double res=0.0;
    for(j=0; j<10; j++){
        res += A[j+20]*pow(K2p[j*2]   /K2[j*2],2.0) 
        + A[j+30]*   (K2p[j*2]   /K2[j*2])
        *(K2p[j*2+1] /K2[j*2+1])*2.0
        + A[j+40]*pow(K2p[j*2+1] /K2[j*2+1],2.0)
        + A[j]       *K2pp[j*2]  /K2[j*2] 
        + A[j+10]    *K2pp[j*2+1]/K2[j*2+1];
    }
    return res;
}
double getHess2(double* A, double* K2, double* K2p, double* K2d, double* K2pd){
    long j;
    double res=0.0;
    for(j=0; j<10; j++){
        res += A[j+20]*(K2p[j*2]   /K2[j*2])  *(K2d[j*2]   /K2[j*2])
        + A[j+30]*(K2p[j*2]   /K2[j*2])  *(K2d[j*2+1] /K2[j*2+1])
        + A[j+30]*(K2p[j*2+1] /K2[j*2+1])*(K2d[j*2]   /K2[j*2])
        + A[j+40]*(K2p[j*2+1] /K2[j*2+1])*(K2d[j*2+1] /K2[j*2+1])
        
        + A[j]       *K2pd[j*2]  /K2[j*2] 
        + A[j+10]    *K2pd[j*2+1]/K2[j*2+1];
    }
    return res;
}

double getHess3(double* A, double* K2, double* K2p){
    long j;
    double res=0.0;
    for(j=0; j<10; j++){
        res += (A[j]+A[j+20]+A[j+30])*K2p[j*2]/K2[j*2] + (A[j+10]+A[j+30]+A[j+40])*K2p[j*2+1]/K2[j*2+1];
    }
    return res;
}






void getK(double pi, double delta, double phi, double* K, double* w, double* K2){
    // A allele
    K[0] = (1.0-pi)*(1.0-delta)*(1.0-phi);
    K[1] = (1.0-pi)*(delta)*(1.0-phi);
    K[2] = (pi)*(1.0-delta)*(1.0-phi);
    K[3] = (pi)*(delta)*(1.0-phi);
    // B allele
    K[4] = (1.0-pi)*(delta)*(phi);
    K[5] = (1.0-pi)*(1.0-delta)*(phi);
    K[6] = (pi)*(delta)*(phi);
    K[7] = (pi)*(1.0-delta)*(phi);
    expand(dipc, dipx, 10, K, w, K2);
}

void getK0s(double pi, double* K0){// K0p = K'/K  K0pp = K''/K
    double* K0p = K0+3;
    double* K0pp= K0+6;
    K0[0]   =  2.0*(1.0-pi);                 K0[1] = 1.0;   K0[2]   = 2.0*pi;
    K0p[0]  = -2.0*(1.0-pi)*pi;              K0p[1] = 0.0;  K0p[2]  = 2.0*(1.0-pi)*pi;
    K0pp[0] = -2.0*(1.0-pi)*pi*(1.0-2.0*pi); K0pp[1] = 0.0; K0pp[2] = 2.0*(1.0-pi)*pi*(1.0-2.0*pi);
    //K0p[0]  = -pi;              K0p[1] = 0.0;  K0p[2]  = (1.0-pi);
    //K0pp[0] = -pi*(1.0-2.0*pi); K0pp[1] = 0.0; K0pp[2] = (1.0-pi)*(1.0-2.0*pi);
}

void getK2s(double pi, double delta, double phi, double* ktmp, double* w, double* K2){
    double* K2p   = K2+20;
    double* K2d   = K2+40;
    double* K2h   = K2+60;
    double* K2pp  = K2+80;
    double* K2pd  = K2+100;
    double* K2ph  = K2+120;
    double* K2dd  = K2+140;
    double* K2dh  = K2+160;
    double* K2hh  = K2+180;
    getK(pi,delta,phi,ktmp,w,K2);
    getdKdPi(pi,delta,phi,ktmp,w,K2p);
    getd2KdPi2(pi,delta,phi,ktmp,w,K2pp);
    getdKdDelta(pi,delta,phi,ktmp,w,K2d);
    getd2KdDelta2(pi,delta,phi,ktmp,w,K2dd);
    getdKdPhi(pi,delta,phi,ktmp,w,K2h);
    getd2KdPhi2(pi,delta,phi,ktmp,w,K2hh);
    getd2KdPidDelta(pi,delta,phi,ktmp,w,K2pd);
    getd2KdPidPhi(pi,delta,phi,ktmp,w,K2ph);
    getd2KdDeltadPhi(pi,delta,phi,ktmp,w,K2dh);
}

void getdKdPi(double pi, double delta, double phi, double* K, double* w, double* K2){
    // A allele
    K[0] = -pi*(1.0-pi)*(1.0-delta)*(1.0-phi);
    K[1] = -pi*(1.0-pi)*(delta)*(1.0-phi);
    K[2] = (pi)*(1.0-pi)*(1.0-delta)*(1.0-phi);
    K[3] = (pi)*(1.0-pi)*(delta)*(1.0-phi);
    // B allele
    K[4] = -pi*(1.0-pi)*(delta)*(phi);
    K[5] = -pi*(1.0-pi)*(1.0-delta)*(phi);
    K[6] = (pi)*(1.0-pi)*(delta)*(phi);
    K[7] = (pi)*(1.0-pi)*(1.0-delta)*(phi);
    expand(dipc, dipx, 10, K, w, K2);
}

void getd2KdPidDelta(double pi, double delta, double phi, double* K, double* w, double* K2){
    // A allele
    K[0] = pi*(1.0-pi)*delta*(1.0-delta)*(1.0-phi);
    K[1] = -pi*(1.0-pi)*delta*(1.0-delta)*(1.0-phi);
    K[2] = -(pi)*(1.0-pi)*delta*(1.0-delta)*(1.0-phi);
    K[3] = (pi)*(1.0-pi)*delta*(1.0-delta)*(1.0-phi);
    // B allele
    K[4] = -pi*(1.0-pi)*delta*(1.0-delta)*(phi);
    K[5] = pi*(1.0-pi)*delta*(1.0-delta)*(phi);
    K[6] = (pi)*(1.0-pi)*delta*(1.0-delta)*(phi);
    K[7] = -(pi)*(1.0-pi)*delta*(1.0-delta)*(phi);
    expand(dipc, dipx, 10, K, w, K2);
}
void getd2KdPidPhi(double pi, double delta, double phi, double* K, double* w, double* K2){
    // A allele
    K[0] = pi*(1.0-pi)*(1.0-delta)*phi*(1.0-phi);
    K[1] = pi*(1.0-pi)*(delta)*phi*(1.0-phi);
    K[2] = -(pi)*(1.0-pi)*(1.0-delta)*phi*(1.0-phi);
    K[3] = -(pi)*(1.0-pi)*(delta)*phi*(1.0-phi);
    // B allele
    K[4] = -pi*(1.0-pi)*(delta)*(phi)*(1.0-phi);
    K[5] = -pi*(1.0-pi)*(1.0-delta)*(phi)*(1.0-phi);
    K[6] = (pi)*(1.0-pi)*(delta)*(phi)*(1.0-phi);
    K[7] = (pi)*(1.0-pi)*(1.0-delta)*(phi)*(1.0-phi);
    expand(dipc, dipx, 10, K, w, K2);
}


void getd2KdPi2(double pi, double delta, double phi, double* K, double* w, double* K2){
    // A allele
    K[0] = -pi*(1.0-pi)*(1.0-2.0*pi)*(1.0-delta)*(1.0-phi);
    K[1] = -pi*(1.0-pi)*(1.0-2.0*pi)*(delta)*(1.0-phi);
    K[2] = (pi)*(1.0-pi)*(1.0-2.0*pi)*(1.0-delta)*(1.0-phi);
    K[3] = (pi)*(1.0-pi)*(1.0-2.0*pi)*(delta)*(1.0-phi);
    // B allele
    K[4] = -pi*(1.0-pi)*(1.0-2.0*pi)*(delta)*(phi);
    K[5] = -pi*(1.0-pi)*(1.0-2.0*pi)*(1.0-delta)*(phi);
    K[6] = (pi)*(1.0-pi)*(1.0-2.0*pi)*(delta)*(phi);
    K[7] = (pi)*(1.0-pi)*(1.0-2.0*pi)*(1.0-delta)*(phi);
    expand(dipc, dipx, 10, K, w, K2);
}



void getdKdDelta(double pi, double delta, double phi, double* K, double* w, double* K2){
    // A allele
    K[0] = -(1.0-pi)*(delta)*(1.0-delta)*(1.0-phi);
    K[1] =  (1.0-pi)*(delta)*(1.0-delta)*(1.0-phi);
    K[2] = -(pi)    *(delta)*(1.0-delta)*(1.0-phi);
    K[3] =  (pi)    *(delta)*(1.0-delta)*(1.0-phi);
    // B allele
    K[4] =  (1.0-pi)*(delta)*(1.0-delta)*(phi);
    K[5] = -(1.0-pi)*(delta)*(1.0-delta)*(phi);
    K[6] =  (pi)    *(delta)*(1.0-delta)*(phi);
    K[7] = -(pi)    *(delta)*(1.0-delta)*(phi);
    expand(dipc, dipx, 10, K, w, K2);
}
void getd2KdDeltadPhi(double pi, double delta, double phi, double* K, double* w, double* K2){
    // A allele
    K[0] =  (1.0-pi)*(delta)*(1.0-delta)*phi*(1.0-phi);
    K[1] = -(1.0-pi)*(delta)*(1.0-delta)*phi*(1.0-phi);
    K[2] =  (pi)    *(delta)*(1.0-delta)*phi*(1.0-phi);
    K[3] = -(pi)    *(delta)*(1.0-delta)*phi*(1.0-phi);
    // B allele
    K[4] =  (1.0-pi)*(delta)*(1.0-delta)*(phi)*(1.0-phi);
    K[5] = -(1.0-pi)*(delta)*(1.0-delta)*(phi)*(1.0-phi);
    K[6] =  (pi)    *(delta)*(1.0-delta)*(phi)*(1.0-phi);
    K[7] = -(pi)    *(delta)*(1.0-delta)*(phi)*(1.0-phi);
    expand(dipc, dipx, 10, K, w, K2);
}
void getd2KdDelta2(double pi, double delta, double phi, double* K, double* w, double* K2){
    // A allele
    K[0] = -(1.0-pi)*(delta)*(1.0-delta)*(1.0-2.0*delta)*(1.0-phi);
    K[1] =  (1.0-pi)*(delta)*(1.0-delta)*(1.0-2.0*delta)*(1.0-phi);
    K[2] = -(pi)    *(delta)*(1.0-delta)*(1.0-2.0*delta)*(1.0-phi);
    K[3] =  (pi)    *(delta)*(1.0-delta)*(1.0-2.0*delta)*(1.0-phi);
    // B allele
    K[4] =  (1.0-pi)*(delta)*(1.0-delta)*(1.0-2.0*delta)*(phi);
    K[5] = -(1.0-pi)*(delta)*(1.0-delta)*(1.0-2.0*delta)*(phi);
    K[6] =  (pi)    *(delta)*(1.0-delta)*(1.0-2.0*delta)*(phi);
    K[7] = -(pi)    *(delta)*(1.0-delta)*(1.0-2.0*delta)*(phi);
    expand(dipc, dipx, 10, K, w, K2);
}


void getdKdPhi(double pi, double delta, double phi, double* K, double* w, double* K2){
    // A allele
    K[0] = -(1.0-pi)*(1.0-delta)*(phi)*(1.0-phi);
    K[1] = -(1.0-pi)*(delta)    *(phi)*(1.0-phi);
    K[2] = -(pi)    *(1.0-delta)*(phi)*(1.0-phi);
    K[3] = -(pi)    *(delta)    *(phi)*(1.0-phi);
    // B allele
    K[4] = (1.0-pi)*(delta)    *(phi)*(1.0-phi);
    K[5] = (1.0-pi)*(1.0-delta)*(phi)*(1.0-phi);
    K[6] = (pi)    *(delta)    *(phi)*(1.0-phi);
    K[7] = (pi)    *(1.0-delta)*(phi)*(1.0-phi);
    expand(dipc, dipx, 10, K, w, K2);
}
void getd2KdPhi2(double pi, double delta, double phi, double* K, double* w, double* K2){
    // A allele
    K[0] = -(1.0-pi)*(1.0-delta)*(phi)*(1.0-phi)*(1.0-2.0*phi);
    K[1] = -(1.0-pi)*(delta)    *(phi)*(1.0-phi)*(1.0-2.0*phi);
    K[2] = -(pi)    *(1.0-delta)*(phi)*(1.0-phi)*(1.0-2.0*phi);
    K[3] = -(pi)    *(delta)    *(phi)*(1.0-phi)*(1.0-2.0*phi);
    // B allele
    K[4] = (1.0-pi)*(delta)    *(phi)*(1.0-phi)*(1.0-2.0*phi);
    K[5] = (1.0-pi)*(1.0-delta)*(phi)*(1.0-phi)*(1.0-2.0*phi);
    K[6] = (pi)    *(delta)    *(phi)*(1.0-phi)*(1.0-2.0*phi);
    K[7] = (pi)    *(1.0-delta)*(phi)*(1.0-phi)*(1.0-2.0*phi);
    expand(dipc, dipx, 10, K, w, K2);
}








double dip1[6];
double dip2[40];

void init_nbem(){
    dip1[0]=2.0;
    dip1[1]=1.0;
    dip1[2]=0.0;
    dip1[3]=0.0;
    dip1[4]=1.0;
    dip1[5]=2.0;
    
    long l, i, j;
    clear1(dip2,40);
    l=0;
    for(i=0; i<4; i++){
        for(j=i; j<4; j++){
            dip2[l+10*i]++;
            dip2[l+10*j]++;
            l++;
        }
    }
}

//K2={A1 B1 A2 B2 ... An Bn}
void expand(double* s1, double* s2, long N, double* K, double* w, double* K2){
    long i, j;
    for(i=0; i<8; i++){K[i] *= (2.0*oasr);}// 
    //for(i=0; i<8; i++){K[i] *= 2.0;}// total becomes 1 for diploid not 2!
    clear1(K2, 2*N);
    for(i=0; i<N; i++){
        for(j=0; j<2; j++){
            if(      s1[i*2+j]<=0.5 && s2[i*2+j]<=0.5){ K2[i*2]+=K[0]*w[j]; K2[i*2+1]+=K[4]*w[j];
            }else if(s1[i*2+j]<=0.5 && s2[i*2+j]> 0.5){ K2[i*2]+=K[1]*w[j]; K2[i*2+1]+=K[5]*w[j];
            }else if(s1[i*2+j]> 0.5 && s2[i*2+j]<=0.5){ K2[i*2]+=K[2]*w[j]; K2[i*2+1]+=K[6]*w[j];
            }else if(s1[i*2+j]> 0.5 && s2[i*2+j]> 0.5){ K2[i*2]+=K[3]*w[j]; K2[i*2+1]+=K[7]*w[j]; }
        }
    }
}

void expandForSun(double* s1, long N, double* K, double* K2){
    long i,j;
    clear1(K2, 2*N);
    for(i=0; i<N; i++){
        for(j=0; j<2; j++){
            if(s1[i*2+j]<=0.5 && j==0){K2[i*2]+=K[0]; K2[i*2+1]+=K[4];
            }else if(s1[i*2+j]<=0.5 && j==1){K2[i*2]+=K[1]; K2[i*2+1]+=K[5];
            }else if(s1[i*2+j]>0.5 && j==0){K2[i*2]+=K[2]; K2[i*2+1]+=K[6];
            }else if(s1[i*2+j]>0.5 && j==1){K2[i*2]+=K[3]; K2[i*2+1]+=K[7];}
        }
    }
}

void printMapGen(double* z, long N, long L){
    long i, j, mapg;
    double pp;
    if(L==3){
        for(i=0; i<N; i++){
            pp=0.0;
            for(j=0; j<3; j++){
                if(z[i*3+j] > pp){
                    mapg=j;
                    pp = z[i*3+j];
                }
            }
            fprintf(stderr, "%ld ", mapg);
        }
        fprintf(stderr, "\n");
    }else{
        for(i=0; i<N; i++){
            double p0 = z[i*10+0]+z[i*10+3]+z[i*10+7];
            double p1 = z[i*10+1]+z[i*10+4]+z[i*10+5]+z[i*10+8];
            double p2 = z[i*10+2]+z[i*10+6]+z[i*10+9];
            if(p0>=p1 && p0>=p2){
                fprintf(stderr, "0 ");
            }else if(p1>=p0 && p1>=p2){
                fprintf(stderr, "1 ");
            }else if(p2>=p0 && p2>=p1){
                fprintf(stderr, "2 ");
            }else{
                fprintf(stderr, "N ");
            }
        }
        fprintf(stderr, "\n");
    }
}




double vsum(double* x, double* y, long N){
    long i;
    double res=0.0;
    for(i=0; i<N; i++){
        res += x[i]*y[i];
    }
    return res;
}
double vsum2(double* x, double* y, long N){
    long i;
    double res=0.0;
    for(i=0; i<N; i++){
        res += x[i]*y[i]*y[i];
    }
    return res;
}

double vsum3(double* x, double* y, double* z, long N){
    long i;
    double res=0.0;
    for(i=0; i<N; i++){
        res += x[i]*y[i]*z[i];
    }
    return res;
}



double vsum3old(double* x, double* y, long J, long M){
    long j,m;
    double res=0.0;
    for(j=0; j<J; j++){
        res += x[j*(M+1)+M]*y[j*M+0]*y[j*M+1];
        for(m=0; m<M; m++){
            res += x[j*(M+1)+m]*y[j*M+m]*y[j*M+m];
        }
    }
    return res;
}


double vsum4(double* x, double* y1, double* y2, long J, long M){
    long j,m;
    double res=0.0;
    for(j=0; j<J; j++){
        res += 0.5*x[j*(M+1)+M]*y1[j*M+0]*y2[j*M+1];
        res += 0.5*x[j*(M+1)+M]*y2[j*M+0]*y1[j*M+1];
        for(m=0; m<M; m++){
            res += x[j*(M+1)+m]*y1[j*M+m]*y2[j*M+m];
        }
    }
    return res;
}

double vsum5(double* x, double* y, long N){
    long i;
    double res=0.0;
    for(i=0; i<N; i++){
        res += (x[i*2]+x[i*2+1])*y[i];
    }
    return res;
}



double getCE(double* A, double* H1, double* h1, long j0){
    long k;
    double res=0.0;
    for(k=jg[j0]; k<jg[j0]+kg[j0]; k++){
        if(h1[j0]>0.0){
            res += (A[k]+A[k+10])*H1[k]/h1[j0];
        }
    }
    return res;
}
double getCE2(double* K2, double* K2p, double* A, double* H1, double* h1, long j0){
    long k;
    double res=0.0;
    for(k=jg[j0]; k<jg[j0]+kg[j0]; k++){
        if(h1[j0]>0.0){
            res += (A[k]*K2p[k*2]/K2[k*2]+A[k+10]*K2p[k*2+1]/K2[k*2+1])*H1[k]/h1[j0];
        }
    }
    return res;
}
double getCV(double* A, double* H1, double* h1, long j0){
    long k;
    double res=0.0;
    double mm = getCE(A, H1, h1, j0);
    for(k=jg[j0]; k<jg[j0]+kg[j0]; k++){
        if(h1[j0]>0.0){
            res += pow(A[k]+A[k+10]-mm, 2.0)*H1[k]/h1[j0];
        }
    }
    return res;
}
double getCV2(double* K2, double* K2p, double* A, double* H1, double* h1, long j0){
    long k;
    double res=0.0;
    double mm = getCE2(K2, K2p, A, H1, h1, j0);
    for(k=jg[j0]; k<jg[j0]+kg[j0]; k++){
        if(h1[j0]>0.0){
            res += pow(A[k]*K2p[k*2]/K2[k*2]+A[k+10]*K2p[k*2+1]/K2[k*2+1]-mm, 2.0)*H1[k]/h1[j0];
        }
    }
    return res;
}

double getCCV(double* K2, double* K2p, double* A, double* H1, double* h1, long j0){
    long k;
    double res=0.0;
    double mm = getCE(A, H1, h1, j0);
    double mmp = getCE2(K2, K2p, A, H1, h1, j0);
    for(k=jg[j0]; k<jg[j0]+kg[j0]; k++){
        if(h1[j0]>0.0){
            res += (A[k]+A[k+10]-mm)*(A[k]*K2p[k*2]/K2[k*2]+A[k+10]*K2p[k*2+1]/K2[k*2+1]-mmp)*H1[k]/h1[j0];
        }
    }
    return res;
}
double getCCV2(double* K2, double* K2p, double* K2d, double* A, double* H1, double* h1, long j0){
    long k;
    double res=0.0;
    double mmp = getCE2(K2, K2p, A, H1, h1, j0);
    double mmd = getCE2(K2, K2d, A, H1, h1, j0);
    for(k=jg[j0]; k<jg[j0]+kg[j0]; k++){
        if(h1[j0]>0.0){
            res += (A[k]*K2p[k*2]/K2[k*2]+A[k+10]*K2p[k*2+1]/K2[k*2+1]-mmp) * (A[k]*K2d[k*2]/K2[k*2]+A[k+10]*K2d[k*2+1]/K2[k*2+1]-mmd)*H1[k]/h1[j0];
        }
    }
    return res;
}

// work : length 45
// M = 2
void getInformation(double* y, double* Y, double* h0, double* H0, double* h1, double* H1, double* H2, double* ki, double* dki, double* ki2, double* K0, double* K2, double* K, double* km, long Lx, long N0, long J0, long J, double beta, double th, double pi, double delta, double phi, double asNonasRatio, double* work, double* hess, integer* ipiv, double* a, double* A){
    
    long i, i0, j, l, i1;
    double lmij, thij, denom, dlAB, thijA, thijB, kij, k_p, k_d, k_h;
    
    double mu=exp(beta);
    
    double* K2p   = K2+20;
    double* K2d   = K2+40;
    double* K2h   = K2+60;
    
    double* K0p = K0+3;
    
    double dLx = (double)Lx;
    if(Lx==0){dLx=1.0;}
    
    double* ESp;
    double* ESd;
    double* ESh;
    double* ESb;
    double* ESt;
    double* VSp;
    double* VSd;
    double* VSh;
    double* VSt;
    
    double* VSpd;
    double* VSph;
    double* VSpt;
    double* VSdh;
    double* VSdt;
    double* VSht;
    
    ESp = work;
    ESd = work+3;
    ESh = work+6;
    ESb = work+9;
    ESt = work+12;
    
    VSp = work+15;
    VSd = work+18;
    VSh = work+21;
    VSt = work+24;
    VSpd = work+27;
    VSph = work+30;
    VSpt = work+33;
    VSdh = work+36;
    VSdt = work+39;
    VSht = work+42;
    
    
    // as counts    
    
    getK0s(pi, K0);
    
    
    clear1(a, J0*5);
    clear1(A, 5*J);
    
    mu=exp(beta);
    for(i0=0; i0<N0; i0++){
        
        getK2s(pi,delta,phi,K,ki2+i0*2,K2);
        
        for(j=0; j<J0; j++){
            VSp[j]=VSd[j]=VSh[j]=VSt[j]=VSpd[j]=VSph[j]=VSpt[j]=VSdh[j]=VSdt[j]=VSht[j]=0.0;
            if(h0[i0*J0+j]>0.0){
                if(Lx>0){
                    kij = ((1.0-asNonasRatio)*K0[j]  + asNonasRatio/oasr/dLx*vsum5(K2  +jg[j]*2, H2+i0*10+jg[j], kg[j])/h0[i0*J0+j]);
                    k_p = ((1.0-asNonasRatio)*K0p[j] + asNonasRatio/oasr/dLx*vsum5(K2p +jg[j]*2, H2+i0*10+jg[j], kg[j])/h0[i0*J0+j])/kij;
                    k_d = (                            asNonasRatio/oasr/dLx*vsum5(K2d +jg[j]*2, H2+i0*10+jg[j], kg[j])/h0[i0*J0+j])/kij;
                    k_h = (                            asNonasRatio/oasr/dLx*vsum5(K2h +jg[j]*2, H2+i0*10+jg[j], kg[j])/h0[i0*J0+j])/kij;
                }else{
                    kij = K0[j];
                    k_p = K0[j]/kij;
                    k_d = k_h = 0.0;
                }
                thij = th*ki[i0]*kij;
                lmij = mu*ki[i0]*dki[i0]*kij;
                
                denom = 1.0 + lmij/thij;
                a[j]      = (y[i0]-lmij)/denom;
                a[j+J0]   = thij * (digamma( thij + y[i0]) - digamma( thij) + log(thij) - log( lmij + thij) + (lmij - y[i0])/(lmij + thij) );
                
                ESp[j] = (a[j]+a[j+J0])*k_p;
                ESd[j] = (a[j]+a[j+J0])*k_d;
                ESh[j] = (a[j]+a[j+J0])*k_h;
                ESb[j] = a[j];
                ESt[j] = a[j+J0];
            }
        }
        
        
        for(l=0; l<Lx; l++){
            i=i0+l*N0;
            if(randomize>0){i1=i+N0;}else{i1=i0;}
            for(j=0; j<J; j++){
                thijA = th*ki[i1]*K2[j*2]*km[i/N0];
                thijB = th*ki[i1]*K2[j*2+1]*km[i/N0];
                
                dlAB  = digamma(thijA+thijB)  - digamma(thijA+thijB + Y[i*2]+Y[i*2+1]) ;
                
                A[j]     = thijA * (digamma(thijA + Y[i*2])    - digamma(thijA)  + dlAB );
                A[j+J]   = thijB * (digamma(thijB + Y[i*2+1])  - digamma(thijB)  + dlAB );
            }
            for(j=0; j<J0; j++){
                if(h1[i0*3+j]>0.0){
                    ESp[j] += getCE2(K2, K2p, A, H1 + i0*J + l*N0*J, h1+i0*3, j);
                    ESd[j] += getCE2(K2, K2d, A, H1 + i0*J + l*N0*J, h1+i0*3, j);
                    ESh[j] += getCE2(K2, K2h, A, H1 + i0*J + l*N0*J, h1+i0*3, j);
                    ESt[j] += getCE(          A, H1 + i0*J + l*N0*J, h1+i0*3, j);
                    VSp[j] += getCV2(K2, K2p, A, H1 + i0*J + l*N0*J, h1+i0*3, j);
                    VSd[j] += getCV2(K2, K2d, A, H1 + i0*J + l*N0*J, h1+i0*3, j);
                    VSh[j] += getCV2(K2, K2h, A, H1 + i0*J + l*N0*J, h1+i0*3, j);
                    VSt[j] += getCV(          A, H1 + i0*J + l*N0*J, h1+i0*3, j);
                    
                    VSpd[j] += getCCV2(K2, K2p, K2d, A, H1 + i0*J + l*N0*J, h1+i0*3, j);
                    VSph[j] += getCCV2(K2, K2p, K2h, A, H1 + i0*J + l*N0*J, h1+i0*3, j);
                    VSdh[j] += getCCV2(K2, K2d, K2h, A, H1 + i0*J + l*N0*J, h1+i0*3, j);
                    
                    VSpt[j] += getCCV(K2, K2p, A, H1 + i0*J + l*N0*J, h1+i0*3, j);
                    VSdt[j] += getCCV(K2, K2d, A, H1 + i0*J + l*N0*J, h1+i0*3, j);
                    VSht[j] += getCCV(K2, K2h, A, H1 + i0*J + l*N0*J, h1+i0*3, j);
                }
            }
        }
        
        hess[0] -= vsum(h1+i0*3, VSp, 3)  + vsum2(h1+i0*3, ESp, 3);
        hess[1] -= vsum(h1+i0*3, VSpd, 3) + vsum3(h1+i0*3, ESp, ESd, 3);
        hess[2] -= vsum(h1+i0*3, VSph, 3) + vsum3(h1+i0*3, ESp, ESh, 3);
        hess[3] -=                          vsum3(h1+i0*3, ESp, ESb, 3);
        hess[4] -= vsum(h1+i0*3, VSpt, 3) + vsum3(h1+i0*3, ESp, ESt, 3);
        
        hess[11] -= vsum(h1+i0*3, VSd, 3)  + vsum2(h1+i0*3, ESd, 3);
        hess[12] -= vsum(h1+i0*3, VSdh, 3) + vsum3(h1+i0*3, ESd, ESh, 3);
        hess[13] -=                          vsum3(h1+i0*3, ESd, ESb, 3);
        hess[14] -= vsum(h1+i0*3, VSdt, 3) + vsum3(h1+i0*3, ESd, ESt, 3);
        
        hess[22]-= vsum(h1+i0*3, VSh, 3)  + vsum2(h1+i0*3, ESh, 3);
        hess[23]-=                          vsum3(h1+i0*3, ESh, ESb, 3);
        hess[24]-= vsum(h1+i0*3, VSht, 3) + vsum3(h1+i0*3, ESh, ESt, 3);
        
        hess[33]-=                          vsum2(h1+i0*3, ESb, 3);
        hess[34]-=                          vsum3(h1+i0*3, ESb, ESt, 3);
        
        hess[44]-= vsum(h1+i0*3, VSt, 3)  + vsum2(h1+i0*3, ESt, 3);
        
        /*
         clear1(work, 50);
         for(j=0; j<J; j++){
         work[j]    += (a[ig[j]] + a[ig[j]+J0])*K0p[ig[j]]/K0[ig[j]] + A[j]*K2p[j*2]/K2[j*2] + A[j+J]*K2p[j*2+1]/K2[j*2+1];
         work[j+10] += A[j]*K2d[j*2]/K2[j*2] + A[j+J]*K2d[j*2+1]/K2[j*2+1];
         work[j+20] += A[j]*K2h[j*2]/K2[j*2] + A[j+J]*K2h[j*2+1]/K2[j*2+1];
         work[j+30] += a[ig[j]];
         work[j+40] += a[ig[j]+J0] + A[j] + A[j+J];
         }
         for(j=0; j<5; j++){wshift(work+j*10, H1+i*J, J);}
         */
        //H1[i*J+j] * 
    }
    
    
    inverseLapack(hess, 5, ipiv, hess+25);
}


double getIA2(double* gl, long N){
    long i;
    double IA=0.0;
    double tot=0.0;
    double th=0.0;
    for(i=0; i<N; i++){
        th += ( (gl[i*10+2]+gl[i*10+6]+gl[i*10+9]) + (gl[i*10+1]+gl[i*10+4]+gl[i*10+5]+gl[i*10+8])/2.0 );
        tot+= 1.0;
    }
    th /= tot;
    if(th<=0.0 || th>=1.0){
        return 1.0;
    }
    for(i=0; i<N; i++){
        IA += (4.0*(gl[i*10+2]+gl[i*10+6]+gl[i*10+9])+(gl[i*10+1]+gl[i*10+4]+gl[i*10+5]+gl[i*10+8]) - pow(2.0*(gl[i*10+2]+gl[i*10+6]+gl[i*10+9])+(gl[i*10+1]+gl[i*10+4]+gl[i*10+5]+gl[i*10+8]),2.0))/2.0/th/(1.0-th);
    }
    IA /= tot;
    return 1.0-IA;
}

double getIA1(double* gl, long N){
    long i;
    double IA=0.0;
    double tot=0.0;
    double th=0.0;
    for(i=0; i<N; i++){if(gl[i*3]>=0.0){
        th += (gl[i*3+2]+gl[i*3+1]/2.0);
        tot+= 1.0;
    }}
    th /= tot;
    if(th<=0.0 || th>=1.0){
        return 1.0;
    }
    for(i=0; i<N; i++){if(gl[i*3]>=0.0){
        IA += (4.0*gl[i*3+2]+gl[i*3+1] - pow(2.0*gl[i*3+2]+gl[i*3+1],2.0))/2.0/th/(1.0-th);
    }}
    IA /= tot;
    return 1.0-IA;
}

double getCov(double* h0, double* h1, long N0){
    long j;
    double m0=0.0, m1=0.0, v0=0.0, v1=0.0, v01=0.0;
    for(j=0; j<N0; j++){
        m0 += 2.0*h0[j*3+2] + h0[j*3+1];
        m1 += 2.0*h1[j*3+2] + h1[j*3+1];
    }
    m0 /= (double)N0;
    m1 /= (double)N0;
    for(j=0; j<N0; j++){
        v0 += 4.0*h0[j*3+2] + h0[j*3+1];
        v1 += 4.0*h1[j*3+2] + h1[j*3+1];
        v01+= 4.0*h1[j*3+2]*h0[j*3+2] + 2.0 * (h1[j*3+2]*h0[j*3+1] + h1[j*3+1]*h0[j*3+2]) + h1[j*3+1]*h0[j*3+1];
    }
    v0 /= (double)N0;
    v1 /= (double)N0;
    v01/= (double)N0;
    
    return pow(v01-m0*m1, 2.0)/(v0-m0*m0)/(v1-m1*m1);
}

double getCov2(double* h0, double* h1, long N0){
    long j;
    double m0=0.0, m1=0.0, v0=0.0, v1=0.0, v01=0.0;
    for(j=0; j<N0; j++){
        m0 += 2.0*getAA(h0+j*10) + getAR(h0+j*10);
        m1 += 2.0*getAA(h1+j*10) + getAR(h1+j*10);
    }
    m0 /= (double)N0;
    m1 /= (double)N0;
    for(j=0; j<N0; j++){
        v0 += 4.0*getAA(h0+j*10) + getAR(h0+j*10);
        v1 += 4.0*getAA(h1+j*10) + getAR(h1+j*10);
        v01+= 4.0*getAA(h1+j*10)*getAA(h0+j*10) + 2.0 * (getAA(h1+j*10)*getAR(h0+j*10) + getAR(h1+j*10)*getAA(h0+j*10)) + getAR(h1+j*10)*getAR(h0+j*10);
    }
    v0 /= (double)N0;
    v1 /= (double)N0;
    v01/= (double)N0;
    
    return pow(v01-m0*m1, 2.0)/(v0-m0*m0)/(v1-m1*m1);
}


double getAA(double* gl){
    return gl[2]+gl[6]+gl[9];
}
double getAR(double* gl){
    return gl[1]+gl[4]+gl[5]+gl[8];
}


void gencall(double* H1, double* Z, double* Zx, long N, long L, double* exon){
    long i, l, il;
    double ratio1, ratio2;
    double* z;
    double* zx;
    long ll=0;
    for(l=0;l<L;l++){
        z=Z+N*2*l;
        if(fabs(exon[l])>1.5){
            zx = Zx + N*2*ll;
            for(i=0; i<N; i++){
                if((1.0-z[i*2])*z[i*2+1] + z[i*2]*(1.0-z[i*2+1])>0.0){
                    ratio2 = (1.0-z[i*2])*z[i*2+1]/((1.0-z[i*2])*z[i*2+1] + z[i*2]*(1.0-z[i*2+1])); // p(0|1)/(p(0|1)+p(1|0))
                    ratio1 = 1.0-ratio2;
                }else{
                    ratio1 = ratio2 = 0.0;
                }
                il = ll*N*10 + i*10;
                zx[i*2]   = ratio1*(H1[il+1]+H1[il+4]+H1[il+5]+H1[il+8]) + (H1[il+2] + H1[il+6] + H1[il+9]);
                zx[i*2+1] = ratio2*(H1[il+1]+H1[il+4]+H1[il+5]+H1[il+8]) + (H1[il+2] + H1[il+6] + H1[il+9]);
            }
            ll++;
        }
    }
}

void gencallAlt(double* H1, double* h1, double* Z, double* Zx, long N, long L, long Lx, double* exon, long csnp){
    double* work; work=(double*)calloc(16,sizeof(double));
    
    long i, l, il;
    //double* z;
    double* zx;
    double* zc;
    double* zc0;
    double tmp;
    zc0= Z+N*2*csnp;
    zc = Zx+N*2*Lx;
    for(i=0; i<N; i++){
        H3toZ(h1+3*i, zc+2*i);
        if((zc0[2*i]-zc0[2*i+1])*(zc[2*i]-zc[2*i+1])<0.0){
            tmp = zc[2*i];
            zc[2*i]=zc[2*i+1];
            zc[2*i+1]=tmp;
        }
    }
    long ll=0;
    //double p01, p10, q01, q10, q11;
    for(l=0;l<L;l++){
        //z=Z+N*2*l;
        if(fabs(exon[l])>1.5){
            zx = Zx + N*2*ll;
            for(i=0; i<N; i++){
                il = ll*N*10 + i*10;
                
                H3toZx(H1+il, zc+i*2, zx+i*2, work);
                
            }
            ll++;
        }
    }
    free(work);
}





void H3toZx(double* H1, double* zc, double* zx, double* work){
    sum2one(H1,10);
    
    int itr;
    //double p01 = (1.0-zc[0])*zc[1];
    //double p10 = zc[0]*(1.0-zc[1]);
    double p01, p10, q01, q10;
    
    //double p00 = H1[0]+H1[1]+H1[2];
    double p11 = H1[7]+H1[8]+H1[9];
    double phet= H1[3]+H1[4]+H1[5]+H1[6];
    
    double q00 = H1[0]+H1[3]+H1[7];
    double q11 = H1[2]+H1[6]+H1[9];
    double qhet= H1[1]+H1[4]+H1[5]+H1[8];
    
    double ph, de;
    ph = 0.51;
    de = 0.51;
    double* Z;
    Z=work;
    
    // exceptions 
    if(q11>=1.0){
        zx[0]=zx[1]=1.0; return;
    }else if(q00>=1.0){
        zx[0]=zx[1]=0.0; return;
    }
    double lkhd;
    double lkhd0 = -1.0e31;
    for(itr=0; itr<100; itr++){
        
        //Z[0]=p00*q00;
        Z[1]=H1[1]*(1.0-de);// p00*q01
        Z[2]=H1[1]*de;      // p00*q10
        //Z[3]=p00*q11;
        
        Z[4] =H1[3]*(1.0-ph);                                      // p01*q00;
        Z[5] =H1[3]*ph;                                            // p10*q00;
        if((1.0-ph)*de + ph*(1.0-de)>0.0){
            Z[6] =H1[4]*(1.0-ph)*de      /((1.0-ph)*de + ph*(1.0-de)); // p01*q10;
            Z[7] =H1[4]*ph*(1.0-de)      /((1.0-ph)*de + ph*(1.0-de)); // p10*q01;
        }else{Z[6]=Z[7]=0.0;}
        if((1.0-ph)*(1.0-de)+ph*de>0.0){
            Z[8] =H1[5]*(1.0-ph)*(1.0-de)/((1.0-ph)*(1.0-de)+ph*de);   // p01*q01;
            Z[9] =H1[5]*     ph *     de /((1.0-ph)*(1.0-de)+ph*de);   // p10*q10;
        }else{Z[8]=Z[9]=0.0;}
        Z[10]=H1[6]*(1.0-ph);                                      // p01*q11;
        Z[11]=H1[6]*ph;                                            // p10*q11;
        
        //Z[12]=p11*q11;
        Z[13]=H1[8]*(1.0-de);
        Z[14]=H1[8]*de;
        //Z[15]=p11*q11;
        
        ph = (Z[5]+Z[7]+Z[9]+Z[11])/(Z[4]+Z[6]+Z[8]+Z[10] + Z[5]+Z[7]+Z[9]+Z[11]);
        de = (Z[2]+Z[6]+Z[9]+Z[14])/(Z[1]+Z[7]+Z[8]+Z[13] + Z[2]+Z[6]+Z[9]+Z[14]);
        
        p01=phet*(1.0-ph);
        p10=phet*(ph);
        q01=qhet*(1.0-de);
        q10=qhet*(de);
        
        lkhd =  de*(1.0-ph)+(1.0-de)*ph>0.0 ? H1[4]*log(de*(1.0-ph)+(1.0-de)*ph):0.0;
        lkhd += (1.0-de)*(1.0-ph)+de*ph>0.0 ? H1[5]*log((1.0-de)*(1.0-ph)+de*ph):0.0;
        if(fabs(lkhd-lkhd0)<1e-8){break;}else{lkhd0=lkhd;}
        
    }
    
    if((p01+p11>p10+p11 && zc[0]>zc[1]) || (p01+p11<p10+p11 && zc[0]<zc[1])){
        zx[1] = q10 + q11;
        zx[0] = q01 + q11;
    }else{
        zx[0] = q10 + q11;
        zx[1] = q01 + q11;
    }
}


void H3toZanalitic(double* h3, double* z){
    double a=h3[0]/(h3[0]+h3[1]+h3[2]);
    double b=h3[1]/(h3[0]+h3[1]+h3[2]);
    z[0]=(2.0*a + b - sqrt(-4.0*a + 4.0*pow(a,2.0) + 4.0*a*b + pow(b,2.0)))/2.0;
    z[1]=(2.0*a + b + sqrt(-4.0*a + 4.0*pow(a,2.0) + 4.0*a*b + pow(b,2.0)))/2.0;
}

void H3toZ(double* h3, double* z){
    int i;
    double tot=asum(h3,3);
    if(tot<=0.0){z[0]=z[1]=0.0; return;}
    h3[0]/=tot;
    h3[1]/=tot;
    h3[2]/=tot;
    if(h3[2]+h3[1]<=0.0){z[0]=z[1]=0.0; return;}
    if(h3[0]+h3[1]<=0.0){z[0]=z[1]=1.0; return;}
    if(h3[0]+h3[2]<=0.0){z[0]=0.0; z[1]=1.0; return;}
    if(h3[0]<0.0){z[0]=z[1]=0.5;return;}
    
    double phi, phi0, p, q;
    if(h3[2]>h3[0]){
        p=(h3[2]+h3[1]/2.0)*1.001;
        q=(h3[2]+h3[1]/2.0)*0.999;
        phi0 = phi = p*(1.-q)/(p*(1.-q)+q*(1.-p));
        for(i=0; i<10000; i++){
            p = h3[2]+h3[1]*phi;
            q = h3[2]+h3[1]*(1.-phi);
            phi=p*(1.-q)/(p*(1.-q)+q*(1.-p));
            if(fabs(phi-phi0)<1e-10){break;}else{phi0=phi;}
        }
    }else{
        p=(h3[0]+h3[1]/2.0)*0.999;
        q=(h3[0]+h3[1]/2.0)*1.001;
        phi0 = phi = p*(1.-q)/(p*(1.-q)+q*(1.-p));
        for(i=0; i<10000; i++){
            p = h3[0]+h3[1]*phi;
            q = h3[0]+h3[1]*(1.-phi);
            phi=p*(1.-q)/(p*(1.-q)+q*(1.-p));
            if(fabs(phi-phi0)<1e-10){break;}else{phi0=phi;}
        }
        p=1.0-p;
        q=1.0-q;
    }
    if(p>q){double tmp=p; p=q; q=tmp;}
    z[0]=p;
    z[1]=q;
}

void H3toZx2(double* H1, double* zc, double* zx, double* work){
    sum2one(H1,10);
    
    //double p01 = (1.0-zc[0])*zc[1];
    //double p10 = zc[0]*(1.0-zc[1]);
    double p01, p10, q01, q10;
    
    //double p00 = H1[0]+H1[1]+H1[2];
    //double p11 = H1[7]+H1[8]+H1[9];
    //double phet= H1[3]+H1[4]+H1[5]+H1[6];
    
    double q00 = H1[0]+H1[3]+H1[7];
    double q11 = H1[2]+H1[6]+H1[9];
    double qhet= H1[1]+H1[4]+H1[5]+H1[8];
    
    double qs[3]={q00,qhet,q11};
    H3toZ(qs, zx);
    
    //double* Z;
    //Z=work;
    
    // exceptions 
    if(q11>=1.0){
        zx[0]=zx[1]=1.0; return;
    }else if(q00>=1.0){
        zx[0]=zx[1]=0.0; return;
    }
    
    double lkhd1, lkhd2;
    
    p01=(1.0-zc[0])*zc[1];
    p10=zc[0]*(1.0-zc[1]);
    q01=(1.0-zx[0])*zx[1];
    q10=zx[0]*(1.0-zx[1]);
    
    lkhd1 = H1[4]*log(p01*q10+p10*q01) + H1[5]*log(p01*q01+p10*q10) ;
    
    p01=(1.0-zc[0])*zc[1];
    p10=zc[0]*(1.0-zc[1]);
    q01=zx[0]*(1.0-zx[1]);
    q10=(1.0-zx[0])*zx[1];
    
    lkhd2 = H1[4]*log(p01*q10+p10*q01) + H1[5]*log(p01*q01+p10*q10) ;
    
    
    if(lkhd1<lkhd2){
        double tmp=zx[0];
        zx[0] = zx[1];
        zx[1] = tmp;
    }
}



void randomPerm(double* y, double* Y, double* Z, double* ki, double* ki2, double* exon, int L, int N, int m){
    int i, l;
    double* Zr;
    Zr = Z+N*(L+m+1)*2;
    double* dki;
    dki = ki+N+m*N;
    double* yt;  yt = (double*)calloc(N*3, sizeof(double));
    int* rord;   rord = (int*)calloc(N, sizeof(int));
    int* rord12; rord12 = (int*)calloc(2, sizeof(int));
    for(i=0; i<N; i++){rord[i]=i;}
    rord12[0]=0; rord12[1]=1;
    //Zr = (double*)calloc(N*2*m, sizeof(double));
    
    // as permutation
    int ll=0;
    for(l=0;l<L;l++){
        if(fabs(exon[l])>1.5){// in a feature
            //for(i=0; i<N; i++){rord[i]=i;}
            getRandomOrder((int)N, rord);
            for(i=0; i<N; i++){
                //getRandomOrder(2, rord12);
                if(rord12[0]==0){
                    Zr[N*2*ll+i*2+0] = Z[N*2*l+rord[i]*2+0];
                    Zr[N*2*ll+i*2+1] = Z[N*2*l+rord[i]*2+1];
                }else{
                    Zr[N*2*ll+i*2+0] = 1.0 - Z[N*2*l+rord[i]*2+0];
                    Zr[N*2*ll+i*2+1] = 1.0 - Z[N*2*l+rord[i]*2+1];
                }
                yt[i*2+0]      = Y[ll*2*N+rord[i]*2+rord12[0]];
                yt[i*2+1]      = Y[ll*2*N+rord[i]*2+rord12[1]];
                ki[(ll+1)*N+i] = ki[rord[i]];
                ki2[(ll+1)*N*2+i*2+0] = ki2[rord[i]*2+rord12[0]];
                ki2[(ll+1)*N*2+i*2+1] = ki2[rord[i]*2+rord12[1]];
            }
            cblas_dcopy(N*2, yt, 1, Y+ll*2*N, 1);
            ll++;
        }
    }
    // pop permutation
    //for(i=0; i<N; i++){rord[i]=i;}
    getRandomOrder((int)N, rord);
    for(i=0; i<N; i++){
        yt[i]     =  y[rord[i]];
        yt[i+N]   =  ki[rord[i]];
        yt[i+N*2] =  dki[rord[i]];
    }
    cblas_dcopy(N, yt,     1, y,   1);
    cblas_dcopy(N, yt+N,   1, ki,  1);
    cblas_dcopy(N, yt+N*2, 1, dki, 1);
    
    // vector free
    free(yt);
    free(rord);
    free(rord12);
}

