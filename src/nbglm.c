#include "util.h"
#include "nbem.h"

#define EPSBETA 1.0e-5
#define EPSTHETA 1.0e-8
#define EPSLKHD 1.0e-5


void ForwardSolveUnused(double* R, long p, double* y, double* x, long ldx){
    // To solve R^T x = y
    // ldy : leading dim of y
    long i, j;
	clear1(x, p);
	if(R[0]>0.0){ x[0] = y[0]/R[0];}
    
	double tm = 0.0;
    for(i=1; i<p; i++){
        tm = 0.0;
        for(j=0; j<i; j++){
            tm += R[p*i+j]*x[j];
        }
	if(R[p*i+i]>0.0){x[i] = (y[i]-tm)/R[p*i+i];}
    }
}



void BackSolve(double* R, long p, double* y, double* x){
    // To solve R x = y
    
    long i, j;
	clear1(x, p);
	//if(R[p*p-1]>0.0){ 
	x[p-1] = y[p-1]/R[p*p-1];
	//}
    
	double tm = 0.0;
    for(i=p-2; i>=0; i--){
        tm = 0.0;
        for(j=p-1; j>i; j--){
            tm += R[p*j+i]*x[j];
        }
		//if(R[p*i+i]>0.0){
		x[i] = (y[i]-tm)/R[p*i+i];
		//}
    }
}



void QR(double* X, long n, long p, double* R){
    // QR decomp using Gram-Schmidt orthogonalization
    // X is replaced by Q
    
	long i, j;
    
	clear1(R, p*p);
	
	for(i=0; i<p; i++){
		// R3
		for(j=0; j<i; j++){
			R[i*p+j] = cblas_ddot(n, X+n*j, 1, X+n*i, 1);
			cblas_daxpy(n, -R[i*p+j], X+n*j, 1, X+n*i, 1);
		}
		R[i*p+i] = cblas_dnrm2(n, X+n*i, 1);
		if(R[i*p+i]>0.0){
			cblas_dscal(n, 1.0/R[i*p+i], X+n*i, 1);
		}
	}
}


double getLkhdGlm(double* y, double* X, double* ki, double* w, long N0, long P, double* lambda, double* beta, double th){
    double lkhd=0.0;
    double lmij, thij;
    int i;
    for(i=0; i<N0; i++){
        lmij = lambda[i];
#ifdef NOK
        thij = th;
#else
        thij = th*ki[i];
#endif
        lkhd += w[i] * (lgamma(thij+y[i]) - (y[i]+thij)*log(lmij+thij) - lgamma(thij) + thij*log(thij) - lgamma(y[i]+1.0) + y[i]*log(lmij));
        
        //lgamma(thij+y[i0]) - lgamma(y[i0]+1.0) - lgamma(thij) + y[i0]*log(lmij) + (thij)*log(thij) - (y[i0]+thij)*log(lmij+thij);
        
    }
    
    lkhd += ((double)P)*log(th)/2.0 - ((double)P)*log(2.0*M_PI)/2.0 - ((double)P)*log(sigma2)/2.0 - th*pow(beta[0]-beta0,2.0)/sigma2/2.0;
    for(i=1; i<P; i++){lkhd += - th*pow(beta[i],2.0)/sigma2/2.0;}
    
    lkhd += kappa*log(omega/2.0)/2.0 + (kappa/2.0-1.0)*log(th) - (omega/2.0)*th - lgamma(kappa/2.0);
    
    
    
    //lkhd -= log(s)/2.0;
    ///lkhd += (0.5 + kappa/2.0-1.0)*log(th) - 0.5*log(2.0*M_PI);
    //lkhd += kappa/2.0*log(omega/2.0) - omega*th/2.0 - lgamma(kappa/2.0);
    //lkhd -= th*pow(beta[0]-beta0,2.0)/sigma2/2.0;
    
    
    return lkhd;
}



double olm(double* y, double* X, double* ki, double* w, long N0, long P, double* beta, double* work){
    double* Q;
    double* R;
    double* z;
    double* gamma;
    double dN=0.0;
    clear1(work, N0*P+P*P+N0+P);
    Q     = work;
    R     = work+N0*P;
    z     = work+N0*P+P*P;
    gamma = work+N0*P+P*P+N0;
    long i, j;
    for(i=0; i<N0; i++){
        z[i]=log((y[i]+1.0)/ki[i])*w[i];
        for(j=0; j<P; j++){
            Q[i+j*N0] = X[i+j*N0]*w[i];
        }
        dN += w[i];
    }
    QR(Q, N0, P, R);
    cblas_dgemv(CblasColMajor, CblasTrans, N0, P, 1.0, Q, N0, z, 1, 0.0, gamma, 1);
    BackSolve(R, P, gamma, beta);
    double v=0.0;
    for(i=0; i<N0; i++){
        for(j=0; j<P; j++){
            z[i] -= X[i+j*N0]*beta[j]*w[i];
        }
        v += z[i]*z[i]/(dN-(double)P);
    }
    return v;
}


// M : num of alleles not num of xSNP
// J : num of genotype statt
double nbglm(double* y, double* X, double* ki, double* dki, double* w, long N0, long P, double* beta, double* ptheta, int fixTheta, double* work) {
    
    long i, j, itr, itr_beta, itr_theta;
    double stepsize = 1.0;
    ///double A0=(double)(N0*M0*J0);
    double A0=(double)(N0);
    
    if(verbose3>0){fprintf(stderr, "N0=%ld P=%ld beta=%lf th=%lf s=%lf kappa=%lf omega=%lf\n", N0, P, beta[0], ptheta[0], sigma2, kappa, omega);}
    
    double th, th1;
    th = ptheta[0];
    // init lkhd
    double lkhd, lkhd0;
    lkhd0 = lkhd = -1.0/0.0;
    
    double lkhd_all = -1.0/0.0;
    
    //sigma2=10000.;
    
    // grad and hess
    double gt;
    double ht;
    //double gb;
    //double hb;
    double lmij, thij, thijA, thijB, d2lijdthAdthB;
    
    
    
    double* Xt;
    double* Q;
    double* R;
    double* z;
    double* v;
    double* lambda;
    double* eta;
    double* gbeta;
    double* gamma;
    double* beta1;
    double* dbeta;
    
    
    //P=5;
    
    long conv=0;
    double stuck_b;
    double stuck_t;
    
    
    clear1(work, 2*(N0+P)*P + P*P + 2*(N0+P) + 2*N0 + 4*P);
     
    Xt      = work;//(double*)calloc(m*p, sizeof(double));
    Q       = work  +(N0+P)*P;// (double*)calloc(m*p, sizeof(double));
    R       = work+2*(N0+P)*P;//(double*)calloc(p*p, sizeof(double));
    z       = work+2*(N0+P)*P+P*P;//(double*)calloc(m, sizeof(double));
    v       = work+2*(N0+P)*P+P*P  +(N0+P);//(double*)calloc(m, sizeof(double));
    lambda  = work+2*(N0+P)*P+P*P+2*(N0+P);//(double*)calloc(n, sizeof(double));
    eta     = work+2*(N0+P)*P+P*P+2*(N0+P)  +N0;//(double*)calloc(n, sizeof(double));
    gbeta   = work+2*(N0+P)*P+P*P+2*(N0+P)+2*N0;//(double*)calloc(p, sizeof(double));
    gamma   = work+2*(N0+P)*P+P*P+2*(N0+P)+2*N0+P;//(double*)calloc(p, sizeof(double));
    
    beta1   = work+2*(N0+P)*P+P*P+2*(N0+P)+2*N0+P*2;
    dbeta   = work+2*(N0+P)*P+P*P+2*(N0+P)+2*N0+P*3;
    
    for(j=0; j<P; j++){cblas_dcopy(N0, X+j*N0, 1, Xt+j*(N0+P), 1);}
    for(j=0; j<P; j++){Xt[N0+j+j*(N0+P)] = 1.0/sqrt(sigma2);}
    
    if(verbose3>0){fprintf(stderr, "beta0: %lf\ninit beta: ", beta0); print(beta, P);}
    
    int jump_t=0;
    double xi=1e-20;
    for(itr=0; itr<1000; itr++){
        xi=1.0;
        clear1(gbeta, P);
        for(itr_beta=0; itr_beta<10000; itr_beta++){            
            
            cblas_dgemv(CblasColMajor, CblasNoTrans, N0, P, 1.0, X, N0, beta, 1, 0.0, eta, 1);
            for(i=0;  i<N0;   i++){lambda[i] = exp(eta[i]+log(ki[i]));  v[i] = sqrt(lambda[i]/(1.0+lambda[i]/th/ki[i]))*w[i];}// lambda v
            for(i=N0; i<N0+P; i++){v[i] = sqrt(th);}// v
            cblas_dcopy((N0+P)*P, Xt, 1, Q, 1);// Q
            for(i=0; i<(N0+P); i++){
                if(i<N0){ z[i] = (eta[i]+(y[i]-lambda[i])/lambda[i])*v[i]*w[i]; }else if(i==N0){z[i] = beta0/sqrt(sigma2)*v[i]; }else{ z[i] = 0.0; }// z
                for(j=0; j<P; j++){
                    Q[i+j*(N0+P)] *= v[i];// Q
                }
            }
            QR(Q, N0+P, P, R);
            //printM(R, 2,2);
            cblas_dgemv(CblasColMajor, CblasTrans, N0+P, P, 1.0, Q, N0+P, z, 1, 0.0, gamma, 1);
            //fprintf(stderr, "\ngamma: "); print(gamma, P);
            BackSolve(R, P, gamma, dbeta);
            cblas_daxpy(P, -1.0, beta, 1, dbeta, 1);// delta beta
            double alpha=2.0;
            stuck_b=1.0;
            //xi = cblas_ddot(P, gbeta, 1, dbeta, 1);
            for(j=0; j<20; j++){
                cblas_dcopy(P, beta, 1, beta1, 1);// beta -> beta1
                cblas_daxpy(P, alpha, dbeta, 1, beta1, 1);// beta1 <- beta1 + alpha * dbeta
                cblas_dgemv(CblasColMajor, CblasNoTrans, N0, P, 1.0, X, N0, beta1, 1, 0.0, eta, 1);
                for(i=0; i<N0; i++){lambda[i] = exp(eta[i]+log(ki[i]));}
                lkhd=getLkhdGlm(y, X, ki, w, N0, P, lambda, beta1, th);
                //if(verbose3>1){fprintf(stderr, "       lkhd0=%lf lkhd1=%lf %ld beta1: ", lkhd0, lkhd, j); print(beta1, P);}
                if(verbose3>10){fprintf(stderr, "Beta  [%ld]: lkhd=%lf beta=%lf\n", j, lkhd, beta1[0]);}
                
                if(lkhd0<lkhd+xi){lkhd0=lkhd; cblas_dcopy(P, beta1, 1, beta, 1); stuck_b=0.0; break;}else{alpha/=2.0;}
            }
            //cblas_dcopy(P, beta1, 1, beta, 1); 
            //lkhd0=lkhd;
            
            //fprintf(stderr, "\nbeta: "); print(beta, P);
            
            // new lambda
            cblas_dgemv(CblasColMajor, CblasNoTrans, N0, P, 1.0, X, N0, beta, 1, 0.0, eta, 1);
            for(i=0; i<N0; i++){lambda[i] = exp(eta[i]+log(ki[i]));}
            
            // grad for beta
            for(i=0; i<N0; i++){z[i] = (y[i]-lambda[i])/(1.0+lambda[i]/th/ki[i])*w[i];}
            cblas_dgemv(CblasColMajor, CblasTrans, N0, P, 1.0, X, N0, z, 1, 0.0, gbeta, 1);
            gbeta[0] -= th*(beta[0]-beta0)/sigma2;
            for(j=1; j<P; j++){
                gbeta[j] -= th*beta[j]/sigma2;
            }
            
            
            //fprintf(stderr, "\ngbeta: "); print(gbeta, P);
            
            if(asum(gbeta, P)<EPSBETA){
                stuck_b=0.0;
                if(verbose3>2){fprintf(stderr, "beta=%lf [beta END]\n", beta[0]);}
                break;
            }else{
                xi /= 2.0;
                if(verbose3>2){fprintf(stderr, "       lkhd0=%lf lkhd1=%lf i=%ld stuck=%lf grad=%lf beta1: ", lkhd0, lkhd, itr_beta, stuck_b, asum(gbeta, P)); print(beta1, 3);}
            }
            
            if(stuck_b>0.5){break;}
            
        }
        
        if(itr==0 && verbose3>0){fprintf(stderr, "Init likelihood:  %lf gradb=%lf\n", lkhd0, asum(gbeta, P));}
        
        
        xi=1.0;
        if(fixTheta==0){for(itr_theta=0; itr_theta<1000; itr_theta++){
            gt=ht=0.0;
            for(i=0; i<N0; i++){
#ifdef NOK
                thij = th;
#else
                thij = th*ki[i];
#endif
                lmij = lambda[i];
                gt += w[i] * thij/th*(          digamma(thij + y[i]) -  digamma(thij) + log(thij) - log( lmij + thij) + (lmij - y[i])/(lmij + thij) );
                ht += w[i] * thij*thij/th/th*( trigamma(thij + y[i]) - trigamma(thij) + 1.0/thij  - 1.0/(lmij + thij) - (lmij - y[i])/pow(lmij + thij, 2.0) );
            }
            
#ifdef LDGAMMA
            gt += (kappa/2.0)/th - omega/2.0;
            ht += -(kappa/2.0)/th/th;
#else
            gt += (kappa/2.0 - 1.0)/th - omega/2.0;
            ht += -(kappa/2.0-1.0)/th/th;
#endif
            gt += ((double)P)/2.0/th;
            gt += -pow(beta[0]-beta0,2.0)/sigma2/2.0;
            for(i=1; i<P; i++){ gt -= beta[i]*beta[i]/sigma2/2.0; }
            ht += -((double)P)/2.0/th/th;
            
            gt = gt*th;
            ht = ht*th*th+gt;
            
            // foo
            double alpha=1.0;
            stuck_t=1.0;
            for(i=0; i<20; i++){
            	th1 = exp( log(th) + alpha*sign1(gt)*fabs(gt/ht) );
                if(th1<1e-10){th1=1e-10;}else if(th1>1e10){th1=1e10;}
                lkhd=getLkhdGlm(y, X, ki, w, N0, P, lambda, beta, th1);
                if(verbose3>3){fprintf(stderr, "Theta [%ld]: lkhd=%lf th=%lf\n", i, lkhd, th1);}
                if(lkhd0<lkhd+xi){lkhd0=lkhd; th=th1; stuck_t=0.0; break;}else{alpha/=2.0;}
            }
            
            
            /*if(th1<1.0e-3){convt=1;
                if(jump_t==0){th=theta0;jump_t++; }else{th=1.0e-3; gt=gt; }
            }else if(th1>1000.){convt=1;
                if(jump_t==0){th=theta0;jump_t++; }else{th=1000.; gt=gt;}
            }else{th=th1; convt=0;}*/
            
            if(isnan(th)>0){fprintf(stderr, "th becomes nan in nbem step=%lf", sign1(gt)*fabs(gt/ht)); th=theta0; break;}
            
            
            if(fabs(gt/ht)<EPSTHETA){
                stuck_t=0.0;
                if(verbose3>2){fprintf(stderr,"   theta=%lf g=%lf h=%lf ss=%lf[theta END]\n", th, gt, ht, alpha);}
                break;
            }else{
                xi /= 2.0;
                if(verbose3>2)fprintf(stderr, "   itrth=%ld lkhd=%lf theta=%lf g=%lf h=%lf ss=%lf stuck=%lf\n",  itr_theta, lkhd0, th, gt/ht, ht, alpha, stuck_t);
            }
            
            if(stuck_t>0.5){break;}
            
        }}else{gt=0.0;}
        
        
        
        
        //Terminal condition
        if(asum(gbeta, P)<EPSBETA && fabs(gt/ht)<EPSTHETA && fabs(lkhd_all-lkhd0)<EPSLKHD){
            conv=1;
            if(verbose3>1){fprintf(stderr," itrall=%ld theta=%lf gt=%lf stt=%lf beta=%lf stb=%lf[glm END]\n", itr, th, gt, gt/ht, beta[0], asum(gbeta, P));}
            break;
        }else{
            conv=0;
            if(asum(gbeta, P)>1e-8){conv+=2;}
            if(fabs(gt/ht)>1e-8){conv+=4;}
            if(fabs(lkhd_all-lkhd0)>1e-8){conv+=8;}
            if(verbose3>1){fprintf(stderr," itrall=%ld theta=%lf gt=%lf stt=%lf beta=%lf stb=%lf\n", itr, th, gt, gt/ht, beta[0], asum(gbeta, P));}
            lkhd_all=lkhd0;
        }
        
        if(stuck_t*stuck_b>0.5){break;}
    }
    
    
    //################ posterior mean
    cblas_dgemv(CblasColMajor, CblasNoTrans, N0, P, 1.0, X, N0, beta, 1, 0.0, eta, 1);
    for(i=0;  i<N0;   i++){lambda[i] = exp(eta[i]+log(ki[i]));  v[i] = sqrt(lambda[i]/(1.0+lambda[i]/th/ki[i]))*w[i];}// lambda v
    for(i=N0; i<N0+P; i++){v[i] = sqrt(th);}// v
    cblas_dcopy((N0+P)*P, Xt, 1, Q, 1);// Q <- Xt
    for(i=0; i<(N0+P); i++){
        for(j=0; j<P; j++){
            Q[i+j*(N0+P)] *= v[i];// Q
        }
    }
    QR(Q, N0+P, P, R);
    for(i=0; i<N0; i++){// lambda is replaced by x^T H^-1 x/2
        lambda[i]=0.0;
        for(j=0; j<P; j++){
            lambda[i]+=pow(Q[i+j*(N0+P)]/v[i],2.0)/2.0;
        }
    }
    //################## pm end
    
    
    cblas_dgemv(CblasColMajor, CblasNoTrans, N0, P-1, 1.0, X+N0, N0, beta+1, 1, 0.0, eta, 1);
    for(i=0; i<N0; i++){dki[i] = exp(eta[i]+log(ki[i]) + 0.*lambda[i]);}// new!
    //printSep(ki, N0, '\n');
    double sdeta = sd(eta, N0);
    double meta = mean(eta, N0);
    double sdki = lsd(dki, N0);
    double mki = lmean(dki, N0);
    
    if(verbose3>0)fprintf(stderr, "mean eta=%lf sd eta=%lf sdki=%lf mki=%lf\n", meta, sdeta, sdki, mki);
    
    //fprintf(stderr, "Trimmed ki region: (%lf, %lf)\n", exp(-sdki*3.0+lmean(ki,N0)), exp(sdki*3.0+lmean(ki,N0)));
    //for(i=0; i<N0; i++){kii = exp(eta[i]+log(ki[i])); if(fabs(log(kii/ki[i]))>sdeta*3.0+meta){fprintf(stderr, "%ld %lf %lf\n", i, kii, ki[i]); ki[i]=1.0; w[i]=0.0;}else{ki[i]=kii;} }
    for(i=0; i<N0; i++){
        //if(ki[i]>11.0){w[i]=0.0;}
        if(fabs(log(dki[i]))>sdki*5+mki){
            if(verbose3>0){fprintf(stderr, "*%ld\t%lf\t%lf\t%lf\n", i, dki[i], dki[i]*exp(beta[0]), y[i]/dki[i]);} 
            //dki[i]=1.0; w[i]=0.0;// masking strange samples
        }else{
            if(verbose3>1){fprintf(stderr, " %ld\t%lf\t%lf\t%lf\n", i, dki[i], dki[i]*exp(beta[0]), y[i]/dki[i]);}
        }
        dki[i] /= ki[i];
    }
    //if(verbose3>0){fprintf(stderr, "Ki [%lf]: ", asum(ki, N0)); printSep(ki, N0, '\n');}
    
    if(verbose3>0) fprintf(stderr, "convstatus=%ld %lf stuck=%lf", conv, asum(w, N0), stuck_b+stuck_t);
    
    if(verbose3>0){ fprintf(stderr, "Final likelihood: %lf\nFinal theta: %lf\nFinal beta: ", lkhd0, th); printSep(beta, P, '\n'); fprintf(stderr, "Gradient th: %lf\nGradient be: %lf\n\n", gt, asum(gbeta, P));}
    *ptheta = th;
    clear1(work, (N0+P)*P + P*P + 2*(N0+P) + 2*N0 + 4*P);
    //sigma2=10000.0;
    return lkhd0;
}



