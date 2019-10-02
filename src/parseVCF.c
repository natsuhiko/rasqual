#include "parseVCF.h"
#include "util.h"
char* format=NULL;
char* infostr=NULL;
//VCF_info vinfo;// = NULL;//{VT_OTHER, 0.0, 0.0};
char* buf=NULL;
char* cell=NULL;
int ret=0;
int i0=0;
int* formatID;

int uniqHap(int* dip, int N, int L){
        int l=0,j=0,k=0;
        int U=1;
        int flag;
        for(j=1; j<N*2; j++){
                for(k=0; k<U; k++){
                        flag=0;
                        for(l=0; l<L; l++){
                                if(dip[k+l*2*N]==dip[j+l*2*N]){flag++;}
                        }
                        if(flag==L){break;}
                }
                if(k==U){
                        for(l=0; l<L; l++){dip[U+l*2*N]=dip[j+l*2*N];}
                        U++;
                }
        }
        return U;
}

double getRsq(int* dip1, int* dip2, long N){
        double m1, m2, v1, v2, cv, dN;
        m1=m2=v1=v2=cv=0.0;
	dN = (double)N;
        int i;
        for(i=0; i<N; i++){
		if(dip2[i*2]+dip2[i*2+1]==1){
			m2++; v2++;
		}else if(dip2[i*2]+dip2[i*2+1]==2){
			m2+=2.0; v2+=4.0;
		}
                if((dip1[i*2]+dip1[i*2+1])==1){
                	m1++; v1++;
			if(dip2[i*2]+dip2[i*2+1]==1){
				cv++;
			}else if(dip2[i*2]+dip2[i*2+1]==2){
				cv+=2.0;
			}
		}else if((dip1[i*2]+dip1[i*2+1])==2){
			m1+=2.0; v1+=4.0;
			if(dip2[i*2]+dip2[i*2+1]==1){
				cv+=2.0;
			}else if(dip2[i*2]+dip2[i*2+1]==2){
				cv+=4.0;
			}
		}
        }
	m1 /= dN;
	m2 /= dN;
	v1 /= dN;
	v2 /= dN;
	cv /= dN;
	v1 -= m1*m1;
	v2 -= m2*m2;
	cv -= m1*m2;
        double r1= cv/sqrt(v1)/sqrt(v2);
//if(isnan(r1)>0){for(i=0;i<N;i++){fprintf(stderr, "%d,%d\n", dip2[i*2], dip2[i*2+1]);};}
//if(isnan(r1)>0){for(i=0;i<N;i++){fprintf(stderr, "%d,%d\n", dip1[i*2], dip1[i*2+1]);};}
	if(isnan(r1)>0){ fprintf(stderr, "%lf %lf %lf %lf %lf\n", m1, m2, v1 ,v2, cv); return 0.0;}else{return r1;}
}

double getRsqHap(int* dip1, int* dip2, long N){
	double a,b,c,d;
	a=b=c=d=0.0;
	int i;
	for(i=0; i<N*2; i++){
		if(dip1[i]==0&&dip2[i]==0){
			a++;
		}else 	if(dip1[i]==1&&dip2[i]==1){
			d++;
		}else 	if(dip1[i]==1&&dip2[i]==0){
			c++;
		}else 	if(dip1[i]==0&&dip2[i]==1){
			b++;
		}
	}
	double rsq=(a*d-b*c)/sqrt((a+b)*(a+c)*(b+d)*(c+d));
	return rsq;
}

double getAF(double* gen, long N){
	long i;
	double res=0.0;
	double denom=0.0;
	for(i=0; i<N; i++){if(gen[i]>=0.0){
		res += gen[i];
		denom += 1.0;
	}}
	res/=(denom*2);
	return res;
	//if(res>0.5){return 1.0-res;}else{return res;}
}
double getWAF(double* gen, double* w, long N){
	long i;
	double res=0.0;
	double denom=0.0;
	for(i=0; i<N; i++){if(gen[i]>=0.0){
		res += gen[i]*w[i];
		denom += w[i];
	}}
	res/=(denom*2);
	return res;
	//if(res>0.5){return 1.0-res;}else{return res;}
}
double getIAfromAP(double* ap, double* gl, long N){
	int i;
	for(i=0; i<N; i++){
		gl[i*3]   = (1.0-ap[i*2])*(1.0-ap[i*2+1]);
		gl[i*3+1] = (1.0-ap[i*2])*ap[i*2+1] + ap[i*2]*(1.0-ap[i*2+1]);
		gl[i*3+2] = ap[i*2]*ap[i*2+1];
	}
	return getIA(gl, N);
}
double getWIAfromAP(double* ap, double* gl, double* w, long N){
	int i;
	for(i=0; i<N; i++){
		gl[i*3]   = (1.0-ap[i*2])*(1.0-ap[i*2+1]);
		gl[i*3+1] = (1.0-ap[i*2])*ap[i*2+1] + ap[i*2]*(1.0-ap[i*2+1]);
		gl[i*3+2] = ap[i*2]*ap[i*2+1];
	}
	return getWIA(gl, w, N);
}

double getIA(double* gl, long N){
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
double getWIA(double* gl, double* w, long N){
	long i;
	double IA=0.0;
	double tot=0.0;
	double th=0.0;
	for(i=0; i<N; i++){if(gl[i*3]>=0.0){
		th += w[i]*(gl[i*3+2]+gl[i*3+1]/2.0);
		tot+= w[i];
	}}
	th /= tot;
	if(th<=0.0 || th>=1.0){
		return 1.0;
	}
	for(i=0; i<N; i++){if(gl[i*3]>=0.0){
		IA += w[i]*(4.0*gl[i*3+2]+gl[i*3+1] - pow(2.0*gl[i*3+2]+gl[i*3+1],2.0))/2.0/th/(1.0-th);
	}}
	IA /= tot;
	return 1.0-IA;
}

double getWHWEfromAP(double* ap, int* dip, double* w, long N){
	int i;
	for(i=0; i<N; i++){
		dip[i*2]   = ap[i*2]>0.5   ? 1 : 0;
		dip[i*2+1] = ap[i*2+1]>0.5 ? 1 : 0;
	}
	return getWHWE(dip, w, N);
}
double getHWEfromAP(double* ap, int* dip, long N){
	int i;
	for(i=0; i<N; i++){
		dip[i*2]   = ap[i*2]>0.5   ? 1 : 0;
		dip[i*2+1] = ap[i*2+1]>0.5 ? 1 : 0;
	}
	return getHWE(dip, N);
}

double getHWE(int* dip, long N){
	long i;
	double g0, g1, g2;
	g0=g1=g2=0.0;
	double tot=0.0;
	for(i=0; i<N; i++){
		if((dip[i*2]+dip[i*2+1])==0){
			g0 += 1.0;tot += 1.0;
		}else if((dip[i*2]+dip[i*2+1])==1){
			g1 += 1.0;tot += 1.0;
		}else if((dip[i*2]+dip[i*2+1])==2){
			g2 += 1.0;tot += 1.0;
		}
	}
	double p = (g1+g2*2)/(2*tot);
	double q = 1.0-p;
	double f = pow(g2-p*p*tot,2.0)/p/p/tot + pow(g1-2.0*p*q*tot,2.0)/2.0/p/q/tot + pow(g0-q*q*tot,2.0)/q/q/tot;
	return f;
}
double getWHWE(int* dip, double* w, long N){
	long i;
	double g0, g1, g2;
	g0=g1=g2=0.0;
	double tot=0.0;
	for(i=0; i<N; i++){
		if((dip[i*2]+dip[i*2+1])==0){
			g0 += w[i];tot += w[i];
		}else if((dip[i*2]+dip[i*2+1])==1){
			g1 += w[i];tot += w[i];
		}else if((dip[i*2]+dip[i*2+1])==2){
			g2 += w[i];tot += w[i];
		}
	}
	double p = (g1+g2*2)/(2*tot);
	double q = 1.0-p;
	double f = pow(g2-p*p*tot,2.0)/p/p/tot + pow(g1-2.0*p*q*tot,2.0)/2.0/p/q/tot + pow(g0-q*q*tot,2.0)/q/q/tot;
	return f;
}
double getWCR(double* gen, double* w, long N){
	long i;
	double cr=0.0;
	double tot=0.0;
	for(i=0; i<N; i++){
		if(gen[i]>=0.0){
			cr += w[i];
		}
		tot += w[i];
	}
	return cr/tot;
}
double getCR(double* gen, long N){
	long i;
	double cr=0.0;
	double tot=0.0;
	for(i=0; i<N; i++){
		if(gen[i]>=0.0){
			cr += 1.0;
		}
		tot += 1.0;
	}
	return cr/tot;
}

void gl2ap(double* g, double* h, int* dip){
	int i;
	if(g[0]<0.0){h[0]=h[1]=0.5;return;}
	//double gtot=g[0]+g[1]+g[2];
	//g[0]/=gtot;
	//g[1]/=gtot;
	//g[2]/=gtot;
	double phi, phi0;
	double p=(g[2]+g[1]/2.0)*1.001;
	double q=(g[2]+g[1]/2.0)*0.999;
	phi0 = phi = p*(1.-q)/(p*(1.-q)+q*(1.-p));
	for(i=0; i<10000; i++){
		p = g[2]+g[1]*phi;
		q = g[2]+g[1]*(1.-phi);
		phi=p*(1.-q)/(p*(1.-q)+q*(1.-p));
		if(fabs(phi-phi0)<1e-10){break;}else{phi0=phi;}
	}
	if(p>q){double tmp=p; p=q; q=tmp;}
	if(dip[0]==0){
		h[0]=p;
		h[1]=q;
	}else{
		h[0]=q;
		h[1]=p;
	}
	if(dip[0]+dip[1]!=1){
		h[0] = h[1] = (h[0]+h[1])/2.0;
	}
}

int parseCell(char* cell, int* dip, double* gl, double* ap, long* ase, double* dose, int* formatID){
	int i, ii0=0, k=0;
	for(i=0; i<strlen(cell)+1; i++){
		if(cell[i]==':' || cell[i]=='\0'){
            //if(cell[i]==':'){fprintf(stderr, "%s ", cell+ii0);}else{fprintf(stderr, "%s\n", cell+ii0);}
			if(i-ii0>0){
				if(formatID[k]==FORMAT_GT){
					if(cell[ii0]=='.'){
						dip[0]=dip[1]=9;
					}else{
						if(cell[ii0+1]=='|'){
							sscanf(cell+ii0, "%d|%d", dip, dip+1);
						}else if(cell[ii0+1]=='/'){
							sscanf(cell+ii0, "%d/%d", dip, dip+1);
						}
					}
				}else if(formatID[k]==FORMAT_GL || formatID[k]==FORMAT_PP){
					sscanf(cell+ii0, "%lf,%lf,%lf", gl, gl+1, gl+2);
					if(formatID[k]==FORMAT_PP){
						if(gl[0]<0.0 || gl[1]<0.0 || gl[2]<0.0 || gl[0]>1.0 || gl[1]>1.0 || gl[2]>1.0){fprintf(stderr, "Posterior Probability (%lf, %lf, %lf) is not in the appropriate range [0,1]...\n", gl[0], gl[1], gl[2]); int j; for(j=0; j<3; j++){if(gl[j]>1.0){gl[j]=1.0;}}  }
						if(gl[0]<1.0e-5){gl[0]=1.0e-5;}
						if(gl[1]<1.0e-5){gl[1]=1.0e-5;}
						if(gl[2]<1.0e-5){gl[2]=1.0e-5;}
					}else if(formatID[k]==FORMAT_GL){
						if(LOG10>0){
							if(gl[0]>0.0 || gl[1]>0.0 || gl[2]>0.0){fprintf(stderr, "Genotype likelihood (%lf, %lf, %lf) is not in the appropriate range [-Inf,0]...\n", gl[0], gl[1], gl[2]);}
							gl[0] = pow(10.0,gl[0]);
							gl[1] = pow(10.0,gl[1]);
							gl[2] = pow(10.0,gl[2]);
						}
					}
					double tot=gl[0]+gl[1]+gl[2];
					gl[0]/=tot; gl[1]/=tot; gl[2]/=tot;
					//if(tot>2.0){gl[0]=gl[1]=gl[2]=-1.0; dip[0]=dip[1]=-1;}
				}else if(formatID[k]==FORMAT_AP){
					sscanf(cell+ii0, "%lf,%lf", ap, ap+1);
					//if(ap[0]>=0.0){ap[0]=log10(1.0-1e-5);}
					//if(ap[1]>=0.0){ap[1]=log10(1.0-1e-5);}
					if(LOG10>0){
						if(ap[0]<=-5.0){ ap[0]=0.0; }else{ ap[0] = pow(10.0,ap[0]); }
						if(ap[1]<=-5.0){ ap[1]=0.0; }else{ ap[1] = pow(10.0,ap[1]); }
					}
				}else if(formatID[k]==FORMAT_AS){
					sscanf(cell+ii0, "%ld,%ld", ase, ase+1);
				}else if(formatID[k]==FORMAT_BF){
					sscanf(cell+ii0, "%lf,%lf", ap, ap+1);
				}else if(formatID[k]==FORMAT_DS){
					sscanf(cell+ii0, "%lf", dose);
				}else if(formatID[k]==FORMAT_AP_ORIG){
					//sscanf(cell+ii0, "%lf,%lf", ap, ap+1);
				}else if(formatID[k]==FORMAT_OTHER){
			
				}
			}else{dip[0]=dip[1]=9;}
			k++;
			ii0=i+1;
		}
	}
	return 0;
}

void dose2ap(int* dip, double* ap, double ds){
        if((dip[0]==0&&dip[1]==0) || (dip[0]==1&&dip[1]==1)){
               ap[0]=ap[1]=ds/2.0;
        }else if(dip[0]==0&&dip[1]==1){
               if(ds>1.0){ap[0]=ds-1.0; ap[1]=1.0;}else{ap[0]=0.0; ap[1]=ds;}
        }else if(dip[0]==1&&dip[1]==0){
               if(ds>1.0){ap[1]=ds-1.0; ap[0]=1.0;}else{ap[1]=0.0; ap[0]=ds;}
        }else{
               ap[0]=ap[1]=0.5;
	}
}
void gt2ap(int* dip, double* ap){
	ap[0] = dip[0]==0 ? 0.0 : 1.0;
	ap[1] = dip[1]==0 ? 0.0 : 1.0;
}

int main2(int argc, char** argv){
	int i;
	if(argc==1){fprintf(stderr, "Usage:\n  ./parseVCF <gzipped file> [-N] [-header]\n");return -1;}
	
	for(i=0; i<argc; i++){if(strcmp(argv[i],"-header")==0){return parseHeader(argc, argv);}}
	long N=0;
	for(i=0; i<argc-1; i++){if(strcmp(argv[i],"-N")==0){N=atoi(argv[i+1]);}}
	fprintf(stderr, "Sample size=%ld\n", N);
	if(N==0){fprintf(stderr, "Sample size is inappropriate\n"); return -1;}
	
	//formatID = (int*)calloc(100, sizeof(int));
	
	//buf=(char*)calloc(BUFS,sizeof(char));
	//cell=(char*)calloc(CELLS,sizeof(char));
	
	if(init()==0){fprintf(stderr, "momory allocation failuer\n"); return 0;}
	
	int* dip;
	double* gen;
	double* gl;
	double* ap;
	//double* ap2;
	dip=(int*)calloc(N*2, sizeof(int));
	gen=(double*)calloc(N, sizeof(double));
	gl =(double*)calloc(N*3, sizeof(double));
	ap =(double*)calloc(N*2, sizeof(double));
	//ap2 =(double*)calloc(N*2, sizeof(double));
	
	//while(parseLine(dip, gen, gl, ap, N)>0){
	//	fprintf(stderr, "\n\n");
	//}
	return 1;
}

void Log10(){
	if(LOG10==1){LOG10=0;}else{LOG10=1;}
}

int init(){
	LOG10=1;
	//VCF_info vinfo = {VT_OTHER, 0.0, 0.0};
	format = (char*)calloc(1000,sizeof(char));
	infostr = (char*)calloc(100,sizeof(char));
	buf=(char*)calloc(BUFS,sizeof(char));
	cell=(char*)calloc(CELLS,sizeof(char));
	formatID = (int*)calloc(100, sizeof(int));
	if(infostr==NULL || format==NULL || buf==NULL || cell==NULL){return 0;}else{return 1;}
}

int Fread_parseVCF(char* buf, FILE *fp){
	i0=0;
	return fread(buf, sizeof(char), BUFS, fp);
}

int existFlag(int* formatID, int nfield, int flag){ 
	int j;
	for(j=0; j<nfield; j++){if(formatID[j]==flag){return 1;}}
	return 0;
}

// gen will be < 0 and ap=c(-1,-1) gl=(-1,-1,-1) if log10(gl) was (0,0,0)

// GZFILE


// STDIN
int parseLine(char* chr, long* pos, char* rs, char** al, VCF_info* vinfo, int* dip, double* dose, double* gl, double* ap, long* ase, long N, int genType){
	return parseLine0(chr, pos, rs, al, vinfo, dip, dose, gl, ap, ase, N, genType, stdin);
}


int parseLine0(char* chr, long* pos, char* rs, char** al, VCF_info* vinfo, int* dip, double* dose, double* gl, double* ap, long* ase, long N, int genType, FILE *fp){
	int i=0,j=0,l=0;
	int nfield=0;
	int warningFlag=0;
	int ncol=0;
	
	do{
		for(i=i0; i<ret; i++){
			if(buf[i]=='\t' || buf[i]=='\n'){
				cell[l]='\0';
				if(ncol==0){
					strcpy(chr, cell);
				}else if(ncol==1){
					sscanf(cell, "%ld", pos);
				}else if(ncol==2){
					strcpy(rs, cell);
				}else if(ncol==3){
					al[0]=(char*)calloc(strlen(cell)+1, sizeof(char)); strcpy(al[0], cell);
				}else if(ncol==4){
					al[1]=(char*)calloc(strlen(cell)+1, sizeof(char)); strcpy(al[1], cell);
				}else if(ncol==7){
					parseInfo(cell, vinfo);
				}else if(ncol==8){
					nfield = parseFormat(cell, formatID);
				}else if(ncol<5){
				
				}else if(ncol>8){
					parseCell(cell, dip+(ncol-9)*2, gl+(ncol-9)*3, ap+(ncol-9)*2, ase+(ncol-9)*2, dose+(ncol-9), formatID);
					//estimating allelic prob
					if(existFlag(formatID, nfield, FORMAT_AP)==0){// no AP
						if(existFlag(formatID, nfield, FORMAT_GT)>0){// GT
							if(existFlag(formatID, nfield, FORMAT_GL)>0 || existFlag(formatID, nfield, FORMAT_PP)>0){// gen likelihood
								gl2ap(gl+(ncol-9)*3, ap+(ncol-9)*2, dip+(ncol-9)*2);
							}else if(existFlag(formatID, nfield, FORMAT_DS)>0){// dose
								dose2ap(dip+(ncol-9)*2, ap+(ncol-9)*2, dose[ncol-9]);
							}else{
								gt2ap(dip+(ncol-9)*2, ap+(ncol-9)*2);
							}
//gt2ap(dip+(ncol-9)*2, ap+(ncol-9)*2);
							dose[ncol-9] = (double)(dip[(ncol-9)*2]+dip[(ncol-9)*2+1]);
							if(dip[(ncol-9)*2]==9||dip[(ncol-9)*2+1]==9){dose[ncol-9] = 9;}
						}else{// allele freq
							dose[ncol-9] = (vinfo->AF)*2.0;
							ap[(ncol-9)*2] = ap[(ncol-9)*2+1] = vinfo->AF;
						}
					}else{
						dose[ncol-9] = ap[(ncol-9)*2]+ap[(ncol-9)*2+1];
					}
				}
				l=0;
				ncol++;
				warningFlag=0;
				if(buf[i]=='\n'){i0=i+1; return 1;}
			}else{
				if(l<CELLS-1){
					cell[l++]=buf[i];
				}else{
					if(warningFlag==0){
						fprintf(stderr, "cell size exceed %d %s\n", CELLS, cell); warningFlag=1; 
					}; 
					if(ncol>8){return -1;} 
				}
			}
		}
	}while((ret=Fread_parseVCF(buf, fp))>0);
	return 0;
}

int startWith(const char *pre, const char *str)
{
    return strncmp(pre, str, strlen(pre));
}


int parseInfo(char* str, VCF_info* vinfo){
	int i;
	int nfield=1;
	for(i=0; i<strlen(str); i++){if(str[i]==';'){nfield++;}}
	int ns=0;
	vinfo->RSQ = -1.0;
	for(i=0; i<nfield; i++){
		sscanf(str+ns, "%[^;];", infostr);
		if(strcmp(infostr,"VT=SNP")==0){
			vinfo->VT=VT_SNP;
		}else if(startWith("RSQ=",infostr)==0){
			sscanf(infostr, "RSQ=%lf", &(vinfo->RSQ));
        }else if(startWith("AF=",infostr)==0){
            sscanf(infostr, "AF=%lf", &(vinfo->AF));
        }else if(startWith("DR2=",infostr)==0){
            sscanf(infostr, "DR2=%lf", &(vinfo->RSQ));
		}else if(startWith("IMP2=",infostr)==0){
			sscanf(infostr, "IMP2=%lf,%lf,%lf", &(vinfo->AF), &(vinfo->RSQ), &(vinfo->CER));
		}
		ns += strlen(infostr)+1;
	}
	return 0;
}

int parseFormat(char* str, int* formatID){
	int i;
	int nfield=1;
	for(i=0; i<strlen(str); i++){if(str[i]==':'){nfield++;}}
	
	int ns=0;
	for(i=0; i<nfield; i++){
		if(i<nfield-1){
			sscanf(str+ns, "%[^:]:", format);
		}else{
			sscanf(str+ns, "%s", format);
		}
		ns += strlen(format)+1;
//fprintf(stderr, "%d %s ", ns, format);		
		if(strcmp(format,"GT")==0){
			formatID[i]=FORMAT_GT;
		}else if(strcmp(format,"GL")==0){
			formatID[i]=FORMAT_GL;
		}else if(strcmp(format,"AP")==0){
			formatID[i]=FORMAT_AP;
		}else if(strcmp(format,"GP")==0){
			formatID[i]=FORMAT_PP;
		}else if(strcmp(format,"PP")==0){
			formatID[i]=FORMAT_PP;
		}else if(strcmp(format,"AS")==0){
			formatID[i]=FORMAT_AS;
		}else if(strcmp(format,"RD")==0){
			formatID[i]=FORMAT_RD;
		}else if(strcmp(format,"BF")==0){
			formatID[i]=FORMAT_BF;
		}else if(strcmp(format,"DS")==0){
			formatID[i]=FORMAT_DS;
		}else if(strcmp(format,"AP2")==0){
			formatID[i]=FORMAT_AP_ORIG;
		}else{
			formatID[i]=FORMAT_OTHER;
		}
	}
	return nfield;
}

int parseHap(char* hap);

int parseGL(char* gl);



int parseHeader(int argc, char** argv){ 
	int i;
	gzFile f1=NULL;
	if((f1 = gzopen(argv[1], "rb6f"))==NULL){fprintf(stderr, "gzread failed. no file!\n"); return -1;};

	if(f1==NULL){fprintf(stderr, "gzread failed. no file!\n"); return -1;}
	
	long Nrow=0, Nhead=0;
	long Ncol=0; 
	long l=0;
	int ret;
	char* buf;
	int flagHeader=0;
	int flagBreak=0;
	long skip=0;
	buf=(char*)calloc(1000,sizeof(char));
	while((ret = gzread(f1, buf, sizeof(buf))) != 0 && ret!=-1){
		for(i=0;i<ret/sizeof(char);i++){
			if(buf[i]=='\n'){
				//fprintf(stderr, "nrow=%ld, ncol=%ld\n", Ncol, Nrow);
				Nrow++;
				l=0;
				if(flagHeader==1){Ncol++; flagBreak=1; break;}else{Ncol=0; flagHeader=0;}
			}else if(buf[i]==' ' || buf[i]=='\t'){
				if(flagHeader<2){
					Ncol++;
				}
			}else{
				if(l==0){
					if(buf[i]=='#'){
						flagHeader=1;
						Nhead++;
					}
				}else if(l==1){
					if(buf[i]=='#'){
						flagHeader=2;
					}
				}
				l++;
			}
			if(flagHeader==2){skip++;}
		}
		if(flagBreak==1){break;}
	}
	//fprintf(stderr, "Num of Variants = %ld\n", Nrow-Nhead);
	gzseek(f1, (skip+Nhead*2-2)*sizeof(char), SEEK_SET);
	//for(i=0; i<skip+Nhead; i++){gzread(f1, buf, 1);};
	char** header;
	header = (char**)calloc(Ncol, sizeof(char*));
	for(i=0; i<Ncol; i++){header[i]=(char*)calloc(100, sizeof(char));}
	flagBreak=0;flagHeader=0;Nrow=Ncol=0;l=0;
	while((ret = gzread(f1, buf, sizeof(buf))) != 0 && ret!=-1){
		for(i=0;i<ret/sizeof(char);i++){
			if(buf[i]=='\n'){
				header[Ncol][l]='\0';
				Ncol++;
				flagBreak=1;break;
			}else if(buf[i]==' ' || buf[i]=='\t'){
				header[Ncol++][l]='\0';
				l=0;
			}else{
				header[Ncol][l++]=buf[i];
			}
		}
		if(flagBreak==1){break;}
	}
	for(i=0; i<Ncol; i++){fprintf(stderr, "%s ", header[i]);}
	fprintf(stderr, "\n");
	fprintf(stdout, "%ld", Ncol);
	return Ncol;
}


