#include "sort.h"
#include "util.h"
#include <time.h>

typedef struct{
    double val;
    int ord;
}ORDER;

typedef struct{
    long val;
    int ord;
}ORDER_long;

double getBH(double* y, double* q, int n){
	int i;
	double en=0.0;
	for(i=0; i<n; i++){
		if(y[i]>0.0){en++;}
	}
	double* ran;
	ran = (double*)calloc(n, sizeof(double));
	getRank(y, n, ran);
	for(i=0; i<n; i++){
		if(y[i]>0.0){
			q[i] = (log(pchisq(y[i], 1.0)) - log(ran[i]) + log(en)) / log(10.0);
		}else{q[i]=0.0;}
	}
	free(ran);
	int* ord;
	ord = (int*)calloc(n, sizeof(int));
	getOrder(y, n, ord);
	double minBH=0.0;
	if( y[ord[n-1]]>0.0 ){minBH=min(0.0, q[ord[n-1]]);}
	for(i=n-1; i>0; i--){
		if(y[ord[i-1]]>0.0){if( q[ord[i-1]]>minBH ){ q[ord[i-1]]=minBH; }else{ minBH = q[ord[i-1]]; } }
	}
	return minBH;
}

double getMinBH(double* y, int n){
	int i;
	double en=0.0;
	for(i=0; i<n; i++){
		if(y[i]>0.0){en++;}
	}
	double* ran;
	ran = (double*)calloc(n, sizeof(double));
	getRank(y, n, ran);
	double minBH = 0.0;
	double minBH1;
	for(i=0; i<n; i++){
		if(y[i]>0.0){
			minBH1 = log(pchisq(y[i], 1.0)) - log(ran[i]) + log(en);
            //printf("%lf %lf log(p)=%lf\n", y[i], minBH1, log(pchisq(y[i], 1.0)));
			if(minBH1 < minBH){
				minBH = minBH1;
			}
		}
	}
	free(ran);
	return minBH/log(10.0);
}

void getOrder(double* y, int n, int* ord){
    ORDER* array;
    int i;
    array=(ORDER*)calloc(n, sizeof(ORDER));
    for(i=0; i<n; i++){
        array[i].val=-y[i];
        array[i].ord=i;
    }
    qsort(array, n, sizeof(ORDER), compare_ORDER);
    for(i=0; i<n; i++){
		ord[i] = array[i].ord;
    }
    free(array);
}


void getRank(double* y, int n, double* ran){
    ORDER* array;
    int i;
    array=(ORDER*)calloc(n, sizeof(ORDER));
    for(i=0; i<n; i++){
        array[i].val=-y[i];
        array[i].ord=i;
    }
    qsort(array, n, sizeof(ORDER), compare_ORDER);
    for(i=0; i<n; i++){
        array[i].val=(double)array[i].ord;
        array[i].ord=i;
    }
    qsort(array, n, sizeof(ORDER), compare_ORDER);
    for(i=0; i<n; i++){ran[i]=(double)array[i].ord+1.0;}
    free(array);
}

int sortAndUniqLong(long* x, int n){
	qsort(x, n, sizeof(long), compare_long);
	int i, l=0;
	for(i=1; i<n; i++){
		if(x[l]==x[i]){
			x[l]=x[i];
		}else{
			x[++l]=x[i];
		}
	}
	return l+1;
}

void randomise(double* y, int n){
	double* y1;
	int* ord;
	int i;
	ord=(int*)calloc(n, sizeof(int));
	for(i=0; i<n; i++){ord[i]=i;}
	y1=(double*)calloc(n, sizeof(double));
	getRandomOrder(n, ord);
	for(i=0; i<n; i++){
		y1[i]=y[ord[i]];
	}
	for(i=0; i<n; i++){
		y[i]=y1[i];
	}
	free(y1);
	free(ord);
}
void randomise2(double* y, int n, int k, int* batch){
	// y: k x n matrix (col major)
	// k: number of exons
	// n: samples size
	// batch: n dim vector of [0,nb) values
	double* y1;
	int* ord;
	int i;
	int nb=0;
	for(i=0; i<n; i++){if(batch[i]>nb){nb=batch[i];}}; nb++;
	ord=(int*)calloc(n, sizeof(int));
	for(i=0; i<n; i++){ord[i]=i;}
	y1=(double*)calloc(n*k, sizeof(double));
	getRandomOrder2(n, ord, nb, batch);
	for(i=0; i<n; i++){
		//cblas_dcopy(k, y+ord[i]*k, 1, y1+i*k, 1);
	}
	//cblas_dcopy(k*n, y1, 1, y, 1);
	free(y1);
	free(ord);
}
void randomize3(double*y, double* Y, int n, int k){
	// y : n x 1
	// Y : n x k (row major?)
	double* y1;
	double* Y1;
	int i;
	int* ord;
	ord=(int*)calloc(n, sizeof(int));
	for(i=0; i<n; i++){ord[i]=i;}
	Y1=(double*)calloc(n*k, sizeof(double));
	y1=(double*)calloc(n, sizeof(double));
	getRandomOrder(n, ord);
	for(i=0; i<n; i++){
		y1[i]=y[ord[i]];
		//cblas_dcopy(k, Y+ord[i]*k, 1, Y1+i*k, 1);
	}
	//cblas_dcopy(n,   y1, 1, y, 1);
	//cblas_dcopy(k*n, Y1, 1, Y, 1);
	free(ord);
	free(y1);
	free(Y1);
}


void randomise4(double*y, double* Y, double* ki, double* X, int n, int k, int p){
	// y : n x 1
	// Y : 2 x n x k 
	// k : n x 1
	// X : n x p
	double* y1;
	double* Y1;
	double* X1;
	double* ki1;
	int i;
	int* ord;
	ord=(int*)calloc(n, sizeof(int));
	for(i=0; i<n; i++){ord[i]=i;}
	Y1=(double*)calloc(2*n*k, sizeof(double));
	X1=(double*)calloc(n*p, sizeof(double));
	y1=(double*)calloc(n, sizeof(double));
	ki1=(double*)calloc(n, sizeof(double));
	getRandomOrder(n, ord);
	for(i=0; i<n; i++){
		y1[i]=y[ord[i]];
		ki1[i]=ki[ord[i]];
		//cblas_dcopy(k, Y+ord[i]*2,   2*n, Y1+i*2,   2*n);
		//cblas_dcopy(k, Y+ord[i]*2+1, 2*n, Y1+i*2+1, 2*n);
		//cblas_dcopy(p, X+ord[i]    ,   n, X1+i    ,   n);
	}
	//cblas_dcopy(n,     y1, 1, y, 1);
	//cblas_dcopy(n,    ki1, 1, ki,1);
    //fprintf(stderr, "randomized!\n!");
	//cblas_dcopy(2*n*k, Y1, 1, Y, 1);
    //	cblas_dcopy(n*p,   X1, 1, X, 1);
	free(ord);
	free(y1);
	free(Y1);
	free(X1);
	free(ki1);
}




void getRandomOrder2(int n, int* ord, int nb, int* batch){
	int i, j, l;
	int* ord2;
	ord2=(int*)calloc(n, sizeof(int));
	for(i=0; i<nb; i++){
		l=0;
		for(j=0; j<n; j++){
			if(batch[j]==i){
				ord2[l++]=ord[j];
			}
		}
		getRandomOrder(l, ord2);
		l=0;
		for(j=0; j<n; j++){
			if(batch[j]==i){
				ord[j]=ord2[l++];
			}
		}
	}
	free(ord2);
}

void roc(double* x, int* y, int n){
	ORDER* array;
	int i;
	array=(ORDER*)calloc(n, sizeof(ORDER));
	for(i=0; i<n; i++){
		array[i].val=x[i];
		array[i].ord=y[i];
	}
	qsort(array, n, sizeof(ORDER), compare_ORDER);
	for(i=0; i<n; i++){
		x[i]=array[i].val;
		y[i]=array[i].ord;
	}
}

void getRandomOrder(int n, int* ord){
	ORDER* array;
	int i;
	//srand((unsigned)(time(NULL)+getpid())); //srand( (unsigned int)time( NULL ) );
	array=(ORDER*)calloc(n, sizeof(ORDER));
	for(i=0; i<n; i++){
		array[i].val=rand();
		array[i].ord=ord[i];
		//fprintf(stderr, "%d %lf\n", array[i].ord, array[i].val);
	}
	qsort(array, n, sizeof(ORDER), compare_ORDER);
	for(i=0; i<n; i++){ord[i]=array[i].ord;}
	free(array);
}

int compare_ORDER( const void *c1, const void *c2 )
{
    ORDER test1 = *(ORDER *)c1;
    ORDER test2 = *(ORDER *)c2;
    
    double tmp1 = test1.val;   /* b を基準とする */
    double tmp2 = test2.val;
    
    if(tmp1-tmp2>0){return 1;}else if(tmp1-tmp2<0){return -1;}else{ return 0;}
}

int compare_ORDER_long( const void *c1, const void *c2 )
{
    ORDER_long test1 = *(ORDER_long *)c1;
    ORDER_long test2 = *(ORDER_long *)c2;
    
    long tmp1 = test1.val;   /* b を基準とする */
    long tmp2 = test2.val;
    
    if(tmp1-tmp2>0){return 1;}else if(tmp1-tmp2<0){return -1;}else{ return 0;}
}

int compare_long( const void *c1, const void *c2 )
{       
    long tmp1 = *(long*)c1;
    long tmp2 = *(long*)c2;
    
    if(tmp1-tmp2>0){return 1;}else if(tmp1-tmp2<0){return -1;}else{ return 0;}
}

int compare_doubles (const void *a, const void *b){
    double temp = *(double *)a - *(double *)b;
    if (temp > 0)
        return 1;
    else if (temp < 0)
        return -1;
    else
        return 0;
}
/*
 int main_sort2(void){
 int i;
 double x[] = {1.,2.,3.,4.,5.,6.};
 int batch[] = {0,1,0,1,0,2};
 randomise2(x, 6, 1, batch);
 for(i=0; i<6; i++){printf("%lf ",x[i]);}
 printf("\n");
 }
 
 int mainsort(void){
 long x[]={1,10,2,4,3,5,6,8,7,9};
 int i,l;
 l=sortAndUniqLong(x,10);
 for(i=0; i<l; i++){printf("%ld ", x[i]);}
 }
 */




