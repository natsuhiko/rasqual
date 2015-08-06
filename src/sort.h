#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <zlib.h>
//#include <cblas.h>

double getMinBH(double* y, int n);
double getBH(double* y, double* q, int n);
void getRank(double* y, int n, double* ran);
void getOrder(double* y, int n, int* ord);

void randomise2(double* y, int n, int k, int* batch);
void randomise3(double* y, double* Y, int n, int k);
void randomise4(double* y, double* Y, double* ki, double* X, int n, int k, int p);
void getRandomOrder2(int n, int* ord, int nb, int* batch);
int sortAndUniqLong(long* x, int n);
int compare_long(const void *a, const void *b);
int compare_doubles (const void *a, const void *b);
int compare_ORDER( const void *c1, const void *c2 );
void getRandomOrder(int n, int* array);
void randomise(double* y, int n);
