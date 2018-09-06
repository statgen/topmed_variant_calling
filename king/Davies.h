#ifndef __DAVIES_h__
#define __DAVIES_h__
//Davies R.B., Algorithm AS 155: The Distribution of a Linear Combination of chi-2 Random Variables,
//Journal of the Royal Statistical Society. Series C (Applied Statistics), 29(3), p. 323-333, (1980)
//void  Davies(double* lb1, double* nc1, int* n1, int *r1, double *sigma, double *c1, int *lim1, double *acc, double* trace, int* ifault, double *res);
#include <stdio.h>
double Davies(double c1, double* lb1, int r1, int *n1=NULL, double *nc1=NULL, double sigma=0, int lim1=10000, double acc=0.0001);

#endif

