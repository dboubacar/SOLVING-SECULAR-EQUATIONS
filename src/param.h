#ifndef DEF_PARAM
#define DEF_PARAM
#include "param.h"
#include <sys/time.h>


typedef struct{
  int DIM;
  double rho;
  double *zeta;
  double *delta;
}PARAM;

PARAM * init_param();
double g(PARAM *p,double x,int k);
double f(PARAM *p,double x);
double f_prime(PARAM *p,double x);
double f_seconde(PARAM *p,double x);
double h(PARAM *p,double x,int k);
void calcul_a_b_c(PARAM *p,double *a,double *b,double *c,double y,int k);
void assert();







#endif
