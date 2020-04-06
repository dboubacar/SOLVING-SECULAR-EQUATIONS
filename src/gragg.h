#ifndef DEF_GRAGG
#define DEF_GRAGG
#include "param.h"


double initial_monotone(PARAM *p,int k);
double initial_non_monotone(PARAM *p);
double y_suivant(double y_prec,double a,double b,double c);
int stop_monotone(double *y,double *y1,double *y2,int i);
int stop_non_monotone(PARAM* p,double y,int k);
double zero_gragg(PARAM * p,int k);








#endif
