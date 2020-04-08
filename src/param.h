#ifndef DEF_PARAM
#define DEF_PARAM


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
void calcul_a(PARAM *p,double *a,double y,int k);
void calcul_b(PARAM *p,double *b,double y,int k);
void calcul_c(PARAM *p,double *c,double y,int k);
void free_param(PARAM *p);

/*voir array_free(struct Array *array_ptr)
{
    free(array_ptr->data);
    array_ptr->data = NULL;    // précaution contre double array_free
    array_ptr->size = 0;
}*/



void assert();







#endif
