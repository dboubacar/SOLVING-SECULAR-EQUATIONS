#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "param.h"

void assert(){
  printf("Erreur a la ligne %d du fichier %s\n", __LINE__, __FILE__);
}
///dddddddddddddddddddddddddddbonnnnnnnnnnnnnnnnnnnnnnnnnnnnnnjour


PARAM * init_param(){

  PARAM *p=(PARAM*)malloc(sizeof(PARAM));
  p->zeta= (double *)malloc(sizeof(double)*p->DIM);
  p->delta= (double *)malloc(sizeof(double)*p->DIM);
  p->DIM=4;
  p->delta[0]=1; p->delta[1]=2;
  p->delta[2]=3 ;p->delta[3]=4;


  p->zeta[0] =0.0058 ;  p->zeta[1]=0.5773;  p->zeta[2]=0.5773;  p->zeta[3] =0.5773;
  p->rho  = pow(10, -2);


/*  p->zeta[0] =0.408044 ;  p->zeta[1]=0.373222;  p->zeta[2]=0.663589;  p->zeta[3] =0.50384;
  p->delta[0]=1.31585 ;p->delta[1]=1.44085;
  p->delta[2]=2.40323 ;p->delta[3]=3.71694;
  p->rho = -0.125;*/


  return p;
}

//la fonction g



double g(PARAM *p,double x,int k){
  double som=0.0;
  for (int j=1 ; j<p->DIM; j++){
     if ((j != k) && (j != k+1)){
       som += (pow(p->zeta[j],2)) / (p->delta[j]-x) ;
     }
  }
  return (p->rho +som);
}


double f(PARAM *p,double x){
  double som=0.0;
  for (int i=0 ; i<p->DIM; i++){
     som += (pow(p->zeta[i],2)) / (p->delta[i]-x) ;
  }
  return 1 + p->rho *som;
}



double f_prime(PARAM *p,double x){
  double som=0.0;
  for (int i=0 ; i<p->DIM; i++){
     som += (pow(p->zeta[i],2)) / pow((p->delta[i]-x),2);

  }
  return p->rho *som;
}

double f_seconde(PARAM *p,double x){
  double som=0.0;
  for (int i=0 ; i<p->DIM; i++){
     som += (2*pow(p->zeta[i],2)) / pow((p->delta[i]-x),3) ;
  }
  return p->rho *som;

}



double h(PARAM *p,double x,int k){
  double d1=p->delta[k]-x;
  double d2=p->delta[k+1]-x;
  if(d1==0 && d2==0) assert();
  return ((pow(p->zeta[k],2)) /d1)+((pow(p->zeta[k+1],2))/d2) ;
}

void calcul_a_b_c(PARAM *p,double *a,double *b,double *c,double y,int k){
  double DELTA_k = p->delta[k] - y;
  double DELTA_k1 = p->delta[k+1] - y;
  double val_f=f(p,y);
  double val_fp=f_prime(p,y);
  double val_fs=f_seconde(p,y);
  //printf("f ok %f\n",val_fp);
  *a=(DELTA_k+DELTA_k1)*val_f-DELTA_k*DELTA_k1*val_fp;
  *b= (DELTA_k *DELTA_k1 )*val_f;
  *c=val_f -(DELTA_k+DELTA_k1)*val_fp+DELTA_k*DELTA_k1*val_fs/2;
}
