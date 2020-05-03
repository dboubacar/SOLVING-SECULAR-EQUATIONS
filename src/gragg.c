/*
Résolution de l’équation séculaire
  Par: Boubacar Diallo  et
       Jeffrey  Kebey

Dans le cadre du Pojet SFPN

gragg.c
------

Les fonctions utilisees dans la methode de Gragg
pour la resolution de l'equation seculaire f(x)=0.
*/

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include "gragg.h"
#include "param.h"



/*calcul le lambda initial sur l'interval monotone*/
double initial_monotone(const PARAM *p,int k){
  double DELTA,a,b,c;
  double arg=(p->delta[k]+p->delta[k+1])/2;
	DELTA=p->delta[k+1] - p->delta[k];
  c=fun_g(p,arg,k);
  if(fun_f(p,arg)>=0){
    a=c*DELTA+pow(p->zeta[k],2)+pow(p->zeta[k+1],2);
    b=pow(p->zeta[k],2)*DELTA;
  }else{
    a=-c*DELTA+pow(p->zeta[k],2)+pow(p->zeta[k+1],2);
    b=pow(p->zeta[k+1],2)*DELTA;
  }
  return a<=0?(a - sqrt(fabs(a*a - (4*b*c))))/(2 * c)+p->delta[k]:(2 * b)/(a + sqrt(fabs(a*a - (4*b*c))))+p->delta[k];

}

/*calcul le lambda initial sur l'interval non monotone*/
double initial_non_monotone(const PARAM *p){
  double DELTA,a,b,c,arg,delta_n_plus1;
  int n=p->DIM-1;
  delta_n_plus1=deltaNPlus1(p);
  arg=(p->delta[n]+delta_n_plus1)/2;
	DELTA=p->delta[n] - p->delta[n-1];
  c=fun_g(p,arg,n-1);
  a=-c*DELTA+pow(p->zeta[n-1],2)+pow(p->zeta[n],2);
  b=-pow(p->zeta[n],2)*DELTA;
  if(fun_f(p,arg)<=0){
    if(c<=-fun_h(p,delta_n_plus1,n-1)){
      return delta_n_plus1;
    }else{
      return a>=0?(a + sqrt(fabs(a*a - (4*b*c))))/(2 * c)+ p->delta[n]:(2 * b) / (a - sqrt(fabs(a*a - (4*b*c))))  +p->delta[n];
    }
  }else{
    return a>=0?(a + sqrt(fabs(a*a - (4*b*c))))/(2 * c)+ p->delta[n]:(2 * b) / (a - sqrt(fabs(a*a - (4*b*c))))  +p->delta[n];
  }

}

/*Le lambda suivant egal au lamba precedent(y)
ajouté à une valeure d'approximation(n)*/
double y_next(const PARAM *p, double y,int k){
  double DELTA_K,DELTA_KPLUS1,a,b,c,n;
  double val_f=fun_f(p,y);
  double val_fp=fun_fprime(p,y);
  double val_fs=fun_fseconde(p,y);
  DELTA_K=p->delta[k] - y;
  DELTA_KPLUS1=k<(p->DIM-1)?(p->delta[k+1]-y):deltaNPlus1(p);
  a=(DELTA_K+DELTA_KPLUS1)*val_f-DELTA_K*DELTA_KPLUS1*val_fp;
  b= (DELTA_K *DELTA_KPLUS1 )*val_f;
  c=val_f -(DELTA_K+DELTA_KPLUS1)*val_fp+DELTA_K*DELTA_KPLUS1*val_fs/2;
  n=a>0?(2*b / (a + sqrt(fabs(a*a - 4*b*c)))):((a - sqrt(fabs(a*a - 4*b*c))) / 2*c);
  return y+n;
}

/*Determine le critere d'arret sur l'interval monotone*/
bool stop_monotone(double *y,double *y1,double *y2,int j){
  *y2=*y1;*y1=*y;
  if(j<2){
    return 0;
  }else{
    return (*y1 - *y2) * (*y-*y1)<=0;
  }
}

/*Determine le critere d'arret sur l'interval non monotone*/
bool stop_non_monotone(const PARAM* p,double y,int k){
  double DELTAK,tho,e,r1,r2,som=0;
  double delta_k ,delta_kplus1;
  delta_k=p->delta[k];
  delta_kplus1=k<(p->DIM-1)?p->delta[k+1]:deltaNPlus1(p);
  DELTAK=fabs(y-delta_k) < fabs(y-delta_kplus1)?delta_k:delta_kplus1;
  tho=y-(DELTAK);
  for(int j=0; j<=k; j++) {
		som+=(k-j+6) * fabs((p->zeta[j]*p->zeta[j])/(p->delta[j] -y));
  }
  for(int j=k; j<=p->DIM-1; j++) {
      som+=(j-k+5) * fabs((p->zeta[j]*p->zeta[j])/(p->delta[j] -y));
	}
  e=2*p->rho+som+fabs(fun_f(p,(DELTAK+tho)));
  r1=e*MACHEPS+MACHEPS*fabs(tho)* fabs(fun_fprime(p,DELTAK-tho));
  r2=fabs(fun_f(p,DELTAK+tho));
  return infOuEgale(r2,r1);
}

/*Trouve le zero du gragg sur l'interavle
  delta_k et delta_kplus1*/
Secular * lambda_gragg(const PARAM * p,int k){
  double y,y1=0,y2=0;
  bool stop=false;
  Secular *gr=malloc(sizeof(Secular));
  gr->nbIter=1;
  y=k<(p->DIM-1)?initial_monotone(p,k):initial_non_monotone(p);
  while(!stop){
    if((p->rho > 0 && k != p->DIM-1) || (p->rho < 0 && k != 0)) {
      stop=stop_monotone(&y,&y1,&y2,gr->nbIter);
    }else{
      stop=stop_non_monotone(p,y,k);
    }
    y=y_next(p,y,k);
    gr->nbIter++;
    gr->lambda=y;

  }
  return gr;
}

 /* La fonction gragg */
 Secular * gragg(const PARAM *p){
   Secular *tab=(Secular*)malloc(p->DIM*sizeof(Secular));
   Secular *tmp;
   for(int i=0;i<p->DIM;i++){
     tmp=lambda_gragg(p,i);
     tab[i].nbIter=tmp->nbIter;
     tab[i].lambda=tmp->lambda;
     free(tmp);
   }
   return tab;

 }
