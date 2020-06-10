/*
Résolution de l’équation séculaire
  Par: Boubacar Diallo  et
       Jeffrey  Kebey

Dans le cadre du Pojet SFPN

hybrid.c
------

Les fonctions utilisees dans la methode de Hybrid
pour la resolution de l'equation seculaire f(x)=0.
*/

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include "hybrid.h"
#include "param.h"
#include "gragg.h"


/* La methode Middle Way sur les deux et trois poles*/
double middleWay(const PARAM *p,double y,int k,bool useF){
  double DELTA_K,DELTA_KM,a,b,c,appr,valf,valfprim;
  int n=p->DIM-1;
  DELTA_K=p->delta[k] - y;
  DELTA_KM=k<n?(p->delta[k+1]-y):(p->delta[n-1] - y);
  valf=useF?fun_f(p,y):fun_fm(p,y,k);
  valfprim=useF?fun_fprime(p,y):fun_fmprime(p,y,k);
  a=(DELTA_K+DELTA_KM)*valf-DELTA_K*DELTA_KM*valfprim;
  b= (DELTA_K *DELTA_KM )*valf;
  if(k<n){
    c = valf - DELTA_K*fun_psiprime(p,y,k) - DELTA_KM*fun_phiprime(p,y,k);
    appr= a <= 0 ? (a - sqrt(fabs(a*a - 4*b*c))) / 2*c :  2*b / (a + sqrt(fabs(a*a - 4*b*c)));
  }else{
    c=valf - DELTA_KM * fun_psiprime(p,y,n-1) - (pow(p->zeta[n],2))/DELTA_K;
		appr= a >= 0 ? (a + sqrt(fabs(a*a - 4*b*c))) / 2*c :  2*b / (a - sqrt(fabs(a*a - 4*b*c)));
  }
  return y+appr;
}


/*La methode Fixed Weight sur les deux et trois poles*/
double fixedWeight(const PARAM *p,double y,int k,bool useF){
  double DELTA_K,DELTA_KPLUS,valf,valfprim,a,b,c,appr;
  double delta_k,delta_kplus;
  int n=p->DIM-1;
  valf=useF?fun_f(p,y):fun_fm(p,y,k);
  valfprim=useF?fun_fprime(p,y):fun_fmprime(p,y,k);
  delta_k=p->delta[k];
  delta_kplus=k<n?p->delta[k+1]:deltaNPlus1(p);
  DELTA_K=delta_k - y;
  DELTA_KPLUS=delta_kplus-y;
  a=(DELTA_K+DELTA_KPLUS)*valf-DELTA_K*DELTA_KPLUS*valfprim;
  b= (DELTA_K *DELTA_KPLUS )*valf;
  bool cas=fabs(y - p->delta[k]) < fabs(p->delta[k+1] - y);
  if(!cas && k<n){
    c = valf - DELTA_K*valfprim - pow((p->zeta[k+1]/DELTA_KPLUS),2)*(delta_kplus-delta_k);
  }else{
    c = valf - DELTA_KPLUS*valfprim- (pow(p->zeta[k],2)/pow(DELTA_K,2))*(delta_k-delta_kplus);
  }
  appr= a <= 0 ? (a - sqrt(fabs(a*a - 4*b*c))) / 2*c :  2*b / (a + sqrt(fabs(a*a - 4*b*c)));
  return y + appr;
}

/*Fait une commutation entre usage à deux  ou trois poles
  en utilisant la methode Middle Way ou Fixed Weight selon la situation
  */
double commute(const PARAM *p,double y,int k ,bool isFixe){
  double ynew=0;
	if(fun_fm(p,y,k)> 0){//usage de deux poles
    ynew=isFixe?fixedWeight(p,y,k,true):middleWay(p,y,k,true);
  }else {//usage de trois poles
    ynew=isFixe?fixedWeight(p,y,k,false):middleWay(p,y,k,false);
  }
  return ynew;
}

/*Trouve le zero du hybrid sur l'interavle
  delta_k et delta_kplus1 en basculant entre les deux poles
  */
Secular * lambda_hybrid(const PARAM * p,int k){
  double y,ynew,fnew,fprec;
  Secular *hyb=malloc(sizeof(Secular));
  bool isFixe=true;
  bool stop=false;
  y=k<(p->DIM-1)?initial_monotone(p,k):initial_non_monotone(p);
  fprec=fun_f(p,y);
  ynew=commute(p,y,k,isFixe);
  fnew=fun_f(p,ynew);
  if(fnew< 0 && fabs(fnew) > 0.1 * fabs(fprec)) {
     isFixe=false;
   }
   //stop=stop_non_monotone(p,ynew,k);
   hyb->nbIter=0;
   hyb->lambda=ynew;
   while(!stop){
     fprec=fnew;
     ynew=commute(p,ynew,k,isFixe);
     fnew=fun_f(p,ynew);
     if(fnew * fprec > 0 && fabs(fnew) > 0.1 * fabs(fprec)){
       isFixe=!isFixe;
     }
     stop=stop_non_monotone(p,ynew,k);
     hyb->nbIter++;
     hyb->lambda=ynew;
     if(hyb->nbIter>=4) break;
   }
   return hyb;
}

/* La fonction hybrid*/
Secular * hybrid(const PARAM *p){
  Secular *tab=(Secular*)malloc(p->DIM*sizeof(Secular));
  fprintf(stderr, "\033[22;36m[HYBRID] Starting iterative solver\n\033[0m");
  fprintf(stderr, "\033[22;36m\r         ----->iteration\n\033[0m");
  #pragma omp parallel for
  for(int i=0;i<p->DIM;i++){
    Secular *tmp=lambda_hybrid(p,i);
    tab[i].nbIter=tmp->nbIter;
    tab[i].lambda=tmp->lambda;
    free(tmp);
  }
  return tab;
}
