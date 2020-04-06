#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>

#include "gragg.h"
#include "param.h"
/**
 * y0 initial sur 0<=k<n-2
 */
double initial_monotone(PARAM *p,int k){
  double DELTA,a,b,c;
  double arg=(p->delta[k]+p->delta[k+1])/2;
	DELTA=p->delta[k+1] - p->delta[k];
  c=g(p,arg,k);
  if(f(p,arg)>=0){
    a=c*DELTA+pow(p->zeta[k],2)+pow(p->zeta[k+1],2);
    b=pow(p->zeta[k],2)*DELTA;
  }else{
    a=-c*DELTA+pow(p->zeta[k],2)+pow(p->zeta[k+1],2);
    b=pow(p->zeta[k+1],2)*DELTA;
  }
  if(a<=0){
    return (a - sqrt(fabs(a*a - (4*b*c))))/(2 * c)+p->delta[k];
  }else{
    return (2 * b)/(a + sqrt(fabs(a*a - (4*b*c))))+p->delta[k];
  }
}
/**
 * y0 initial sur k=n-1
 */
double initial_non_monotone(PARAM *p){
  double DELTA,a,b,c,somme;
  int n=p->DIM-1;
  double arg=(p->delta[n]+p->delta[n-1])/2;
	DELTA=p->delta[n] - p->delta[n-1];

  c=g(p,arg,n);

  if(f(p,arg)<=0){
    if(c<=-h(p,p->delta[n],n))   {
      somme=0;
      for(int i=0; i<n; i++) {
        somme+=pow(p->zeta[i],2);
      }
      return p->zeta[n] + somme/p->rho;
    }else{
      a=-c*DELTA+pow(p->zeta[n-1],2)+pow(p->zeta[n],2);
      b=-pow(p->zeta[n],2)*DELTA;
      if(a>=0){
        return (a + sqrt(fabs(a*a - (4*b*c))))/(2 * c)+ p->delta[n];
      }
      else{
        return (2 * b) / (a - sqrt(fabs(a*a - (4*b*c))))  +p->delta[n];
      }
    }
  }else{
    a=-c*DELTA+pow(p->zeta[n-1],2)+pow(p->zeta[n],2);
    b=-pow(p->zeta[n],2)*DELTA;
    if(a>=0){
      return (a + sqrt(fabs(a*a - (4*b*c))))/(2 * c)+ p->delta[n];
    }
    else{
      return (2 * b) / (a - sqrt(fabs(a*a - (4*b*c))))  +p->delta[n];
    }
  }

}

/**
 * le y suivant sachant y_prec precedent
 */
 double y_suivant(double y_prec,double a,double b,double c){
   double apr=a<=0?(a - sqrt(a*a - (4*b*c)))/(2 * c):(2 * b)/(a + sqrt(a*a - (4*b*c)));
   //printf("approch:%f\n",a);

   return y_prec+apr;
 }
 /**
  * critere d'arret de l'interval non monotone
  */

int stop_non_monotone(PARAM* p,double y,int k){
  double t,som1=0,som2=0,e,r;
  int K;
  double em=pow(2, -52);
  t=y-(p->delta[k]);
	K=fabs(y-p->delta[k]) < fabs(y-p->delta[k+1])?k:k+1;

  for(int j=0; j<k; j++) {
		som1+=(k-j+6) * fabs(pow(p->zeta[j],2)/(p->delta[j] -y));
  }
  for(int j=k; j<p->DIM; j++) {
      som2+=(j-k+5) * fabs(pow(p->zeta[j],2)/(p->delta[j] -y));
	}
  e=2*p->rho+som1+som2+fabs( f(p,(p->delta[K]+t)));
  //printf("le e %f\n",e);

  r=e*em+em*fabs(t)* f_prime(p,p->delta[k]-t);
  double r1=fabs(f(p,p->delta[k]+t));

  return r1<=r;
}
/**
 * critere d'arret de l'interval monotone
 */
int stop_monotone(double *y,double *y1,double *y2,int j){
  *y2=*y1;*y1=*y;
  if(j<2){
    return 0;
  }else{
    return (*y1 - *y2) * (*y-*y1)<=0;
  }
}

/**
 * trouver le zero de gragg
 */

double zero_gragg(PARAM * p,int k){
  double a, b, c,y,y1,y2;
  int i=0,stop=0;
  if(k<p->DIM-1){y=initial_monotone(p,k);}
  else{ y=initial_non_monotone(p);}
  //printf("y0%f\n",y);
  //printf("y0  %f",y);

  while(!stop && i<10){

    calcul_a_b_c(p,&a,&b,&c,y,k);
    /*if(i==0){
      printf("a:\n b:%f\n c %f\n",a,b,c);
    }*/
    if((p->rho > 0 && k != p->DIM-1) || (p->rho < 0 && k != 0)) {
      stop=stop_monotone(&y,&y1,&y2,i);
    }else{
      stop=stop_non_monotone(p,y,k);
    }
    y=y_suivant(y,a,b,c);
    i++;
    //printf("y%d :%f\n",i,y);

  }
  return y;
}
