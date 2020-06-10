/*
Résolution de l’équation séculaire
  Par: Boubacar Diallo  et
       Jeffrey  Kebey

Dans le cadre du Pojet SFPN

param.c
------

Ces fonctions  implementent les fonctions mathematiques utilisees dans
la methode de Gragg et hybride pour la resolution l'equation seculaire.
*/


#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "param.h"



/*
  Definie les parametres de la fonction f(x)
*/
PARAM * init_param(int dimension){
  double norm,som=0.0;
  PARAM *p=malloc(sizeof(PARAM));
  double beta[5]={1,pow(10, -1),pow(10, -3), pow(10, -6), pow(10, -10)};
  srand(time(NULL));
  int betai= rand()%5;
  p->DIM=dimension;
  p->zeta= malloc(sizeof(double)*p->DIM);
  p->delta= malloc(sizeof(double)*p->DIM);
  for(int i=0; i<dimension; i++) {
    p->delta[i]=(rand()/(double)RAND_MAX)+beta[betai]+i;
    p->zeta[i]=(rand()/(double)RAND_MAX)+beta[betai]+1;
    som += pow(p->zeta[i],2);
  }
  norm=sqrt(som);


  for(int i=0;i<dimension;i++){
    p->zeta[i]=p->zeta[i]/norm;
  }
  p->rho=fmod((double)rand(),0.986602)+0.1;
  return p;
}

/*la fonction f(x)*/
double fun_f(const PARAM *p,double x){
  double som=0.0;
  for (int i=0 ; i<p->DIM; i++){
     som += (p->zeta[i]*p->zeta[i])/(p->delta[i]-x) ;
  }
  return 1 + p->rho *som;
}

/*la derivée premiére de la fonction f(x)*/
double fun_fprime(const PARAM *p,double x){
  double som=0.0;
  for (int i=0 ; i<p->DIM; i++){
     som += pow(p->zeta[i],2)/ pow((p->delta[i]-x),2);
  }
  return p->rho *som;
}

/*la derivée seconde de la fonction f(x)*/
double fun_fseconde(const PARAM *p,double x){
  double som=0.0;
  for (int i=0 ; i<p->DIM; i++){
     som += (2*pow(p->zeta[i],2)) / pow((p->delta[i]-x),3) ;
  }
  return p->rho *som;

}

/*la fonction g(x)*/
double fun_g(const PARAM *p,double x,int k){
  double som=0.0;
  for (int j=1 ; j<p->DIM; j++){
     if ((j != k) && (j != k+1)){
       som += (p->zeta[j]*p->zeta[j]) / (p->delta[j]-x) ;
     }
  }
  return (p->rho +som);
}

/*la fonction h(x)*/
double fun_h(const PARAM *p,double x,int k){
  double d1=p->delta[k]-x;
  double d2=p->delta[k+1]-x;
  return ((pow(p->zeta[k],2)) /d1)+((pow(p->zeta[k+1],2))/d2) ;
}

/*la derivée premiére de la fonction psi(x)*/
double fun_psiprime(const PARAM *p,double x,int k){
  double som=0.0;
  for (int i=0 ; i<=k; i++){
     som += (pow(p->zeta[i],2)) / pow((p->delta[i]-x),2);

  }
  return p->rho *som;
}

/*la derivée premiére de la fonction phi(x)*/
double fun_phiprime(const PARAM *p,double x,int k){
  double som=0.0;
  for (int i=k+1; i<p->DIM; i++){
     som += (pow(p->zeta[i],2)) / pow((p->delta[i]-x),2);

  }
  return p->rho *som;
}

/*la fonction fm(x) utilisee dans la menthode hybride pour les 3 poles*/
double fun_fm(const PARAM *p,double x,int k){
  double som=0.0;
  for (int i=0 ; i<p->DIM; i++){
    if(i!=k){
     som += (pow(p->zeta[i],2)) / (p->delta[i]-x) ;
   }
  }
  return 1+p->rho*som;
}

/*la derivée premiére de la fonction phi(x)*/
double fun_fmprime(const PARAM *p,double x,int k){
  double som=0.0;
  for (int i=0 ; i<p->DIM; i++){
    if(i!=k){
     som += (pow(p->zeta[i],2)) / pow(p->delta[i]-x,2) ;
   }
  }
  return p->rho*som;
}

/* Seconds (wall-clock time) since an arbitrary point in the past */
double wtime()
{
	struct timeval ts;
	gettimeofday(&ts, NULL);
	return (double)ts.tv_sec + ts.tv_usec / 1e6;
}

/*cette fonction permet de resourde le probleme d'egalite
en double precision avec une erreur relative*/
bool infOuEgale(double a, double b){
   if(a==b)   return 1;
   double aa= fabs(a);
   double bb= fabs(b);
   double absError=aa-bb;

   return   absError <= EPSILON   ||   ( absError / (aa>bb? aa:bb) )  <=  EPSILON;
}

/*calcul le dernier delta sur l'interval non monotone*/
double deltaNPlus1(const PARAM *p){
  double somme=0;
  int n=p->DIM-1;
  for(int i=0; i<=n; i++) {
    somme+=pow(p->zeta[i],2);
  }
  return p->delta[n]+(somme/p->rho);

}

/*affiche les resultats des methodes*/
void print_secular(Secular *secular,int n,char titre[],double time){
  int max=0,total=0;
//  printf("\033[22;36mSOLUTION:(lambda,iteration)\n\033[0m");
 printf("\033[22;36m%s\033[0m",titre);
  if(n<50)
  printf("\n[");
  for(int i=0;i<n;i++){
    if(n<30){
      if(i<n-1){
        printf("(%f,%d) ",secular[i].lambda,secular[i].nbIter);
      }else{
        printf("(%f,%d)]",secular[i].lambda,secular[i].nbIter);

      }
    }
    total+=secular[i].nbIter;
    if(secular[i].nbIter>max){
      max=secular[i].nbIter;
    }
  }




  printf("\033[22;32m\nTotal iterations:\033[0m");
  printf("%d  ",total);
  printf("\033[22;32mAverage iterations:\033[0m");
  printf("%d  ",total/n);
  printf("\033[22;32mLargest iterations:\033[0m");
  printf("%d  ",max);
  printf("\033[22;32mTime:\033[0m");
  printf("%g  ",time);

  printf("\n");
}

/*une petite aide pour utiliser solv_secular*/
void help(){
  printf("USAGE: solv_secular [method] [size]\n"
          "Solving Secular Equations Stably and Efficiently\n"
          "method:\n"
          "       0:Gragg's\n"
          "       1:Hybrid\n"
          "       2:Gragg's & Hybrid\n"
          "size:  size of matrix\n"
          "Exemple: solv_secular 0 4\n");
  }
void free_data(const PARAM *p,Secular *s1,Secular *s2){
  free(s1);
  free(s2);
  free((double *)p->zeta);
  free((double *)p->delta);
  free((PARAM *)p);
}
