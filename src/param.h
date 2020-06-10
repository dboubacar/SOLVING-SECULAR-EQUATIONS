/*
Résolution de l’équation séculaire
  Par: Boubacar Diallo  et
       Jeffrey  Kebey

Dans le cadre du Pojet SFPN

gragg.h
------

Contient les prototypes des fonctions de param.c
*/


#ifndef DEF_PARAM
#define DEF_PARAM

#include <stdbool.h>
#include <time.h>
#include <sys/time.h>
#include <omp.h>

#define  EPSILON  1e-15//erreur abolue
#define  MACHEPS pow(2,-52) //unité d'arrondi

//structure pour les parametres de la fonction f(x)
typedef struct{
  int DIM;
  double rho;
  double *zeta;
  double *delta;
}PARAM;
/*Structure pour le resultat de deux methode*/
typedef struct{
  int nbIter;//nombre d'iteration
  double lambda;//le lamba trouvé par la methode
}Secular;

PARAM * init_param(int dimension);
double fun_f(const PARAM *p,double x);
double fun_fprime(const PARAM *p,double x);
double fun_fseconde(const PARAM *p,double x);
double fun_g(const PARAM *p,double x,int k);
double fun_h(const PARAM *p,double x,int k);
double fun_psiprime(const PARAM *p,double x,int k);
double fun_phiprime(const PARAM *p,double x,int k);
double fun_fm(const PARAM *p,double x,int k);
double fun_fmprime(const PARAM *p,double x,int k);
bool infOuEgale(double a, double b);
double deltaNPlus1(const PARAM *p);
void print_secular(Secular *secular,int n,char titre[],double time);
void help();
double wtime();
void free_data(const PARAM *p,Secular *s1,Secular *s2);



#endif
