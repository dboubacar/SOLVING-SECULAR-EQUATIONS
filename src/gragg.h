/*
Résolution de l’équation séculaire
  Par: Boubacar Diallo  et
       Jeffrey  Kebey

Dans le cadre du Pojet SFPN

gragg.h
------

Contient les prototypes des fonctions de gragg.c
*/


#ifndef DEF_GRAGG
#define DEF_GRAGG
#include "param.h"


double initial_monotone(const PARAM *p,int k);
double initial_non_monotone(const PARAM *p);
double y_next(const PARAM *p, double y,int k);
bool stop_non_monotone(const PARAM* p,double y,int k);
bool stop_monotone(double *y,double *y1,double *y2,int i);
Secular * lambda_gragg(const PARAM * p,int k);
Secular * gragg(const PARAM *p);




#endif
