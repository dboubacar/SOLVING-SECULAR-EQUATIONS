/*
Résolution de l’équation séculaire
  Par: Boubacar Diallo  et
       Jeffrey  Kebey

Dans le cadre du Pojet SFPN

hybrid.h
------

Contient les prototypes des fonctions de hybrid.c
*/


#ifndef DEF_HYBRID
#define DEF_HYBRID
#include "param.h"


double middleWay(const PARAM *p,double y,int k,bool useF);
double fixedWeight(const PARAM *p,double y,int k,bool useF);

double commute(const PARAM *p,double y,int k ,bool isFixe);

Secular * lambda_hybrid(const PARAM * p,int k);

Secular * hybrid(const PARAM *p);










#endif
