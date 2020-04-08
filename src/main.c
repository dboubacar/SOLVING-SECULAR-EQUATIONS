#include <stdio.h>
#include <stdlib.h>
#include "gragg.h"
#include "param.h"

int main()
{

    PARAM * p= init_param();
    //double x=1.33282;
    //double r=gragg_zero(p,1);
    //double r=f_seconde(p,x);
    /*double a, b, c;
    double y=1.764171;*/


  //  calcul_a_b_c(p,&a,&b,&c,y,1);
    //double yn=y_suivant(y,a,b,c);*/
    double* resultats= malloc(4*sizeof(double));

    /*resultats[0]=5;
    resultats[1]=2;
    resultats[2]=3;
    resultats[3]=1;*/

    /*rsol0 :1.000000
    sol1 :2.003316
    sol2 :3.003333
    sol3 :0.000000*/
    //double resultats[4];



    gragg(p,resultats);
    /*resultats[0]=zero_gragg(p,0);
    resultats[1]=zero_gragg(p,1);
    resultats[2]=zero_gragg(p,2);
    resultats[3]=zero_gragg(p,3);*/
    //printf("GRAGG y nouveau %f\n",resultats[0]);
   printf("GRAG \nsol0 :%f\nsol1 :%f\nsol2 :%f\nsol3 :%f\n",resultats[0],resultats[1],resultats[2],resultats[3]);
   free_param(p);

    //printf("GRAG \nsol0 :%f\nsol1 :%f\nsol2 :%f\nsol3 :%f\n",sol0,sol1,sol2,sol3);
    return 0;
}
