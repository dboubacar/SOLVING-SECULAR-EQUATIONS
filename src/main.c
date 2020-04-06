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
    double y=1.764171;

    calcul_a_b_c(p,&a,&b,&c,y,1);
    double yn=y_suivant(y,a,b,c);*/
    double sol0=zero_gragg(p,0);
    double sol1=zero_gragg(p,1);
    double sol2=zero_gragg(p,2);
    double sol3=zero_gragg(p,3);
    //printf("GRAGG y nouveau %f\n",r);

    printf("GRAG \nsol0 :%f\nsol1 :%f\nsol2 :%f\nsol3 :%f\n",sol0,sol1,sol2,sol3);
    return 0;
}
