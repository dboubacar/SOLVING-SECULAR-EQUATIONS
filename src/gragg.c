#include <stdio.h>
#include <stdlib.h>

#define TAILLE_D 4

double f_y(double d1, double d2 , double y, double c, double s, double S){
  return c + (s/(d1 - y)) + (S/(d2 - y));
}

double f1_y(double d1, double d2 , double y, double s, double S){
  return  (s/((d1 - y)*(d1 - y))) + (S/((d2 - y)*(d2 - y)));
}

double f2_y(double d1, double d2 , double y, double s, double S){
  return 2*((s/((d1 - y)*(d1 - y)*(d1 - y))) + (S/((d2 - y)*(d2 - y)*(d2 - y))));
}


double w_0(double tabd[], double tabu[], double d1, double d2, double y, int k ){
  double  sw0 = 0;
  for( i = 0; i < TAILLE_D; i++){
    if ( (i != k) && (i != (k+1))){
      sw0 +=  (tabu[i]*tabu[i])* ((tabd[i] - d2)/((d1 - y)*(d1 - y)*(d1 - y)));
    }
  }
  sw0 = (tabu[k]*tabu[k]) + (((d1 - y)*(d1 - y)*(d1 - y))/(d1 - d2)) * sw0;
  return sw0;
}

double w_1(double tabd[], double tabu[], double d1, double d2, double y, int k ){
  double  sw1 = 0;
  for( i = 0; i < TAILLE_D; i++){
    if ( (i != k) && (i != (k+1))){
      sw1 +=  (tabu[i]*tabu[i])* ((tabd[i] - d1)/((tabd[i] - y)*(tabd[i] - y)*(tabd[i] - y)));
    }
  }
  sw1 = (tabu[k+1]*tabu[k+1]) + (((d2 - y)*(d2 - y)*(d2 - y)))/(d2 - d1)) * sw1;
  return sw1;
}



int main(int argc, char const *argv[]) {
  double tabd[TAILLE_D] = {1,2,3,4};
  double ui[TAILLE_D]  = {2,2,2,2};
  double sol[TAILLE_D];

  double a, b, c, appr, rho;

for( int j =0, j < TAILLE_D - 1, j++){
    k = j;
    while(true){

    }

  }

  return 0;
}
