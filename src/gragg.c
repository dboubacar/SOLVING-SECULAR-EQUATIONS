#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define TAILLE_D 4

int main(int argc, char const *argv[]) {
  printf("akdjigeuz");
  /* on choisit arbitrairement beta = 10^-3*/
  double a, y,  b, c,t ,appr, rho, del, h, h_prime, h_second, s, s_1, c_1, f, f_prime;
  double del_k, del_k1;
  double em=pow(2,-52);
  double beta = pow(10,-3);
  double tab[TAILLE_D];
  tab[0] = 1;
  tab[1] = 2 - beta;
  tab[2] = 2 + beta;
  tab[3] = pow((double) 3,(double) (1/3));
  double ui[TAILLE_D]  = {1,2,3,4};
  rho = sqrt(pow(2,2) + pow(beta,2) + pow(beta,2) + pow(2,2));
  rho  = pow((double) rho, -2);
  printf("rho eest de valeur : %f", rho);

  double sol[TAILLE_D];
  printf("debut de la boucle");

  for( int k =0; k < TAILLE_D - 1;k++){
    /* Calcul de du y initial */
    double z = 0;
    /* Savoir duquel delta on est le plus proche */
    if ( k == TAILLE_D - 1){
      /*On choisit arbitrairement un delta k+1 */
      double df_k1 = tab[k] + 1 ;
      y =  5;
      del_k = tab[k] - y;
      del_k1 = tab[k+1] - y;

    }else{



    for (int i = 0; i<TAILLE_D; i++){
      z += pow(ui[i],2)/(tab[i] - ((tab[k]+tab[k+1])/2));
    }
    z = 1 + rho * z;
    double g = 0;
    for (int i = 0; i<TAILLE_D; i++){
      if ((i != k) && (i != k+1)){
        g += pow(ui[i],2)/(tab[i] - ((tab[k] - tab[k+1])/2));
      }
    }
    g =  rho * g;
    if(z >= 0 ){
      a = g * (tab[k] - tab[k+1]) + (pow(ui[k],2) + pow(ui[k+1],2));
      b = pow(ui[k], 2) * (tab[k] - tab[k+1]);
      del = tab[k];

    }else{
      a = -g * (tab[k] - tab[k+1]) + (pow(ui[k],2) + pow(ui[k+1],2));
      b = - pow(ui[k+1], 2) * (tab[k] - tab[k+1]);
      del = tab[k+1];
    }
    if (a <= 0 ){
      t = (a - sqrt(a*a - (4*b*g)))/(2 * g);
    }else{
      t = (2 * b)/(a + sqrt(a*a - (4*b*g)));
    }


    /* y esy le lambda initial */
    y = t + del;
    printf("la valeur de y est %f",y);
    del_k = tab[k] - y;
    del_k1 = tab[k+1] - y;
    }
    /* Iintialisation de la boucle */
    h = 0;
    h_prime = 0;
    h_second = 0;
    int stop=1;
    while(stop){

      for (int i = 0; i<TAILLE_D; i++){
          h +=  pow(ui[i],2)/(tab[i] -y);
          h_prime +=  (rho * pow(ui[i],2))/(pow((tab[i] -y), 2));
          h_second +=  ((2 * rho * pow(ui[i],2) * tab[i]) - (2 * rho * pow(ui[i],2) * y) )/(pow((tab[i] -y), 4));
        }
      h = 1 + rho * h;
      s = ((pow(del_k,3) * del_k1)/(del_k - del_k1)) * ((h_prime/del_k1) - (h_second/2));
      s_1 = ((pow(del_k,3) * del_k1)/(del_k1 - del_k)) * ((h_prime/del_k) - (h_second/2));
      c_1 = h - (del_k - del_k1 ) * h_prime + del_k * del_k1 * (h_second/2);
      //Calcul du f(y),f'(y) et f"(y)
      double f_y=c_1+(s/(tab[k]-y))+(s_1/(tab[k+1]-y));
      double f1_y=s/(pow(tab[k]-y,2)) + s_1/(pow(tab[k+1]-y,2));
      double f2_y=2*(s/(pow(tab[k]-y,3)) + s_1/(pow(tab[k+1]-y,3)));
      //Calcul a1,b2,c2 pour l'approximation
      double a1=(del_k+del_k1)*f_y-del_k*del_k1*f1_y;
      double b1=(del_k+del_k1)*f_y;


      //Calcul de l'approximation

      appr=2*b1/(a1+sqrt(a1*a1-4*b1*c_1));

      //critere d'arret
      double f_dto=c_1+(s/(tab[k]-(tab[k]+t)))+(s_1/(tab[k+1]-(tab[k]+t)));
      double f1_dto=s/(pow(tab[k]-(tab[k]-t),2)) + s_1/s/(pow(tab[k+1]-(tab[k]-t),2));

      double s1=0,s2=0;
      for (int i = 0; i<TAILLE_D; i++){
        double tm=abs(pow(ui[i],2)/(tab[i]-y));
        s1 += (k-i+6)*tm;
        if(i>=k+1){
          s2 +=(i-k+5)*tm;
        }
      }
      double e=2*rho+s1+s2+f_dto;
      stop=abs(f_dto)<=(e*em+em*abs(t)*f1_dto)?0:1;
      y=y+appr;

    }
    sol[k]=y;


  }
  for (int i = 0; i<TAILLE_D; i++){
    printf("Solution de %d est : %f \n",i,sol[i]);
  }

  return 0;
}
