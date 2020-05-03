#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "gragg.h"
#include "param.h"
#include "hybrid.h"
#include <math.h>


double my_gettimeofday(){
	struct timeval tmp_time;
	gettimeofday(&tmp_time, NULL);
	return tmp_time.tv_sec + (tmp_time.tv_usec * 1.0e-6L);
}

int main(int argc, char **argv)
{
  Secular *secu1=NULL;
  Secular *secu2=NULL;
  int choix,sizeMatrice=4;
  double time1,time2,debut;
  char titre[200];
  if(argc == 3) {
    choix=atoi(argv[1]);
    sizeMatrice=atoi(argv[2]);
    if((choix<0||choix>2)||sizeMatrice<2){
      help();
      exit(1);
    }
  }
  else {
    help();
    exit(1);
  }

  const PARAM * p=init_param(sizeMatrice);
  if(choix==0){
    printf("\033[22;36m*************SOLVING SECULAR EQUATIONS WITH GRAGG'S SCHEME\n\033[0m");
    debut=my_gettimeofday();
    secu1=gragg(p);
    time1=my_gettimeofday()-debut;
  }else if (choix==1){
    printf("\033[22;36m*************SOLVING SECULAR EQUATIONS WITH HYBRID SCHEME\n\033[0m");
    debut=my_gettimeofday();
    secu1=hybrid(p);
    time1=my_gettimeofday()-debut;
  }else{
    printf("\033[22;36m*********SOLVING SECULAR EQUATIONS WITH GRAGG'S & HYBRID\n\033[0m");
    debut=my_gettimeofday();
    secu1=gragg(p);
    time1=my_gettimeofday()-debut;
    debut=my_gettimeofday();
    secu2=hybrid(p);
    time2=my_gettimeofday()-debut;

  }

  printf("\033[22;36mDATA:\n\033[0m");
  printf("\033[22;32mzeta:\033[0m");
  printf("[");
  for(int i=0;i<p->DIM;i++){
    if(i<p->DIM-2){
      printf("%f ",p->zeta[i]);
    }else{
      printf("%f",p->zeta[i]);
    }
  }
  printf("]\n");
  printf("\033[22;32mdelta:\033[0m");
  printf("[");
  for(int i=0;i<p->DIM;i++){
    if(i<p->DIM-2){
      printf("%f ",p->delta[i]);
    }else{
      printf("%f",p->delta[i]);
    }
  }
  printf("]\n");
  printf("\033[22;32mrho:\033[0m");
  printf("%f\n",p->rho);

  if(choix==2){
    strcpy(titre, "SOLUTION WITH GRAGG'S:(lambda,iteration):");
    print_secular(secu1,sizeMatrice,titre,time1);
    strcpy(titre, "SOLUTION WITH HYBRID:(lambda,iteration):");
    print_secular(secu2,sizeMatrice,titre,time2);
  }else{
    strcpy(titre, "SOLUTION:(lambda,iteration):");
    print_secular(secu1,sizeMatrice,titre,time1);
  }
	//PARAM *tmp=p;
	free_data(p,secu1,secu2);
  return 0;
}
