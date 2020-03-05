#include <stdlib.h>
#include <stdio.h>
#define TAILLE_Z 10
#define TAB 4

int main(int argc, char *argv[])
{


/* Calculer les valeurs propres de D et A? puis les stocker dans un tableau, ici on simulera des valeurs */

	double delta[4] = {1,2,3,4};
	double v[4] = {1,2,3,4};

/*Initialisation */
	double  y,w0,w1,p,al;
	
/* Choisir une valeur sur la fonction, p est initialisé arbitrairement très petit */

	p = -10000;
	al = 1;
/* Choix des dk et dk+1 , stockage de leure indice */
	double dk, dka;
	int k,ka; /* ka = k+1:*/
	
/* début de la boucle de caclcul du zéro */
	for( int j =0, j < TAILLE_Z, j++}{
		k = j;
		if (j == TAILLE_Z - 1){ /*Dans ce cas, on choisit un nombre arbitrairement */
			ka = 11;
		}else{
			ka = j+1;
		}
		double x = 2; /* Normalement, la valeur de x dépendra du choix de l'intervalle */
		while( p < x){
			double sy = 0;
			/* Boucle de la somme pour y */
			for(int i=0; i<TAB; i ++){
				if ( (i != k) && (i != ka){
					sy +=  (v[i]**2/(delta[i] - x)) * ((delta[i] - delta[k])/(delta[i] - x)) * ((delta[i] - delta[ka])/(delta[i] - x));
				}
			}
			y = 3*x - ( al + delta[k] + delta[ka] ) + sy;
		
			double sw0 = 0;
			/* Boucle de la somme pour w0 */
			for(int i=0; i<TAB; i ++){
				if ( (i != k) && (i != ka){
					sw0 +=  (v[i]**2/(delta[i] - x)**2)  * ((delta[i] - delta[ka])/(delta[i] - x));
				}
			}
			w0 = v[k]**2 + ((delta[k] - x)**3/(delta[k] - delta[ka]))*(1 + sw0);
			
			/* Boucle de la somme pour w1 */
			for(int i=0; i<TAB; i ++){
				if ( (i != k) && (i != ka){
					sw1 +=  (v[i]**2/(delta[i] - x)**2)  * ((delta[i] - delta[k])/(delta[i] - x));
				}
			}
			w1 = v[ka]**2 + ((x - delta[ka])**3/(delta[k] - delta[ka]))*(1 + sw0);
			
			p = y + (w0/(delta[k] - x)) + (w1/(delta[ka] - x));
			
			x = p;
		}
		
		printf("La %d eme valeur de lambda est %d", j, x);
	}
	
	return 0;


}
