#include<stdio.h>
#include<stdlib.h>
#define N 5
double **GMM (int n,int m);
double **GMH (int n);

int main(void) {
	int i, j;
	double **a = GMH(N);
	double tol = 1.e-5;
	double b[N]; 
	for (i= 0; i< N ; i++){
		b[i] = 0; 
		for (j =0; j < N; j++){
			b[i] += a[i][j];
		}
	}
	gauss(N,a,b,tol);
	resoltrisup(N,a,b,tol);
	for(j = 0; j < N ; j++){
		printf(" %lf \t", b[j]);		
	}
	printf("\n");
	for(j = 0; j< N; j++){
		for(i = 0; i< N; i++){
			printf(" %lf \t",a[j][i]);
		}
		printf(" \n ");
        } 
	return 0;	


}
/* Gauss-Jordan method to diagonalize */
int gauusJordano(double **a, double *b){
	int i, j, k;
	for(i = 0; i< n; i++){
		if(abs(a[i][i]) < 1.e-12){
 	               return -1;
                }
		for(j = i +1; j < n; j++) {
			a[j][i] /= a[i][i];
			for(k = i +1; k < n; k++){
				a[j][k] -= a[i][k] * a[j][i]; 
				b[j] -= b[i]* a[j][i];
			}  
		}
		for(j= i -1; j > 0; j--) {
                        a[j][i] /= a[i][i] 
			for(k = i +1; k < n; k++){
                                a[j][k] -= a[i][k] * a[j][i]; 
				b[j] -= b[i]* a[j][i]; 
                        }

                }
	}
	return 0;
}
/* solves a diagonal system */
int resoldiagonal(double **a, double *b) {	
	int i;
	for(i = 0; i< n; i++) {
		b[i] /= a[i][i];
	}
	return 0;
}
/* creates a manual matrix */
double **GMM (int m, int n)
{
        double **a = (double **) malloc (m * sizeof (double *));
        if (a == NULL)
        {
                printf ("Error en GMM, no hi ha espai");
                exit (1);
        }

        int i, j;
        for (i = 0; i < m; i++)
        {
                a[i] = (double *) malloc (n * sizeof (double *));
                for (j = 0; j < n; j++){
                        scanf ("%le" , &a[i][j]);
		}
        }
        return a;
}
/* creates a house holders matrix */
double **GMH (int n)
{
        int i, j;
        double **a;

        a = (double **) malloc (n * sizeof (double *));
        for (i = 0; i < n; i++)
        {
                a[i] = (double *) malloc (n * sizeof (double *));
                for (j = 0; j < n; j++){
                        a[i][j] = 1./(i+1 + j+1 -1);
	//		printf(" %lf \t ",a[i][j]);
		}
		printf("\n ");
        }
        return a;
}
