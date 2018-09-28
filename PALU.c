#include<stdio.h>
#include<stdlib.h>
#define N 5
int resoltrisup(int n, double **a, double *b, double tol);
int gauss (int n, double **a, double *b, double tol);
int gauss2 (int n, double **a, int *p, double tol);
int palu(int n, double **a, int *p, double tol);
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
			printf(" %la \t",a[j][i]);
		}
		printf(" \n ");
        } 
	return 0;	


}
int resoltrisup (int n, double **a, double *b, double tol){
	int i,j;
	for(i = n-1; i>= 0; i--) {
		for(j=i+1; j<= n-1; j++){
			b[i] -= a[i][j]*b[j];
		}
//		if (abs(a[i][i])< tol){
//			return 1;
//		}
//		else{	
			b[i]= b[i]/ a[i][i];
//		}
	}
	return 0;
}
/*simple pivoting */
int gauss (int n, double **a, double *b, double tol){
	int i,k,j;
	for(k=0;k<=n-1;k++){
		for(i=k+1;i<=n-1;i++){
//			if ((a[k][k] != 0) && (abs(a[k][k])> tol)){ 
				a[i][k] =a[k][k]/ a[i][k]; 
//			}
//			else { return 1;	}
			for(j=k+1;j<=n-1;j++){
				a[i][j]-= a[i][k]*a[k][j];
			}
			 b[i] -= a[i][k]*b[k];
		}
	}
	return 0;

}
/* relative partial pivot */
int gauss2 (int n, double **a, int *p, double tol){

	int i,k,j,l,tmp,tmp1;
	double *tp;
	for(k=0;k<=n-1;k++){
		tmp = k;
		for(i=k+1;i<=n-1;i++){
			for(l=k;l<=n-1;l++){
				if(abs(a[l][k]) > abs(a[tmp][k])){	
					tmp = l;
				}
			}
			if(tmp!=k){ 
				tp = a[k];
				a[k]=a[tmp];
				a[tmp]=tp;
                                tmp1= p[k]; 
			 	p[k] = p[tmp];
				p[tmp]=tmp1;
			}
			if(abs(a[k][k]) > tol){
				a[i][k]=a[i][k] / a[k][k]; 
			}
			else { return 0;  }	
			for(j=k+1;j<=n-1;j++){
				a[i][j]-= a[i][k]*a[k][j];
			}
		}
	}
	return 1;
}
/* factorización PA=LU */
int palu(int n, double **a, int *p, double tol){
        int var=gauss2(n,a,p,tol);
        if (var!=0) { return 1; }
        else{ return 0; }
        // fila i de pa es fila p(i) de a
	// L 1 en la diagonal y factores por debajo 
}

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


