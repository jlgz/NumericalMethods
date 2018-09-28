#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#define M_PI  3.14159265358979323846


void difdiv(int n,  double *x, double *f);
double horner(double z,int n, double *x, double *c);
void  equidistantes(double a,double b,double *x, int n);
double fun( double x);
void writefun(double *x, double *f,int n);
int main(void) {
	int n =4,i;
	double x[n+1],f[n+1];  
	equidistantes(-1,1,x,n);
	for(i=0;i<=n;i++){
		f[i]=fun(x[i]);
	}
	difdiv(n,x,f);
	writefun(x,f,n); 
	return 0;
}

void difdiv(int n,  double *x, double *f){
	int i,j,k;
	for( k = 0; k< n; k ++){
		for(i = n; i> k; i--) {
			f[i] =( f[i] - f[i-1])/ (x[i] - x[i-k-1]);	
		}	
	}
}
double horner(double z,int n, double *x, double *c){ 
	double r = c[n];
	int i;
	for(i =n -1; i>= 0; i --){
		r = r*(z -x[i]) + c[i];
	}
	return r;
}
// a<b
void  equidistantes(double a,double b,double *x, int n){
	int i;
	double dif = (b- a) / n;
	for(i = 0; i <=  n; i ++) {
		x[i] = a + dif* i;
	} 
}
/*
void Chebyshev(double a, double b, double *x) {
	for(i =0; i<n; i++) {
		x[i] = cos(M_PI*(2*i +1) / (2*(n+1)));
	}
}
*/
double fun(double x){
	return 1/(1+25*pow(x,2));
}
void writefun(double *x, double *f,int n)	{
	FILE *fp;
	int i;
	double tmp;
	fp= fopen("fichero.txt", "w");
	for (i = 0; i<= n;i++){
		tmp = horner(x[i],n,x,f);
		fprintf(fp,"%lf \t %lf \n" , x[i],tmp);
	}
	fclose(fp);
}
/* orden n= 4
 k = 0
 f[x2,x3] = f[3] - f[2] / x[3] - x[2] = f[3]
 f[2] = f[2] - f[1] / x[2] - x[1] = f[x1,x2]
 ...............
k = 1
f[3] = f[x2,x3] - f[x1,x2] / x[3] - x[1] = f[x1,x2,x3]
f[2] = f[x1,x2] - f[x0,x1] / x[2] -x[0] = f[x0,x1,x2]
k = 2
f[3] = f[x1,x2,x3] - f[x0,x1,x2] / x[3] - x[0] = f[x0,x1,x2,x3]
*/
