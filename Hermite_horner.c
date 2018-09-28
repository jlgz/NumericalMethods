#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#define M_PI  3.14159265358979323846
void difdivherm( int m, double *x,double *f);
void horner2( int m, double z, double *x, double *f,double p[2]);
void  equidistantes2(double a,double b,double *x, int n,int m);
void writefun(double *x, double *f,int n);
double fun(double x);
double dfun(double x);
int main(void) {
	int n = 10;
	int m = 2*n +1;
	int i;
	double p[2];
	double f[m+1];
	double x[m+1];
	equidistantes2(-1,1,x,n,m);
	for( i=0; i< m ; i++) {
		f[i] = fun(x[i]);
		f[i + 1] = dfun(x[i]);
		i++;
	} 	
	difdivherm(m,x,f);
	horner2(m,1,x,f,p);
	writefun(x,f,m);
	return 0;
}
/* hermite interpolation, x axes, f function and derivative values*/
void difdivherm( int m, double *x,double *f){
	int i,j,k;
	for(i = m; i> 0; i--) {
		if ( x[i] != x[i-1] ) {
                f[i] =( f[i] - f[i-2])/ (x[i] - x[i-1]); 
   		 }
	}
	for( k = 1; k< m; k ++){
		for(i = m; i> k; i--) {
			f[i] =( f[i] - f[i-1])/ (x[i] - x[i-k-1]);	
		}	
	}
}
/* f = cos(x)
f(0) = 1 , df(0) = 0, f(pi) = -1 , df(pi) = 0
x = [0,0,pi,pi] f = [1,0,-1,0]

k = 0: orden 2:
x[3] == x[2] ---> f[3] = f[3] = df(pi) = f[pi,pi] 
!= ---> f[2] = (f(pi) - f(0)) / (pi - 0) = -2 / pi = f[0,pi]
== ----> f[1] = df(0) = 0= f[0,0] 
k = 1: orden 3: 
f[3] = (f[pi,pi] - f[0,pi]) / (pi - 0) = f[0,pi,pi]
f[2] = (f[0,pi] - f[0,0] ) / pi- 0 = f[0,0,pi] 
k=2
f[3] = f[0,0,pi,pi]*/
/* horner evaluation method, p[0] = f(z), p[1] = df(z) */
void horner2( int m, double z, double *x, double *f,double p[2]){
	int j;
	p[0] = f[m];
	p[1] = 0;
	for(j = m-1; j>= 0; j--) {
		
		p[1] =p[1] *(z- x[j]) + p[0];
		p[0] = p[0] * (z - x[j]) + f[j];

	}
}
/*creates equidistant axes */
void  equidistantes2(double a,double b,double *x, int n,int m){
	int i,h=0;
	double dif = (b- a) / n;
	for(i = 0; i <=  n; i ++) {
		x[i+h] = a + dif* h;
		x[i+1+h]= x[i+h];
		h ++;
	} 
}
void writefun(double *x, double *f,int n)	{
	FILE *fp;
	int i;
	double tmp[2];
	fp= fopen("fichero.txt", "w");
	for (i = 0; i<= n;i++){
		horner2(n,x[i],x,f,tmp);
		fprintf(fp,"%lf \t %lf \t" , x[i],tmp[0]);
		i++;
		fprintf(fp,"%lf \t %lf \n" , x[i],tmp[1]);
	}
	fclose(fp);
}
double fun(double x){
	return 1/(1+25*pow(x,2));
}
double dfun(double x){
	return 50*x / pow(1+25*pow(x,2),2);
}
	
