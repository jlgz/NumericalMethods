#include<stdio.h>
#include<stdlib.h>
#include<math.h>
double d_f_z(double z);
double newton_1var(double zo,double tol);
double f_z(double z);
int newton_2var(double *xy,double tol, double *g_f_1, double *g_f_2,double xo,double yo,double ho,int j);
double f(double x, double y);
double g(double x, double y,double xo,double yo,double ho);
int grad_f1( double x, double y,double *g_f_1);
int grad_f2(double x, double y, double xo,double yo,double *g_f_2);
int inv_dif(double *g_f_1, double *g_f_2,double **inv);


int main(void) {
	FILE *fp;
	double zo = 0,xo,yo;
	double *xy = (double *)malloc(sizeof(double)*2);
	double *g_f_1 = (double *)malloc(sizeof(double)*2);
	double *g_f_2 = (double *)malloc(sizeof(double)*2);
	int i;
	double tmp[2],tmp2[2];
 	fp= fopen("fichero.txt", "w");
	zo = newton_1var(zo,10.e-15); // corte con el exe y
	zo = sqrt(zo); //y = z² 
	xy[0] = 0;
	xy[1] = zo; 
	tmp2[0] = 0; //tangente anterior
	tmp2[1] = 0;
	for(i=0;i<= 100000; i++){
		xo = xy[0];
		yo = xy[1];
		grad_f1(xy[0],xy[1],g_f_1);
		tmp[0] = -g_f_1[1] /sqrt(pow(g_f_1[0],2) + pow(g_f_1[1],2)); //calculo del vector tangente unitario
		tmp[1] = g_f_1[0] / sqrt(pow(g_f_1[0],2) + pow(g_f_1[1],2));
		if( (tmp[0]*tmp2[0] + tmp[1]*tmp2[1]) < 0){ //se comprueba el tangente anterior
			tmp[0] *= -1;
			tmp[1] *= -1;
		}
		tmp[0] *= 0.01;
		tmp[1] *= 0.01;
		/*printf("i %d :( %lf , %lf ) \n",i,xy[0],xy[1]);
		printf( "tangente: %lf %lf \n",tmp[0],tmp[1]);
		printf("f = %.15lf \n ",f(xy[0],xy[1]));*/
		xy[0] += tmp[0]; //primera aproximacion
		xy[1] += tmp[1];
		newton_2var(xy,10.e-15,g_f_1,g_f_2,xo,yo,0.01, i); // newton..
		tmp2[0] = tmp[0]; //se guarda el tangentente
		tmp2[1] = tmp[1];
		fprintf(fp,"%.15lf \t %.15lf \n" , xy[0],xy[1]);
		if( pow(xy[0],2) + pow(xy[1]- zo,2) -pow(0.01,2) < 10.e-15 && i >2) {break;} // distancia 0.01 del primero
		
	}
	free(xy);
	free(g_f_1);
	free(g_f_2);
	fclose(fp);
	return 0;
}
/*
 f(x,y)
*/
double f(double x, double y) {
	return (3*pow(x, 2) +3*pow(y,2) -1)*(pow(x,2)+pow(y,2) -5)*(pow(x,2)  + pow(y,2) - 3*x +2)  +1;
}
/*
 g(x,y) = (x-xo)² + (y-yo)² - h²
*/
double g(double x, double y,double xo,double yo,double ho){
	return pow(x - xo,2) + pow(y - yo,2) - pow(ho,2);
}
/*
f(0,z) z= y²
*/
double d_f_z(double z){
	return 9*pow(z,2)-20*z -27;
}
/*
 df(0,z) z=y²
*/
double f_z(double z){
	return 3*pow(z,3) -10 *pow(z,2) -27*z + 11;
}
/*
 gradiente de f(x,y)
*/
int grad_f1( double x, double y,double *g_f_1) {
	g_f_1[0] = (6*x*(pow(x,2) + pow(y,2) -5) + 2*x*(3*pow(x,2) + 3*pow(y,2) -1))*(pow(x,2) + pow(y,2) - 3 * x + 2)+(3*pow(x,2) + 3*pow(y,2) -1)*(pow(x,2) + pow(y,2)  -5)*(2*x - 3);
	g_f_1[1] =  (6*y*(pow(x,2) + pow(y,2) -5) + 2*y*(3*pow(x,2) + 3*pow(y,2) -1))*(pow(x,2) + pow(y,2) - 3 * x + 2)+(3*pow(x,2) + 3*pow(y,2) -1)*(pow(x,2) + pow(y,2)  -5)*2*y;
	return 0;  
}
/*
 gradiente de g(x,y) = (x- xo)² + (y- yo)² - h²
*/
int grad_f2(double x, double y, double xo,double yo,double *g_f_2){
	g_f_2[0] = 2*(x-xo);
	g_f_2[1] = 2*(y-yo);
	return 0;
}
/*
inversa de F(x,y)  = (f(x,y),g(x,y))
*/
int inv_dif(double *g_f_1, double *g_f_2, double **inv) {
	double delta = g_f_1[0] * g_f_2[1] - g_f_2[0]*g_f_1[1];
	inv[0][0] = g_f_2[1] / delta;
	inv[0][1] = -g_f_1[1] / delta;
	inv[1][0] = -g_f_2[0] / delta;
	inv[1][1] = g_f_1[0] / delta;
	return 0;
}
/*
newton de f(0,z) = 0 z= y²
*/
double newton_1var(double zo,double tol){
	int i;
	for(i = 0; i <= 100; i++) {
		zo = zo - f_z(zo) /d_f_z(zo);
		if( fabs(f_z(zo)) < tol) {
			break;
		}
	}
	return zo;
}
/*
newton de F(x,y) = ( 0, 0 )
		g(x,y,h) = ( 0, 0 )
*/
int newton_2var(double *xy,double tol, double *g_f_1, double *g_f_2,double xo,double yo,double ho,int j){
	double f_,g_;
	int i;
	double **inv = (double **)malloc(sizeof(double *) * 2);
	inv[0] = (double *) malloc( sizeof(double) *2);
	inv[1] = (double *) malloc( sizeof(double) *2);
	for(i = 0; i<= 1000;i++) {
		grad_f1(xy[0],xy[1],g_f_1); //gradientes
		grad_f2(xy[0],xy[1],xo,yo,g_f_2);
		//if (j == 21){ 
			//printf("xy: %lf %lf \n",xy[0],xy[1]);
			//printf("gradiente1 %lf , %lf \n",g_f_1[0],g_f_1[1]);
			//printf("gradiente2 %lf , %lf \n",g_f_2[0],g_f_2[1]);
		//}		
		inv_dif(g_f_1,g_f_2,inv); //inversa
		f_ = f(xy[0],xy[1]);
		g_ = g(xy[0],xy[1],xo,yo,ho);
		if (fabs(f_) < tol && fabs(g_) < tol){ break;} 
		xy[0] = xy[0] - inv[0][0]*f_ - inv[0][1]*g_;
		xy[1] = xy[1] - inv[1][0]*f_ - inv[1][1]*g_;	
	}
	free(inv[0]);
	free(inv[1]);
	free(inv);
	return 0;	
}
