#include <stdio.h>
#include <stdlib.h>
#include <math.h>

double extrapolar(int h,int j,double f0h,double f0h2,double f0h4,double f0h8);
double dif(int centro,double h, double *f);
double integ(double h,double *f,int ai, int bi,double a, double b);

int main(int argc, char *argv[]) {
	double f[17] = {0.00,3.11,5.35,6.79,7.60,7.91,8.01,7.98,7.93,7.85,7.74,7.60,7.43,7.22,6.97,6.72,6.45};
	double f0h = dif(8,24,f);
	double f0h2= dif(8,12,f);
	double f0h4 = dif(8,6,f);
	double f0h8 = dif(8,3,f); // dif central para h / 8
	printf("diferencia central en 3 etapas de extrapolacion %lf \n" ,extrapolar(1,3,f0h,f0h2,f0h4,f0h8));
	f0h = integ(24,f,0,16,0,48);
	f0h2 = integ(12,f,0,16,0,48);
	f0h4 = integ(6,f,0,16,0,48);
	f0h8 = integ(3,f,0,16,0,48);
	printf("intregral por trapecios con 3 etapas de extrapolacion %lf \n" ,extrapolar(1,3,f0h,f0h2,f0h4,f0h8));
	
	return 0;
}
/* 3-stage Richardson-Romberg extrapolation */
double extrapolar(int h,int j,double f0h,double f0h2,double f0h4,double f0h8){
	if(j != 1){
		return (pow(4,j)*extrapolar(h+1,j-1,f0h,f0h2,f0h4,f0h8) - extrapolar(h,j-1,f0h,f0h2,f0h4,f0h8))/ (pow(4,j)-1);
    }
    else{
    	if (h == 1){
    		return (4*f0h2 - f0h)/ 3;
    	}
		if (h == 2){
    		return (4*f0h4 - f0h2)/3;	
		}
		if (h = 3) {
			return (4*f0h8 - f0h4)/3;
		}
	    
	}
}
/* central diference */
double dif(int centro,double h, double *f){
	int paso =  h / 3; // ha de ser entero (la tabla ha de cumplir el formato del paso a dividir por potencias de 2) (va de 3 en 3)
	return (f[centro + paso] - f[centro - paso]) / (2*h);
}
/* integral by trapezoids bi = b /3 */
double integ(double h,double *f,int ai, int bi,double a, double b){
	int paso = h / 3; //entero si coincide el formato de dividir por 2
	int xi = ((b-a)/h) - 1;
	int i;
	double sum = 0;
	for(i = 1; i<= xi; i++){
		sum += f[ai+i*paso];
	}
	return h*(0.5* f[ai] + sum + 0.5*f[bi]);
	
}
