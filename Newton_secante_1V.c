#include <stdio.h>
#include <stdlib.h>
#include <math.h>
/* run this program using the console pauser or add your own getch, system("pause") or input loop */
double fun(double x);
double der(double x);
int secant(double x0,double x1,double prec,int it);
int newton(double xn,double prec, int it);
//precision n %.Nlf
int main(int argc, char *argv[]) {
	double a,b,x0,x1;
	int m,i;
	printf("Entra los extremos del intervalo  (a< b) \n");
	scanf("%lf %lf",&a,&b);
	printf("Entra numero de abcisas \n");
	scanf("%d",&m);
	for(i = 0; i<= m ; i ++){
		double abcisa = a+i*((b-a)/m);
		printf("z, f(z): %lf %lf \n",abcisa,fun(abcisa));
	}
	printf("Entra las 2 mejores aproximaciones de 0  \n");
	scanf(" %lf %lf",&x0,&x1);
	secant(x0,x1,10e-4,1);
	newton((x0-x1)/2,10e-14,1);
	return 0;
}
double fun(double x){
	return 8*pow(x,3) - 44*pow(x,2) + 70*x -25;
}
double der(double x){
	return 24*pow(x,2)-88*x +70;
}
int secant(double x0, double x1,double prec, int it){
	double xn;
	if((fabs(x1 - x0)< prec )&& (fabs(fun(x1)) < prec)) {
		printf("cero encontrado por metodo secante: %lf \n ",x1);
		return 0;
	}
	else{
		if(it== 100){
			return -1;
		}
		else{
			xn = x1 - fun(x1)*((x1-x0)/(fun(x1)- fun(x0)));
			return secant(x1,xn,prec,it+1);
		}
	}
}
int newton(double xn,double prec, int it){
	if(fabs(fun(xn)) < prec) {
		printf("cero encontrado por metodo newton/raphson: %lf \n",xn);
		return 0;
	}
	else{
	
		if(it== 10){
			return -1;
		}
		else{
			xn = xn - (fun(xn)/der(xn));
			return newton(xn,prec,it+1);
		}
	}
}

