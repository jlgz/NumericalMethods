#include<stdio.h>
#include<stdlib.h>
#include<math.h>

// funciones principales
int jacobi_14(double *x1,double *x2,int n,int itmax);
int gauss_seidel_14(double *x1,double *x2, int n,int itmax);
int sor_14(double *x1,double *x2,double w,int n,int itmax);
int steepest_descent_14(double *x1, double w,int n,int itmax);

//funciones auxiliares
double error_14(double *x,int n);
double alpha_15_14(double *xg,double *ptp, int n);
int gradiente_15_14(double *x1,double *xg, int n);
double error_14_2(double *x1, double *x2,int n);

int main() {
	//tamano de la matriz ejercicio 14
	int n = 1000000;
	//se reserva memoria para las soluciones de los metodos
	double *xj = (double*) malloc(n* sizeof(double ));
	double *xj2 = (double*) malloc(n* sizeof(double )); 
	double *xgs = (double *) malloc(n* sizeof(double));
	double *xgs2 = (double *) malloc(n* sizeof(double));	
	double *xsor = (double *) malloc(n* sizeof(double));
	double *xsor2 = (double *) malloc(n* sizeof(double));
	double *xsd = (double *) malloc(n* sizeof(double));

	//inicializaciones
	int i,it,itmin= 100;
	double wmin = 0.1;
	double w;
	for(i=0;i<n;i++){
		xsd[i] = 0;			
		xj[i] = 0;
		xgs[i]= 0;
		xsor[i] = 0;
	}
	//resultados metodos (numero de iteraciones)
	printf("finalizado jacobi en %d iteraciones \n",jacobi_14(xj,xj2,n,100));
	int itgs = gauss_seidel_14(xgs,xgs2,n,100);		
	printf("finalizado gauss seidel en %d iteraciones \n", itgs);
	printf("finalizado steepest descent en %d iteraciones \n", steepest_descent_14(xsd,1,n,100));
	printf(" calculando.. un buen w para SOR... \n");	
	//busqueda de un buen w en SOR
	for(w = 0.3; w< 2; w+= 0.4){ //w = 1.09 30 iteraciones
		it = sor_14(xsor,xsor2,w,n,100);
		if(it < itmin){
			itmin = it;
			wmin = w;
		}
		for(i = 0; i<n;i++){ xsor[i] = 0;}
	} 
	printf("finalizado SOR en %d iteraciones(%d iteraciones menos que gauss seidel) con w = %lf \n",itmin,itgs - itmin,wmin); 
	//busqueda de un buen w en steepest descent
	printf(" calculando.. un buen w para steepest descent... \n");
	itmin = 140;
	for(w = 0.5; w< 2; w+= 0.4){ 
		for(i = 0; i<n;i++){ xsd[i] = 0;}
		it = steepest_descent_14(xsd,w,n,100); // w = 0.9 36 iteraciones
		if(it < itmin){
			itmin = it;
			wmin = w;
		}
	}
	printf("finalizado steepest descent en %d iteraciones con w = %lf \n",itmin,wmin);
	//liberar memoria
	free(xj);
	free(xj2);
	free(xgs);
	free(xgs2);
	free(xsor);
	free(xsor2);
	free(xsd);
	return 0;
}
//soluciona el sistema del ejercicio 14 por gauss-seidel y retorna el numero de iteraciones.
int gauss_seidel_14(double *x1,double *x2,int n,int itmax){
	int i,it;
	double *tmp;
	for(it= 0;it<itmax;it++) {	
		x2[0] = 1/3.0 *(-x1[2] - x1[n-2]+ 1.0/n);
		x2[1] = 1/3.0 *(-x1[3] - x1[n-1]+ 2.0/n);
		for(i = 2;i<n-2; i++){
			x2[i] = 1/3.0 *(-x2[i-2] -x1[i+2] + ((double)(i+1))/n);
		}
		x2[n-2] = 1/3.0 *(-x2[0] -x2[n-4] + ((double)(n-1))/n);
		x2[n-1] = 1/3.0 *(-x2[1] -x2[n-3] + 1);
		if(((2.0/3)/(1 - 2.0/3))*error_14_2(x1,x2,n)< pow(10,-12)){break;} //condicion de parada ||Xk - Xk-1|| *||B||/(1-||B||)<10^-12    (||Bgs||inf<||BJ||inf=2/3)
		tmp = x2;
		x2 = x1;
		x1 = tmp;
	}
	return it;
}
// retorna la propagacion del error de la solucion x en el sistema del ejercicio 14 (||Ax -b||).
double error_14(double *x,int n){
	int i;
	double max;
	max = fabs(3*x[0] + x[2] + x[n-2] - 1.0/n);
	double tmp;
	tmp = fabs(3*x[1] + x[3] + x[n-1] - 2.0/n);	
	if(max < tmp) {
        	max = tmp;
    }	
	for(i =2; i< n-2;i++){
		tmp = fabs(x[i -2] + 3* x[i] + x[i+2] -((double) (i +1) )/n);
		if(max < tmp) {
			max = tmp; 
		}
	}
	for(i = n-2; i < n; i++) {
		tmp = fabs(x[i-n+2] + x[i-2] + 3*x[i] - ((double)(i+1)) /n);
		if(max < tmp) {
                        max = tmp; 
        	}
	}
	return max;
}
//norma infinito de la diferencia de 2 vectores iteraciones consecutivos
double error_14_2(double *x1,double *x2,int n) {
	double max, tmp;
	int i;
	max = fabs(x1[0] - x2[0]);
	for(i = 1; i<n; i++) {
		tmp = fabs(x1[i] - x2[i]);
		if(tmp> max){ max = tmp;}
	}
	return max;
}
//soluciona  el sistema del ejercicio 14 por el metodo de jacobi y retorna el numero de iteraciones.
int jacobi_14(double *x1,double *x2,int n, int itmax) {
	int i,it;
	double *tmp;
	for(it=0;it<itmax;it++){
		x2[0] = 1/3.0 *(-x1[2] - x1[n-2]+ 1.0/n);
		x2[1] = 1/3.0 *(-x1[3] - x1[n-1]+ 2.0/n);
		for(i = 2;i<n-2; i++){
			x2[i] = 1/3.0 *(-x1[i-2] -x1[i+2] + ((double)(i+1))/n);
		}
		x2[n-2] = 1/3.0 *(-x1[0] -x1[n-4] + ((double)(n-1))/n);
		x2[n-1] = 1/3.0 *(-x1[1] -x1[n-3] + 1);
		if(((2.0/3)/(1 - 2.0/3))*error_14_2(x1,x2,n)< pow(10,-12)){break;} // (||Bgs||inf<||BJ||inf=2/3)||Xk - Xk-1|| *||B||/(1-||B||)<10^-12
		tmp = x2;
		x2 = x1;
		x1 = tmp;
	}
	return it;
}
//soluciona el sistema del ejercicio 14 por el metodo SOR con parametro w y retorna el numero de iteraciones. La solucion queda guardada en x.
int sor_14(double *x1,double *x2,double w,int n,int itmax){
	int i,it;
	double *tmp;
	for(it = 0; it < itmax; it++){               // ||Bw|| < |1 - w| + w2/3 
		x2[0] = x1[0] + w/3.0 *(-3*x1[0] -x1[2] - x1[n-2]+ 1.0/n);
		x2[1] = x1[1] + w/3.0 *(-3*x1[1]-x1[3] - x1[n-1]+ 2.0/n);
		for(i = 2;i<n-2; i++){
			x2[i] = x1[i]+ w/3.0 *(-3*x1[i]-x2[i-2] -x1[i+2] + ((double)(i+1))/n);
		}
		x2[n-2] = x1[n-2] + w/3.0 *(-3*x1[n-2]-x2[0] -x2[n-4] + ((double)(n-1))/n);
		x2[n-1] = x1[n-1] + w/3.0 *(-3*x1[n-1] -x2[1] -x2[n-3] + 1);
		if (((2.0/3)/(1 - 2.0/3))*error_14_2(x1,x2,n)< pow(10,-12)){ break;} //en los casos en que SOR converge mas rapido que G-S misma condicion de parada
		tmp = x2;
		x2 = x1;
		x1 = tmp;
	}
	return it;
}
// soluciona el sistema del ejercicio 14 por el metodo steepest descent con factor de ralajacion w y retorna numero de iteraciones. Q(x) = 1/2 xt A x - bt x -> grad(Q(x)) = Ax -b
int steepest_descent_14(double *x1,double w, int n, int itmax){
	int i,it;
	double *ptp = malloc(sizeof(double));
	double alphk;
	double *xg = (double *) malloc(sizeof(double) *n);
	for(it=0;it<itmax;it++){
		gradiente_15_14(x1,xg,n);
		alphk = w*alpha_15_14(xg,ptp,n);
		for(i = 0; i<n; i++){
			x1[i] = x1[i] - alphk* xg[i];
		}
		if(pow(ptp[0],0.5) <pow(10,-12)){break;}
	}
	free(ptp);
	free(xg);
	return it;
}
//retorna el alpha sub-k del metodo steepest descent en el paso k+1. xg es el vector gradiente de Q(x) evaluado en la solucion k-esima.
double alpha_15_14(double *xg, double *ptp, int n) {
	int i;
	double sum = 0;
	double sum2 = 0;
	for(i = 0; i <2; i++) {
		sum += xg[i]*xg[i];
		sum2 += xg[i] * (3*xg[i] + xg[n-2 + i] + xg[i+2]);
	}
	for (i = 2; i<n-2; i++){
		sum += xg[i]*xg[i];
		sum2 += xg[i] * (3*xg[i] + xg[i-2] + xg[i+2]);
	}
	for(i = n-2; i < n; i++) {
		sum += xg[i]*xg[i];
		sum2 += xg[i] * (3*xg[i] + xg[i - n + 2] + xg[i-2]);
	}
	ptp[0] = sum;
	return sum / sum2;
}
// evalua el vector gradiente de Q(x) en x1. La evalucion queda guardada en xg.
int gradiente_15_14(double *x1,double *xg, int n){
	int i;
	xg[0] = 3 * x1[0] + x1[2] + x1[n-2] - 1.0/n;
	xg[1] = 3 * x1[1] + x1[3] + x1[n-1] - 2.0/n;
	for (i = 2; i< n-2; i++){
		xg[i] = x1[i-2] + 3*x1[i] + x1[i+2] - (double)(i+1)/n;
	}
	xg[n-2] = x1[0] + x1[n-4]+ 3*x1[n-2] - (double)(n-1)/n;
	xg[n-1] = x1[1] + x1[n-3]+ 3*x1[n-1] - 1;
	return 0;
}
