#include<stdio.h>
#include<stdlib.h>
int main(void){
	return 0;
}
double normavect(int n,double *b){
	double tmp=0;
	int i;
	for(i=0;i<=n-1;i++){
		if(abs(b[i]) > tmp){
			tmp = abs(b[i]);
		}
	}
	return tmp;
}
double normamatrix(int n,double **a){
	double sum=0,tmp=0;
	int i,j;
	for(i=0;i<=n-1;i++){
		for(j=0;j<=n-1;j++){
			sum+=abs(a[i][j]);
		}
		if (sum>tmp){ tmp =sum; }
		sum=0;
	}
	return tmp;
} 
