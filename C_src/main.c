#include "headers.h"

int main(int argc, char **argv){
	int i, n, maxiter;
	double *lb, *ub, *xott, *xin;
	double fbest, eps;

	setdim(&n);

	lb    = (double *)malloc(n*sizeof(double));
	ub    = (double *)malloc(n*sizeof(double));
	xott  = (double *)malloc(n*sizeof(double));
	xin   = (double *)malloc(n*sizeof(double));

	setbounds(n,lb,ub);

	startp(n,xin);

	maxiter = 1000*n;
	eps     = 1.e-6;
	sdbox(n,lb,ub,xin,funct,xott,&fbest,maxiter,eps);

	printf("fbest=%f\n",fbest);
	for(i=0;i<n;i++) {
		printf("xbest[%d]=%f\n",i,xott[i]);
	}

	return(0);
}


