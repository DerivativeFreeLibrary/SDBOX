#include <stdlib.h>
#include <stdio.h>
#include <math.h>

// functions declarations

void setdim(int *n);
void setbounds(int n, double* lb, double* ub);
void startp(int n, double* xin);
double funct(int n, double* x);

void sdbox(int n, 
	   double *lb, 
	   double *ub, 
	   double *xin,
	   double funct(int, double*),
	   double *xott, 
	   double *fbest, 
	   int maxiter, 
	   double eps);

void linesearch(int i, 
		int n, 
		double *x, 
		double *z, 
		double *d, 
		double alpha_max,
		double *alpha, 
                double funct(int, double*),
		double gamma, 
		double f, 
		double fz, 
		double *f_ext, 
		int *nf);

int arresto(int n, 
	    double *alpha_tilde, 
	    double tol, 
	    int kfstop, 
	    double *fstop, 
	    double fstoptol);


