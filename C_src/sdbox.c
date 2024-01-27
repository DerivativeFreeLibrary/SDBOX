/*============================================================================================
!    SDBOX - C implementation of a Derivative-Free algorithm for bound 
!    constrained optimization problems 
!    Copyright (C) 2011  G.Liuzzi, A. Risi
!
!    This program is free software: you can redistribute it and/or modify
!    it under the terms of the GNU General Public License as published by
!    the Free Software Foundation, either version 3 of the License, or
!    (at your option) any later version.
!
!    This program is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!    GNU General Public License for more details.
!
!    You should have received a copy of the GNU General Public License
!    along with this program.  If not, see <http://www.gnu.org/licenses/>.
!
!    S. Lucidi, M. Sciandrone. A Derivative-Free Algorithm for Bound Constrained Optimization, 
!    Computational Optimization and Applications, 21(2): 119-142 (2002)
!    DOI: 10.1023/A:1013735414984
!
==============================================================================================
*/

#include "headers.h"

void sdbox(int n, double *lb, double *ub, double *xin, 
	   double funct(int, double*),
	   double *xott, double *fbest, 
		   int maxiter, double tol){
	int i, k, ifstop, kfstop;
	int nf, fall;
	int flag_dir;
	double *x, *d, *z, *alpha_tilde, *fstop;
	double gamma, theta, maxalpha;
	double f, fz, f_ext, alpha_max, alpha, fstoptol;
    	
	x			 = (double *)malloc(n  *sizeof(double));
	d			 = (double *)malloc(n  *sizeof(double));
	z			 = (double *)malloc(n  *sizeof(double));
	alpha_tilde  = (double *)malloc(n  *sizeof(double));
	fstop		 = (double *)malloc((n+1)*sizeof(double));

	nf = 0;
	for(i=0;i<n;i++){
		alpha_tilde[i] = 1.0;
		d[i]		   = 1.0;
		x[i]		   = xin[i];
	}
	gamma    = 1.e-6;	
	theta    = 0.5;
	fstoptol = 0.01*tol;
	f        = funct(n,x);
	nf       = nf+1;
	f_ext    = f;

	
	// CICLO PRINCIPALE
	k             = 1;
	ifstop        = 0;
	kfstop        = 1;
	fstop[ifstop] = f;
	while (k <=  maxiter){
		maxalpha = alpha_tilde[0];
		for(i=1;i<n;i++) if(maxalpha<alpha_tilde[i]) maxalpha = alpha_tilde[i];

		printf("obj.f: %f max_alpha: %f k: %d nf: %d\n", f, maxalpha,k,nf);

		for(i=0;i<n;i++) z[i] = x[i];
        if(arresto(n,alpha_tilde,tol,kfstop,fstop,fstoptol)>0) break;
        fall = 0;
		for(i=0;i<n;i++){
			flag_dir = 1;
			//--------------
			//STEP 1--> di
			//--------------
			f = f_ext;
			if(d[i] > 0.0){ 
				alpha_max = ub[i]-x[i];
			}else{
				alpha_max = x[i]-lb[i];
			}
			alpha = (alpha_tilde[i] < alpha_max ? alpha_tilde[i] : alpha_max);
			z[i] = x[i]+alpha*d[i];
			fz = funct(n,z);
			nf = nf+1;
			if((alpha > 0.0) && (fz <= f-gamma*pow(alpha,2.0))){
				//--------------
				//STEP 3(1)
				//--------------
			  linesearch(i,n,x,z,d,alpha_max,&alpha,funct,gamma,f,fz,&f_ext,&nf);
				alpha_tilde[i] = (alpha < 1.0 ? alpha : 1.0);
				flag_dir = 0;
			}
			//--------------
			//STEP 2--> -di
			//--------------
			if(flag_dir){
				if(d[i] > 0.0){
					alpha_max = x[i]-lb[i];
				}else{
					alpha_max = ub[i]-x[i];
				}
				alpha = (alpha_tilde[i] < alpha_max ? alpha_tilde[i] : alpha_max);
				z[i] = x[i]-alpha*d[i];
				fz = funct(n,z);
				nf = nf+1;
				if((alpha > 0.0) && (fz <= f-gamma*pow(alpha,2.0))){
					//--------------
					//STEP 3(2)
					//--------------
					d[i] = -d[i];
					linesearch(i,n,x,z,d,alpha_max,&alpha,funct,gamma,f,fz,&f_ext,&nf);
					alpha_tilde[i] = (alpha < 1.0 ? alpha : 1.0);
				}else{
					alpha = 0.0;
					alpha_tilde[i] = 
						theta*(alpha_tilde[i] < alpha_max ? alpha_tilde[i] : alpha_max);
					}
			}
			//--------------
			//STEP 4
			//--------------
			x[i] = x[i] + alpha*d[i];
			z[i] = x[i];
			f = f_ext;
			if(alpha > 0.0){
				ifstop = ifstop + 1;
				if(ifstop <= n){
					if(kfstop < ifstop) kfstop = ifstop;
				}else{
					ifstop = 0;
				}
				fstop[ifstop] = f;
			}else{
				fall = fall + 1;
			}
		}
		k = k+1;
	}
	
	maxalpha = alpha_tilde[0];
	for(i=1;i<n;i++) if(maxalpha<alpha_tilde[i]) maxalpha = alpha_tilde[i];
	for(i=0;i<n;i++) xott[i] = x[i];
	*fbest = f;

	printf("obj.f: %f max_alpha: %f k: %d nf: %d\n", f, maxalpha,k,nf);

	return;
}

int arresto(int n, double *alpha_tilde, double tol, int kfstop, double *fstop, 
			double fstoptol){
	int i;
	double fm, sq, maxalpha;
	
	maxalpha = alpha_tilde[0];
	for(i=1;i<n;i++) maxalpha = ((maxalpha < alpha_tilde[i]) ? alpha_tilde[i] : maxalpha);

	if(maxalpha <= tol){
		return(1);
	}else{
		if(kfstop == n){
			fm = 0.0;
			for(i=0;i<=n;i++) fm = fm + fstop[i];
			fm = fm / (double)(n+1);
			sq = 0.0;
			for(i=0;i<=n;i++) sq = sq + pow((fstop[i] - fm),2.0);
			sq = sqrt(sq / (double)(n+1));
			if(sq <= fstoptol){
				return(1);
			}
		}
	}
	
	return(0);
}

void linesearch(int i, int n, double *x, double *z, double *d, double alpha_max,
		double *alpha, double funct(int, double*),
		double gamma, double f, double fz, double *f_ext, int *nf){
	double delta, alpha_tilde;
	
	//INIZIALIZZAZIONI
	delta = 0.5;
	*f_ext = fz;
	alpha_tilde = ( alpha_max < (*alpha)/delta ? alpha_max : (*alpha)/delta );
	z[i] = x[i]+alpha_tilde*d[i]; 
	fz = funct(n,z);
	*nf = *nf+1;
	while((alpha_tilde < alpha_max ) && (fz<=f-gamma*pow(alpha_tilde,2.0))){
		*f_ext = fz;
		*alpha = alpha_tilde;
		alpha_tilde = ( alpha_max < (*alpha)/delta ? alpha_max : (*alpha)/delta );
		z[i] = x[i]+alpha_tilde*d[i];
		fz = funct(n,z);
		*nf = *nf+1;
	}
	return ;
}
