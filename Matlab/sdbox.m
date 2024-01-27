function [xout,fout,nf,tcpu]=sdbox(x,lb,ub,options,f)
%function [xout,fout,nf,tcpu]=sdbox(x,lb,ub,options,f,model,y,weights,varargin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% sdbox is a local minimization algorithm
%INPUT 
% f = function to be minimized
% lb: lower bound, it is a vector Nx1(the algorithm needs a lower and upper bound, they cannot be infinity) 
% ub: upper bound, it is a vector Nx1
% options is a struct to set:
%   tol = tolerance in the stopping criterion. The default is 1.0d-4
%   maxiter = maximum allowed number of iterations. The default is 10000.
%   verbose = it can be 0 or 1, if it is 1 the default is 1. If it is zero the algorithm does not print 
%   anything during the execution. If verbose=1 then tha algorithm prints at each iteration the iteration 
%   number, the current number of function evaluations, the current best function value and the current 
%   maximum stepsize
% If options is empty, then the default value are used.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% OUTPUT:
% xout best point found
% fout objective function value in xout
% nf number of function evaluations required to find xout
% tcpu cpu time 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ((isempty(options))),
    
    tol=1.d-4;
    maxiter=3000;
    verbose =1; 
    maxfeval=3000;
else
    
    tol=options.tol;
    maxiter=options.maxiter;
    verbose =options.verbose;
    maxfeval=options.maxfeval;
end
tc=cputime;

    n=length(x);
	nf = 0;
	alpha_tilde = ones(n,1);
	alpha0 =1;
    d = ones(n,1);
	gamma = 1.d-6;	
	theta = 0.5;
    fu=f(x);
	nf = nf+1;
	f_ext = fu;
    
    
    if (verbose ==1) 
        disp(['iteration     ', 'nf    ' , 'f    ', 'alpha_max']);
    end

	for k=1:maxiter
		z = x;
		if((max(alpha_tilde) <= tol) | (nf > maxfeval)) 
           xout=x;
           fout=fu;
           break;
        end
		for i=1:n
			flag_dir = true;
			fu = f_ext;
			if(d(i) > 0.d0) 
				alpha_max = ub(i)-x(i);
			else
				alpha_max = x(i)-lb(i);
			end 
			alpha0 = min(alpha_tilde(i),alpha_max);
			z(i) = x(i)+alpha0*d(i);
            
			fz=f(z);
			nf = nf+1;
			if((alpha0 > 0.d0)&(fz <= fu-gamma*alpha0^2)) 
                %keyboard;
                [alpha0,f_ext,nf]=linesearch(i,x,z,d,alpha_max,gamma,fu,fz,nf, alpha0,f);
                alpha_tilde(i) = alpha0;
				flag_dir = false;
			end 
			if(flag_dir) 
				if(d(i) > 0.d0) 
                    alpha_max = x(i)-lb(i);
				else
					alpha_max = ub(i)-x(i);
				end 
				alpha0 = min(alpha_tilde(i),alpha_max);
				z(i) = x(i)-alpha0*d(i);
				fz=f(z);
				nf = nf+1;
				if((alpha0 > 0.d0)&(fz <= fu-gamma*alpha0^2)) 
					d(i) = -d(i);
                    [alpha0,f_ext,nf]=linesearch(i,x,z,d,alpha_max,gamma,fu,fz,nf, alpha0,f);
   					alpha_tilde(i) = alpha0;
				else
					alpha0 = 0.d0;
					alpha_tilde(i) = theta*min(alpha_tilde(i),alpha_max);
				end 
			end 
			x(i) = x(i) + alpha0*d(i);
			z(i) = x(i);
			fu = f_ext;
        end 
        if (verbose ==1) 
            %disp([num2str(k),'  ', num2str(nf), '  ' , num2str(fu), '  ' , num2str(max(alpha_tilde))]);
            fprintf(' %5d  %6d  %13.6e  %13.6e\n',k,nf,fu,max(alpha_tilde));
        end
    end
    tcpu = cputime-tc;

end

function [alpha0,f_ext,nf]=linesearch(i,x,z,d,alpha_max,gamma,fu,fz,nf, alpha0,f)
	delta = 0.5d0;
	f_ext = fz;
	alpha_tilde = min(alpha_max,alpha0/delta);
	z(i) = x(i)+alpha_tilde*d(i); 
    fz=f(z);
	nf = nf+1;
	while((alpha_tilde < alpha_max )&(fz<=fu-gamma*alpha_tilde^2))
		f_ext = fz;
		alpha0 = alpha_tilde;
		alpha_tilde = min(alpha_max,alpha0/delta);
		z(i) = x(i)+alpha_tilde*d(i);
		fz=f(z);
		nf = nf+1;
    end 
    return
end    

