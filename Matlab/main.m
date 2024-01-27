clear all;
close all;
clc;

tol = 1.d-4;
maxit = 100000000;
lb=[-10;-10;-10;-10];
ub=[ 10; 10; 10; 10];

% random sequence initialization
rng(137885)

x = (lb+ub)./2;

options = struct('tol',1.0d-6,'maxiter',1000,'maxfeval',1000,'verbose',1);

[pout,fout,nf,tcpu]=sdbox(x,lb,ub,options,@powell);
