function [s,niter,E_out,time]=Smoothing_conTime(b,M,Mt,D,Dt,lambda,opts,m,n)
%  Created on Feb 1st, 2013
%  Author: Zhao Tan
%  1. This code implemented the analysis compressive senings using FISTA.
%  2. The convergence is depended on the rho and accuracy requirement.
%  3. The idea on continuation is used here to speedup the algorithm.
%  4. It has the accuracy as normal CVX method with no regularization and
%  also the same computing speed when problem dimension is small. r the
%  ratio to change the acc and rho for every iteration.
%  5. More optimization can be implemented.
% Inputs:
% -b  the measurement
% -M  the measuring operator or matrix
% -Mt adjoint operator of M
% -D  the sensing operator or matrix
% -Dt the ajoint operator of D
% -lambda tradeoff between l2 and l1 norm
% -opts the options for this algorithm
% -opts.C number of continuation
% -opts.r changing rate of continuation
% -opts.acc the accuracy we want to obtain
% -opts.rhof the final value of rho
% -opts.maxiter the maximum iteration number
% Outputs:
% -s the reconstructed signal
% -niter number of iteratios used to achieve the desired accuracy
time=[];
E_out=[];
optsin.print=opts.print;

if isa(M,'function_handle')
   x_ref=Mt(b);
else
   x_ref = M'*b;
end
X_iter=opts.xplug;

if isa(D,'function_handle')
    Dx_ref = D(x_ref);
else
    Dx_ref = D*x_ref;
end

mu0 =0.9*max(abs(Dx_ref));
Gamma = (opts.muf/mu0)^(1/opts.C);
optsin.mu= mu0;
Gammat= (opts.accf/0.1)^(1/opts.C);
optsin.acc = 0.1;

maxiter=opts.maxiter;%max iteration number overall
optsin.maxiter=opts.innerstep;%max iteration number for one continuation
optsin.x0=opts.x0;

niter=0;

for c=1:opts.C
    if c==1
         [X_iter,iter_num,E,time_out]=SmoothingMethodTime(b,M,Mt,D,Dt,lambda,X_iter,optsin,m,n);
    else
         [X_iter,iter_num,E,time_out]=SmoothingMethodTime(b,M,Mt,D,Dt,lambda,X_iter,optsin,m,n);
    end
    E_out=[E_out,E];
    niter=niter+iter_num;
    optsin.acc=optsin.acc*Gammat;
    optsin.mu=optsin.mu*Gamma;
    time=[time,time_out];
    if niter>maxiter
        break;
    end
end
s=X_iter;