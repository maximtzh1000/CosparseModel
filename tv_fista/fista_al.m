function [X_out,niter]=fista_al(b,M,Mt,D,D_temp,lambda,pars)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%This function implements FISTA for solving the linear inverse problem with 
% the total variation regularizer and either reflexive or periodic boundary
% conditions
%
% Based on the paper
% Amir Beck and Marc Teboulle, "Fast Gradient-Based Algorithms for Constrained
% Total Variation Image Denoising and Deblurring Problems"
% -----------------------------------------------------------------------
% Copyright (2008): Amir Beck and Marc Teboulle
% 
% FISTA is distributed under the terms of 
% the GNU General Public License 2.0.
% 
% Permission to use, copy, modify, and distribute this software for
% any purpose without fee is hereby granted, provided that this entire
% notice is included in all copies of any software which is or includes
% a copy or modification of this software and in all copies of the
% supporting documentation for such software.
% This software is being provided "as is", without any express or
% implied warranty.  In particular, the authors do not make any
% representation or warranty of any kind concerning the merchantability
% of this software or its fitness for any particular purpose."
% ----------------------------------------------------------------------
% INPUT
%
% b............................. The observed image vector
% lambda ...................... Regularization parameter
% M............................ measurement operator
% Mt .......................... adjoint operator
% D ........................... sensing operator
% Dt .......................... adjoint sensing operator
% pars.................................Parameters structure
% pars.MAXITER ..................... maximum number of iterations
%                                                      (Default=100)
% pars.denoiseiter .......... number of iterations of the denoising inner
%                                                      problem (default=10)
% pars.acc stopping criteria
% OUTPUT
% 
% X_out ......................... Solution of the problem

[m,n]=size(D_temp);
MAXITER=pars.MAXITER;
denoiseiter=1000;
parsin.tv='l1';
parsin.epsilon=10^-6;
parsin.MAXITER=denoiseiter;
parsin.print=0;
acc=pars.acc;

%The Lipschitz constant of the gradient of ||M(X)-b||^2
if isa(M,'function_handle')
   MtM=@(x) Mt(M(x));
   normM=my_normest(MtM,MtM,n);
else
   normM=normest(M'*M);
end
    
L=4*normM;

% initialization
X_iter=zeros(n,1);
Y=X_iter;
t_new=1;

niter=0;
value_old=inf;

for i=1:MAXITER
    % store the old value of the iterate and the t-constant
    X_old=X_iter;
    t_old=t_new;
    % gradient step
    Y=Y-2/L*Mt(M(Y)-b);
     
    %invoking the denoising procedure 
    if (i==1)
        [X_iter,iter,fun_denoise,P]=denoise_bound_init(reshape(Y,sqrt(n),sqrt(n)),2*lambda/L,-inf,inf,[],parsin);
    else
        [X_iter,iter,fun_denoise,P]=denoise_bound_init(reshape(Y,sqrt(n),sqrt(n)),2*lambda/L,-inf,inf,P,parsin);
    end
    X_iter=X_iter(:);
    niter=iter+niter;
    % Compute the total variation and the function value and store it in
    % the function values vector fun_all if exists.
    if isa(M,'function_handle')
         value_L2=(M(X_iter)-b)'*(M(X_iter)-b);
    else
         value_L2=(M*X_iter-b)'*(M*X_iter-b);
    end
    value=value_L2+lambda*norm(D(X_iter),1);
    if(abs(value_old-value)<acc)
       break;
    end
    
    %updating t and Y
    t_new=(1+sqrt(1+4*t_old^2))/2;
    Y=X_iter+(t_old-1)/t_new*(X_iter-X_old);
    value_old=value;
end

X_out=X_iter;