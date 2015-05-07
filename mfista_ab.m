function [s,niter,f_out,E_out]=mfista_ab(b,M,Mt,D,Dt,lambda,pars,m,n)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This code is based on Amir Beck's 2009 paper about total variation
% Then it is revised on the author to fit the analysis compressive sensing
% case. We use a monotone version here
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
% niter.......................... Number of iterations 

%to do list: it's not monotone when A is not identitiy matrix
global normM;

MAXITER=pars.MAXITER;
maxiter=pars.maxiter;
parsin.epsilon=10^-5;
parsin.MAXITER=pars.denoiseiter;
parsin.print=pars.print;
parsin.x0=pars.x0;
acc=pars.acc;
f_out=[];
E_out=[];

%The Lipschitz constant of the gradient of ||M(X)-b||^2
    
L=4*normM;

% initialization
X_iter=pars.xplug;
Y=X_iter;
t_new=1;

niter=0;
value_old=inf;
OK=0;

for i=1:MAXITER
    % store the old value of the iterate and the t-constant
    X_old=X_iter;
    t_old=t_new;
    % gradient step
    if isa(M,'function_handle')
        Y=Y-2/L*Mt(M(Y)-b);
    else
        Y=Y-2/L*Mt*(M*Y-b);
    end
     
    %invoking the denoising procedure 
    if (i==1)
        [Z_iter,iter,P,f_val,E]=fista_ab_init(b,Y,lambda,2*lambda/L,M,Mt,D,Dt,[],m,parsin);
    else
        [Z_iter,iter,P,f_val,E]=fista_ab_init(b,Y,lambda,2*lambda/L,M,Mt,D,Dt,P,m,parsin);
    end
    E_out=[E_out,E];
    niter=iter+niter;
    f_out=[f_out;f_val];
    if niter>maxiter
        break;
    end
    % Compute the total variation and the function value and store it in
    % the function values vector fun_all if exists.
    if isa(M,'function_handle')
         value_L2=(M(Z_iter)-b)'*(M(Z_iter)-b);
    else
         value_L2=(M*Z_iter-b)'*(M*Z_iter-b);
    end
    if isa(D,'function_handle')
        value=value_L2+lambda*norm(D(Z_iter),1);
    else
        value=value_L2+lambda*norm(D*Z_iter,1);
    end
    if value<value_old
        X_iter=Z_iter;
        %value=value_old;
    else
        X_iter=X_old;
    end
    if(abs(value_old-value)<acc&&OK==5)
       break;
    elseif(abs(value_old-value)<acc)
        OK=OK+1;
    else
       OK=0;
    end   
    %updating t and Y
    t_new=(1+sqrt(1+4*t_old^2))/2;
    Y=X_iter+(t_old/t_new)*(Z_iter-X_iter)+(t_old-1)/t_new*(X_iter-X_old);
    value_old=min(value,value_old);
end

s=X_iter;