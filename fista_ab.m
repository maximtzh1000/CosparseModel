function [X_out,niter]=fista_ab(b,M,Mt,D,Dt,lambda,pars,m,n)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This code is based on Amir Beck's 2009 paper about total variation
% Then it is revised on the author to fit the analysis compressive sensing
% case
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

MAXITER=pars.MAXITER;
denoiseiter=1000;
parsin.epsilon=10^-6;
parsin.MAXITER=denoiseiter;
parsin.print=pars.print;
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
    if isa(M,'function_handle')
        Y=Y-2/L*Mt(M(Y)-b);
    else
        Y=Y-2/L*Mt*(M*Y-b);
    end
     
    %invoking the denoising procedure 
    if (i==1)
        [X_iter,iter,P]=fista_ab_init(Y,2*lambda/L,D,Dt,[],m,parsin);
    else
        [X_iter,iter,P]=fista_ab_init(Y,2*lambda/L,D,Dt,P,m,parsin);
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
    if isa(D,'function_handle')
        value=value_L2+lambda*norm(D(X_iter),1);
    else
        value=value_L2+lambda*norm(D*X_iter,1);
    end
    if(abs(value_old-value)<acc)
       break;
    end
    
    %updating t and Y
    t_new=(1+sqrt(1+4*t_old^2))/2;
    Y=X_iter+(t_old-1)/t_new*(X_iter-X_old);
    value_old=value;
end

X_out=X_iter;