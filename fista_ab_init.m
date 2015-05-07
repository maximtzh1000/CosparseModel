function [X_den,iter,P,f_out,E_out]=fista_ab_init(b0,b,lambda0,lambda,M,Mt,D,Dt,P_init,m,pars)
%This function implements the FISTA method for TV denoising problems. This function is used inside the deblurring
% procedure since it uses a warm-start strategy
%
% INPUT
% b0 .............................measurements from outter iteration
% b ..............................an observed noisy image.
% lambda ........................ parameter
% pars.................................parameters structure
% pars.MAXITER ..................... maximum number of iterations
%                                                      (Default=100)
% pars.epsilon ..................... tolerance for relative error used in
%                                                       the stopping criteria (Default=1e-4)
% pars.print ..........................  1 if a report on the iterations is
%                                                       given, 0 if the  report is silenced
% pars.tv .................................. type of total variation
%                                                      penatly.  'iso' for isotropic (default)
%                                                      and 'l1' for nonisotropic
%  
% OUTPUT
% X_den ........................... The solution of the problem 
%                                            min{||X-Xobs||^2+2*lambda*TV(X)}
% iter .............................  Number of iterations required to get
%                                            an optimal solution (up to a tolerance)
% fun_all ......................   An array containing all the function
%                                             values obtained during the
%                                             iterations

% Assigning parameres according to pars and/or default values
global normD;

flag=exist('pars');
f_out=[];
if (flag&isfield(pars,'MAXITER'))
    MAXITER=pars.MAXITER;
else
    MAXITER=100;
end
if (flag&isfield(pars,'epsilon'))
    epsilon=pars.epsilon;
else
    epsilon=1e-4;
end
if(flag&isfield(pars,'print'))
    prnt=pars.print;
else
    prnt=1;
end

clear P
clear R
if(isempty(P_init))
    P=zeros(m,1);
    R=zeros(m,1);
else
    P=P_init;
    R=P_init;
end

n=length(b);
tk=1;
tkp1=1;
count=0;
i=0;

G=zeros(n,1);
fval=inf;

L=2*lambda^2*normD;
E_out=[];

while((i<MAXITER)&(count<5))
    fold=fval;
    %%%%%%%%%
    % updating the iteration counter
    i=i+1;
    %%%%%%%%%
    % Storing the old value of the current solution
    Gold=G;
    %%%%%%%%%%
    %Computing the gradient of the objective function
    Pold=P;
    tk=tkp1;
    if isa(D,'function_handle')
        G=b-lambda*Dt(R);
        Q=-2*lambda*D(G);
    else
        G=b-lambda*Dt*R;
        Q=-2*lambda*D*G;
    end

    %%%%%%%%%%
    % Taking a step towards minus of the gradient
    P=R-Q/L;
    %%%%%%%%%%
    % Peforming the projection step
    P=P./(max(abs(P),1));
           
    %%%%%%%%%%
    %Updating R and t
    tkp1=(1+sqrt(1+4*tk^2))/2;
    
    R=P+(tk-1)/(tkp1)*(P-Pold);
    
    re=norm(G-Gold)/norm(G);
    if (re<epsilon)
        count=count+1;
    else
        count=0;
    end
    if isa(D,'function_handle')
        C=b-lambda*Dt(P);
        f_primal=0.5*(M(G)-b0)'*(M(G)-b0)+lambda0*norm(D(G),1);
    else
        C=b-lambda*Dt*P;
        f_primal=0.5*(M*G-b0)'*(M*G-b0)+lambda0*norm(D*G,1);
    end
    fval=norm(C)^2;
    f_out=[f_out;f_primal];
    e=norm(G-pars.x0)/norm(pars.x0);
    E_out=[E_out,e];
    if(prnt)
        fprintf('iter= %5d value = %10.10f %10.10f',i,f_primal,norm(G-Gold)/norm(G));
        if (fval>fold)
            fprintf('  *\n');
        else
            fprintf('   \n');
        end
    end
end
X_den=G;
iter=i;