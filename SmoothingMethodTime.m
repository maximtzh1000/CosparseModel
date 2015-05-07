function [x,niter,E_out,time]=SmoothingMethodTime(b,M,Mt,D,Dt,lambda,x_init,opts,m,n)
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
 
maxiter=opts.maxiter;
acc=opts.acc;
mu=opts.mu;
print=opts.print;
fmean = realmin/10;
 
if(isempty(x_init))
    if isa(M,'function_handle')
       x_old=Mt(b);
    else
        x_old=M'*b;
    end
else
    x_old=x_init;
end
 
v=x_old;
niter=0;
f_out=[];
E_out=[];
 
if isa(M,'function_handle')
    Mtb=Mt(b);
    MtM=@(x) Mt(M(x));
else
    MtM=M'*M;
    Mtb=M'*b;
end
 
if isa(D,'function_handle')
    DtD=@(x) Dt(D(x));
else
    DtD=D'*D;
end
 
t_old=1;
if (isa(D,'function_handle')==1)&& (isa(M,'function_handle')==1)
        L=my_normest(MtM,MtM,n,10^-6,50)+lambda*my_normest(DtD,DtD,n,10^-6,50)/mu;%2 instead of 1 to give us more window
else
        L=norm(M)^2+lambda*norm(D)^2/mu;
end
    
iter_num=0;
OK=0;
val_old=inf;
time=[];
    while(iter_num<maxiter)
        t=(1+sqrt(1+4*t_old^2))/2;
        if isa(D,'function_handle')
            uk = D(v);
        else
            uk = D*v;
        end    
        uk = uk./max(mu,abs(uk));
    
        if isa(Dt,'function_handle')
           df = Dt(uk);
        else
           df = D'*uk;
        end
        if (isa(D,'function_handle')==1)&& (isa(M,'function_handle')==1)
            x=v-(MtM(v)-Mtb+lambda*df)/L;           
        else
            x=v-(MtM*v-Mtb+lambda*df)/L;
        end
        z=x;
        
        if isa(M,'function_handle')
           value_L2=0.5*(M(x)-b)'*(M(x)-b);
        else
           value_L2=0.5*(M*x-b)'*(M*x-b);
        end
        if isa(D,'function_handle')
           temp=D(x);
           f_val=value_L2+lambda*norm(D(x),1);
        else
           temp=D*x;
           f_val=value_L2+lambda*norm(D*x,1);
        end
%         
%         %calculate the Huber function
%         ind= find(abs(temp)<mu);
%         f_temp=zeros(m,1);
%         f_temp(ind)=temp(ind).^2/(2*mu);
%         ind=setdiff(1:m,ind);
%         f_temp(ind)=abs(temp(ind))-mu/2;
%         val=value_L2+lambda*sum(f_temp);
%         
%         if (val>val_old)
%             x=x_old;
%         else
%             val_old=val;
%         end
        
        %v=x+(t_old-1)/t*(x-x_old);
        v=x+(t_old/t)*(z-x)+(t_old-1)/t*(x-x_old);
        qp = abs(f_val - mean(fmean))/mean(fmean);%stop test
        if qp <= acc && OK; break;end
        if qp <= acc && ~OK; OK=1; end
        fmean = [f_val,fmean];
 
        if (length(fmean) > 10) fmean = fmean(1:10);end      
        iter_num=iter_num+1;
        if(print)
            fprintf('iter= %5d value = %10.10f\n',iter_num,norm(x_old-x));
        end
        t_old=t;
        x_old=x;
%         f_out=[f_out;f_val];
        e=norm(x_old-opts.x0)/norm(opts.x0);
        E_out=[E_out,e];
        time=[time,cputime];
    end
    niter=niter+iter_num;