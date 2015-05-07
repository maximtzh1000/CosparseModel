function [v,niter,E_out,time]=LDP(b,M,Mt,D,Dt,lambda,x_init,opts,m,n)
%  Created on Oct 26th, 2013
%  Author: Zhao Tan
%  1. This code implemented the analysis compressive senings using conjuguate gradient method.
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
alpha=0.05;
beta=0.6;
 
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

iter_num=0;
OK=0;
    
df=gradient_val(v);

if (isa(D,'function_handle')==1)&& (isa(M,'function_handle')==1)
   g=MtM(v)-Mtb+lambda*df;           
else
   g=MtM*v-Mtb+lambda*df;
end
dv=-g;
f0=func_val(v);

time=[];
while(iter_num<maxiter)
     t=1;
     while (func_val(v+t*dv)>f0+alpha*t*real(g'*dv))
         t=beta*t;
     end
     f0=func_val(v+t*dv);
     v_old=v;
     v=v+t*dv;
     
        if isa(M,'function_handle')
           squareterm=0.5*(M(v)-b)'*(M(v)-b);
        else
           squareterm=0.5*(M*v-b)'*(M*v-b);
        end
     
        if isa(D,'function_handle')
           f_val=squareterm+lambda*norm(D(v),1);
        else
           f_val=squareterm+lambda*norm(D*v,1);
        end
        
        g_old=g;
                
        df=gradient_val(v);
        
        if (isa(D,'function_handle')==1)&& (isa(M,'function_handle')==1)
            g=MtM(v)-Mtb+lambda*df;           
        else
            g=MtM*v-Mtb+lambda*df;
        end
        
        gamma=norm(g,2)^2/norm(g_old,2)^2;
        
        dv=-g+gamma*dv;
               
        qp = abs(f_val - mean(fmean))/mean(fmean);%stop test
        if qp <= acc && OK; break;end
        if qp <= acc && ~OK; OK=1; end
        fmean = [f_val,fmean];
 
        if (length(fmean) > 10) fmean = fmean(1:10);end      
        iter_num=iter_num+1;
        if(print)
             fprintf('iter= %5d value = %10.10f\n',iter_num,norm(v_old-v));
        end
        
        f_out=[f_out;f_val];
        e=norm(v-opts.x0)/norm(opts.x0);
        E_out=[E_out,e];
        time=[time,cputime];
end
    niter=niter+iter_num;

    
    function [funcval]=func_val(x)
        if isa(D,'function_handle')
          l1temp=D(x);
        else
          l1temp=D*x;
        end
        
        if isa(M,'function_handle')
           value_L2=0.5*(M(x)-b)'*(M(x)-b);
        else
           value_L2=0.5*(M*x-b)'*(M*x-b);
        end
        
        %calculate the Huber function
        ind= find(abs(l1temp)<mu);
        f_temp=zeros(m,1);
        f_temp(ind)=l1temp(ind).^2/(2*mu);
        ind=setdiff(1:m,ind);
        f_temp(ind)=abs(l1temp(ind))-mu/2;
        funcval=value_L2+lambda*sum(f_temp);        
    end

    function [grad]=gradient_val(x)
        if isa(D,'function_handle')
           uk = D(x);
        else
           uk = D*x;
        end    
           uk = uk./max(mu,abs(uk));
    
        if isa(Dt,'function_handle')
           grad = Dt(uk);
        else
           grad = D'*uk;
        end
        
    end
end