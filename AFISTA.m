function [x,z,niter,f_out,E_out]=AFISTA(b,M,Mt,D,Dt,lambda,x_init,z_init,opts,m,n)
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
rho=opts.rho;
print=opts.print;
fmean = realmin/10;

if(isempty(z_init))
    if isa(M,'function_handle')
       x_old=Mt(b);
       z_old=D(x_old);
    else
        x_old=M'*b;
        z_old=D*x_old;
    end
else
    z_old=z_init;
    x_old=x_init;
end

y=z_old;
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
% if (isa(D,'function_handle')==1)&& (isa(M,'function_handle')==1)
%         tMatrix=@(x) Bigmatrix(x,M,Mt,D,Dt,rho,m,n);
%         L=1.0*my_normest(tMatrix,tMatrix,m+n,10^-6,50);%2 instead of 1 to give us more window
% else
%         tMatrix=[MtM+DtD*rho,-Dt*rho;-D*rho,rho*eye(m)];
%         L=normest(tMatrix);
% end
if (isa(D,'function_handle')==1)&& (isa(M,'function_handle')==1)
        L=my_normest(MtM,MtM,n,10^-6,50)+(my_normest(DtD,DtD,n,10^-6,50)+1)*rho;%2 instead of 1 to give us more window
else
        L=norm(M)^2+(1+norm(D)^2)*rho;
end

iter_num=0;
OK=0;
val_old=inf;
    while(iter_num<maxiter)
        t=(1+sqrt(1+4*t_old^2))/2;
        if (isa(D,'function_handle')==1)&& (isa(M,'function_handle')==1)
            x=v-(MtM(v)+rho*DtD(v)-(Mtb+rho*Dt(y)))/L;
            y_tmp=y-rho*(y-D(x))/L;            
        else
            x=v-(MtM*v+rho*DtD*v-(Mtb+rho*Dt*y))/L;
            y_tmp=y-rho*(y-D*x)/L; 
        end

        ind=find(abs(y_tmp)>lambda/L);
        z=zeros(m,1);
        z(ind)=(abs(y_tmp(ind))-lambda/L).*sign(y_tmp(ind));
        
        x0=x;
        z0=z;
        
        if isa(M,'function_handle')
           value_L2=0.5*(M(x)-b)'*(M(x)-b);
        else
           value_L2=0.5*(M*x-b)'*(M*x-b);
        end
        if isa(D,'function_handle')
           f_val=value_L2+lambda*norm(D(x),1);
           value=value_L2+lambda*norm(z,1)+0.5*rho*(D(x)-z)'*(D(x)-z);
        else
           f_val=value_L2+lambda*norm(D*x,1);
           value=value_L2+lambda*norm(z,1)+0.5*rho*(D*x-z)'*(D*x-z); 
        end
        
        if (value>val_old)
            x=x_old;
            z=z_old;
        else
            val_old=value;
        end
        
        v=x+(t_old/t)*(x0-x)+(t_old-1)/t*(x-x_old);
        y=z+(t_old/t)*(z0-z)+(t_old-1)/t*(z-z_old);
        
        qp = abs(f_val - mean(fmean))/mean(fmean);%stop test
        if qp <= acc && OK; break;end
        if qp <= acc && ~OK; OK=1; end
        fmean = [f_val,fmean];

        if (length(fmean) > 10) fmean = fmean(1:10);end

%         if (abs(value_old-value)<acc)
%             break;
%         end
        
        iter_num=iter_num+1;
        if(print)
            fprintf('iter= %5d value = %10.10f %10.10f\n',iter_num,value,norm(x_old-x)+norm(z_old-z));
        end
        z_old=z;
        t_old=t;
        x_old=x;
        %value_old=value;
        f_out=[f_out;f_val];
        e=norm(x_old-opts.x0)/norm(opts.x0);
        E_out=[E_out,e];
    end
    niter=niter+iter_num;