function [x,z,niter,f_out,E_out]=AFISTA_old(b,M,Mt,D,Dt,lambda,x_init,z_init,opts,m,n)
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

if(isempty(z_init))
    z_old=zeros(m,1);
    x_old=zeros(n,1);
else
    z_old=z_init;
    x_old=x_init;
end
y=z_old;
v=x_old;
value_old=inf;
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
        tMatrix=@(x) Bigmatrix(x,M,Mt,D,Dt,rho,m,n);
        L=4*my_normest(tMatrix,tMatrix,m+n);%4 instead of 2 to give us more window
else
        tMatrix=[MtM/rho+DtD,Dt;D,eye(m)];
        %TTT=@(x) matrix(x,M,Mt,D,Dt,rho,m,n);
        L=4*normest(tMatrix);
        %L2=4*my_normest(TTT,TTT,n+m);
end
    
iter_num=0;
OK=0;
    while(iter_num<maxiter)
        t=(1+sqrt(1+4*t_old^2))/2;
        if (isa(D,'function_handle')==1)&& (isa(M,'function_handle')==1)
            x=v-2*(MtM(v)/rho+DtD(v)-(Mtb/rho+Dt(y)))/L;
            y_tmp=y-2*(y-D(v))/L;            
        elseif (isa(D,'function_handle')==0)&& (isa(M,'function_handle')==1)
            x=v-2*t0*(MtM(v)+rho*DtD*v-(Mtb+rho*Dt*y));
            y_tmp=y-2*t1*(y-D*v);
        elseif (isa(D,'function_handle')==1)&& (isa(M,'function_handle')==0)
            x=v-2*t0*(MtM*v+rho*DtD(v)-(Mtb+rho*Dt(y)));
            y_tmp=y-2*t1*(y-D(v));
        else
            x=v-2*(MtM*v/rho+DtD*v-(Mtb/rho+Dt*y))/L;
            y_tmp=y-2*(y-D*v)/L; 
        end

        ind=find(abs(y_tmp)>lambda/(rho*L));
        z=zeros(m,1);
        z(ind)=(abs(y_tmp(ind))-lambda/(L*rho)).*sign(y_tmp(ind));
        y=z+(t_old-1)/t*(z-z_old);
        v=x+(t_old-1)/t*(x-x_old);
        if isa(M,'function_handle')
           value_L2=(M(x)-b)'*(M(x)-b);
        else
           value_L2=(M*x-b)'*(M*x-b);
        end
        if isa(D,'function_handle')
           f_val=value_L2+lambda*norm(D(x),1);
           value=value_L2+lambda*norm(z,1)+rho*(D(x)-z)'*(D(x)-z);
        else
           f_val=value_L2+lambda*norm(D*x,1);
           value=value_L2+lambda*norm(z,1)+rho*(D*x-z)'*(D*x-z); 
        end
        if(abs(value_old-value)<acc&&OK==5)
            break;
        elseif(abs(value_old-value)<acc)
            OK=OK+1;
        else
            OK=0;
        end
        
        iter_num=iter_num+1;
        if(print)
            fprintf('iter= %5d value = %10.10f %10.10f\n',iter_num,value,norm(x_old-x)+norm(z_old-z));
        end
        z_old=z;
        t_old=t;
        x_old=x;
        value_old=value;
        f_out=[f_out;f_val];
        e=norm(x_old-opts.x0)/norm(opts.x0);
        E_out=[E_out,e];
    end
    niter=niter+iter_num;