function [u,v,niter]=Decomp_dual(b,M,Mt,D,Dt,lambda,u_old,v_old,opts,m,n)

maxiter=opts.maxiter;
acc=opts.acc;
rho=opts.rho;
print=opts.print;
fmean = realmin/10;

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

if(isempty(u_old))
    u_old=b;
    v_old=zeros(m,1);
end

x=u_old;
y=v_old;
niter=0;

t_old=1;
if (isa(D,'function_handle')==1)&& (isa(M,'function_handle')==1)
        tMatrix=@(x) Bigmatrix(x,M,Mt,D,Dt,rho,m,n);
        L=1.0*my_normest(tMatrix,tMatrix,m+n,10^-6,50);%2 instead of 1 to give us more window
else
        tMatrix=[M*Mt*rho+eye(size(M,1)),rho*M*Dt;rho*D*Mt,rho*D*Dt];
        L=normest(tMatrix);
end

iter_num=0;
OK=0;
    while(iter_num<maxiter)
        t=(1+sqrt(1+4*t_old^2))/2;
            u=x-((x-b)+rho*M*(M'*x+D'*y))/L;
            v=y-rho*D*(M'*x+D'*y)/L; 

        ind=find(abs(v)>lambda);
        
        x=u+(t_old-1)/t*(u-u_old);
        y=v+(t_old-1)/t*(v-v_old);
        f_val=0.5*norm(b-u)^2+rho/2*norm(M'*u+D'*v)^2;
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
            fprintf('iter= %5d value = %10.10f %10.10f\n',iter_num,f_val,norm(u_old-u)+norm(v_old-v));
        end
        u_old=u;
        t_old=t;
        v_old=v;
        %value_old=value;
    end
    niter=niter+iter_num;
