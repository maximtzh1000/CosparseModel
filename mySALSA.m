function [x,niter,E,time]=mySALSA(b,M,Mt,D,Dt,lambda,opts,m,n)

maxiter=opts.maxiter;
acc=opts.acc;
print=opts.print;
fmean = realmin/10;

addpath('tv_fista');
mu=5*10^-3;

MtM=@(x) Mt(M(x));
Mtb=Mt(b);
invLS = @(x) (x - (1/(1+mu))*MtM(x) )/mu;

iter_num=0;
OK=0;
E=[];

time=[];
%initilization 
x_old=Mt(b);
d=zeros(length(x_old),1);
pars.MAXITER=40;
pars.tv='l1';
pars.print=0;

while(iter_num<maxiter)
    v=denoise_bound(reshape(x_old-d,sqrt(n),sqrt(n)),lambda/mu,-Inf,Inf,pars);
    v=v(:);
    r = Mtb + mu*(v+d);
    x = invLS(r);
    d = d + (v - x);
    
    f_val=0.5*(M(x)-b)'*(M(x)-b)+lambda*norm(D(x),1);
    
    qp = abs(f_val - mean(fmean))/mean(fmean);%stop test
    if qp <= acc && OK; break;end
    if qp <= acc && ~OK; OK=1; end
    fmean = [f_val,fmean];
 
    if (length(fmean) > 10) fmean = fmean(1:10);end      
    iter_num=iter_num+1;
    if(print)
        fprintf('iter= %5d value = %10.10f, %10.10f\n',iter_num,f_val,norm(x_old-x));
    end
    x_old=x;
    e=norm(x_old-opts.x0)/norm(opts.x0);
    E=[E,e];
    time=[time,cputime];
    
end
niter=iter_num;