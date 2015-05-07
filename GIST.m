function [x,niter,E,time]=GIST(b,M,Mt,D,Dt,lambda,opts,m,n)
% this function implement the method proposed in Loris 11 paper

maxiter=opts.maxiter;
acc=opts.acc;
print=opts.print;
fmean = realmin/10;


x_old=Mt(b);

w=zeros(m,1);

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
    DDt=@(x) D(Dt(x));
else
    DDt=D*D';
end

tau=my_normest(MtM,MtM,n,10^-6,50);
tau=1.50/tau;

delta=my_normest(DDt,DDt,m,10^-6,50);
delta=0.90/delta;

iter_num=0;
OK=0;
E=[];

time=[];
while(iter_num<maxiter)
    g=x_old+tau*Mtb-tau*MtM(x_old);
    w=w+delta/tau*D(g-tau*Dt(w));
    ind=find(abs(w)>lambda);
    w(ind)=lambda*w(ind)./abs(w(ind));
    x=g-tau*Dt(w);
    f_val=0.5*(M(x)-b)'*(M(x)-b)+lambda*norm(D(x),1);
    
    qp = abs(f_val - mean(fmean))/mean(fmean);%stop test
    if qp <= acc && OK; break;end
    if qp <= acc && ~OK; OK=1; end
    fmean = [f_val,fmean];
 
    if (length(fmean) > 10) fmean = fmean(1:10);end      
    iter_num=iter_num+1;
    if(print)
        fprintf('iter= %5d value = %10.10f\n',iter_num,norm(x_old-x));
    end
    x_old=x;
    e=norm(x_old-opts.x0)/norm(opts.x0);
    E=[E,e];
    time=[time,cputime];
end
niter=niter+iter_num;