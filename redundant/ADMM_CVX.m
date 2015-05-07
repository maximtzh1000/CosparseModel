function [s]=ADMM_CVX(b,M,D,lambda,x_init,rho,maxiter,acc)

fmean = realmin/10;
iter_num=0;
OK=0;
z=D*x_init;
[m,n]=size(D);
d=zeros(m,1);

while (iter_num<maxiter)
    cvx_begin quiet
              variable x(n)
              minimize (0.5*sum(square_abs(b-(M*x)))+0.5*rho*sum(square_abs((D*x)-z-d)))
    cvx_end
    cvx_begin quiet
              variable z(m)
              minimize (lambda*norm(z,1)+0.5*rho*sum(square_abs((D*x)-z-d)))
    cvx_end
    d=d-(D*x-z);
    f_val=lambda*norm(D*x,1)+0.5*norm(M*x-b,2)^2
    qp = abs(f_val - mean(fmean))/mean(fmean);%stop test
    if qp <= acc && OK; break;end
    if qp <= acc && ~OK; OK=1; end
    fmean = [f_val,fmean];
 
    if (length(fmean) > 10) fmean = fmean(1:10);end      
    iter_num=iter_num+1;
end
s=x;