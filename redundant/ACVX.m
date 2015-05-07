function [s]=ACVX(b,M,D,lambda)

[m,n]= size(D);

cvx_begin quiet
     variables x(n)
     %minimize(sum(square_abs(b-(M*x)))+(lambda*norm(z,1))+(rho*sum(square_abs(z-(D*x)))))
     minimize(sum(square_abs(b-(M*x)))+(lambda*norm(D*x,1)))
cvx_end

s=x;