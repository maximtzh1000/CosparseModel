function [y]=matrix(x,M,Mt,D,Dt,rho,m,n)

y=zeros(m+n,1);
y(1:n)=Mt*M*x(1:n)/rho+Dt*D*x(1:n)+Dt*x(n+1:n+m);
y(n+1:n+m)=D*x(1:n)+x(n+1:n+m);