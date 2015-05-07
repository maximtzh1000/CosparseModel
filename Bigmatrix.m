function [y]=Bigmatrix(x,M,Mt,D,Dt,rho,m,n)

y=zeros(m+n,1);
y(1:n)=Mt(M(x(1:n)))+rho*Dt(D(x(1:n)))-rho*Dt(x(n+1:n+m));
y(n+1:n+m)=-rho*D(x(1:n))+rho*x(n+1:n+m);