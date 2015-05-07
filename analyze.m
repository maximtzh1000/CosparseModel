function [E]=analyze(X,x0)

n=size(X,2);
E=zeros(n,1);

for i=1:n
    x=X(:,i);
    E(i,1)=norm(x-x0)/norm(x0);
end