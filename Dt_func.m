function X=Dt_func(P,m,n)

m1=m-1;
n1=n;
m2=m;
n2=n-1;

P1=P(1:m1*n1);
P1=reshape(P1,m1,n1);

P2=P(m1*n1+1:m1*n1+m2*n2);
P2=reshape(P2,m2,n2);

X=zeros(m,n);
X(1:m-1,:)=P1;
X(:,1:n-1)=X(:,1:n-1)+P2;
X(2:m,:)=X(2:m,:)-P1;
X(:,2:n)=X(:,2:n)-P2;

X=X(:);