function [P_par,P_ort]=orthogonalspace(A)

[d,n]=size(A); %d is the dimension and n is the number of basis in A
k=rank(A);
[Q,R]=qr(A);
P_par=Q(:,1:k)*Q(:,1:k)';%projection matrix for parrallel space
if k<d
    P_ort=Q(:,k+1:d)*Q(:,k+1:d)';
else
    P_ort=zeros(d,d);
end
