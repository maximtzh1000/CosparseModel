function [D]=normalize(D)
% normalize the matrix D by each column

[~,n]=size(D);
for i=1:n
    D(:,i)=D(:,i)/norm(D(:,i));
end