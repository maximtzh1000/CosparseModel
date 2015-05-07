function y = Dt_oper(x, m, n)
%2D difference operation for signal x

yc=reshape(x,m,n);
d1=diff(yc);
d2=diff(yc');
y=[d1(:);d2(:)];