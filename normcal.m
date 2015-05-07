function normcal(M,Mt,n,D,Dt)

global normM;
global normD;

if isa(M,'function_handle')
   MtM=@(x) Mt(M(x));
   normM=my_normest(MtM,MtM,n,10^-6,50);
else
   normM=normest(M'*M);
end

if isa(D,'function_handle')
   DtD=@(x) Dt(D(x));
   normD=my_normest(DtD,DtD,n,10^-6,50);
else
   normD=normest(D'*D);
end