%this is the montecarlo simulation to test the performance of smoothing and
%decomposition method
% Author: Zhao Tan Date: April 28th, 2013

sigma=0.05:0.05:1;
r=0.05:0.05:1;
M=50; %number of monto carlo runs
E_decomp=zeros(length(sigma),length(r));
E_smooth=zeros(length(sigma),length(r));
for i=1:length(sigma)
    disp(['the value of sigma is ',num2str(sigma(i))]);
    for j=1:length(r)
        disp(['the value of r is ',num2str(r(j))]);
        for k=1:M
            disp(['The number of monte carlo is ',num2str(k)]);
            [error_decomp,error_smooth]=randomtest(sigma(i),r(j));
            E_decomp(i,j)=E_decomp(i,j)+error_decomp;
            E_smooth(i,j)=E_smooth(i,j)+error_smooth;
        end
        E_decomp(i,j)=E_decomp(i,j)/M;
        E_smooth(i,j)=E_smooth(i,j)/M;
    end
end
save('E_decomp');
save('E_smooth');