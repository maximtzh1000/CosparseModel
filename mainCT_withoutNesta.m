clear;clc;close all;
addpath('/Users/zhaotan/Dropbox/CS_Analysis Model/Toyexample/NESTA_v1.1');
randomseed = RandStream('mcg16807','Seed',1);
RandStream.setGlobalStream(randomseed);
global normM;
global normD;

P = 256.0*phantom(256);
[N] = size(P,1);
x=P(:);

figure(1)
imagesc(P);
colormap(gray)
axis image
axis equal
axis off
title(sprintf('Original image (%d x %d = %d pixels)',N,N,N*N))

x=x/256;

p=2*(N-1)*N;%length of operator D
N0=N*N;%length of x
sigma=1*10^-3;

pro=15; %number of projections

[mask_temp,Mh,mi,mhi] = LineMask(pro,N);
mask = fftshift(mask_temp);
A = @(x)  masked_FFT(x,mask,N,N);
At = @(x) (masked_FFT_t(x,mask));
AtA = @(x) (ifft2c(mask.*fft2c(x))) ;

D=@(x) D_func(x,N,N);
Dt=@(p) Dt_func(p,N,N);

normcal(A,At,N0,D,Dt);

b=A(x);
L=length(b);
noise=sigma*randn(L,1);
b=b+noise;
normcal(A,At,N0,D,Dt);
SNR=20*log10(norm(b)/norm(noise));
lambda=0.001;

%% get the ground true by running the NESTA algorithm for 10000 iterations
tic
La=normM;
opts=[];
opts.mu=10^-4;
opts.TolVar=10^-20;
opts.U=D;
opts.Ut=Dt;
opts.normU=normD;
opts.maxiter=10000;
opts.print=1;
opts.x0=x;
opts.xplug=At(b);
%X_init=M'*b;

%[x0,niter,residuals,outputData,E_out,opts] = MyCore_Nesterov_UP(A,At,b,lambda,La,opts.mu,opts);
x0=x;

%% begin the comparsion
mu=[10^-2,10^-3,10^-4,];
l=length(mu);
maxiter=5000;

F_decomp=zeros(maxiter,l);
F_smooth=zeros(maxiter,l);
E_decomp=zeros(maxiter,l);
E_smooth=zeros(maxiter,l);

for i=1:l
    [f_smooth,f_decomp,e_smooth,e_decomp]=shepplogantest(x0,b,A,At,D,Dt,lambda,p,N0,mu(i),maxiter);
    F_decomp(:,i)=f_decomp';
    F_smooth(:,i)=f_smooth';
    E_decomp(:,i)=e_decomp;
    E_smooth(:,i)=e_smooth;
end

%% plot the results
fntsz = 12; lwdth = 1; %display parameter
range=1:100:5000;
figure(1);
plot(range,F_decomp(range,1),'-bs','linewidth',lwdth);
hold on;
plot(range,F_smooth(range,1),'--b*','linewidth',lwdth);
plot(range,F_decomp(range,2),'-gx','linewidth',lwdth);
plot(range,F_smooth(range,2),'--g^','linewidth',lwdth);
plot(range,F_decomp(range,3),'-r<','linewidth',lwdth);
plot(range,F_smooth(range,3),'--ro','linewidth',lwdth);
xlim([1,5000]);
h=legend('DFISTA \rho=10^2\lambda','SFISTA \mu=10^{-2}\lambda^{-1}','DFISTA \rho=10^3\lambda','SFISTA \mu=10^{-3}\lambda^{-1}','DFISTA \rho=10^4\lambda','SFISTA \mu=10^{-4}\lambda^{-1}',fntsz);
set(h,'FontSize',fntsz);
xlabel('Iteration Number','fontsize',fntsz);
ylabel('Objective Function F(x)','fontsize',fntsz);

figure(2);
range=1:100:5000;
semilogy(range,E_decomp(range,1),'-bs','linewidth',lwdth);
hold on;
semilogy(range,E_smooth(range,1),'--b*','linewidth',lwdth);
semilogy(range,E_decomp(range,2),'-gx','linewidth',lwdth);
semilogy(range,E_smooth(range,2),'--g^','linewidth',lwdth);
semilogy(range,E_decomp(range,3),'-r<','linewidth',lwdth);
semilogy(range,E_smooth(range,3),'--ro','linewidth',lwdth);
xlim([1,5000]);
ylim([10^-3,10]);
h=legend('DFISTA \rho=10^2\lambda','SFISTA \mu=10^{-2}\lambda^{-1}','DFISTA \rho=10^3\lambda','SFISTA \mu=10^{-3}\lambda^{-1}','DFISTA \rho=10^4\lambda','SFISTA \mu=10^{-4}\lambda^{-1}',fntsz);
set(h,'FontSize',fntsz);
xlabel('Iteration Number','fontsize',fntsz);
ylabel('Error E(x)','fontsize',fntsz);