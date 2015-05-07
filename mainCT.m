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

pro=12; %number of projections
[mhi]=RandMask_rect(pro,N,N);
A=@(x) A_fhp_rect(x,mhi,N,N);
At=@(y) At_fhp_rect(y,mhi,N,N);

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

[x0,niter,residuals,outputData,E_out,opts] = MyCore_Nesterov_UP(A,At,b,lambda,La,opts.mu,opts);

%% begin the comparsion
mu=[10^-3,10^-4];
l=length(mu);
maxiter=5000;

F_decomp=zeros(maxiter,l);
F_smooth=zeros(maxiter,l);
F_nesta=zeros(maxiter,l);
E_decomp=zeros(maxiter,l);
E_smooth=zeros(maxiter,l);
E_nesta=zeros(maxiter,l);
for i=1:2
    [f_smooth,f_decomp,f_nesta,e_smooth,e_decomp,e_nesta]=shepplogantest(x0,b,A,At,D,Dt,lambda,p,N0,mu(i),maxiter);
    F_decomp(:,i)=f_decomp';
    F_smooth(:,i)=f_smooth';
    F_nesta(:,i)=f_nesta';
    E_decomp(:,i)=e_decomp;
    E_smooth(:,i)=e_smooth;
    E_nesta(:,i)=e_nesta;
end

%% plot the results
fntsz = 20; lwdth = 1; %display parameter

figure(1);
plot(F_decomp(1:100,1),'linewidth',lwdth);
hold on;
plot(F_smooth(1:100,1),'g','linewidth',lwdth);
plot(F_nesta(2:101,1),'r','linewidth',lwdth);
plot(F_decomp(1:100,2),'--','linewidth',lwdth);
plot(F_smooth(1:100,2),'--g','linewidth',lwdth);
plot(F_nesta(2:101,2),'--r','linewidth',lwdth);
% plot(F_decomp(1:100,3),'r','linewidth',lwdth);
% plot(F_smooth(1:100,3),'o-r','linewidth',lwdth);
% plot(F_nesta(2:101,3),'--r','linewidth',lwdth);
xlim([1,100]);
legend('Decomposition \rho=10^3\lambda','Smoothing \mu=10^{-3}\lambda^{-1}','NESTA \mu=10^{-3}\lambda^{-1}','Decomposition \rho=10^4\lambda','Smoothing \mu=10^{-4}\lambda^{-1}','NESTA \mu=10^{-4}\lambda^{-1}',fntsz);
xlabel('Iteration Number','fontsize',fntsz);
ylabel('Objective Function F(x)','fontsize',fntsz);

figure(2);
semilogy(E_decomp(:,1),'linewidth',lwdth);
hold on;
semilogy(E_smooth(:,1),'g','linewidth',lwdth);
semilogy(E_nesta(:,1),'r','linewidth',lwdth);
semilogy(E_decomp(:,2),'--','linewidth',lwdth);
semilogy(E_smooth(:,2),'--g','linewidth',lwdth);
semilogy(E_nesta(:,2),'--r','linewidth',lwdth);
% plot(F_decomp(1:100,3),'r','linewidth',lwdth);
% plot(F_smooth(1:100,3),'o-r','linewidth',lwdth);
% plot(F_nesta(2:101,3),'--r','linewidth',lwdth);
xlim([1,5000]);
legend('Decomposition \rho=10^3\lambda','Smoothing \mu=10^{-3}\lambda^{-1}','NESTA \mu=10^{-3}\lambda^{-1}','Decomposition \rho=10^4\lambda','Smoothing \mu=10^{-4}\lambda^{-1}','NESTA \mu=10^{-4}\lambda^{-1}',fntsz);
xlabel('Iteration Number','fontsize',fntsz);
ylabel('Relative Error E(x)','fontsize',fntsz);