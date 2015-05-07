%This M file compares the Conjugate method with the SFISTA we proposed

clear;clc;close all;
addpath('/Users/zhaotan/Dropbox/CS_Analysis Model/Toyexample/NESTA_v1.1');
randomseed = RandStream('mcg16807','Seed',2);
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

b=A(x);
L=length(b);
noise=sigma*randn(L,1);
b=b+noise;
x_org=x;
normcal(A,At,N0,D,Dt);

SNR=20*log10(norm(b)/norm(noise));
%lambda = sigma*sqrt(2*log(p));
lambda=0.001;
mu=10^-4;
rho=lambda/mu;
%rho=10^6*lambda*sigma;
acc=10^-30;
maxiter=3000;

%% This part is developed based on Conjugate gradient method
opts=[];
opts.mu=mu;
opts.acc=acc;
opts.maxiter=300;
opts.print=1;
opts.x0=x_org;
opts.xplug=At(b);

tic
time0=cputime;
[x_s,n,E_s,time_s]=LDP(b,A,At,D,Dt,lambda,[],opts,p,N0);
time_s=time_s-time0;
disp(['Time for LDP is']);
toc
disp(['Number of iteration for LDP is ',num2str(n)]);

% %% GIST method
% tic
% opts=[];
% opts.acc=10^-10;
% opts.maxiter=maxiter;
% opts.print=1;
% opts.x0=x_org;
% opts.xplug=At(b);
% 
% t0=cputime;
% [x_ist,n,E_ist,time_ist]=GIST(b,A,At,D,Dt,lambda,opts,p,N0);
% time_ist=time_ist-t0;
% 
% disp(['Time for generalized IST is']);
% toc
% disp(['Number of iteration for GIST is ',num2str(n)]);

%% This part is developed based on Smoothing method
opts=[];
opts.mu=mu;
opts.acc=acc;
opts.maxiter=maxiter;
opts.print=1;
opts.x0=x_org;
opts.xplug=At(b);

tic
time0=cputime;
[x_fista,n,E_fista,time_fista]=SmoothingMethodTime(b,A,At,D,Dt,lambda,[],opts,p,N0);
time_fista=time_fista-time0;
disp(['Time for Smoothing is']);
toc
disp(['Number of iteration for Smoothing is ',num2str(n)]);

% %% This part is developed based on Decomposition method
% tic
% opts=[];
% opts.rho=rho;
% opts.acc=acc;
% opts.maxiter=maxiter;
% opts.print=1;
% opts.x0=x_org;
% opts.xplug=At(b);
% 
% time0=cputime;
% [s,~,n,E_decomp,time_decomp]=AFISTATime(b,A,At,D,Dt,lambda,[],[],opts,p,N0);
% time_decomp=time_decomp-time0;
% disp(['Time for decomposition is']);
% toc
% disp(['Number of iteration for decomposition is ',num2str(n)]);
%% plot results
figure(2)
imagesc(reshape(real(x_s)*256,N,N));
colormap(gray)
axis image
axis equal
axis off
title(sprintf('reconstructed image using LDP (%d x %d = %d pixels)',N,N,N*N))

figure(3)
imagesc(reshape(real(x_fista)*256,N,N));
colormap(gray)
axis image
axis equal
axis off
title(sprintf('reconstructed image using Smoothing (%d x %d = %d pixels)',N,N,N*N))

error=norm(x_fista-x)/norm(x);
disp(['The error of Smoothing is ',num2str(error)]);
error_xs=norm(x_s-x)/norm(x);
disp(['The error of LDP method is ',num2str(error_xs)]);
% 
% error_s=norm(s-x)/norm(x);
% disp(['The error of LDP method is ',num2str(error_s)]);

%% plot the function value 
% figure(5)
fntsz = 14; lwdth = 2; %display parameter

% analyze error with iteration
figure(6)
% range=1:100:maxiter;
% Range=1:10:300;
semilogy(time_s,E_s,'--r','linewidth',lwdth);
hold on;
%semilogy(time_ist(range),E_ist(range),'-go','linewidth',lwdth);
%semilogy(time_decomp(range),E_decomp(range),'-md','linewidth',lwdth);
semilogy(time_fista,E_fista,'linewidth',lwdth);
xlim([0,250]);
h_legend=legend('CGD \mu=10^{-4}\lambda^{-1}','SFISTA \mu=10^{-4}\lambda^{-1}',fntsz);
set(h_legend,'FontSize',fntsz);
xlabel('CPU Time(sec.)','fontsize',fntsz);
ylabel('Error E(x)','fontsize',fntsz);