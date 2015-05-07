%This demo use several images to show the 
%algorithms for analysis compressive sensing algorithm
%The M is taken as 2D discrete fourier transforms
%The D matrix is taken as 2D difference operator
%Samples are taken in a CT fashion

%To do list: use 2D fourier transformation to form the measurement matrix

clear;clc;close all;
%addpath('/Users/zhaotan/Dropbox/CS_Analysis Model/Toyexample/NESTA_v1.1');
randomseed = RandStream('mcg16807','Seed',3);
RandStream.setGlobalStream(randomseed);

fntsz = 14; lwdth = 1; %display parameter

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

b=A(x);
L=length(b);
noise=sigma*randn(L,1);
b=b+noise;

normcal(A,At,N0,D,Dt);

SNR=20*log10(norm(b)/norm(noise));
%lambda = sigma*sqrt(2*log(p));
lambda=0.004;
mu=10^-5;
rho=lambda/mu;
acc=10^-7;
maxiter=3000;

%% This part is based on the method developed by Amir Beck in 2009
% tic
% pars.acc=acc;
% pars.MAXITER=130;
% pars.maxiter=1300;
% pars.denoiseiter=5;
% pars.print=1;
% pars.x0=x;
% pars.xplug=At(b);
% [x_den,n_den,f_den,E_den]=mfista_ab(b,A,At,D,Dt,lambda,pars,p,N0);
% disp(['Time for dual method is']);
% toc
% disp(['Number of iteration for dual problem is ',num2str(n_den)]);

%% This part is developed based on Smoothing Method
tic
opts=[];
opts.muf=mu;
opts.innerstep=10000;
opts.accf=acc;
opts.r=0.5;
opts.C=5;
opts.maxiter=maxiter;
opts.print=1;
opts.x0=x;
opts.xplug=At(b);

[x_s,n,f_s,E_s]=Smoothing_con(b,A,At,D,Dt,lambda,opts,p,N0);
disp(['Time for Smoothing with continuation is']);
toc
disp(['Number of iteration for Smoothing is ',num2str(n)]);

%% Decomposition with continuation
tic
opts=[];
opts.rhof=rho;
opts.innerstep=1000;
opts.accf=acc;
opts.r=0.5;
opts.C=5;
opts.maxiter=maxiter;
opts.print=1;
opts.x0=x;
opts.xplug=At(b);

[s,n,f_afista,E_afista]=AFISTA_con(b,A,At,D,Dt,lambda,opts,p,N0);
disp(['Time for AFISTA with continuation is']);
toc
disp(['Number of iteration for FISTA is ',num2str(n)]);

%% plot results
figure(2)
imagesc(reshape(x_s*256,N,N));
colormap(gray)
axis image
axis equal
axis off
title(sprintf('reconstructed image using Smoothing Method (%d x %d = %d pixels)',N,N,N*N))

figure(3)
imagesc(reshape(s*256,N,N));
colormap(gray)
axis image
axis equal
axis off
title(sprintf('reconstructed image using Decomposition Method (%d x %d = %d pixels)',N,N,N*N))

error=norm(s-x)/norm(x);
disp(['The error of decomposition with continuation is ',num2str(error)]);
error_xk=norm(x_s-x)/norm(x);
disp(['The error of smoothing with continuation is ',num2str(error_xk)]);

%% plot the function value 
figure(5)
plot(f_afista,'linewidth',lwdth);
hold on;
plot(f_s,'r','linewidth',lwdth);
legend('fval decomposition','fval smoothing','fontsize',fntsz);
xlabel('Iteration Number');
ylabel('Objective Function F(x)');

%% analyze error with iteration
figure(6)
plot(E_afista,'linewidth',lwdth);
hold on;
plot(E_s,'r','linewidth',lwdth);
legend('fval decomposition','fval smoothing','fontsize',fntsz);
