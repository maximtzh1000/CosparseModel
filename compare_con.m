% this file is used to compare the Smoothing method with and without
% continuation
% Used in the final version of the paper

clear;clc;close all;
%addpath('/Users/zhaotan/Dropbox/CS_Analysis Model/Toyexample/NESTA_v1.1');
addpath('MRIinit');
randomseed = RandStream('mcg16807','Seed',1);
RandStream.setGlobalStream(randomseed);
global normM;
global normD;

% P = imread('short12.jpg');
% P=double(P);
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
%lambda = sigma*sqrt(2*log(p));
lambda=0.001;
mu=10^-5;
rho=lambda/mu;
acc=10^-10;
ratio=0.1;

%% get the ground true of x0 using NESTA
tic
La=normM;
opts=[];
opts.mu=mu;
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

%% method with continuation using smoothing method
tic
opts=[];
opts.muf=mu;
opts.innerstep=10000;
opts.accf=acc;
opts.r=ratio;
opts.C=5;
opts.maxiter=1000;
opts.print=1;
opts.x0=x0;
opts.xplug=At(b);

t0=cputime;
[x_con,n,E_con,time_con]=Smoothing_conTime(b,A,At,D,Dt,lambda,opts,p,N0);
time_con=time_con-t0;
disp(['Time for Smoothing with continuation is']);
toc
disp(['Number of iteration for Smoothing is ',num2str(n)]);

maxiter=n;

%% method with ADMM
tic
opts=[];
opts.acc=10^-10;
opts.maxiter=250;
opts.print=1;
opts.x0=x0;
opts.xplug=At(b);

t0=cputime;
[x_admm,n,E_admm,time_admm]=mySALSA(b,A,At,D,Dt,lambda,opts,p,N0);
time_admm=time_admm-t0;
disp(['Time for ADMM method is']);
toc
disp(['Number of iteration for ADMM method is ',num2str(n)]);

%% Try with the generalized IST method 
tic
opts=[];
opts.acc=10^-10;
opts.maxiter=2000;
opts.print=1;
opts.x0=x0;
opts.xplug=At(b);

t0=cputime;
[x_ist,n,E_ist,time_ist]=GIST(b,A,At,D,Dt,lambda,opts,p,N0);
time_ist=time_ist-t0;

disp(['Time for generalized IST is']);
toc
disp(['Number of iteration for GIST is ',num2str(n)]);


%% Method of smoothing without continuation
tic
opts=[];
opts.mu=10^-4;
opts.acc=10^-10;
opts.maxiter=2000;
opts.print=1;
opts.x0=x0;
opts.xplug=At(b);

t0=cputime;
[s,n,E,t]=SmoothingMethodTime(b,A,At,D,Dt,lambda,[],opts,p,N0);
t=t-t0;
disp(['Time for smoothing is']);
toc
disp(['Number of iteration for smoothing is ',num2str(n)]);

%% Ploting the results
fntsz = 14; lwdth = 1; %display parameter
figure(2)
imagesc(reshape(real(x_con)*256,N,N));
colormap(gray)
axis image
axis equal
axis off
%title(sprintf('reconstructed image using Smoothing Method with Cont.(%d x %d = %d pixels)',N,N,N*N))

figure(3)
imagesc(reshape(real(s)*256,N,N));
colormap(gray)
axis image
axis equal
axis off
%title(sprintf('reconstructed image without continuation (%d x %d = %d pixels)',N,N,N*N))

error_con=norm(x_con-x)/norm(x);
disp(['The error of smoothing with continuation is ',num2str(error_con)]);
error=norm(s-x)/norm(x);
disp(['The error of smoothing is ',num2str(error)]);

%% plot the function value 
figure(4)
lwdth=2;
%range=1:30:maxiter;
semilogy(time_ist,E_ist,'-.g','linewidth',lwdth);
hold on
semilogy(t,E,':r','linewidth',lwdth);
semilogy(time_con,E_con,'-b','linewidth',lwdth);
semilogy(time_admm,E_admm,'--m','linewidth',lwdth);
h=legend('GIST','SFISTA Without Cont','SFISTA With Cont.','SALSA',fntsz);
set(h,'FontSize',fntsz);
xlabel('CPU Time (sec.)','fontsize',fntsz);
ylabel('Error E(x)','fontsize',fntsz);
xlim([0,200]);
% figure(4)
% plot(f,'linewidth',lwdth);
% hold on;
% plot(f_con,'r','linewidth',lwdth);
% xlim([1,737]);
% legend('Without continuation','With continuation',fntsz);
% xlabel('Iteration Number','fontsize',fntsz);
% ylabel('Objective Function F(x)','fontsize',fntsz);
