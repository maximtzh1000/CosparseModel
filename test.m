%This demo use several images to show the 
%algorithms for analysis compressive sensing algorithm
%The M is taken as 2D discrete fourier transforms
%The D matrix is taken as 2D difference operator
%Samples are taken in a MRI fashion

clear;clc;close all;
addpath('/Users/zhaotan/Dropbox/CS_Analysis Model/Toyexample/NESTA_v1.1');
randomseed = RandStream('mcg16807','Seed',2);
RandStream.setGlobalStream(randomseed);
global normM;
global normD;

fntsz = 14; lwdth = 1; %display parameter

%P = double(imread('pics/MRI_Coronal_Brain.jpg'));
P = double(imread('pics/3DMR_Chest.jpg'));
%P = double(imread('pics/3DMR_Renal_Arteries.jpg'));
%P = double(imread('pics/cameraman.tif'));
%P = double(imread('pics/lena.jpg'));
%P = 256.0*phantom(256);
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
sigma=1*10^-5;
d = 2;
[mhi] = RandMask(double(N/d),double(N/d),N,N);

A=@(x) A_fhp_rect(x,mhi,N,N);
At=@(y) At_fhp_rect(y,mhi,N,N);

D=@(x) D_func(x,N,N);
Dt=@(p) Dt_func(p,N,N);

b=A(x);
L=length(b);
noise=sigma*randn(L,1);
b=b+noise;

x_org=x;
normcal(A,At,N0,D,Dt);

SNR=20*log10(norm(b)/norm(noise));
lambda = sigma*sqrt(2*log(p));
%lambda=0.04;
mu=10^-5;
rho=lambda/mu;
%rho=10^6*lambda*sigma;
acc=10^-20;
maxiter=1000;

%% This part is developed based on Smoothing method
tic
opts=[];
opts.mu=mu;
opts.acc=acc;
opts.maxiter=maxiter;
opts.print=1;
opts.x0=x_org;
opts.xplug=At(b);

[x_s,n,f_s,E_s]=SmoothingMethod(b,A,At,D,Dt,lambda,[],opts,p,N0);
disp(['Time for smoothing is']);
toc
disp(['Number of iteration for smoothing is ',num2str(n)]);

%% AFISTA with continuation
tic
opts=[];
opts.rho=rho;
opts.acc=acc;
opts.maxiter=maxiter;
opts.print=1;
opts.x0=x_org;
opts.xplug=At(b);

[s,~,n,f_afista,E_afista]=AFISTA(b,A,At,D,Dt,lambda,[],[],opts,p,N0);
disp(['Time for AFISTA with continuation is']);
toc
disp(['Number of iteration for FISTA is ',num2str(n)]);

%% This part is developed based on NESTA method
tic
La=normM;
opts=[];
opts.mu=mu;
opts.TolVar=acc;
opts.U=D;
opts.Ut=Dt;
opts.normU=normD;
opts.maxiter=maxiter;
opts.print=1;
opts.x0=x_org;
opts.xplug=At(b);
%X_init=M'*b;

[xk,niter,residuals,outputData,E_out,opts] = MyCore_Nesterov_UP(A,At,b,lambda,La,mu,opts);
f_nesta=residuals(:,2);
disp(['Time for Smoothing method is']);
toc
disp(['Number of iteration for NESTA method is ',num2str(n)]);

%% plot results
figure(2)
imagesc(reshape(x_s*256,N,N));
colormap(gray)
axis image
axis equal
axis off
title(sprintf('reconstructed image using Smoothing (%d x %d = %d pixels)',N,N,N*N))

figure(3)
imagesc(reshape(s*256,N,N));
colormap(gray)
axis image
axis equal
axis off
title(sprintf('reconstructed image using Decomposition (%d x %d = %d pixels)',N,N,N*N))

error=norm(s-x)/norm(x);
disp(['The error of decomposition is ',num2str(error)]);
error_xs=norm(x_s-x)/norm(x);
disp(['The error of smoothing method is ',num2str(error_xs)]);
error_xk=norm(xk-x)/norm(x);
disp(['The error of NESTA method is ',num2str(error_xk)]);

%% plot the function value 
figure(5)
fntsz = 20; lwdth = 1; %display parameter
plot(f_afista(1:maxiter),'linewidth',lwdth);
hold on;
plot(f_s(1:maxiter),'r','linewidth',lwdth);
plot(f_nesta(1:maxiter),'g','linewidth',lwdth);
legend('fval decomposition','fval smoothing','fval nesta',fntsz);
xlabel('Iteration Number');
ylabel('Objective Function F(x)');

% analyze error with iteration
figure(6)
plot(E_afista,'linewidth',lwdth);
hold on;
plot(E_s,'r','linewidth',lwdth);
legend('Error decomposition','Error smoothing');