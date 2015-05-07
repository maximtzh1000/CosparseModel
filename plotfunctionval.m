%this is used to plot different function values for every iteration
%This demo use several images to show the 
%algorithms for analysis compressive sensing algorithm
%The M is taken as 2D discrete fourier transforms
%The D matrix is taken as 2D difference operator
%Samples are taken in a MRI fashion

%To do list: use 2D fourier transformation to form the measurement matrix

clear;clc;close all;
addpath('/Users/zhaotan/Dropbox/CS_Analysis Model/Toyexample/NESTA_v1.1');
randomseed = RandStream('mcg16807','Seed',2);
RandStream.setGlobalStream(randomseed);

fntsz = 14; lwdth = 1; %display parameter

%P = double(imread('pics/MRI_Coronal_Brain.jpg'));
%P = double(imread('pics/3DMR_Chest.jpg'));
%P = double(imread('pics/3DMR_Renal_Arteries.jpg'));
%P = double(imread('pics/cameraman.tif'));
%P = double(imread('pics/lena.jpg'));
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

pro=16; %number of projections
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
lambda = sigma*sqrt(2*log(p));
rho=10^7*sigma^2;
acc=10^-6;
maxiter=100;

%% This part is based on the method developed by Amir Beck in 2009
tic
pars.acc=acc;
pars.MAXITER=10000;
pars.maxiter=maxiter;
pars.print=1;
[x_den,n_den,f_den]=mfista_ab(b,A,At,D,Dt,lambda,pars,p,N0);
disp(['Time for dual method is']);
toc
disp(['Number of iteration for dual problem is ',num2str(n_den)]);

%% This part is developed based on NESTA paper using ALASSO with continuation
tic
muf = 10^-5;
La=my_normest(A,At,N*N,10^-4)^2;
opts = [];
opts.U = D;
opts.Ut = Dt;
opts.Verbose = 1;
opts.MaxIntIter=5;
opts.TolVar=10^-20;
opts.stopTest=1;
opts.maxiter=maxiter;

[xk,niter,residuals,outputData,opts] =Core_Nesterov_UP(A,At,b,lambda/2,La,muf,opts);
f_nesta=2*residuals(:,2);
%[xk] = NESTA(M,D,b,lambda,mu);
disp(['Time for NESTA with continuation is']);
toc
disp(['Number of iteration for NESTA is ',num2str(niter)]);

%% AFISTA with continuation
tic
opts=[];
opts.rho=rho;
opts.acc=acc;
opts.r=0.5;
opts.C=5;
opts.maxiter=maxiter;
opts.print=1;
[s,n,f_afista]=AFISTA(b,A,At,D,Dt,lambda,opts,p,N0);
disp(['Time for AFISTA with continuation is']);
toc
disp(['Number of iteration for FISTA is ',num2str(n)]);

%% plot results
% figure(2)
% imagesc(reshape(xk*256,N,N));
% colormap(gray)
% axis image
% axis equal
% axis off
% title(sprintf('reconstructed image using NESTA (%d x %d = %d pixels)',N,N,N*N))
% 
% figure(3)
% imagesc(reshape(s*256,N,N));
% colormap(gray)
% axis image
% axis equal
% axis off
% title(sprintf('reconstructed image using FISTA (%d x %d = %d pixels)',N,N,N*N))
% 
% figure(4)
% imagesc(reshape(x_den*256,N,N));
% colormap(gray)
% axis image
% axis equal
% axis off
% title(sprintf('reconstructed image using Amir Beck method (%d x %d = %d pixels)',N,N,N*N))

% error=norm(s-x)/norm(x);
% disp(['The error of AFISTA with continuation is ',num2str(error)]);
% error_xk=norm(xk-x)/norm(x);
% disp(['The error of NESTA with continuation is ',num2str(error_xk)]);
% error_den=norm(x_den-x)/norm(x);
% disp(['The error of Amir Beck is ',num2str(error_den)]);

%% plot the function value 
figure(5)
plot(f_afista(1:maxiter),'linewidth',lwdth);
hold on;
plot(f_nesta(1:maxiter),'r','linewidth',lwdth);
plot(f_den(1:maxiter),'g','linewidth',lwdth);
legend('fval fista','fval nesta','fval mfista','fontsize',fntsz);