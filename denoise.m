%This demo use several images to show the 
%algorithms for denoising. But l1 norm might not be a good choice here
%The M is taken as 2D discrete fourier transforms
%The D matrix is taken as 2D difference operator

%To do list: use 2D fourier transformation to form the measurement matrix

clear;clc;close all;
addpath('/Users/zhaotan/Dropbox/CS_Analysis Model/Toyexample/NESTA_v1.1');

%P = double(imread('pics/MRI_Coronal_Brain.jpg'));
%P = double(imread('pics/3DMR_Chest.jpg'));
%P = double(imread('pics/3DMR_Renal_Arteries.jpg'));
P = double(imread('pics/cameraman.tif'));
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
sigma=2*10^-2;

A=@(x) A_eye(x);
At=@(y) A_eye(y);

D=@(x) D_func(x,N,N);
Dt=@(p) Dt_func(p,N,N);

b=A(x);
L=length(b);
noise=sigma*randn(L,1);
b=b+noise;

SNR=20*log10(norm(b)/norm(noise));
lambda = sigma*sqrt(2*log(p));
%lambda=0.02;
rho=10^8*2*sigma^2;
acc=10^-6;

figure(2)
imagesc(reshape(b*256,N,N));
colormap(gray)
axis image
axis equal
axis off
title(sprintf('The noisy image (%d x %d = %d pixels)',N,N,N*N))

%% This part is based on the method developed by Amir Beck in 2009
% tic
% pars.acc=acc;
% pars.MAXITER=100;
% pars.print=0;
% [x_den,n_den]=mfista_ab(b,A,At,D,Dt,lambda,pars,p,N0);
% disp(['Time for dual method is']);
% toc
% disp(['Number of iteration for dual problem is ',num2str(n_den)]);
% 
% %% This part is developed based on NESTA paper using ALASSO with continuation
% tic
% muf = 10^-8;
% La=my_normest(A,At,N*N,10^-4)^2;
% opts = [];
% opts.U = D;
% opts.Ut = Dt;
% opts.Verbose = 0;
% opts.MaxIntIter=5;
% opts.TolVar=10^-10;
% opts.stopTest=1;
% [xk,niter,residuals,outputData,opts] =NESTA_UP(A,At,b,lambda,La,muf,opts);
% %[xk] = NESTA(M,D,b,lambda,mu);
% disp(['Time for NESTA with continuation is']);
% toc
% disp(['Number of iteration for NESTA is ',num2str(niter)]);

%% AFISTA with continuation
tic
opts=[];
opts.rhof=rho;
opts.accf=acc;
opts.r=0.5;
opts.C=5;
opts.maxiter=10000;
[s,n]=AFISTA_con(b,A,At,D,Dt,lambda,opts,p,N0);
disp(['Time for AFISTA with continuation is']);
toc
disp(['Number of iteration for FISTA is ',num2str(n)]);

%% plot results
figure(3)
imagesc(reshape(xk*256,N,N));
colormap(gray)
axis image
axis equal
axis off
title(sprintf('reconstructed image using NESTA (%d x %d = %d pixels)',N,N,N*N))

figure(4)
imagesc(reshape(s*256,N,N));
colormap(gray)
axis image
axis equal
axis off
title(sprintf('reconstructed image using FISTA (%d x %d = %d pixels)',N,N,N*N))

figure(5)
imagesc(reshape(x_den*256,N,N));
colormap(gray)
axis image
axis equal
axis off
title(sprintf('reconstructed image using Amir Beck method (%d x %d = %d pixels)',N,N,N*N))

error=norm(s-x)/norm(x);
disp(['The error of AFISTA with continuation is ',num2str(error)]);
error_xk=norm(xk-x)/norm(x);
disp(['The error of NESTA with continuation is ',num2str(error_xk)]);
error_den=norm(x_den-x)/norm(x);
disp(['The error of Amir Beck is ',num2str(error_den)]);
error_noise=norm(b-x)/norm(x);
disp(['The power of orginal noise is ',num2str(error_noise)]);