%try to solve the min norm(Dx,1)+norm(b-Mx,2)^2 in a method looks like
%FISTA for synthesis model
% This code is used for analysis compressve sensing model with data
% generated from Gaussian distribution. FISTA, FISTA method with
% continuation, FISTA method with projection and CVX method are implemented
% and AFISTA_con has the best result. 
% Improvements: 1.NESTA method 2. Change the numerical results with other
% models.

%To do list: Understand what is tight frame and generate random tight frame

clear;close all;clc;
addpath('/Users/zhaotan/Dropbox/CS_Analysis Model/Toyexample/NESTA_v1.1');
randomseed = RandStream('mcg16807','Seed',4);
RandStream.setGlobalStream(randomseed);

fntsz = 14; lwdth = 1; %display parameter

sigma=0.9;
r=0.10;
N=120;%dimension of the signal
L=floor(sigma*N);%number of measurements
N0=144;%number of dictionary
D=randn(N0,N)*1; %special case
%D_tmp=normalize(D');
[D_tmp,~]=qr(D);
D=D_tmp(:,1:N);%get the tight frame
M=randn(L,N);%measurement matrix
%[Q,R]=qr(M');
%M=Q(:,1:L)';
x=10*ones(N,1);
K=N-floor(r*L);%number of rows in D which are not orthogonal to x
ind=randsample(N0,K);%choose the index of these rows
[P_par,P_ort]=orthogonalspace(D(ind,:)');
x=P_ort*x;
x_org=x;
normcal(M,M',N,D,D');

Sigma=1*10^-5;
noise=Sigma*randn(L,1);
b=M*x+noise;%received signal
b=M*x;%received signal

SNR=20*log10(norm(b)/norm(noise));
lambda = Sigma*sqrt(2*log(N0));
%lambda=0.001;
mu=1*10^-9;
rho=lambda/mu;
acc=10^-30;
maxiter=20000;

%% This part is developed based on Smoothing Method
tic
opts=[];
opts.muf=mu;
opts.innerstep=10000;
opts.accf=acc;
opts.r=0.5;
opts.C=5;
opts.maxiter=maxiter;
opts.print=0;
opts.x0=x;
opts.xplug=zeros(N,1);

[x_s,n,f_s,E_s]=Smoothing_con(b,M,M',D,D',lambda,opts,N0,N);
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
opts.print=0;
opts.x0=x;
opts.xplug=zeros(N,1);

[s,n,f_afista,E_afista]=AFISTA_con(b,M,M',D,D',lambda,opts,N0,N);
disp(['Time for AFISTA with continuation is']);
toc
disp(['Number of iteration for FISTA is ',num2str(n)]);

%% ADMM method using CVX
tic
x_init=zeros(N,1);
rho=0.01;
[s_admm]=ADMM_CVX(b,M,D,lambda,x_init,rho,50,acc);
disp(['Time for ADMM is']);
toc
%% calculate error
error=norm(s-x)/norm(x);
disp(['The error of Decomposition method is ',num2str(error)]);
error_xs=norm(x_s-x)/norm(x);
disp(['The error of Smoothing method is ',num2str(error_xs)]);
error_admm=norm(s_admm-x)/norm(x);
disp(['The error of ADMM method is ',num2str(error_admm)]);

%% analyze error with iteration
figure(1)
semilogy(E_afista,'linewidth',lwdth);
hold on;
semilogy(E_s,'r','linewidth',lwdth);
legend('fval decomposition','fval smoothing','fontsize',fntsz);
