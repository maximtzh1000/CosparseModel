%try to solve the min lambda*norm(Dx,1)+1/2norm(b-Mx,2)^2 in a method with
%decomposition and smoothing method
% This code is used for analysis compressve sensing model with data
% generated from Gaussian distribution.

clear;close all;clc;
%addpath('/Users/zhaotan/Dropbox/CS_Analysis Model/Toyexample/NESTA_v1.1');
randomseed = RandStream('mcg16807','Seed',1);
RandStream.setGlobalStream(randomseed);

fntsz = 14; lwdth = 1; %display parameter

sigma=0.9;
r=0.5;
N=120;%dimension of the signal
L=floor(sigma*N);%number of measurements
N0=144;%number of dictionary
D=randn(N0,N)*1; 
%D_tmp=normalize(D');
[D_tmp,noused]=qr(D);
D=D_tmp(:,1:N);%get the tight frame
M=randn(L,N);%measurement matrix
x=10*ones(N,1);
K=N-floor(r*L);%number of rows in D which are not orthogonal to x
ind=randsample(N0,K);%choose the index of these rows
[P_par,P_ort]=orthogonalspace(D(ind,:)');
x=P_ort*x;
x_org=x;
normcal(M,M',N,D,D');

%Sigma=1*10^-5;
%noise=Sigma*randn(L,1);
%b=M*x+noise;%received signal
b=M*x;%received signal

%SNR=20*log10(norm(b)/norm(noise));
%lambda = Sigma*sqrt(2*log(N0));
lambda=0.004;
mu=1*10^-5;
rho=lambda/mu;
acc=10^-30;
maxiter=3000;

%% This part is developed based on Smoothing method
tic
opts=[];
opts.mu=mu;
opts.acc=acc;
opts.maxiter=maxiter;
opts.print=1;
opts.x0=x_org;
X_init=zeros(N,1);
%X_init=M'*b;

[x_s,n,f_s,E_s]=SmoothingMethod(b,M,M',D,D',lambda,X_init,opts,N0,N);
disp(['Time for Smoothing method is']);
toc
disp(['Number of iteration for Smoothing method is ',num2str(n)]);

%% This part is developed based on decomposition method
tic
opts=[];
opts.rho=rho;
opts.acc=acc;
opts.maxiter=maxiter;
opts.print=1;
opts.x0=x_org;
X_init=zeros(N,1);
%X_init=M'*b;
%Z_init=D*X_init;
Z_init=zeros(N0,1);

[s,notused,n,f_afista,E_afista]=AFISTA(b,M,M',D,D',lambda,X_init,Z_init,opts,N0,N);
disp(['Time for Decomposition method is']);
toc
disp(['Number of iteration for Decomposition method is ',num2str(n)]);

%% This part is developed based on NESTA method
tic
La=norm(M)^2;
opts=[];
opts.mu=mu;
opts.TolVar=acc;
opts.U=D;
opts.Ut=D';
opts.maxiter=maxiter;
opts.print=1;
opts.x0=x_org;
opts.xplug=zeros(N,1);
%X_init=M'*b;

[xk,niter,residuals,outputData,E_nesta,opts] = MyCore_Nesterov_UP(M,M',b,lambda,La,mu,opts);
f_nesta=residuals(:,2);
disp(['Time for Smoothing method is']);
toc
disp(['Number of iteration for NESTA method is ',num2str(n)]);

%% calculate error
error=norm(s-x)/norm(x);
disp(['The error of Decomposition method is ',num2str(error)]);
error_xs=norm(x_s-x)/norm(x);
disp(['The error of Smoothing method is ',num2str(error_xs)]);
error_xk=norm(xk-x)/norm(x);
disp(['The error of NESTA method is ',num2str(error_xs)]);

%% plot the function value 
figure(1)
plot(f_afista,'linewidth',lwdth);
hold on;
plot(f_s,'r','linewidth',lwdth);
plot(f_nesta,'g','linewidth',lwdth);
legend('fval decomposition','fval smoothing','fval nesta','fontsize',fntsz);

% analyze error with iteration
figure(2)
plot(E_afista,'linewidth',lwdth);
hold on;
plot(E_s,'r','linewidth',lwdth);
plot(E_nesta,'g','linewidth',lwdth);
legend('Error decomposition','Error smoothing','Error nesta');
