function [error_decomp,error_smooth]=randomtest(sigma,r)
%try to solve the min lambda*norm(Dx,1)+1/2norm(b-Mx,2)^2 in a method with
%decomposition and smoothing method
% This code is used for analysis compressve sensing model with data
% generated from Gaussian distribution.
% Input:
% two parameters controls number of measurement and also the level of
% cosparsity
% Output:
% the reconstruction error

%randomseed = RandStream('mcg16807','Seed',1);
%RandStream.setGlobalStream(randomseed);
N=120;%dimension of the signal
L=floor(sigma*N);%number of measurements
N0=144;%number of dictionary
D=randn(N0,N)*1; %special case
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
mu=1*10^-3;
rho=lambda/mu;
acc=10^-20;
maxiter=3000;

%% This part is developed based on Smoothing method
opts=[];
opts.mu=mu;
opts.acc=acc;
opts.maxiter=maxiter;
opts.print=0;
opts.x0=x_org;
X_init=zeros(N,1);
%X_init=M'*b;

[x_s,n,f_s,E_s]=SmoothingMethod(b,M,M',D,D',lambda,X_init,opts,N0,N);

%% This part is developed based on decomposition method
opts=[];
opts.rho=rho;
opts.acc=acc;
opts.maxiter=maxiter;
opts.print=0;
opts.x0=x_org;
X_init=zeros(N,1);
%X_init=M'*b;
%Z_init=D*X_init;
Z_init=zeros(N0,1);

[s,notused,n,f_afista,E_afista]=AFISTA(b,M,M',D,D',lambda,X_init,Z_init,opts,N0,N);

%% calculate error
error_decomp=norm(s-x)/norm(x);
error_smooth=norm(x_s-x)/norm(x);
