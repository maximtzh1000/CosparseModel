%this matlab file tests the alternating project idea of solving generalized
%LASSO problem using the dual form

clear;close all;clc;
%addpath('/Users/zhaotan/Dropbox/CS_Analysis Model/Toyexample/NESTA_v1.1');
randomseed = RandStream('mcg16807','Seed',1);
RandStream.setGlobalStream(randomseed);

fntsz = 14; lwdth = 1; %display parameter

sigma=0.9;
r=0.9;
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
ind1=randsample(N0,K);%choose the index of these rows
[P_par,P_ort]=orthogonalspace(D(ind1,:)');
x=P_ort*x;

%% noiseless case
b=M*x;%received signal
lambda=0.05;
mu=1*10^-8;
rho=lambda/mu;
acc=10^-30;
maxiter=10000;

%% Using FISTA to solve the dual problem
opts=[];
opts.rho=rho;
opts.acc=acc;
opts.maxiter=maxiter;
opts.print=1;
U_init=b;
V_init=zeros(N0,1);

[u,v,niter]=Decomp_dual(b,M,M',D,D',lambda,U_init,V_init,opts,N0,N);

%% from the dual solution get the primal solution
K=N-L; %how many more information we need
[~,ind]=sort(abs(v));
ind=ind(1:K);
D_tmp=D(ind1,:);
b0=[b-u;zeros(length(ind1),1)];
A=[M;D_tmp];
s=A\b0;
norm(s-x)/norm(x)
