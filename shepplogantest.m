function [f_smooth,f_decomp,E_smooth,E_decomp]=shepplogantest(x,b,A,At,D,Dt,lambda,p,N0,mu,maxiter)
%This demo use several images to show the 
%algorithms for analysis compressive sensing algorithm
%The M is taken as 2D discrete fourier transforms
%The D matrix is taken as 2D difference operator
%Samples are taken in a filter back projection fashion
%lambda = sigma*sqrt(2*log(p));
rho=lambda/mu;
%rho=10^6*lambda*sigma;
acc=10^-20;
global normM;
global normD;

%% This part is developed based on Smoothing method
tic
opts=[];
opts.mu=mu;
opts.acc=acc;
opts.maxiter=maxiter;
opts.print=1;
opts.x0=x;
opts.xplug=At(b);

[x_s,n,f_s,E_smooth]=SmoothingMethod(b,A,At,D,Dt,lambda,[],opts,p,N0);
disp(['Time for smoothing is']);
toc
disp(['Number of iteration for smoothing is ',num2str(n)]);

%% This part is developed based on Decomposition method
tic
opts=[];
opts.rho=rho;
opts.acc=acc;
opts.maxiter=maxiter;
opts.print=1;
opts.x0=x;
opts.xplug=At(b);

[s,~,n,f_afista,E_decomp]=AFISTA(b,A,At,D,Dt,lambda,[],[],opts,p,N0);
disp(['Time for decomposition is']);
toc
disp(['Number of iteration for decomposition is ',num2str(n)]);


e_decomp=norm(s-x)/norm(x);
e_smooth=norm(x_s-x)/norm(x);
f_decomp=f_afista(1:maxiter);
f_smooth=f_s(1:maxiter);
