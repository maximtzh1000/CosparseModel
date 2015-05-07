% this is a test file for analysis dictionry learning.
% Author: Zhao Tan, Date: Mar 31th, 2014

clear;close all;clc;
randomseed = RandStream('mcg16807','Seed',1);
RandStream.setGlobalStream(randomseed);

sigma=0.9;
r=0.1;
N=20;%dimension of the signal
L=floor(sigma*N);%number of measurements
N0=40;%number of dictionary
D=randn(N0,N)*1;
x=10*ones(N,1);
K=N-floor(r*L);%number of rows in D which are not orthogonal to x
ind=randsample(N0,K);%choose the index of these rows
[P_par,P_ort]=orthogonalspace(D(ind,:)');
x=P_ort*x;
x_org=x;
