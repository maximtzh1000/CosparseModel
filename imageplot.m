clear; close all;clc;
load E_smooth;
fntsz = 18; lwdth = 1; %display parameter

E_decomp=E_decomp(2:20,2:20);
E_decomp=rot90(E_decomp);
E_smooth=E_smooth(2:20,2:20);
E_smooth=rot90(E_smooth);
portNumberx=0.1:0.1:1;
portNumbery=1:-0.1:0.1;

figure(1);
imagesc(1-E_decomp);
%title('Decomposition method','fontsize',fntsz);
xlabel('\alpha','fontsize',fntsz);
ylabel('\beta','fontsize',fntsz);
colormap('gray');
set(gca,'XTick',1:2:19);
set(gca,'XTicklabel',num2cell(portNumberx));
set(gca,'YTick',1:2:19);
set(gca,'YTicklabel',num2cell(portNumbery));

figure(2);
imagesc(1-E_smooth);
%title('Smoothing method','fontsize',fntsz);
xlabel('\alpha','fontsize',fntsz);
ylabel('\beta','fontsize',fntsz);
colormap('gray');
set(gca,'XTick',1:2:19);
set(gca,'XTicklabel',num2cell(portNumberx));
set(gca,'YTick',1:2:19);
set(gca,'YTicklabel',num2cell(portNumbery));