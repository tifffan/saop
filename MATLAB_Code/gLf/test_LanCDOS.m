clc
clear all;
close all;
randn('seed', 1); 
rand('seed', 1);


%% Load data
N=500;
p=.3;
G=gsp_erdos_renyi(N,p);
%G = gsp_david_sensor_network(500);
G = gsp_compute_fourier_basis(G);


%% CDOS by Lanczos approximation

Mdeg = 50; % Number of Lanczos Steps
nvec =10; % Number of starting vectors 
x=0:.05:G.lmax; % discretized points in the interval [lmin,lMax]
% x=linspace(0, G.lmax, 20);
tic
 LanDOS= LanczosCDOS(G.L, Mdeg, nvec,x);
 lantime=toc
 
 
 %% CDOS by LDLT
tic
cdf_param.num_pts=Mdeg;
G=gsp_spectrum_cdf_approx(G,cdf_param);
LDLtime=toc

LDLDOS=G.spectrum_cdf_approx(x');

%% plot results
 plot(x,LDLDOS,'b','LineWidth',2);
 hold on;
 plot(x,LanDOS,'r','LineWidth',2);
 legend('LDLT','Lanczos');
titl=sprintf('Approximate CDOS Lanczos m=%d steps',Mdeg);
title(titl,'FontSize',16)
 
%  
% %% Inverse CDOS
% xi=0:0.05:1;
% g=chebfun(@(s) G.lmax*LanczosCDOS(G.L, Mdeg, nvec,s),[0,G.lmax],'splitting','on'); 
% invg=inv(g,'splitting','on'); 
% 
% InvLanDOS = invg(xi);
%   
 

