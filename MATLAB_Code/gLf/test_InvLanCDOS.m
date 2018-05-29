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
param.num_pts=20; % for approximating spectral cdf 
G=gsp_spectrum_cdf_approx(G,param); % incorporate alternative ways to compute this (e.g., KPM for DoS)

%% CDOS by Lanczos approximation

Mdeg = 50; % Number of Lanczos Steps
nvec =10; % Number of starting vectors 
x=0:.5:G.lmax; % discretized points in the interval [lmin,lMax]

K=10;


use_Chebypts=1;
if(use_Chebypts)
chebpts=cos((0:K)*pi/K); 
tx=(chebpts+1)/2; 
else
tx=linspace(0, 1, K+1);
end


tic
 InLanDOS= InvLanczosCDOS(G.L, Mdeg, nvec,G.lmax,x',tx');
 lantime=toc
 
 tic
g=chebfun(@(s) G.spectrum_cdf_approx(s),[0,G.lmax],'splitting','on'); 
gi=inv(g,'splitting','on'); % warping function is the inverse of the spectral CDF
InLDLDOS=gi(tx);
 LDLTtime=toc

 plot(tx,InLDLDOS,'b','LineWidth',2);
 hold on;
 plot(tx,InLanDOS,'r','LineWidth',2);
 legend('LDLT','Lanczos');
titl=sprintf('Approximate CDOS Lanczos m=%d steps',Mdeg);
title(titl,'FontSize',16)
 