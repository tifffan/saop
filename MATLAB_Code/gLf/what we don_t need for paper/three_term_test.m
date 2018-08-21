clc
clear all;
close all;
randn('seed', 1); 
rand('seed', 1);

%% Create graph/matrix
graph='comet';
switch graph
    case 'gnp'
        N=500;
        p=.2;
        G=gsp_erdos_renyi(N,p);
    case 'sensor'
        N=500;
        G=gsp_david_sensor_network(N);
    case 'minnesota'
        G=gsp_minnesota(1);
    case 'comet'
        N=500;
        G=gsp_comet(N,50);
    case 'community'
        N=5000;
        G=gsp_community(N);
    case 'net25'
        load('/Users/davidshuman/Dropbox/Current Research Work/MCSFB/Shuni_Thesis/GitHub/mcsfb2018/net25 graph data/net25.mat');
        A=Problem.A;
        A = A - diag(diag(A)); 
        A(4228,6327) = 1;
        A(6327,4228) = 1;
        G=gsp_graph(A);
    otherwise
        error('graph type not recognized');
end
G = gsp_compute_fourier_basis(G);

%% signal to filter
mysig = rand(G.N,1);
h=@(x) 1./(1+x);
smooth_exact=gsp_filter(G,h,mysig);

%% matrix adapted orthogonal polynomials
order=10;
param=struct;
num_vec=20;
num_its=1; % not currently used
Mdeg=30;

%%---Weights using Chebyshev filters and stochastic trace estimator
tic
[ab,absc,weights,Pi]= matrix_adapted_ortho_poly(G,order,num_vec,num_its,param);
c = matrix_adapted_poly_coeff(G, h, absc', weights', Pi, order);
smooth_warped=three_term_recurr_op(G, ab, c, mysig,param);
warped_time=toc


%%-- Weights via. Lanczos CDF
tic
[ab_lan,absc_lan,weights_lan,Pi_lan]= matrix_adapted_ortho_poly_Lan(G,Mdeg,num_vec,order);
c_lan = matrix_adapted_poly_coeff(G, h, absc_lan', weights_lan', Pi_lan, order);
smooth_warped_Lan=three_term_recurr_op(G, ab_lan, c_lan, mysig,param);
Lanwarped_time=toc

% regular shifted Chebyshev filter
G=rmfield(G,'U');
G=rmfield(G,'e');
param.order=order;
tic
smooth_cheb=gsp_filter(G,h,mysig,param);
cheb_time=toc
% plot pdf
% figure;
% plot(absc,weights,'LineWidth',2);

%% plot cdf and compare to ldlt/mono cubic method
figure;
 plot(absc,cumsum(weights),'LineWidth',2);
hold on;
 plot(absc,weights_lan,'r','LineWidth',2);
cdf_param.num_pts=order+1;
tic
G=gsp_spectrum_cdf_approx(G,cdf_param);
CDtime=toc
hold on;
xx=0:.01:G.lmax;
plot(xx,G.spectrum_cdf_approx(xx),'g','LineWidth',2);

%% Errors
error_warped=max(abs(smooth_warped-smooth_exact))
error_Lanwarped=max(abs(smooth_warped_Lan-smooth_exact))
error_cheb=max(abs(smooth_cheb-smooth_exact))

nmse_cheb=sum((smooth_cheb-smooth_exact).^2)/sum(smooth_exact.^2)
nmse_warped=sum((smooth_warped-smooth_exact).^2)/sum(smooth_exact.^2)
nmse_Lanwarped=sum((smooth_warped_Lan-smooth_exact).^2)/sum(smooth_exact.^2)
   
