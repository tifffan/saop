% This code requires the following packages:
% 1) GSPBox, available at https://github.com/epfl-lts2/gspbox 
% 2) Chebfun, available at http://www.chebfun.org/
% 3) CVX, available at http://cvxr.com/cvx/

% It may also be helpful after downloading the GSPBox to go into the 
% 3rdparty/LDL/LDL/MATLAB/ folder and execute ldl_make.m before running
% this code

% If graph legend text overlaps, try remove the subfolder in CVX package
% cvx/lib/narginchk_ from path or rename function narginchk.m 

close all;
clear all;
randn('seed', 18); rand('seed', 18)

% Graph
graph='sensor';
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
        load('/Users/lifan/Desktop/Research/git/spectral-warping/MATLAB_Code/Data/net25.mat');
        %load('/Users/davidshuman/Dropbox/Current Research Work/MCSFB/Shuni_Thesis/GitHub/mcsfb2018/net25 graph data/net25.mat');
        A=Problem.A;
        A = A - diag(diag(A)); 
        A(4228,6327) = 1;
        A(6327,4228) = 1;
        G=gsp_graph(A);
    otherwise
        error('graph type not recognized');
end

%cdf_approx
G=gsp_estimate_lmax(G);
lmax_est=G.lmax;
param.num_pts=50; % for approximating spectral cdf 
param.cdf_method='lanczos'; % can change to 'kpm' or 'lanczos'
param.num_vec=30;
param.order=30;
param.pts=linspace(.1,G.lmax,param.num_pts);
G=spectral_cdf_approx2(G,param); 


% Filters - testing mostly for analytic filters. If we have a
% discontinuity, we need more interpolation points near the discontinuity
filter='inverse';

switch filter
    case 'inverse'
        tau=.5; 
        f=@(x) tau./(tau+x);
    case 'ideal'
        f=@(x) (x<=(G.lmax/2));
    case 'heat'
        tau=1;
        f=@(x) exp(-tau*x);
    otherwise
        error('filter type not recognized');
end

% Poly approx order
K=10;

% Chebyshev
h=chebfun(@(s) f(s),[0,lmax_est],'splitting','on');
c=chebcoeffs(h,K+1); 
c(1)=c(1)*2;

% Legendre
pleg = polyfit(h,K);
plegc=pleg.coeffs;
plegc(1)=plegc(1)*2;

% Regular Lanczos
lanc_param.method='lanczos';
lanc_param.order=K;           % set K = larger constant for testing

% Warped
g=chebfun(@(s) lmax_est*G.spectrum_cdf_approx(s),[0,lmax_est],'splitting','on'); 
gi=inv(g,'splitting','on'); % warping function is the inverse of the spectral CDF

chebpts=cos((0:K)*pi/K)'; 
chebpts_tx=(chebpts+1)*lmax_est/2; 
chebpts_tx=sort(chebpts_tx,'ascend');
chebpts_warped=gi(chebpts_tx);

% Note: since G.spectrum_cdf_approx returns 0 at 0, the first
% interpolation point should always be 0 here. However, there may not be
% any points near the end of the spectrum. So as a hack here, we provide an
% option to force one point to be the estimate of the maximum eigenvalue
% for now

force_max=0;
if force_max
    chebpts_warped(end)=lmax_est;
    chebpts_warped=sort(chebpts_warped,'ascend');
end

dom=domain(0,lmax_est);
p_warped=chebfun.interp1(chebpts_warped,f(chebpts_warped),dom);
p_warped_c=p_warped.coeffs;
p_warped_c(1)=p_warped_c(1)*2;

figure;
Fbar=@(x)lmax_est*G.spectrum_cdf_approx(x);
fbarparam.plot_eigenvalues=0;
gsp_plot_filter(G,Fbar,fbarparam);
hold on;
plot(chebpts_warped,zeros(length(chebpts_warped),1),'xr','LineWidth',...
            4,'MarkerSize',15);
plot(zeros(length(chebpts_warped),1),chebpts_tx,'xr','LineWidth',...
            4,'MarkerSize',15);
%xlim([0,70]);
%ylim([0,70]);
set(gca,'FontSize',24);
%set(gca,'XTick',0:10:70);
%set(gca,'YTick',0:10:70);

% matrix adapted orthogonal polynomials
mop_param=struct;
num_vec=10;
num_its=1; % not currently used

[ab,absc,weights,Pi]= matrix_adapted_ortho_poly(G,K,num_vec,num_its,mop_param);
c_spec_adapted_ortho = matrix_adapted_poly_coeff(G, f, absc', weights', Pi, K);

% Least squares (assumes full knowledge of eigenvalues; just for comparison
% to ideal; not scalable)
G=gsp_compute_fourier_basis(G);
y=f(G.e);
lsc=polyfit(G.e,y,K);

% Compute polynomial approximation values at the actual eigenvalues and 
% the corresponding superior and squared error
y_cheb=gsp_cheby_eval(G.e,c,[0,lmax_est]);
errors_cheb=y-y_cheb;
sup_err_cheb=max(abs(errors_cheb));
se_cheb=sum(errors_cheb.^2);

y_leg=gsp_cheby_eval(G.e,plegc,[0,lmax_est]);
errors_leg=y-y_leg;
sup_err_leg=max(abs(errors_leg));
se_leg=sum(errors_leg.^2);

y_lanczos=G.U'*gsp_filter(G,f,sum(G.U)',lanc_param);   % but approximation is oscillating... fix this
errors_lanczos=y-y_lanczos;
sup_err_lanczos=max(abs(errors_lanczos));
se_lanczos=sum(errors_lanczos.^2);

figure;
plot(G.e, [f(G.e), y_lanczos])

y_cheb_warped=p_warped(G.e);
errors_cheb_warped=y-y_cheb_warped;
sup_err_cheb_warped=max(abs(errors_cheb_warped));
se_cheb_warped=sum(errors_cheb_warped.^2);

y_spec_adapted_ortho=three_term_eval(G,G.e,ab,c_spec_adapted_ortho);
errors_spec_adapted_ortho=y-y_spec_adapted_ortho;
sup_err_spec_adapted_ortho=max(abs(errors_spec_adapted_ortho));
se_spec_adapted_ortho=sum(errors_spec_adapted_ortho.^2);

y_ls=polyval(lsc,G.e);
errors_ls=y-y_ls;
sup_err_ls=max(abs(errors_ls));
se_ls=sum(errors_ls.^2);

sup_err_all = [sup_err_cheb, sup_err_leg, sup_err_lanczos, sup_err_cheb_warped, sup_err_spec_adapted_ortho, sup_err_ls]
se_all = [se_cheb, se_leg, se_lanczos, se_cheb_warped, se_spec_adapted_ortho, se_ls]

% Plots
xmax=max(lmax_est,G.lmax);
delta=xmax/5000;
xx=0:delta:xmax;
xx=xx';
yy=f(xx);

yy_cheb=gsp_cheby_eval(xx,c,[0,lmax_est]);
yy_leg=gsp_cheby_eval(xx,plegc,[0,lmax_est]);
yy_cheb_warped=p_warped(xx);
yy_spec_adapted_ortho=three_term_eval(G,xx,ab,c_spec_adapted_ortho);
yy_ls=polyval(lsc,xx);

% Figure 2: approximations for filter function
figure;
p1=plot(xx,[yy_cheb,pleg(xx),yy_cheb_warped,yy_spec_adapted_ortho,yy_ls,yy],'LineWidth',4);
set(gca,'FontSize',24)
legend(p1,'Chebyshev','Legendre','Warped Interpolation','Spectrum-Adapted Ortho. Poly.','Discrete LS','g','Location','NorthEast');
grid on;
hold on;
xlabel('\lambda');
ylabel('Filter approximations');
plot(G.e,zeros(G.N,1),'xk','LineWidth',2,'MarkerSize',6);

figure;
cc=lines(5);
p1=plot(xx,[yy_cheb_warped,yy],'LineWidth',4);
set(gca,'FontSize',24)
legend(p1,'Warped Interpolation','g','Location','NorthEast');
set(p1, {'color'}, {cc(3,:); cc(5,:)});
grid on;
hold on;
xlabel('\lambda');
plot(G.e,zeros(G.N,1),'xk','LineWidth',2,'MarkerSize',6);
plot(chebpts_warped,f(chebpts_warped),'xr','LineWidth',4,'MarkerSize',15);
%xlim([0,70]);
%set(gca,'XTick',0:10:70);


figure;
p2=semilogy(xx,[abs(yy-yy_cheb),abs(yy-yy_leg),abs(yy-yy_cheb_warped),abs(yy-yy_ls)],'LineWidth',4);
set(gca,'FontSize',24)
hold on;
plot(G.e,ones(G.N,1),'xk','LineWidth',2,'MarkerSize',6);
legend(p2,'Chebyshev','Legendre','Regular Lanczos','Warped Interpolation','Discrete LS','Location','SouthEast');  
xlabel('\lambda');
ylabel('$$|\tilde{g}(\lambda)-g(\lambda)|$$','Interpreter','latex');
grid on;


% Test f(L)b on random signal b
num_tests=5;

ntime_cheb=zeros(num_tests,1);
ntime_leg=zeros(num_tests,1);
ntime_lanczos=zeros(num_tests,1);
ntime_warped=zeros(num_tests,1);
ntime_spec=zeros(num_tests,1);
ntime_ls=zeros(num_tests,1);
ntime_exact=zeros(num_tests,1);

nmse_cheb=zeros(num_tests,1);
nmse_leg=zeros(num_tests,1);
nmse_lanczos=zeros(num_tests,1);
nmse_warped=zeros(num_tests,1);
nmse_spec=zeros(num_tests,1);
nmse_ls=zeros(num_tests,1);


G2=G;
G2.lmax=lmax_est;
G2=rmfield(G2,'U');
G2=rmfield(G2,'e');
lanc_param.method='lanczos';
lanc_param.order=K;


for i=1:num_tests
    b=rand(G.N,1);
    bhat=G.U'*b;
    tic
    gLf_exact=G.U*(f(G.e).*bhat);
    ntime_exact(i)=toc;     % add time to compute fourier basis
    tic
    gLf_ls=G.U*(polyval(lsc,G.e).*bhat);    
    ntime_ls(i)=toc;        % add time to compute fourier basis
    tic
    gLf_cheb=gsp_cheby_op(G2,c,b); 
    ntime_cheb(i)=toc;
    tic
    gLf_leg=gsp_cheby_op(G2,plegc,b); 
    ntime_leg(i)=toc;
    tic
    gLf_lanczos=gsp_filter(G2,f,b,lanc_param);
    ntime_lanczos(i)=toc;
    tic
    gLf_warped=gsp_cheby_op(G2,p_warped_c,b);
    ntime_warped(i)=toc;
    %tic
    %gLf_spec=gsp_cheby_op(G2,c_spec_adapted_ortho,b);
    %ntime_spec(i)=toc;

    nmse_cheb(i)=sum((gLf_cheb-gLf_exact).^2)/sum(gLf_exact.^2);
    nmse_leg(i)=sum((gLf_leg-gLf_exact).^2)/sum(gLf_exact.^2);
    nmse_lanczos(i)=sum((gLf_lanczos-gLf_exact).^2)/sum(gLf_exact.^2);
    nmse_warped(i)=sum((gLf_warped-gLf_exact).^2)/sum(gLf_exact.^2);
    %nmse_spec(i)=sum((gLf_spec-gLf_exact).^2)/sum(gLf_exact.^2);
    nmse_ls(i)=sum((gLf_ls-gLf_exact).^2)/sum(gLf_exact.^2);
    
end

avg_ntime_cheb=mean(ntime_cheb);
avg_ntime_leg=mean(ntime_leg);
avg_ntime_lanczos=mean(ntime_lanczos);
avg_ntime_warped=mean(ntime_warped);
%avg_ntime_spec=mean(ntime_spec);
avg_ntime_ls=mean(ntime_ls);
avg_ntime_exact=mean(ntime_exact);

avg_nmse_cheb=mean(nmse_cheb);
avg_nmse_leg=mean(nmse_leg);
avg_nmse_lanczos=mean(nmse_lanczos);
avg_nmse_warped=mean(nmse_warped);
%avg_nmse_spec=mean(nmse_spec);
avg_nmse_ls=mean(nmse_ls);


%avg_time_all = [avg_ntime_cheb, avg_ntime_leg, avg_ntime_lanczos, avg_ntime_warped, avg_time_spec, avg_ntime_ls, avg_ntime_exact]
%avg_nmse_ls = [avg_nmse_cheb, avg_nmse_leg, avg_nmse_lanczos, avg_nmse_warped, avg_nmse_spec, avg_nmse_ls, avg_nmse_exact]
