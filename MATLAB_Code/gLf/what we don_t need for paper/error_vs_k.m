close all;
clear all;
randn('seed', 18); rand('seed', 18)

% Graph
graph='saylr4';
switch graph
    case 'gnp'
        N=500;
        p=.2;
        G=gsp_erdos_renyi(N,p);
        %G.L=G.L/10;
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
        %G.L=G.L/10;
    case 'net25'
        load('/Users/lfan/Documents/MATLAB/git/spectral-warping/MATLAB_Code/Data/net25.mat');
        %load('/Users/lifan/Desktop/Research/git/spectral-warping/MATLAB_Code/Data/net25.mat');
        %load('/Users/davidshuman/Dropbox/Current Research Work/MCSFB/Shuni_Thesis/GitHub/mcsfb2018/net25 graph data/net25.mat');
        A=Problem.A;
        A = A - diag(diag(A)); 
        A(4228,6327) = 1;
        A(6327,4228) = 1;
        G=gsp_graph(A);
        %G.L=G.L./20;
    case 'saylr4'
        load('/Users/lfan/Documents/MATLAB/git/spectral-warping/MATLAB_Code/Data/saylr4.mat');
        %load('/Users/davidshuman/Dropbox/Current Research Work/...');
        A=Problem.A;
        A = A - diag(diag(A)); 
        G=gsp_graph(A);
        %G.L=G.L/2000;
    case 'gr3030'   % positive definite, G.lmin=0.06
        load('/Users/lfan/Documents/MATLAB/git/spectral-warping/MATLAB_Code/Data/gr_30_30.mat');
        %load('/Users/davidshuman/Dropbox/Current Research Work/...');
        G=struct;
        G.L=Problem.A;
        G.N=size(G.L,1);
    case 'cdde1'   % original data not symmetric
        load('/Users/lfan/Documents/MATLAB/git/spectral-warping/MATLAB_Code/Data/cdde1.mat');
        %load('/Users/davidshuman/Dropbox/Current Research Work/...');
        G=struct;
        G.L=Problem.A;
        G.L=(G.L+G.L')/2;
        G.N=size(G.L,1);
   case 'si2' % positive definite
        load('/Users/lfan/Documents/MATLAB/git/spectral-warping/MATLAB_Code/Data/Si2.mat');
        %load('/Users/davidshuman/Dropbox/Current Research Work/...');
        G=struct;
        G.L=Problem.A;
        G.N=size(G.L,1);
    case 'nd3k'
        load('/Users/lfan/Documents/MATLAB/git/spectral-warping/MATLAB_Code/Data/nd3k.mat');
        %load('/Users/davidshuman/Dropbox/Current Research Work/...');
        G=struct;
        G.N=9000;
        G.L=Problem.A;
        G.L=G.L./20;
    case 'bipt100'
        N=100;
        K=N/4;
        weights=rand(K,K);
        W=[zeros(K,K),weights,zeros(K,K),zeros(K,K);weights',zeros(K,K),weights,zeros(K,K);...
         zeros(K,K),weights',zeros(K,K),weights;zeros(K,K),zeros(K,K),weights',zeros(K,K)];
        num_delete_edge=1000;
        ind1=randi(K,num_delete_edge);
        ind2=randi(K,num_delete_edge)+K;
        for i=1:num_delete_edge
            W(ind1(i),ind2(i))=0;
            W(ind2(i),ind1(i))=0;
        end
        G=gsp_graph(W);
        G.L=G.L/16;
    case 'bipt120'
        N=120;
        K=N/6;
        weights=rand(K,K);
        W=[zeros(K,K),weights,zeros(K,K),zeros(K,K),zeros(K,K),zeros(K,K);weights',zeros(K,K),weights,zeros(K,K),zeros(K,K),zeros(K,K);...
            zeros(K,K),weights',zeros(K,K),weights,zeros(K,K),zeros(K,K);zeros(K,K),zeros(K,K),weights',zeros(K,K),weights,zeros(K,K);...
            zeros(K,K),zeros(K,K),zeros(K,K),weights',zeros(K,K),weights; zeros(K,K),zeros(K,K),zeros(K,K), zeros(K,K),weights',zeros(K,K)];
        num_delete_edge=1000;
        ind1=randi(K,num_delete_edge);
        ind2=randi(K,num_delete_edge)+K;
        for i=1:num_delete_edge
            W(ind1(i),ind2(i))=0;
            W(ind2(i),ind1(i))=0;
        end
        G=gsp_graph(W);
        G.L=G.L/16;
    case 'spd'
        G=struct;
        G.N=512;
        e=0.5*randn(G.N,1)+4;
        G.e=sort(e,'ascend');
        H=hadamard(G.N)/sqrt(G.N);
        G.L=H*diag(G.e)*H';
    otherwise
        error('graph type not recognized');
end

%s=2;
%G.L=G.L+s*speye(G.N);

G2=G;
tic
G2=gsp_compute_fourier_basis(G2);
time_exact_spec=toc;


G=gsp_estimate_lmax(G);
%G.lmin=s;
G.lmin=G2.e(1)*(G2.e(1)>1e-6);
param.cdf_method='kpm'; % can change to 'kpm' or 'lanczos' or 'ldlt'
param.num_pts=10; % for kpm and ldlt
param.num_vec=10; %for lanczos only
param.order=30; %100 works for kpm and ldlt
G=spectral_cdf_approx2(G,param); 
gi=@(s) G.spectrum_inv_cdf_approx((s-G.lmin)/(G.lmax-G.lmin));

% Filters - testing mostly for analytic filters. If we have a
% discontinuity, we need more interpolation points near the discontinuity
filter='heat';

switch filter
    case 'inverse'
        tau=.5; 
        f=@(x) tau./(tau+x);
    case 'low'
        f=@(x) (x<=gi(G.lmax/2));
    case 'mid'
        f=@(x) (x>=gi(G.lmax/3) & x<=gi(G.lmax*2/3));
    case 'heat'
        tau=1;
        f=@(x) exp(-tau*x);
        %f=@(x) exp(-tau*10*x);
    case 'exp'
        tau=1;
        f=@(x) exp(tau*x);
        %f=@(x) exp(tau*10*x);
    case 'log'
        f=@(x) log(x);
    otherwise
        error('filter type not recognized');
end

start_pts_interp='cheb'; % K pts for warped chebyshev interpolation
start_pts_ls='cheb'; % G.N/10 pts for warped LS fitting

%--------------------------------------------------------------------------
% Poly approx order
n_its=13;
poly_orders=[1:n_its]'+2;
all_plots=0; %plot approximated f and errors for every iteration

cheb_sup_err=zeros(n_its,1);
leg_sup_err=zeros(n_its,1);
lanc_sup_err=zeros(n_its,1);
warp_sup_err=zeros(n_its,1);
warp_ls_sup_err=zeros(n_its,1);
spec_sup_err=zeros(n_its,1);
spec_warp_sup_err=zeros(n_its,1);
spec_spline_sup_err=zeros(n_its,1);
spec_wch_sup_err=zeros(n_its,1);
weighted_ls_sup_err=zeros(n_its,1);

cheb_sq_err=zeros(n_its,1);
leg_sq_err=zeros(n_its,1);
lanc_sq_err=zeros(n_its,1);
warp_sq_err=zeros(n_its,1);
warp_ls_sq_err=zeros(n_its,1);
spec_sq_err=zeros(n_its,1);
spec_warp_sq_err=zeros(n_its,1);
spec_spline_sq_err=zeros(n_its,1);
spec_wch_sq_err=zeros(n_its,1);
weighted_ls_sq_err=zeros(n_its,1);

cheb_time=zeros(n_its,1);
leg_time=zeros(n_its,1);
lanc_time=zeros(n_its,1);
warp_time=zeros(n_its,1);
warp_ls_time=zeros(n_its,1);
spec_time=zeros(n_its,1);
spec_warp_time=zeros(n_its,1);
spec_spline_time=zeros(n_its,1);
spec_wch_time=zeros(n_its,1);
weighted_ls_time=zeros(n_its,1);

cheb_anmse=zeros(n_its,1);
leg_anmse=zeros(n_its,1);
lanc_anmse=zeros(n_its,1);
warp_anmse=zeros(n_its,1);
warp_ls_anmse=zeros(n_its,1);
spec_anmse=zeros(n_its,1);
spec_warp_anmse=zeros(n_its,1);
spec_spline_anmse=zeros(n_its,1);
spec_wch_anmse=zeros(n_its,1);
weighted_ls_anmse=zeros(n_its,1);
weighted_ls_2_anmse=zeros(n_its,1);
weighted_ls_spline_anmse=zeros(n_its,1);



for K_ind=1:n_its
K=poly_orders(K_ind);

% Chebyshev
c=gsp_cheby_coeff(G,f,K,1000);

% Legendre
h=chebfun(@(s) f(s),[G.lmin,G.lmax],'splitting','on');
pleg = polyfit(h,K);
plegc=pleg.coeffs;
plegc(1)=plegc(1)*2;

% Lanczos
lanc_param.method='lanczos';
lanc_param.order=K;
x=sum(G2.U')';
[V,H]=lanczos(G.L,K+1,x);
e=eig(H);

% Starting points for interpolation
switch start_pts_interp
    case 'cheb'
        pts_interp=cos((0:K)*pi/K)'; 
        pts_tx_interp=(pts_interp+1)*(G.lmax-G.lmin)/2+G.lmin;
        pts_tx_interp=sort(pts_tx_interp,'ascend');
    case 'even'
        pts_tx_interp=linspace(G.lmin,G.lmax,K+1);
end

pts_warped_interp=gi(pts_tx_interp);
pts_warped_interp(1)=G.lmin;

% Starting points for LS fitting
npts_ls=500;
%npts_ls=max(ceil(G.N/10),K+1);
switch start_pts_ls
    case 'cheb'
        pts_ls=cos((0:npts_ls)*pi/npts_ls)';
        pts_tx_ls=(pts_ls+1)*(G.lmax-G.lmin)/2+G.lmin;
        pts_tx_ls=sort(pts_tx_ls,'ascend');
    case 'even'
        pts_tx_ls=linspace(G.lmin,G.lmax,npts_ls);
    case 'unif'
        pts_tx_ls=(G.lmax-G.lmin)*rand(npts_ls,1)+G.lmin;
        pts_tx_ls=sort(pts_tx_ls,'ascend');
    otherwise
        error('grid type not recognized');
end

pts_warped_ls=gi(pts_tx_ls);
pts_warped_ls(1)=G.lmin;

% Note: since G.spectrum_cdf_approx returns 0 at 0, the first
% interpolation point should always be 0 here. However, there may not be
% any points near the end of the spectrum. So as a hack here, we provide an
% option to force one point to be the estimate of the maximum eigenvalue
% for now

force_max=0;
if force_max
    pts_warped_interp(end)=G.lmax;
    pts_warped_ls(end)=G.lmax;
    pts_warped_interp=sort(pts_warped_interp,'ascend');
    pts_warped_ls=sort(pts_warped_ls,'ascend');
end

% Chebyshev interpolation on warped points
[p_warped,s_warped,mu_warped]=polyfit(pts_warped_interp,f(pts_warped_interp),K);
p_warped_fun=@(x) polyval(p_warped,x,s_warped,mu_warped); 
p_warped_c=gsp_cheby_coeff(G,p_warped_fun,K,1000);

% % Piecewise cubic hermite polynomial interpolation + Chebyshev approximation
% p_warped_pchip=pchip(pts_warped_interp,f(pts_warped_interp));
% p_warped_pchip_fun=@(x) pchip(pts_warped_interp,f(pts_warped_interp),x); 
% p_warped_pchip_interp_c=gsp_cheby_coeff(G,p_warped_pchip_fun,K,1000);

% Many points + LS fitting
[pts_warped_ls,weights_warped_ls]=remove_repeat(pts_warped_ls,0.01);
[p_warped_lsc,s_warped_ls,mu_warped_ls]=weighted_polyfit(pts_warped_ls,f(pts_warped_ls),weights_warped_ls,K);

% Include damping factors against Gibbs oscillations
damping='none';
switch damping
    case 'jackson'
        gamma=ones(K+1,1);
        for k=1:K
            gamma(k+1,1)=((1-k/(K+2))*sin(pi/(K+2))*cos(k*pi/(K+2))...
                    +1/(K+2)*cos(pi/(K+2))*sin(k*pi/(K+2)))/sin(pi/(K+2));
            % Damping factors calculated from eqn(12), M-CSFB; has stronger 
            % damping effect on large degree polynomials than sigma factors
        end
        p_warped_c=p_warped_c.*gamma;
        %p_warped_pchip_interp_c=p_warped_pchip_interp_c.*gamma;
    case  'sigma'
        sigma=ones(K+1,1);
        for k=1:K
            sigma(k+1,1)=sin(k*pi/(K+1))/(k*pi/(K+1));
            % Damping factors calculated from 3.1, Napoli et.al, Efficient 
            % estimation of eigenvalue counts in an interval
        end
        p_warped_c=p_warped_c.*sigma;
        %p_warped_pchip_interp_c=p_warped_pchip_interp_c.*sigma;
    case 'none'
    otherwise
        error('damping type not recognized');
end

% % Figure 1: warped points on approximated CDF
% figure;
% Fbar=@(x)G.lmax*G.spectrum_cdf_approx(x);
% fbarparam.plot_eigenvalues=0;
% gsp_plot_filter(G,Fbar,fbarparam);
% hold on;
% plot(pts_warped_interp,zeros(length(pts_warped_interp),1),'xr','LineWidth',...
%             4,'MarkerSize',15);
% plot(zeros(length(pts_warped_interp),1),pts_tx_interp,'xr','LineWidth',...
%             4,'MarkerSize',15);   
% set(gca,'FontSize',24);


% Matrix/Spectrum adapted orthogonal polynomials
mop_param=struct;
mop_param.grid_order=100;
mop_param.init_poly_order=30;
mop_param.num_vec=10;
mop_param.absc_method='linear'; %linear, even(n-2 linearly spaced + both ends), warp, spline
mop_param.weights_method='pdf'; %count, pdf2(from kpm paper)
[absc,weights]=gen_absc_weights(G,K,mop_param); 
[ab,Pi]=matrix_adapted_ortho_poly(absc,weights,K);
c_spec_adapted_ortho=matrix_adapted_poly_coeff(G, f, absc, weights/sum(weights), Pi, K);

mop_param.grid_order=100;
mop_param.absc_method='warp';
mop_param.weights_method='count';
[absc_warped,weights_warped]=gen_absc_weights(G,K,mop_param);
% [absc_warped,~,ind_warped]=remove_repeat(absc_warped,0.01,weights_warped);
% weights_warped=weights_warped(ind_warped);
[ab_warped,Pi_warped]=matrix_adapted_ortho_poly(absc_warped,weights_warped,K);
c_spec_adapted_ortho_warped=matrix_adapted_poly_coeff(G, f, absc_warped, weights_warped/sum(weights_warped), Pi_warped, K);

mop_param.grid_order=100;
mop_param.absc_method='warp';
mop_param.weights_method='constant';
[absc_wch,weights_wch]=gen_absc_weights(G,K,mop_param);
[absc_wch,weights_wch]=remove_repeat(absc_wch,0.01,weights_wch);
[ab_wch,Pi_wch]=matrix_adapted_ortho_poly(absc_wch,weights_wch*1000,K);
c_spec_adapted_ortho_wch=matrix_adapted_poly_coeff(G, f, absc_wch, weights_wch/sum(weights_wch), Pi_wch, K);

mop_param.grid_order=100;
mop_param.absc_method='spline';
mop_param.weights_method='pdf';
[absc_spline,weights_spline]=gen_absc_weights(G,K,mop_param);
[ab_spline,Pi_spline]=matrix_adapted_ortho_poly(absc_spline,weights_spline,K);
c_spec_adapted_ortho_spline=matrix_adapted_poly_coeff(G, f, absc_spline, weights_spline/sum(weights_spline), Pi_spline, K);

% Weighted LS fitting
mop_param.grid_order=100;
mop_param.absc_method='linear';
mop_param.weights_method='pdf';
[absc_weighted_ls,weights_weighted_ls]=gen_absc_weights(G,K,mop_param);
[p_weighted_lsc,s_weighted_ls,mu_weighted_ls]=weighted_polyfit(absc_weighted_ls,f(absc_weighted_ls),weights_weighted_ls,K);
% mat_absc_weighted_ls=absc_weighted_ls.^(0:K);
% p_weighted_lsc=lscov(mat_absc_weighted_ls,f(absc_weighted_ls),G.spectrum_pdf_approx(absc_weighted_ls));
% p_weighted_lsc=flipud(p_weighted_lsc);

mop_param.grid_order=100;
mop_param.absc_method='linear';
mop_param.weights_method='count';
[absc_linear,weights_linear]=gen_absc_weights(G,K,mop_param);
[p_weighted_lsc_2,s_weighted_ls_2,mu_weighted_ls_2]=weighted_polyfit(absc_linear,f(absc_linear),weights_linear,K);

% Least squares (assumes full knowledge of eigenvalues; just for comparison
% to ideal; not scalable)
y=f(G2.e);
[lsc,s_ls,mu_ls]=polyfit(G2.e,y,K);

% Compute polynomial approximation values at the actual eigenvalues and 
% the corresponding superior and squared error
y_cheb=gsp_cheby_eval(G2.e,c,[G.lmin,G.lmax]);
errors_cheb=y-y_cheb;
sup_err_cheb=max(abs(errors_cheb));
se_cheb=sum(errors_cheb.^2);

y_leg=gsp_cheby_eval(G2.e,plegc,[G.lmin,G.lmax]);
errors_leg=y-y_leg;
sup_err_leg=max(abs(errors_leg));
se_leg=sum(errors_leg.^2);

y_lanczos=G2.U'*gsp_filter(G,f,sum(G2.U')',lanc_param);
errors_lanczos=y-y_lanczos;
sup_err_lanczos=max(abs(errors_lanczos));
se_lanczos=sum(errors_lanczos.^2);

y_cheb_warped=gsp_cheby_eval(G2.e,p_warped_c,[G.lmin,G.lmax]); 
errors_cheb_warped=y-y_cheb_warped;
sup_err_cheb_warped=max(abs(errors_cheb_warped));
se_cheb_warped=sum(errors_cheb_warped.^2);

% y_cheb_warped_pchip_interp=gsp_cheby_eval(G2.e,p_warped_pchip_interp_c,[0,G.lmax]);
% errors_cheb_warped_pchip_interp=y-y_cheb_warped_pchip_interp;
% sup_err_cheb_warped_pchip_interp=max(abs(errors_cheb_warped_pchip_interp));
% se_cheb_warped_pchip_interp=sum(errors_cheb_warped_pchip_interp.^2);

y_warped_ls=polyval(p_warped_lsc,G2.e,s_warped_ls,mu_warped_ls);
errors_warped_ls=y-y_warped_ls;
sup_err_warped_ls=max(abs(errors_warped_ls));
se_warped_ls=sum(errors_warped_ls.^2);

y_spec_adapted_ortho=three_term_eval(G,G2.e,ab,c_spec_adapted_ortho);
errors_spec_adapted_ortho=y-y_spec_adapted_ortho;
sup_err_spec_adapted_ortho=max(abs(errors_spec_adapted_ortho));
se_spec_adapted_ortho=sum(errors_spec_adapted_ortho.^2);

y_spec_adapted_ortho_warped=three_term_eval(G,G2.e,ab_warped,c_spec_adapted_ortho_warped);
errors_spec_adapted_ortho_warped=y-y_spec_adapted_ortho_warped;
sup_err_spec_adapted_ortho_warped=max(abs(errors_spec_adapted_ortho_warped));
se_spec_adapted_ortho_warped=sum(errors_spec_adapted_ortho_warped.^2);

y_spec_adapted_ortho_wch=three_term_eval(G,G2.e,ab_wch,c_spec_adapted_ortho_wch);
errors_spec_adapted_ortho_wch=y-y_spec_adapted_ortho_wch;
sup_err_spec_adapted_ortho_wch=max(abs(errors_spec_adapted_ortho_wch));
se_spec_adapted_ortho_wch=sum(errors_spec_adapted_ortho_wch.^2);

y_spec_adapted_ortho_spline=three_term_eval(G,G2.e,ab_spline,c_spec_adapted_ortho_spline);
errors_spec_adapted_ortho_spline=y-y_spec_adapted_ortho_spline;
sup_err_spec_adapted_ortho_spline=max(abs(errors_spec_adapted_ortho_spline));
se_spec_adapted_ortho_spline=sum(errors_spec_adapted_ortho_spline.^2);

y_weighted_ls=polyval(p_weighted_lsc,G2.e,s_weighted_ls,mu_weighted_ls);
%y_weighted_ls=polyval(p_weighted_lsc,G2.e);
errors_weighted_ls=y-y_weighted_ls;
sup_err_weighted_ls=max(abs(errors_weighted_ls));
se_weighted_ls=sum(errors_weighted_ls.^2);

y_weighted_ls_2=polyval(p_weighted_lsc_2,G2.e,s_weighted_ls_2,mu_weighted_ls_2);
errors_weighted_ls_2=y-y_weighted_ls_2;

y_ls=polyval(lsc,G2.e,s_ls,mu_ls);
errors_ls=y-y_ls;
sup_err_ls=max(abs(errors_ls));
se_ls=sum(errors_ls.^2);

% Method = ["Chebyshev"; "Legendre"; "Regular Lanczos";"Warped Cheby";...
%     "Warped LS";"Spectrum-Adapted Ortho. Poly.";"Spectrum-Adapted w/ Warped Grid";"Spectrum-Adapted w/ Spline";"Spectrum-Adapted w/ Warped Cheby Grid";"Weighted LS";"Discrete LS";"Exact"];
% Sup_Err = [sup_err_cheb;sup_err_leg;sup_err_lanczos;sup_err_cheb_warped;...
%     sup_err_warped_ls;sup_err_spec_adapted_ortho;sup_err_spec_adapted_ortho_warped;sup_err_spec_adapted_ortho_spline;sup_err_spec_adapted_ortho_wch;sup_err_weighted_ls;sup_err_ls;0];
% Sq_Err = [se_cheb;se_leg;se_lanczos;se_cheb_warped;...
%     se_warped_ls;se_spec_adapted_ortho;se_spec_adapted_ortho_warped;se_spec_adapted_ortho_spline;se_spec_adapted_ortho_wch;se_weighted_ls;se_ls;0];

% Plots
if all_plots

xmax=G.lmax;
delta=xmax/5000;
xx=0:delta:xmax;
xx=xx';
yy=f(xx);

yy_cheb=gsp_cheby_eval(xx,c,[G2.e(1),G.lmax]);
yy_leg=gsp_cheby_eval(xx,plegc,[G2.e(1),G.lmax]);
yy_cheb_warped=gsp_cheby_eval(xx,p_warped_c,[G2.e(1),G.lmax]);
%yy_cheb_warped_pchip_interp=gsp_cheby_eval(xx,p_warped_pchip_interp_c,[0,G.lmax]);
yy_cheb_warped_ls=polyval(p_warped_lsc,xx,s_warped_ls,mu_warped_ls);
yy_spec_adapted_ortho=three_term_eval(G,xx,ab,c_spec_adapted_ortho);
yy_spec_adapted_ortho_warped=three_term_eval(G,xx,ab_warped,c_spec_adapted_ortho_warped);
yy_ls=polyval(lsc,xx,s_ls,mu_ls);

% Figure 2: approximations for function g
figure;
colors={[0.308 0.785 0.636]; [0 0.4470 0.7410]; [0.8500 0.3250 0.0980]; [0.9290 0.6940 0.1250]; [0.4660 0.6740 0.1880]; [0.3010 0.7450 0.9330]; [0.4940 0.1840 0.5560]};
plot2=plot(xx,[yy_ls,yy_cheb,pleg(xx),yy_cheb_warped,yy_cheb_warped_ls,yy_spec_adapted_ortho,yy_spec_adapted_ortho_warped],'LineWidth',4);
set(plot2,{'Color'},colors);
set(gca,'FontSize',20);
legend(plot2,'Discrete LS','Chebyshev','Legendre','Warped Chebyshev','Warped LS','Spectrum-Adapted Ortho. Poly.','Spectrum-Adapted w/ Warped Grid','Location','Northwest');  
xlabel('\lambda');
grid on;
hold on;
xlabel('\lambda');
ylabel('Filter approximations');
plot(G2.e, y_lanczos, 'LineWidth',4,'DisplayName','Regular Lanczos','Color',[1 0.6 0.8]);
plot(xx, yy, 'LineWidth',4,'DisplayName','g','Color',[0.9 0.5 0.25]);
plot(G2.e,zeros(G.N,1),'xk','LineWidth',2,'MarkerSize',6,'DisplayName','Eigenvalues');
plot(pts_warped_interp,f(pts_warped_interp),'xr','LineWidth',2,'MarkerSize',10,'DisplayName','Warped Chebyshev Pts');
plot(e,f(e),'xb','LineWidth',2,'MarkerSize',10,'DisplayName','Lanczos Pts');
title(['Graph=',graph,', Filter=',filter,', K=',num2str(K)]);
ylim([0,1]);

% Figure 3: absolute error comparison focused on eigenvalues
figure;
plot3=semilogy(G2.e,[abs(errors_ls),abs(errors_cheb),abs(errors_leg),abs(y-y_lanczos),abs(errors_cheb_warped),abs(errors_warped_ls),abs(errors_spec_adapted_ortho),abs(errors_spec_adapted_ortho_warped)],'-o','LineWidth',1,'MarkerSize',6);
set(gca,'FontSize',24);
morecolors={[0.308 0.785 0.636]; [0 0.4470 0.7410]; [0.8500 0.3250 0.0980]; [1 0.6 0.8]; [0.9290 0.6940 0.1250]; [0.4660 0.6740 0.1880]; [0.3010 0.7450 0.9330]; [0.4940 0.1840 0.5560]};
set(plot3,{'Color'},morecolors);
set(plot3, {'MarkerFaceColor'}, morecolors); 
hold on;
plot(G2.e,ones(G.N,1),'xk','LineWidth',2,'MarkerSize',6);
legend(plot3,'Discrete LS','Chebyshev','Legendre','Lanczos','Warped Chebyshev','Warped LS','Spectrum-Adapted Ortho. Poly.','Spectrum-Adapted w/ Warped Grid','Location','SouthEast');  
xlabel('\lambda');
ylabel('$$|\tilde{g}(\lambda)-g(\lambda)|$$','Interpreter','latex');
title(['G=',graph,', f=',filter,', K=',num2str(K),', LS=',start_pts_ls,', Damp=',damping]);
grid on;
%xlim([0,4.5]);
%ylim([10e-20,10e0]);
end

% Test f(L)b on random signal b
num_tests=50;

ntime_cheb=zeros(num_tests,1);
ntime_leg=zeros(num_tests,1);
ntime_lanczos=zeros(num_tests,1);
ntime_warped=zeros(num_tests,1);
ntime_warped_ls=zeros(num_tests,1);
ntime_spec=zeros(num_tests,1);
ntime_spec_warp=zeros(num_tests,1);
ntime_spec_spline=zeros(num_tests,1);
ntime_spec_wch=zeros(num_tests,1);
ntime_weighted_ls=zeros(num_tests,1);
ntime_weighted_ls_2=zeros(num_tests,1);
ntime_weighted_ls_spline=zeros(num_tests,1);
ntime_ls=zeros(num_tests,1);
ntime_exact=zeros(num_tests,1);

nmse_cheb=zeros(num_tests,1);
nmse_leg=zeros(num_tests,1);
nmse_lanczos=zeros(num_tests,1);
nmse_warped=zeros(num_tests,1);
nmse_warped_ls=zeros(num_tests,1);
nmse_spec=zeros(num_tests,1);
nmse_spec_warp=zeros(num_tests,1);
nmse_spec_spline=zeros(num_tests,1);
nmse_spec_wch=zeros(num_tests,1);
nmse_weighted_ls=zeros(num_tests,1);
nmse_weighted_ls_2=zeros(num_tests,1);
nmse_weighted_ls_spline=zeros(num_tests,1);
nmse_ls=zeros(num_tests,1);

lanc_param.method='lanczos';
lanc_param.order=K;

hpf=@(x) (x>=gi(G.lmax*3/4)); % check if lanczos performs bad when b is
%high pass filtered

for i=1:num_tests
    b=rand(G.N,1);
    %b=gsp_filter(G,hpf,b);
    bhat=G2.U'*b;
    tic
    gLf_exact=G2.U*(f(G2.e).*bhat);
    ntime_exact(i)=toc;
    ntime_exact(i)=ntime_exact(i)+time_exact_spec;
    tic
    gLf_ls=G2.U*(polyval(lsc,G2.e,s_ls,mu_ls).*bhat);
    ntime_ls(i)=toc;
    ntime_ls(i)=ntime_ls(i)+time_exact_spec;
    tic
    gLf_cheb=gsp_cheby_op(G,c,b); 
    ntime_cheb(i)=toc;
    tic
    gLf_leg=gsp_cheby_op(G,plegc,b); 
    ntime_leg(i)=toc;
    tic
    gLf_lanczos=gsp_filter(G,f,b,lanc_param);
    ntime_lanczos(i)=toc;
    tic
    %gLf_warped=poly_op(G,p_warped,b); 
    gLf_warped=G2.U*(polyval(p_warped,G2.e,s_warped,mu_warped).*bhat);
    ntime_warped(i)=toc;
%    tic
%    gLf_warped_pchip_interp=gsp_cheby_op(G,p_warped_pchip_interp_c,b);
%    ntime_spec_warp(i)=toc;
    tic
    %gLf_warped_ls=poly_op(G,p_warped_lsc,b);
    gLf_warped_ls=G2.U*(polyval(p_warped_lsc,G2.e,s_warped_ls,mu_warped_ls).*bhat);
    ntime_warped_ls(i)=toc;
    tic
    gLf_spec=three_term_recurr_op(G,ab,c_spec_adapted_ortho,b);
    ntime_spec(i)=toc;
    tic
    gLf_spec_warp=three_term_recurr_op(G,ab_warped,c_spec_adapted_ortho_warped,b);
    ntime_spec_warp(i)=toc;
    tic
    gLf_spec_spline=three_term_recurr_op(G,ab_spline,c_spec_adapted_ortho_spline,b);
    ntime_spec_spline(i)=toc;
    tic
    gLf_spec_wch=three_term_recurr_op(G,ab_wch,c_spec_adapted_ortho_wch,b);
    ntime_spec_wch(i)=toc;
    tic
    gLf_weighted_ls=G2.U*(polyval(p_weighted_lsc,G2.e,s_weighted_ls,mu_weighted_ls).*bhat);
    %gLf_weighted_ls=G2.U*(polyval(p_weighted_lsc,G2.e).*bhat);
    ntime_weighted_ls(i)=toc;
    tic
    gLf_weighted_ls_2=G2.U*(polyval(p_weighted_lsc_2,G2.e,s_weighted_ls_2,mu_weighted_ls_2).*bhat);
    ntime_weighted_ls_2(i)=toc;
    %tic
    %gLf_weighted_ls_spline=G2.U*(polyval(p_weighted_lsc_spline,G2.e,s_weighted_ls_spline,mu_weighted_ls_spline).*bhat);
    %ntime_weighted_ls_spline(i)=toc;

    nmse_cheb(i)=sum((gLf_cheb-gLf_exact).^2)/sum(gLf_exact.^2);
    nmse_leg(i)=sum((gLf_leg-gLf_exact).^2)/sum(gLf_exact.^2);
    nmse_lanczos(i)=sum((gLf_lanczos-gLf_exact).^2)/sum(gLf_exact.^2);
    nmse_warped(i)=sum((gLf_warped-gLf_exact).^2)/sum(gLf_exact.^2);
    nmse_warped_ls(i)=sum((gLf_warped_ls-gLf_exact).^2)/sum(gLf_exact.^2);
    nmse_spec(i)=sum((gLf_spec-gLf_exact).^2)/sum(gLf_exact.^2);
    nmse_spec_warp(i)=sum((gLf_spec_warp-gLf_exact).^2)/sum(gLf_exact.^2);
    nmse_spec_spline(i)=sum((gLf_spec_spline-gLf_exact).^2)/sum(gLf_exact.^2);
    nmse_spec_wch(i)=sum((gLf_spec_wch-gLf_exact).^2)/sum(gLf_exact.^2);
    nmse_weighted_ls(i)=sum((gLf_weighted_ls-gLf_exact).^2)/sum(gLf_exact.^2);
    nmse_weighted_ls_2(i)=sum((gLf_weighted_ls_2-gLf_exact).^2)/sum(gLf_exact.^2);
    %nmse_weighted_ls_spline(i)=sum((gLf_weighted_ls_spline-gLf_exact).^2)/sum(gLf_exact.^2);
    nmse_ls(i)=sum((gLf_ls-gLf_exact).^2)/sum(gLf_exact.^2);
    
end

avg_ntime_cheb=mean(ntime_cheb);
avg_ntime_leg=mean(ntime_leg);
avg_ntime_lanczos=mean(ntime_lanczos);
avg_ntime_warped=mean(ntime_warped);
avg_ntime_warped_ls=mean(ntime_warped_ls);
avg_ntime_spec=mean(ntime_spec);
avg_ntime_spec_warp=mean(ntime_spec_warp);
avg_ntime_spec_spline=mean(ntime_spec_spline);
avg_ntime_spec_wch=mean(ntime_spec_wch);
avg_ntime_weighted_ls=mean(ntime_weighted_ls);
avg_ntime_weighted_ls_2=mean(ntime_weighted_ls_2);
%avg_ntime_weighted_ls_spline=mean(ntime_weighted_ls_spline);
avg_ntime_ls=mean(ntime_ls);
avg_ntime_exact=mean(ntime_exact);

avg_nmse_cheb=mean(nmse_cheb);
avg_nmse_leg=mean(nmse_leg);
avg_nmse_lanczos=mean(nmse_lanczos);
avg_nmse_warped=mean(nmse_warped);
avg_nmse_warped_ls=mean(nmse_warped_ls);
avg_nmse_spec=mean(nmse_spec);
avg_nmse_spec_warp=mean(nmse_spec_warp);
avg_nmse_spec_spline=mean(nmse_spec_spline);
avg_nmse_spec_wch=mean(nmse_spec_wch);
avg_nmse_weighted_ls=mean(nmse_weighted_ls);
avg_nmse_weighted_ls_2=mean(nmse_weighted_ls_2);
%avg_nmse_weighted_ls_spline=mean(nmse_weighted_ls_spline);
avg_nmse_ls=mean(nmse_ls);

cheb_sup_err(K_ind)=sup_err_cheb;
cheb_sq_err(K_ind)=se_cheb;
cheb_time(K_ind)=avg_ntime_cheb;
cheb_anmse(K_ind)=avg_nmse_cheb;

leg_sup_err(K_ind)=sup_err_leg;
leg_sq_err(K_ind)=se_leg;
leg_time(K_ind)=avg_ntime_leg;
leg_anmse(K_ind)=avg_nmse_leg;

lanc_sup_err(K_ind)=sup_err_lanczos;
lanc_sq_err(K_ind)=se_lanczos;
lanc_time(K_ind)=avg_ntime_lanczos;
lanc_anmse(K_ind)=avg_nmse_lanczos;

warp_sup_err(K_ind)=sup_err_cheb_warped;
warp_sq_err(K_ind)=se_cheb_warped;
warp_time(K_ind)=avg_ntime_warped;
warp_anmse(K_ind)=avg_nmse_warped;

warp_ls_sup_err(K_ind)=sup_err_warped_ls;
warp_ls_sq_err(K_ind)=se_warped_ls;
warp_ls_time(K_ind)=avg_ntime_warped_ls;
warp_ls_anmse(K_ind)=avg_nmse_warped_ls;

spec_sup_err(K_ind)=sup_err_spec_adapted_ortho;
spec_sq_err(K_ind)=se_spec_adapted_ortho;
spec_time(K_ind)=avg_ntime_spec;
spec_anmse(K_ind)=avg_nmse_spec;

spec_warp_sup_err(K_ind)=sup_err_spec_adapted_ortho_warped;
spec_warp_sq_err(K_ind)=se_spec_adapted_ortho_warped;
spec_warp_time(K_ind)=avg_ntime_spec_warp;
spec_warp_anmse(K_ind)=avg_nmse_spec_warp;

spec_spline_sup_err(K_ind)=sup_err_spec_adapted_ortho_spline;
spec_spline_sq_err(K_ind)=se_spec_adapted_ortho_spline;
spec_spline_time(K_ind)=avg_ntime_spec_spline;
spec_spline_anmse(K_ind)=avg_nmse_spec_spline;

spec_wch_sup_err(K_ind)=sup_err_spec_adapted_ortho_wch;
spec_wch_sq_err(K_ind)=se_spec_adapted_ortho_wch;
spec_wch_time(K_ind)=avg_ntime_spec_wch;
spec_wch_anmse(K_ind)=avg_nmse_spec_wch;

weighted_ls_sup_err(K_ind)=sup_err_weighted_ls;
weighted_ls_sq_err(K_ind)=se_weighted_ls;
weighted_ls_time(K_ind)=avg_ntime_weighted_ls;
weighted_ls_anmse(K_ind)=avg_nmse_weighted_ls;

weighted_ls_2_anmse(K_ind)=avg_nmse_weighted_ls_2;
%weighted_ls_spline_anmse(K_ind)=avg_nmse_weighted_ls_spline;

end

parameters = cell2table({graph,param.cdf_method,filter,damping},...
    'VariableNames',{'Graph' 'CDFMethod' 'Filter' 'Damping'})

cheb_table=table(poly_orders,cheb_sup_err,cheb_sq_err,cheb_anmse,cheb_time)
leg_table=table(poly_orders,leg_sup_err,leg_sq_err,leg_anmse,leg_time)
lanc_table=table(poly_orders,lanc_sup_err,lanc_sq_err,lanc_anmse,lanc_time)
warp_table=table(poly_orders,warp_sup_err,warp_sq_err,warp_anmse,warp_time)
warp_ls_table=table(poly_orders,warp_ls_sup_err,warp_ls_sq_err,warp_ls_anmse,warp_ls_time)
spec_table=table(poly_orders,spec_sup_err,spec_sq_err,spec_anmse,spec_time)
spec_warp_table=table(poly_orders,spec_warp_sup_err,spec_warp_sq_err,spec_warp_anmse,spec_warp_time)
spec_spline_table=table(poly_orders,spec_spline_sup_err,spec_spline_sq_err,spec_spline_anmse,spec_spline_time)
spec_wch_table=table(poly_orders,spec_wch_sup_err,spec_wch_sq_err,spec_wch_anmse,spec_wch_time)
weighted_ls_table=table(poly_orders,weighted_ls_sup_err,weighted_ls_sq_err,weighted_ls_anmse,weighted_ls_time)


% Figure 4: convergence comparison - error vs order plot
figure;
plot4=semilogy(poly_orders,[cheb_anmse,lanc_anmse,warp_anmse,warp_ls_anmse,...
    spec_anmse,weighted_ls_anmse],'-o','LineWidth',1,'MarkerSize',10);
set(gca,'FontSize',24);
set(plot4, {'MarkerFaceColor'}, get(plot4,'Color')); 
hold on;
legend(plot4,'Chebyshev','Lanczos','Interpolation','Warped LS',...
    'SAOP','Weighted LS','Location','Northwest');  
xlabel('K');
ylabel('$$||p(A)b-f(A)b||^2/||f(A)b||^2$$','Interpreter','latex');
grid on;
title(['G=',graph,', lmax=',num2str(G.lmax),', f=',filter]);

% figure;
% plot4=semilogy(poly_orders,[cheb_anmse,leg_anmse,lanc_anmse,warp_anmse,warp_ls_anmse,...
%     spec_anmse,spec_warp_anmse,spec_spline_anmse,spec_wch_anmse,weighted_ls_anmse],'-o','LineWidth',1,'MarkerSize',10);
% set(gca,'FontSize',20);
% ninecolors={[0 0.4470 0.7410]; [0.8500 0.3250 0.0980]; [1 0.6 0.8]; [0.9290 0.6940 0.1250];...
%     [0.4660 0.6740 0.1880]; [0.3010 0.7450 0.9330]; [0.4940 0.1840 0.5560];[0.9 0.5 0.25];...
%     [0.5 0.3 0.9]; [0.308 0.785 0.636]};
% set(plot4,{'Color'},ninecolors);
% set(plot4, {'MarkerFaceColor'}, get(plot4,'Color')); 
% hold on;
% legend(plot4,'Chebyshev','Legendre','Lanczos','Interpolation','Warped LS',...
%     'SAOP','SAOP w/ Warped Grid','SAOP w/ Spline',...
%     'SAOP w/ Warped Nodes & Equal Weights','Weighted LS','Location','Northwest');  
% xlabel('K');
% ylabel('$$||p(A)b-f(A)b||^2/||f(A)b||^2$$','Interpreter','latex');
% grid on;
% title(['G=',graph,', lmax=',num2str(G.lmax),', f=',filter]);

% % Figure 5: error vs order plot for weighted LS and spec adapted with the
% % same set of evenly spaced absc and weights from pdf
% figure;
% plot5=semilogy(poly_orders,[weighted_ls_anmse,spec_anmse],'-o','LineWidth',1,'MarkerSize',10);
% set(gca,'FontSize',24);
% set(plot5, {'MarkerFaceColor'}, get(plot5,'Color')); 
% hold on;
% legend(plot5,'Weighted LS pdf','SA even','Location','Northwest');  
% xlabel('K');
% ylabel('$$||p(A)b-f(A)b||^2/||f(A)b||^2$$','Interpreter','latex');
% title(['G=',graph,', lmax=',num2str(G2.lmax),', f=',filter]);
% grid on;


% % Figure 6: absolute errors for weighted LS and spec adapted with the
% % same set of evenly spaced absc and weights from pdf
% figure;
% plot6=semilogy(G2.e,[abs(errors_weighted_ls),abs(errors_spec_adapted_ortho)],'-o','LineWidth',1,'MarkerSize',10);
% set(gca,'FontSize',24);
% set(plot6, {'MarkerFaceColor'}, get(plot6,'Color')); 
% hold on;
% plot(G2.e,ones(G.N,1),'xk','LineWidth',2,'MarkerSize',6);
% legend(plot6,'Weighted LS pdf','SA even','Location','Northwest');  
% xlabel('\lambda');
% ylabel('$$|\tilde{g}(\lambda)-g(\lambda)|$$','Interpreter','latex');
% title(['G=',graph,', lmax=',num2str(G2.lmax),', f=',filter,', K=',num2str(K)]);
% grid on;


% % plot pdf/cdf estimate
% xx=linspace(G.lmin,G.lmax,1000);
% figure; plot(xx,G.spectrum_pdf_approx(xx),'linewidth',4,'displayname','Spectral Density');
% mu=@(t) sum(G2.e<=t)/G.N;
% figure;
% plot(xx,mu(xx),'linewidth',4,'displayname','CDF');hold on;
% plot(xx,G.spectrum_cdf_approx(xx),'linewidth',4,'displayname','CDF estimate');
% legend('show');



