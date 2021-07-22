% Demo for spectrum-adapted polynomial approximation for matrix function 
% f(A), where A is a graph Laplacian matrix

close all;
clear all;

% Graph
N=500;
p=.2;
G=gsp_erdos_renyi(N,p);

% Save a copy of the graph for comparison to exact computation
G2=G;
G2=gsp_compute_fourier_basis(G2);
G2.lmin=min(G2.e);

% Estimate spectral CDF of the graph laplacian matrix
G=gsp_estimate_lmax(G);
G.lmin=(abs(min(G2.e))>1e-6)*min(G2.e);
param.cdf_method='kpm'; % 'kpm', 'lanczos' or 'ldlt'
param.num_pts=10; % parameter for kpm and ldlt
param.num_vec=10; % parameter for lanczos
param.order=30;
G=spectral_cdf_approx2(G,param); 
gi=@(s) G.spectrum_inv_cdf_approx((s-G.lmin)/(G.lmax-G.lmin));

% Filter function
tau=1;
f=@(x) exp(-tau*x);

% Plot the estimated spectral CDF and the actual spectral CDF
xx=linspace(G.lmin,G.lmax,10000);
mu=@(t) sum(G2.e<=t)/G.N;
figure;
plot1=plot(xx,[G.spectrum_cdf_approx(xx)',mu(xx)'],'linewidth',4);hold on;
plot(xx,xx/G.lmax,'k:','LineWidth',1.5);
set(gca,'FontSize',24);
grid on;
legend(plot1,'Estimated Spectral CDF','Actual Spectral CDF','Location','Southeast');  
xlabel('\lambda');
xlim([G.lmin,G.lmax]);

%==========================================================================

% Order for polynomial approximation
K=10;

% Truncated Chebyshev expansion
c=gsp_cheby_coeff([G.lmin,G.lmax],f,K,1000);

% Lanczos
lanc_param.method='lanczos';
lanc_param.order=K;
x=sum(G2.U')';

% Starting points for interpolation
pts_interp=cos((0:K)*pi/K)'; 
pts_tx_interp=(pts_interp+1)*(G.lmax-G.lmin)/2+G.lmin;
pts_tx_interp=sort(pts_tx_interp,'ascend');
pts_warped_interp=gi(pts_tx_interp);
pts_warped_interp(1)=G.lmin;

% Spectrum-adapted interpolation
[p_warped,s_warped,mu_warped]=polyfit(pts_warped_interp,f(pts_warped_interp),K);
p_warped_fun=@(x) polyval(p_warped,x,s_warped,mu_warped); 
p_warped_c=gsp_cheby_coeff([G.lmin,G.lmax],p_warped_fun,K,1000);

% Truncated expansion in spectrum-adapted orthogonal polynomials 
% (equivalent to spectrum-adapted weighted least squares regression of same grid order)
mop_param=struct;
mop_param.grid_order=100;
mop_param.absc_method='linear'; % 'linear', 'even', 'warp', or 'spline'
mop_param.weights_method='pdf';
[absc,weights]=gen_absc_weights(G,K,mop_param); 
[ab,Pi]=matrix_adapted_ortho_poly(absc,weights,K);
c_spec_adapted_ortho=matrix_adapted_poly_coeff(G, f, absc, weights/sum(weights), Pi, K);

% Spectrum-adapted weighted least squares regression
% (equivalent to truncated expansion in spectrum-adapted orthogonal polynomials of same grid order)
mop_param.grid_order=100;
mop_param.absc_method='linear';
mop_param.weights_method='pdf';
[absc_weighted_ls,weights_weighted_ls]=gen_absc_weights(G,K,mop_param);
[p_weighted_lsc,s_weighted_ls,mu_weighted_ls]=weighted_polyfit(absc_weighted_ls,f(absc_weighted_ls),weights_weighted_ls,K);

% Compute polynomial approximation at the eigenvalues and errors
y=f(G2.e);

% Truncated Chebyshev expansion
y_cheb=gsp_cheby_eval(G2.e,c,[G.lmin,G.lmax]);
errors_cheb=y-y_cheb;
sup_err_cheb=max(abs(errors_cheb));
se_cheb=sum(errors_cheb.^2);

% Lanczos
y_lanczos=G2.U'*gsp_filter(G,f,sum(G2.U')',lanc_param);
errors_lanczos=y-y_lanczos;
sup_err_lanczos=max(abs(errors_lanczos));
se_lanczos=sum(errors_lanczos.^2);

% Spectrum-adapted interpolation
y_cheb_warped=gsp_cheby_eval(G2.e,p_warped_c,[G.lmin,G.lmax]); 
errors_cheb_warped=y-y_cheb_warped;
sup_err_cheb_warped=max(abs(errors_cheb_warped));
se_cheb_warped=sum(errors_cheb_warped.^2);

% Truncated expansion in spectrum-adapted orthogonal polynomials 
y_spec_adapted_ortho=three_term_eval(G,G2.e,ab(1:K+2,:),c_spec_adapted_ortho(1:K+1));
errors_spec_adapted_ortho=y-y_spec_adapted_ortho;
sup_err_spec_adapted_ortho=max(abs(errors_spec_adapted_ortho));
se_spec_adapted_ortho=sum(errors_spec_adapted_ortho.^2);

% Spectrum-adapted weighted least squares regression
% (equivalent to truncated expansion in spectrum-adapted orthogonal polynomials of same grid order)
y_weighted_ls=polyval(p_weighted_lsc,G2.e,s_weighted_ls,mu_weighted_ls);
errors_weighted_ls=y-y_weighted_ls;
sup_err_weighted_ls=max(abs(errors_weighted_ls));
se_weighted_ls=sum(errors_weighted_ls.^2);

% Plot absolute errors for different approximation methods
figure;
plot2=semilogy(G2.e,[abs(errors_cheb),abs(errors_lanczos),abs(errors_cheb_warped),abs(errors_weighted_ls)],'-o','LineWidth',1,'MarkerSize',10);
set(gca,'FontSize',24);
set(plot2, {'MarkerFaceColor'}, get(plot2,'Color')); 
hold on;
plot(G2.e,ones(G.N,1),'xk','LineWidth',2,'MarkerSize',6);
legend(plot2,'Chebyshev','Lanczos','Interpolation','Weighted LS','Location','SouthEast');  
xlabel('\lambda');
%ylabel('$$|f(\lambda)-p_{10}(\lambda)|$$','Interpreter','latex');
grid on;

