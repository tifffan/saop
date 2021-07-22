% Figures for cdf estimations and numerical examples of saop paper

close all;
clear all;
randn('seed', 18); rand('seed', 18)

% Graph

graph='bunny';
switch graph
    case 'bunny'
        G=gsp_bunny();
    case 'gnp'
        N=500;
        p=.2;
        G=gsp_erdos_renyi(N,p);
    case 'minnesota'
        G=gsp_minnesota(1);
    case 'net25'
        load('0_data/net25.mat');
        A=Problem.A;
        A = A - diag(diag(A)); 
        A(4228,6327) = 1;
        A(6327,4228) = 1;
        G=gsp_graph(A);
        %G.L=G.L./20;
    case 'cage9'
        load('0_data/cage9.mat');
        G=struct;
        G.L=Problem.A;
        G.L=(G.L+G.L')/2;
        G.N=size(G.L,1);
   case 'si2'
        load('0_data/Si2.mat');
        G=struct;
        G.L=Problem.A;
        G.L=(G.L+G.L')/2;
        G.N=size(G.L,1);
    case 'saylr4'
        load('0_data/saylr4.mat');
        A=Problem.A;
        A = A - diag(diag(A)); 
        G=gsp_graph(A);
   case 'saylr4sc'
        load('0_data/saylr4.mat');
        A=Problem.A;
        A = A - diag(diag(A)); 
        G=gsp_graph(A);
        G.L=G.L/2000;
end

if graph=="si2"
G.L=G.L+0.01*eye(G.N);
end
G2=G;
tic
G2=gsp_compute_fourier_basis(G2);
G2.lmin=min(G2.e);
time_exact_spec=toc;

G=gsp_estimate_lmax(G);
G.lmin=(abs(min(G2.e))>1e-6)*min(G2.e);  %%% Replace with estimate
param.cdf_method='kpm'; % 'kpm' or 'lanczos' or 'ldlt'
param.num_pts=10; % for kpm and ldlt
param.num_vec=10; % for lanczos
param.order=30; % 100 works for kpm and ldlt
G=spectral_cdf_approx2(G,param); 
gi=@(s) G.spectrum_inv_cdf_approx((s-G.lmin)/(G.lmax-G.lmin));

% plot cdf
% xx=linspace(G.lmin,G.lmax,1000);
% mu=@(t) sum(G2.e<=t)/G.N;
% figure;
% plot1=plot(xx,[G.spectrum_cdf_approx(xx)',mu(xx)'],'linewidth',4);hold on;
% plot(xx,xx/G.lmax,'k:','LineWidth',1.5);
% set(gca,'FontSize',24);
% grid on;
% legend(plot1,'Estimated Spectral CDF','Actual Spectral CDF','Location','Southeast');  
% xlabel('\lambda');
% if graph=="si2"
%     ylim([-0.01,1]);
% end
% xlim([G.lmin,G.lmax]);
% savefig(strcat(graph,'_cdf.fig'));
% close all;


% % plot cdf 2
% xx=linspace(G.lmin,G.lmax,1000);
% figure;
% plot2=plot(xx,G.spectrum_cdf_approx(xx),'linewidth',4);hold on;
% set(gca,'FontSize',24);
% grid on; hold on;
% scatter(pts_warped_interp,pts_tx_interp/G.lmax);
% xlabel('\lambda');
% ylabel('Estimated Spectral CDF');
% 
% 
% % plot cdf 2
% xx=linspace(G.lmin,G.lmax,1000);
% figure;
% plot2=plot(xx,p_warped_fun(xx),'linewidth',4);hold on;
% set(gca,'FontSize',24);
% grid on; hold on;
% % scatter(pts_warped_interp,pts_tx_interp/G.lmax);
% xlabel('\lambda');
% ylabel('Estimated function');



%==========================================================================

% filter
% [filts] = spm_parseval_filters(G,4,'unif_itersine');
% f = filts{1}
% filter = "unif_itersine1";

% [filts] = spm_parseval_filters(G,4,'wav_itersine');
% f = filts{1}
% filter = "wav_itersine1";

% [filts] = spm_parseval_filters(G,4,'sa_unif_itersine');
% f = filts{1}
% filter = "sa_unif_itersine1";

% [filts] = spm_parseval_filters(G,4,'sa_wav_itersine');
% f = filts{1}
% filter = "sa_wav_itersine1";

filt_str='unif_itersine';
max_filt_idx=5;
[filts] = spm_parseval_filters(G,max_filt_idx,filt_str);

% for filt_idx=1:max_filt_idx
for filt_idx=1:2:max_filt_idx
f = filts{filt_idx};
filter = strcat(filt_str,num2str(filt_idx));



%==========================================================================



% % plot filter
% xx=linspace(G.lmin,G.lmax,1000);
% figure;
% plot2=plot(xx,f(xx),'linewidth',4);
% set(gca,'FontSize',24);
% grid on;
% xlabel('\lambda');
% ylabel('f(\lambda)');
% if graph=="si2"
%     ylim([-0.01,1]);
% end
% xlim([G.lmin,G.lmax]);
% savefig(strcat('fig_filter_',filter,'.fig'));
% close all;

%==========================================================================

% order for the plot of errors at eigenvalues
K=10;

% Chebyshev
c=gsp_cheby_coeff([G.lmin,G.lmax],f,K,1000);

% Lanczos
lanc_param.method='lanczos';
lanc_param.order=K;
x=sum(G2.U')';

% Starting points for interpolation

interp_points='chebyshev';
switch interp_points
    case 'chebyshev'
        pts_interp=cos((0:K)*pi/K)'; 
        pts_tx_interp=(pts_interp+1)*(G.lmax-G.lmin)/2+G.lmin;
        pts_tx_interp=sort(pts_tx_interp,'ascend');
    case 'uniform'
        pts_tx_interp=linspace(G.lmin,G.lmax,K+1)';
end

pts_warped_interp=gi(pts_tx_interp);
pts_warped_interp(1)=G.lmin;

% Chebyshev interpolation on warped points
[p_warped,s_warped,mu_warped]=polyfit(pts_warped_interp,f(pts_warped_interp),K);
p_warped_fun=@(x) polyval(p_warped,x,s_warped,mu_warped); 
p_warped_c=gsp_cheby_coeff([G.lmin,G.lmax],p_warped_fun,K,1000);

% Matrix/Spectrum adapted orthogonal polynomials
mop_param=struct;
mop_param.grid_order=100;
%mop_param.init_poly_order=30;
%mop_param.num_vec=10;
mop_param.absc_method='linear'; %linear, even(n-2 linearly spaced + both ends), warp, spline
mop_param.weights_method='pdf'; %count, pdf2(from kpm paper)
[absc,weights]=gen_absc_weights(G,K,mop_param); 
[ab,Pi]=matrix_adapted_ortho_poly(absc,weights,K);
c_spec_adapted_ortho=matrix_adapted_poly_coeff(G, f, absc, weights/sum(weights), Pi, K);

% Weighted LS fitting
mop_param.grid_order=100;
mop_param.absc_method='linear';
mop_param.weights_method='pdf';
[absc_weighted_ls,weights_weighted_ls]=gen_absc_weights(G,K,mop_param);
[p_weighted_lsc,s_weighted_ls,mu_weighted_ls]=weighted_polyfit(absc_weighted_ls,f(absc_weighted_ls),weights_weighted_ls,K);

% Compute polynomial approximation values at the actual eigenvalues and 
% the corresponding superior and squared error
y=f(G2.e);

y_cheb=gsp_cheby_eval(G2.e,c,[G.lmin,G.lmax]);
errors_cheb=y-y_cheb;
sup_err_cheb=max(abs(errors_cheb));
se_cheb=sum(errors_cheb.^2);

y_lanczos=G2.U'*gsp_filter(G,f,sum(G2.U')',lanc_param);
errors_lanczos=y-y_lanczos;
sup_err_lanczos=max(abs(errors_lanczos));
se_lanczos=sum(errors_lanczos.^2);

y_cheb_warped=gsp_cheby_eval(G2.e,p_warped_c,[G.lmin,G.lmax]); 
errors_cheb_warped=y-y_cheb_warped;
sup_err_cheb_warped=max(abs(errors_cheb_warped));
se_cheb_warped=sum(errors_cheb_warped.^2);

y_spec_adapted_ortho=three_term_eval(G,G2.e,ab(1:K+2,:),c_spec_adapted_ortho(1:K+1));
errors_spec_adapted_ortho=y-y_spec_adapted_ortho;
sup_err_spec_adapted_ortho=max(abs(errors_spec_adapted_ortho));
se_spec_adapted_ortho=sum(errors_spec_adapted_ortho.^2);

y_weighted_ls=polyval(p_weighted_lsc,G2.e,s_weighted_ls,mu_weighted_ls);
errors_weighted_ls=y-y_weighted_ls;
sup_err_weighted_ls=max(abs(errors_weighted_ls));
se_weighted_ls=sum(errors_weighted_ls.^2);

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
if graph=="si2"
    xlim([-1,50]);
end
% savefig(strcat(graph,'_',filter,'_deg',num2str(K),'.fig'));
close all;

%==========================================================================

% Poly approx order
poly_orders=3:25;
n_its=length(poly_orders);
poly_orders=poly_orders';

cheb_sup_err=zeros(n_its,1);
lanc_sup_err=zeros(n_its,1);
warp_sup_err=zeros(n_its,1);
spec_sup_err=zeros(n_its,1);
weighted_ls_sup_err=zeros(n_its,1);

cheb_sq_err=zeros(n_its,1);
lanc_sq_err=zeros(n_its,1);
warp_sq_err=zeros(n_its,1);
spec_sq_err=zeros(n_its,1);
weighted_ls_sq_err=zeros(n_its,1);

cheb_time=zeros(n_its,1);
lanc_time=zeros(n_its,1);
warp_time=zeros(n_its,1);
spec_time=zeros(n_its,1);
weighted_ls_time=zeros(n_its,1);

cheb_anmse=zeros(n_its,1);
lanc_anmse=zeros(n_its,1);
warp_anmse=zeros(n_its,1);
spec_anmse=zeros(n_its,1);
weighted_ls_anmse=zeros(n_its,1);

%==========================================================================

btype='constant_spectral'

switch btype
    case 'constant'
        num_tests=1;
        BB=ones(G.N,1);
    case 'constant_spectral'
        num_tests=1;
        BBB=(1/sqrt(G.N))*ones(G.N,1);
        BB=G2.U*BBB;
    case 'randn_spectral'
        num_tests=50;
        BBB=(1/sqrt(G.N))*randn(G.N,num_tests);
        BB=G2.U*BBB;
    case 'rand_spectral'
        num_tests=50;
        BBB=(1/sqrt(G.N))*rand(G.N,num_tests);
        BB=G2.U*BBB;
    case 'randn'
        num_tests=50;
        BB=(1/sqrt(G.N))*randn(G.N,num_tests);
    case 'rand'
        num_tests=50;
        BB=(1/sqrt(G.N))*rand(G.N,num_tests);
    case 'low'
        num_tests=50;
        BBB=(1/sqrt(G.N))*rand(G.N,num_tests);
        BB=G2.U*BBB;
        f2=@(x) (x<=gi(G.lmax*1/20)); 
        BB=gsp_filter(G,f2,BB);
    case 'high'
        num_tests=50;
        BBB=(1/sqrt(G.N))*rand(G.N,num_tests);
        BB=G2.U*BBB;
        f2=@(x) (x>=gi(G.lmax*19/20)); 
        BB=gsp_filter(G,f2,BB);
    case 'mid'
        num_tests=50;
        BBB=(1/sqrt(G.N))*rand(G.N,num_tests);
        BB=G2.U*BBB;
        f2=@(x) (x>=gi(G.lmax/3) & x<=gi(G.lmax*2/3));
        BB=gsp_filter(G,f2,BB);
    case 'bunny'
        num_tests=1;
        load('/Users/tiffan/Desktop/SAOP2020/data/pwbunny_signal.mat');
        BB=pwbunny_signal;
    otherwise
        error('unknown b type');
end

% Matrix/Spectrum adapted orthogonal polynomials
mop_param=struct;
mop_param.grid_order=100;  % 100, 200
mop_param.init_poly_order=30; % 30, 50
% mop_param.num_vec=10;
mop_param.absc_method='linear'; %linear, even(n-2 linearly spaced + both ends), warp, spline
mop_param.weights_method='pdf'; %count, pdf2(from kpm paper)
[absc,weights]=gen_absc_weights(G,K,mop_param); 
[ab,Pi]=matrix_adapted_ortho_poly(absc,weights,max(poly_orders));
c_spec_adapted_ortho=matrix_adapted_poly_coeff(G, f, absc, weights/sum(weights), Pi, max(poly_orders));
       
for K_ind=1:n_its
K=poly_orders(K_ind);

% Chebyshev
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

% Chebyshev interpolation on warped points
[p_warped,s_warped,mu_warped]=polyfit(pts_warped_interp,f(pts_warped_interp),K);
p_warped_fun=@(x) polyval(p_warped,x,s_warped,mu_warped); 
p_warped_c=gsp_cheby_coeff([G.lmin,G.lmax],p_warped_fun,K,1000);

% Weighted LS fitting
% mop_param.grid_order=100;
% mop_param.absc_method='linear';
% mop_param.weights_method='pdf';
[absc_weighted_ls,weights_weighted_ls]=gen_absc_weights(G,K,mop_param);
[p_weighted_lsc,s_weighted_ls,mu_weighted_ls]=weighted_polyfit(absc_weighted_ls,f(absc_weighted_ls),weights_weighted_ls,K);


y=f(G2.e);

% Compute polynomial approximation values at the actual eigenvalues and 
% the corresponding superior and squared error
y_cheb=gsp_cheby_eval(G2.e,c,[G.lmin,G.lmax]);
errors_cheb=y-y_cheb;
sup_err_cheb=max(abs(errors_cheb));
se_cheb=sum(errors_cheb.^2);

y_lanczos=G2.U'*gsp_filter(G,f,sum(G2.U')',lanc_param);
errors_lanczos=y-y_lanczos;
sup_err_lanczos=max(abs(errors_lanczos));
se_lanczos=sum(errors_lanczos.^2);

y_cheb_warped=gsp_cheby_eval(G2.e,p_warped_c,[G.lmin,G.lmax]); 
errors_cheb_warped=y-y_cheb_warped;
sup_err_cheb_warped=max(abs(errors_cheb_warped));
se_cheb_warped=sum(errors_cheb_warped.^2);

y_spec_adapted_ortho=three_term_eval(G,G2.e,ab(1:K+2,:),c_spec_adapted_ortho(1:K+1));
errors_spec_adapted_ortho=y-y_spec_adapted_ortho;
sup_err_spec_adapted_ortho=max(abs(errors_spec_adapted_ortho));
se_spec_adapted_ortho=sum(errors_spec_adapted_ortho.^2);

y_weighted_ls=polyval(p_weighted_lsc,G2.e,s_weighted_ls,mu_weighted_ls);
errors_weighted_ls=y-y_weighted_ls;
sup_err_weighted_ls=max(abs(errors_weighted_ls));
se_weighted_ls=sum(errors_weighted_ls.^2);

cheb_sq_err(K_ind)=se_cheb;
lanc_sq_err(K_ind)=se_lanczos;
warp_sq_err(K_ind)=se_cheb_warped;
spec_sq_err(K_ind)=se_spec_adapted_ortho;
weighted_ls_sq_err(K_ind)=se_weighted_ls;

% Test f(L)b on random signal b

ntime_cheb=zeros(num_tests,1);
ntime_lanczos=zeros(num_tests,1);
ntime_warped=zeros(num_tests,1);
ntime_spec=zeros(num_tests,1);
ntime_weighted_ls=zeros(num_tests,1);

nmse_cheb=zeros(num_tests,1);
nmse_lanczos=zeros(num_tests,1);
nmse_warped=zeros(num_tests,1);
nmse_spec=zeros(num_tests,1);
nmse_weighted_ls=zeros(num_tests,1);

lanc_param.method='lanczos';
lanc_param.order=K;

for i=1:num_tests
    b=BB(:,i);
    bhat=G2.U'*b;
    gLf_exact=G2.U*(f(G2.e).*bhat);
    tic
    gLf_cheb=gsp_cheby_op(G,c,b); 
    ntime_cheb(i)=toc;
    tic
    gLf_lanczos=gsp_filter(G,f,b,lanc_param);
    ntime_lanczos(i)=toc;
    tic
    gLf_warped=G2.U*(polyval(p_warped,G2.e,s_warped,mu_warped).*bhat);
    ntime_warped(i)=toc;
    tic
    gLf_spec=three_term_recurr_op(G,ab(1:K+2,:),c_spec_adapted_ortho(1:K+1),b,param);
    ntime_spec(i)=toc;
    tic
    gLf_weighted_ls=G2.U*(polyval(p_weighted_lsc,G2.e,s_weighted_ls,mu_weighted_ls).*bhat);
    ntime_weighted_ls(i)=toc;

    nmse_cheb(i)=sum((gLf_cheb-gLf_exact).^2)/sum(gLf_exact.^2);
    nmse_lanczos(i)=sum((gLf_lanczos-gLf_exact).^2)/sum(gLf_exact.^2);
    nmse_warped(i)=sum((gLf_warped-gLf_exact).^2)/sum(gLf_exact.^2);
    nmse_spec(i)=sum((gLf_spec-gLf_exact).^2)/sum(gLf_exact.^2);
    nmse_weighted_ls(i)=sum((gLf_weighted_ls-gLf_exact).^2)/sum(gLf_exact.^2);
    
end

avg_ntime_cheb=mean(ntime_cheb);
avg_ntime_lanczos=mean(ntime_lanczos);
avg_ntime_warped=mean(ntime_warped);
avg_ntime_spec=mean(ntime_spec);
avg_ntime_weighted_ls=mean(ntime_weighted_ls);

avg_nmse_cheb=mean(nmse_cheb);
avg_nmse_lanczos=mean(nmse_lanczos);
avg_nmse_warped=mean(nmse_warped);
avg_nmse_spec=mean(nmse_spec);
avg_nmse_weighted_ls=mean(nmse_weighted_ls);

cheb_sup_err(K_ind)=sup_err_cheb;
cheb_sq_err(K_ind)=se_cheb;
cheb_time(K_ind)=avg_ntime_cheb;
cheb_anmse(K_ind)=avg_nmse_cheb;

lanc_sup_err(K_ind)=sup_err_lanczos;
lanc_sq_err(K_ind)=se_lanczos;
lanc_time(K_ind)=avg_ntime_lanczos;
lanc_anmse(K_ind)=avg_nmse_lanczos;

warp_sup_err(K_ind)=sup_err_cheb_warped;
warp_sq_err(K_ind)=se_cheb_warped;
warp_time(K_ind)=avg_ntime_warped;
warp_anmse(K_ind)=avg_nmse_warped;

spec_sup_err(K_ind)=sup_err_spec_adapted_ortho;
spec_sq_err(K_ind)=se_spec_adapted_ortho;
spec_time(K_ind)=avg_ntime_spec;
spec_anmse(K_ind)=avg_nmse_spec;

weighted_ls_sup_err(K_ind)=sup_err_weighted_ls;
weighted_ls_sq_err(K_ind)=se_weighted_ls;
weighted_ls_time(K_ind)=avg_ntime_weighted_ls;
weighted_ls_anmse(K_ind)=avg_nmse_weighted_ls;

end


cheb_table=table(poly_orders,cheb_sup_err,cheb_sq_err,cheb_anmse,cheb_time)
lanc_table=table(poly_orders,lanc_sup_err,lanc_sq_err,lanc_anmse,lanc_time)
warp_table=table(poly_orders,warp_sup_err,warp_sq_err,warp_anmse,warp_time)
spec_table=table(poly_orders,spec_sup_err,spec_sq_err,spec_anmse,spec_time)
weighted_ls_table=table(poly_orders,weighted_ls_sup_err,weighted_ls_sq_err,weighted_ls_anmse,weighted_ls_time)

% Figure 4: simulated error vs order plot
% % Use both implementations for SAOP
% figure;
% plot4=semilogy(poly_orders,[cheb_anmse,lanc_anmse,warp_anmse,spec_anmse,...
%    weighted_ls_anmse],'-o','LineWidth',1,'MarkerSize',10);
% set(gca,'FontSize',24);
% set(plot4, {'MarkerFaceColor'}, get(plot4,'Color')); 
% hold on;
% legend(plot4,'Chebyshev','Lanczos','Interpolation','SAOP',...
%    'Weighted LS','Location','Southwest');  
% xlabel('K');
% ylabel('$$||f({\bf A}){\bf b}-p_K({\bf A}){\bf b}||^2/||f({\bf A}){\bf b}||^2$$','Interpreter','latex');
% grid on;
% savefig(strcat(graph,'_',filter,'_',btype,'.fig'));
% close all;

% Use one implementation (Weighted LS) for SAOP
figure;
plot4=semilogy(poly_orders,[cheb_anmse,lanc_anmse,warp_anmse,...
   weighted_ls_anmse],'-o','LineWidth',1,'MarkerSize',10);
set(gca,'FontSize',24);
set(plot4, {'MarkerFaceColor'}, get(plot4,'Color')); 
hold on;
legend(plot4,'Chebyshev','Lanczos','Interpolation',...
   'Weighted LS','Location','Southwest');  
xlabel('K');
% ylabel('$$||f({\bf A}){\bf b}-p_K({\bf A}){\bf b}||^2/||f({\bf A}){\bf b}||^2$$','Interpreter','latex');
grid on;
savefig(strcat(graph,'_',filter,'_',btype,'.fig'));
close all;


% Figure 5: square error vs order plot
% % Use both implementations for SAOP
% figure;
% ysq=sum(f(G2.e).^2);
% plot6=semilogy(poly_orders,[cheb_sq_err/ysq,lanc_sq_err/ysq,warp_sq_err/ysq,spec_sq_err/ysq,...
%     weighted_ls_sq_err/ysq],'-o','LineWidth',1,'MarkerSize',10);
% set(gca,'FontSize',24);
% set(plot6, {'MarkerFaceColor'}, get(plot6,'Color')); 
% hold on;
% legend(plot6,'Chebyshev','Lanczos','Interpolation','SAOP',...
%     'Weighted LS','Location','Southwest');  
% xlabel('K');
% ylabel('$$\sum_{\ell=1}^N \left(f(\lambda_{\ell})-p_K(\lambda_{\ell})\right)^2/\sum_{\ell=1}^N f(\lambda_{\ell})^2$$','Interpreter','latex');
% grid on;
% savefig(strcat(graph,'_',filter,'_everr.fig'));
% close all;

% % Use one implementation (Weighted LS) for SAOP
% figure;
% ysq=sum(f(G2.e).^2);
% plot6=semilogy(poly_orders,[cheb_sq_err/ysq,lanc_sq_err/ysq,warp_sq_err/ysq,spec_sq_err/ysq,...
%     weighted_ls_sq_err/ysq],'-o','LineWidth',1,'MarkerSize',10);
% set(gca,'FontSize',24);
% set(plot6, {'MarkerFaceColor'}, get(plot6,'Color')); 
% hold on;
% legend(plot6,'Chebyshev','Lanczos','Interpolation','SAOP',...
%     'Weighted LS','Location','Southwest');  
% xlabel('K');
% ylabel('$$\sum_{\ell=1}^N \left(f(\lambda_{\ell})-p_K(\lambda_{\ell})\right)^2/\sum_{\ell=1}^N f(\lambda_{\ell})^2$$','Interpreter','latex');
% grid on;
% savefig(strcat(graph,'_',filter,'_everr.fig'));
end