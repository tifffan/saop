close all;
clear all;
randn('seed', 18); rand('seed', 18)

% Graph
graph='gnp';
switch graph
    case 'gnp'
        N=500;
        p=.2;
        G=gsp_erdos_renyi(N,p);
    case 'sensor'
        N=500;
        G=gsp_david_sensor_network(N);
    case 'grassi'
        N=5000;
        G=gsp_sensor(N);
        G = gsp_create_laplacian(G, 'normalized');
    case 'net25'
        load('0_data/net25.mat');
        A=Problem.A;
        A = A - diag(diag(A)); 
        A(4228,6327) = 1;
        A(6327,4228) = 1;
        G=gsp_graph(A);
        n=9520;
   case 'si2'
        load('0_data/Si2.mat');
        G=struct;
        G.L=Problem.A;
        G.L=(G.L+G.L')/2;
        G.N=size(G.L,1);
        N=769;
    case 'cage9'
        load('0_data/CAGE9.mat');
        G=struct;
        G.L=Problem.A;
        G.L=(G.L+G.L')/2;
        G.N=size(G.L,1);
        N=3534;
    case 'saylr4'
        load('0_data/saylr4.mat');
        A=Problem.A;
        A = A - diag(diag(A)); 
        G=gsp_graph(A);
        N=3564;
    case 'saylr4sc'
        load('0_data/saylr4.mat');
        A=Problem.A;
        A = A - diag(diag(A)); 
        G=gsp_graph(A);
        G.L=G.L/2000;
        N=3564;
end

% Add Time Dimension
T=1000;
G = gsp_jtv_graph(G, T);
G.lmin = 0;

% G2 = G with Exact Eigenvalues
G2=G;
tic
G2=gsp_compute_fourier_basis(G2);

svd_time=toc;

disp('graph spectral decomposition computed');
datetime

G2 = gsp_jtv_graph(G2, T);
G2.lmin=min(G2.e);
omega = fftshift(G.jtv.omega);
[ X, Y ] = meshgrid( omega, G2.e );

disp('time-vertex graph built');
datetime
time_exact=toc;

% Graph Spectral CDF and inverse CDF
tic

G=gsp_estimate_lmax(G);
G.lmin=(abs(min(G2.e))>1e-6)*min(G2.e);  %*(G2.e(1)>1e-6); %%% Replace with estimate

time_cheb=toc;
time_lanc=toc;
time_interp=toc;
time_saop=toc;
time_weightedls=toc;

tic

param.cdf_method='kpm'; % can change to 'kpm' or 'lanczos' or 'ldlt'
param.num_pts=10; % for kpm and ldlt
param.num_vec=10; % for lanczos only
param.order=30; % up to 100 works for kpm and ldlt
G=spectral_cdf_approx2(G,param); 
gi=@(s) G.spectrum_inv_cdf_approx((s-G.lmin)/(G.lmax-G.lmin));

time_interp=time_interp+toc;
time_saop=time_saop+toc;
time_weightedls=time_weightedls+toc;


% Filter
filter = "wave";

r = 100;
lcut = G.lmax/2;
wcut = max(G.jtv.omega)/2; % change to a fraction of max 

switch filter
    case "lowpass" % ideal lowpass filter
        hlp = @( lambda,omega ) double(and(abs(lambda)<lcut,abs(omega)<wcut));
    case "approx-lowpass" % a smooth approximation to lowpass filter
        hlp = @(l,w) (1-1./(1+exp(-r*(l-lcut)))).* (1-1./(1+exp(-r*(abs(w)-wcut))));
    case "wave" % wave filter from Grassi's paper
        hlp = @( l,w ) (exp(-abs(2*pi*abs(w)-acos(1-l/G.lmax/2)).^2*r))*sqrt(T)/2;
    case "joint-lowpass" % denoising filter from Grassi's paper
        t1=1;
        t2=1;
        hlp = @(l,w) 1 ./ (1 + t1*l + 2*t2*(1-cos(2*pi*w)));
end

Hlp = hlp( Y,X );

% Plot of filter
y=min(G.jtv.omega):.05:max(G.jtv.omega);
x=0:.05:G.lmax;
figure;surf(y,x,hlp(x',y));
savefig(strcat(graph,'_',filter,'.fig'));
close all;

% Testing signal
btype='constant_spectral';

switch btype
    case 'constant'
        num_tests=1;
        X=ones(G.N,T);
    case 'constant_spectral'
        B=(1/sqrt(G.N*T))*ones(G.N,T);    
        X=(fft((G2.U*B)'))'; %want components of different time frequency and different spectral frequency        
    case 'randn_spectral'
         B=(1/sqrt(G.N*T))*randn(G.N,T);
         X=(fft((G2.U*B)'))';
    case 'rand_spectral'
        B=(1/sqrt(G.N*T))*rand(G.N,T);
        X=(fft((G2.U*B)'))';
    case 'randn'
        X=(1/sqrt(G.N*T))*randn(G.N,T);
    case 'rand'
        X=(1/sqrt(G.N*T))*rand(G.N,T);
    otherwise
        error('unknown b type');
end

exptime_cheb=0;
exptime_lanc=0;
exptime_interp=0;
exptime_saop=0;
exptime_weightedls=0;

tic
Xhat=gsp_tft(G,X);
exptime_exact=svd_time+toc;
exptime_cheb=exptime_cheb+toc;
exptime_lanc=exptime_lanc+toc;
exptime_interp=exptime_interp+toc;
exptime_saop=exptime_saop+toc;
exptime_weightedls=exptime_weightedls+toc;

fAXhat_exact=zeros(G.N,T);
fAXhat_cheb=zeros(G.N,T);
fAXhat_lanc=zeros(G.N,T);
fAXhat_interp=zeros(G.N,T);
fAXhat_saop=zeros(G.N,T);
fAXhat_weightedls=zeros(G.N,T);

% Degrees of polynomial approximation
orders = [1 5:5:30];
%orders = 1:15; % for wave filter?

% Initializing error list for each order: K -> error
err_lowpass_cheb = zeros(numel(orders),1);
err_lowpass_lanc = zeros(numel(orders),1);
err_lowpass_interp = zeros(numel(orders),1);
err_lowpass_saop = zeros(numel(orders),1);
err_lowpass_weightedls = zeros(numel(orders),1);

methods=["Chebyshev","Lanczos","Interpolation","SAOP","WeightedLS"];

error_table=zeros(numel(orders),5);
max_error_table=zeros(numel(orders),5);

for m=1:numel(orders)
    K=orders(m);
    
    fset={};

cheb_coef_mat=zeros(K+1,T);
interp_coef_mat=zeros(K+1,T);
interp_s_set={};
interp_mu_mat=zeros(2,T);
saop_coef_mat=zeros(K+1,T);
weightedls_coef_mat=zeros(K+1,T);
weightedls_s_set={};
weightedls_mu_mat=zeros(2,T);

% Initialize error matrices within a given order: same size as signal
error_mat_cheb=zeros(N,T);
error_mat_lanc=zeros(N,T);
error_mat_interp=zeros(N,T);
error_mat_saop=zeros(N,T);
error_mat_weightedls=zeros(N,T);


% Initialize warped points and discrete measure based on spectral pdf

% 1) Starting points for interpolation
tic

pts_interp=cos((0:K)*pi/K)'; 
pts_tx_interp=(pts_interp+1)*(G.lmax-G.lmin)/2+G.lmin;
pts_tx_interp=sort(pts_tx_interp,'ascend');
pts_warped_interp=gi(pts_tx_interp);
pts_warped_interp(1)=G.lmin;

time_interp=time_interp+toc;


% 2) Matrix/Spectrum adapted orthogonal polynomials
tic

mop_param=struct;
mop_param.grid_order=100;
mop_param.init_poly_order=30;
mop_param.num_vec=10;
mop_param.absc_method='linear'; %linear, even(n-2 linearly spaced + both ends), warp, spline
mop_param.weights_method='pdf'; %count, pdf2(from kpm paper)
[absc,weights]=gen_absc_weights(G,K,mop_param); 
[ab,Pi]=matrix_adapted_ortho_poly(absc,weights,K);

time_saop=time_saop+toc;

% 3) Weighted LS fitting
tic

mop_param.grid_order=100;
mop_param.absc_method='linear';
mop_param.weights_method='pdf';
[absc_weighted_ls,weights_weighted_ls]=gen_absc_weights(G,K,mop_param);

time_weightedls=time_weightedls+toc;

    
for t=1:T
    f=@(lambda) hlp(lambda, omega(t));
    fset{1,t}=f;
    
% For each time frequency, we apply the following methods in 1-D
% Chebyshev
tic

c=gsp_cheby_coeff(G,f,K,K+1); % expect: replicate of Grassi
cheb_coef_mat(:,t)=c;
time_cheb=time_cheb+toc;

% Lanczos
tic

lanc_param.method='lanczos';
lanc_param.order=K;
x=sum(G2.U')';
%[V,H]=lanczos(G.L,K,x);
% e=eig(H);

time_lanc=time_lanc+toc;

% Chebyshev interpolation on warped points
tic

[p_warped,s_warped,mu_warped]=polyfit(pts_warped_interp,f(pts_warped_interp),K);
p_warped_fun=@(x) polyval(p_warped,x,s_warped,mu_warped); 
p_warped_c=gsp_cheby_coeff(G,p_warped_fun,K,1000);
interp_coef_mat(:,t)=p_warped;
interp_s_set{1,t}=s_warped;
interp_mu_mat(:,t)=mu_warped;
time_interp=time_interp+toc;

% Matrix/Spectrum adapted orthogonal polynomials
tic

c_spec_adapted_ortho=matrix_adapted_poly_coeff(G, f, absc, weights/sum(weights), Pi, K);
saop_coef_mat(:,t)=c_spec_adapted_ortho;
time_saop=time_saop+toc;

% Weighted LS fitting
tic

[p_weighted_lsc,s_weighted_ls,mu_weighted_ls]=weighted_polyfit(absc_weighted_ls,f(absc_weighted_ls),weights_weighted_ls,K);
weightedls_coef_mat(:,t)=p_weighted_lsc;
weightedls_s_set{1,t}=s_weighted_ls;
weightedls_mu_mat(:,t)=mu_weighted_ls;
time_weightedls=time_weightedls+toc;
end


for t=1:T
    
    b=Xhat(:,t);
    bhat=G2.U'*b;
    f=fset{1,t};
    
    tic
    fAXhat_exact(:,t)=G2.U*(f(G2.e).*bhat);
    exptime_exact=exptime_exact+toc;
    
    tic
    fAXhat_cheb(:,t)=gsp_cheby_op(G,cheb_coef_mat(:,t),b);
    exptime_cheb=exptime_cheb+toc;

    tic
    fAXhat_lanc(:,t)=gsp_filter(G,f,b,lanc_param);
    exptime_lanc=exptime_lanc+toc;

    tic
    fAXhat_interp(:,t)=G2.U*(polyval(interp_coef_mat(:,t),G2.e,interp_s_set{1,t},interp_mu_mat(:,t)).*bhat);
    exptime_interp=exptime_interp+toc;
    
    tic
    fAXhat_saop(:,t)=three_term_recurr_op(G,ab(1:K+2,:),saop_coef_mat(:,t),b,param);
    exptime_saop=exptime_saop+toc;
    
    tic
    fAXhat_weightedls(:,t)=G2.U*(polyval(weightedls_coef_mat(:,t),G2.e,weightedls_s_set{1,t},weightedls_mu_mat(:,t)).*bhat);
    exptime_weightedls=exptime_weightedls+toc;
    
end

fAX_exact=gsp_itft(G,fAXhat_exact);
fAX_cheb=gsp_itft(G,fAXhat_cheb);
fAX_lanc=gsp_itft(G,fAXhat_lanc);
fAX_interp=gsp_itft(G,fAXhat_interp);
fAX_saop=gsp_itft(G,fAXhat_saop);
fAX_weightedls=gsp_itft(G,fAXhat_weightedls);

error_table(m,:)=[norm((fAX_cheb-fAX_exact),'fro'),norm((fAX_lanc-fAX_exact),'fro'),...
norm((fAX_interp-fAX_exact),'fro'),norm((fAX_saop-fAX_exact),'fro'),norm((fAX_weightedls-fAX_exact),'fro')];

error_table(m,:)=error_table(m,:)./ norm(fAX_exact,'fro');

% max_error_table(m,:)=[max(max(abs(fAX_cheb-fAX_exact))),max(max(abs(fAX_lanc-fAX_exact))),...
%     max(max(abs(fAX_interp-fAX_exact))),max(max(abs(fAX_saop-fAX_exact))),max(max(abs(fAX_weightedls-fAX_exact)))];
% 
% max_error_table(m,:)=max_error_table(m,:)./ max(max(abs(fAX_exact)));

end

% relative error in frobenius norm, both saop and weighted ls displayed
figure;
plot=semilogy(orders,[error_table(:,1),error_table(:,2),error_table(:,3),...
    error_table(:,4), error_table(:,5)],'-o','LineWidth',1,'MarkerSize',10);
set(gca,'FontSize',20);
set(plot, {'MarkerFaceColor'}, get(plot,'Color')); 
hold on;
legend(plot,'Chebyshev','Lanczos','Interpolation','SAOP','Weighted LS','Location','Northwest');  
%ylabel('$$\sum_{\ell=1}^N \left(f(\lambda_{\ell})-p_K(\lambda_{\ell})\right)^2/\sum_{\ell=1}^N f(\lambda_{\ell})^2$$','Interpreter','latex');
grid on;
xlabel('K');
% xlabel('Polynomial Order');
% ylabel('Frobenius Norm of Error');
% title(strcat('Graph=',graph,', Filter=',filter));
savefig(strcat(graph,'_',filter,'_fro_saop.fig'));
close all;

% relative error in frobenius norm, weighted ls displayed
figure;
plot=semilogy(orders,[error_table(:,1),error_table(:,2),error_table(:,3),...
     error_table(:,5)],'-o','LineWidth',1,'MarkerSize',10);
set(gca,'FontSize',20);
set(plot, {'MarkerFaceColor'}, get(plot,'Color')); 
hold on;
legend(plot,'Chebyshev','Lanczos','Interpolation','Weighted LS','Location','Northwest');  
%ylabel('$$\sum_{\ell=1}^N \left(f(\lambda_{\ell})-p_K(\lambda_{\ell})\right)^2/\sum_{\ell=1}^N f(\lambda_{\ell})^2$$','Interpreter','latex');
grid on;
xlabel('K');
% xlabel('Polynomial Order');
% ylabel('Frobenius Norm of Error');
% title(strcat('Graph=',graph,', Filter=',filter));
savefig(strcat(graph,'_',filter,'_fro.fig'));
close all;

% relative error in max norm
% figure;
% plot=semilogy(orders,[max_error_table(:,1),max_error_table(:,2),max_error_table(:,3),...
%     max_error_table(:,4), max_error_table(:,5)],'-o','LineWidth',1,'MarkerSize',10);
% set(gca,'FontSize',20);
% set(plot, {'MarkerFaceColor'}, get(plot,'Color')); 
% hold on;
% legend(plot,'Chebyshev','Lanczos','Interpolation','SAOP','Weighted LS','Location','Northwest');  
% %ylabel('$$\sum_{\ell=1}^N \left(f(\lambda_{\ell})-p_K(\lambda_{\ell})\right)^2/\sum_{\ell=1}^N f(\lambda_{\ell})^2$$','Interpreter','latex');
% grid on;
% xlabel('K');
% % xlabel('Polynomial Order');
% % ylabel('Maximum Error');
% % title(strcat('Graph=',graph,', Filter=',filter));
% savefig(strcat(graph,'_',filter,'_max_saop.fig'));
% close all;