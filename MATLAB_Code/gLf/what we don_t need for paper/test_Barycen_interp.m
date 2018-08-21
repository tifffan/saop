close all;
clear all;
randn('seed', 18); rand('seed', 18);

%%  Graph data
N = 500; % minnesota is larger
graph='gnp';

switch graph
    case 'gnp'
        p=.2;
        G=gsp_erdos_renyi(N,p);
    case 'sensor'
        G=gsp_david_sensor_network(N);
    case 'minnesota'
        G=gsp_minnesota(1);
    case 'comet'
        G=gsp_comet(N,50);
    otherwise
        error('graph type not recognized');
end

G=gsp_estimate_lmax(G);
lmax_est=G.lmax;
param.num_pts=20; % for approximating spectral cdf 
G=gsp_spectrum_cdf_approx(G,param); % incorporate alternative ways to compute this (e.g., KPM for DoS)

%% Filters - testing mostly for analytic filters. If we have a
% discontinuity, we need more interpolation points near the discontinuity
filter='inverse';

switch filter
    case 'inverse'
        tau=1; 
        f=@(x) tau./(tau+x);
    case 'ideal'
        f=@(x) (x<=(G.lmax/2));
    case 'heat'
        tau=1;
        f=@(x) exp(-tau*x);
    otherwise
        error('filter type not recognized');
end


%% Approximations
% Poly approx order
K=10;

% Chebyshev
h=chebfun(@(s) f(s),[0,G.lmax],'splitting','on');
c=chebcoeffs(h,K+1); 
c(1)=c(1)*2;

% % Legendre
% pleg = polyfit(h,K);
% plegc=pleg.coeffs;
% plegc(1)=plegc(1)*2;

% Warped via. LDLT CDF
g=chebfun(@(s) G.spectrum_cdf_approx(s),[0,G.lmax],'splitting','on'); 
gi=inv(g,'splitting','on'); % warping function is the inverse of the spectral CDF


use_Chebypts=1;
if(use_Chebypts)
chebpts=cos((0:K)*pi/K); 
tx=(chebpts+1)/2; 
else
tx=linspace(0, 1, K);
end
    
chebpts_warped=gi(tx);


force_max=1;
if force_max
    chebpts_warped=sort(chebpts_warped,'ascend');
    chebpts_warped(end)=G.lmax;
end

dom=domain(0,G.lmax);
p_warped=chebfun.interp1(chebpts_warped,f(chebpts_warped),dom);
p_warped_c=p_warped.coeffs;
p_warped_c(1)=p_warped_c(1)*2;


% Warped via. Lanczos CDF
Mdeg = 40; % Number of Lanczos Steps
nvec =10; % Number of starting vectors 
x=0:.1:G.lmax; % discretized points in the interval [lmin,lMax]


Lan_cheb_warped=InvLanczosCDOS(G.L, Mdeg, nvec,G.lmax,x',tx');

force_max=1;
if force_max
    Lan_cheb_warped=sort(Lan_cheb_warped,'ascend');
    Lan_cheb_warped(end)=G.lmax;
end

dom=domain(0,G.lmax);
Lan_warped=chebfun.interp1(Lan_cheb_warped,f(Lan_cheb_warped),dom);
Lan_warped_c=Lan_warped.coeffs;
Lan_warped_c(1)=Lan_warped_c(1)*2;

% Least squares (assumes full knowledge of eigenvalues; just for comparison
% to ideal; not scalable)
G=gsp_compute_fourier_basis(G);
y=f(G.e);
lsc=polyfit(G.e,y,K);

%% Compute polynomial approximation values at the actual eigenvalues and 
% the corresponding squared error
y_cheb=gsp_cheby_eval(G.e,c,[0,G.lmax]);
errors_cheb=y-y_cheb;
se_cheb=sum(errors_cheb.^2)

% y_leg=gsp_cheby_eval(G.e,plegc,[0,G.lmax]);
% errors_leg=y-y_leg;
% se_leg=sum(errors_leg.^2)

y_cheb_warped=p_warped(G.e);
errors_cheb_warped=y-y_cheb_warped;
se_cheb_warped=sum(errors_cheb_warped.^2)

y_Lan_warped=Lan_warped(G.e);
errors_Lan_warped=y-y_Lan_warped;
se_Lan_warped=sum(errors_Lan_warped.^2)

y_ls=polyval(lsc,G.e);
errors_ls=y-y_ls;
se_ls=sum(errors_ls.^2)

%% Plots
xmax=max(lmax_est,G.lmax);
delta=xmax/500;
xx=0:delta:xmax;
xx=xx';
yy=f(xx);

yy_cheb=gsp_cheby_eval(xx,c,[0,G.lmax]);
%yy_leg=gsp_cheby_eval(xx,plegc,[0,G.lmax]);
yy_cheb_warped=p_warped(xx);
yy_Lan_warped=Lan_warped(xx);
yy_ls=polyval(lsc,xx);

figure;
p1=plot(xx,[yy,gsp_cheby_eval(xx,c,[0,G.lmax]),yy_cheb_warped,yy_Lan_warped,yy_ls],'LineWidth',2);
set(gca,'FontSize',24)
legend(p1,'f','Chebyshev','Warped Chebyshev','Warped Lan-Cheb','Discrete LS');
grid on;
hold on;
xlabel('\lambda');
ylabel('Filter approximations');
plot(G.e,zeros(G.N,1),'xk','LineWidth',2,'MarkerSize',6);

figure;
p2=semilogy(xx,[abs(yy-yy_cheb),abs(yy-yy_cheb_warped),abs(yy-yy_Lan_warped),abs(yy-yy_ls)],'LineWidth',2);
set(gca,'FontSize',24)
hold on;
plot(G.e,ones(G.N,1),'xk','LineWidth',2,'MarkerSize',6);
legend(p2,'Chebyshev','Warped Chebyshev','Warped Lan-Cheb','Discrete LS','Location','SouthEast');
xlabel('\lambda');
ylabel('$$|\hat{f}(\lambda)-f(\lambda)|$$','Interpreter','latex');
grid on;


% Test f(L)b on random signal b
num_tests=10;
nmse_cheb=zeros(num_tests,1);
%nmse_leg=zeros(num_tests,1);
nmse_warped=zeros(num_tests,1);
nmse_Lanwarped=zeros(num_tests,1);
nmse_ls=zeros(num_tests,1);

for i=1:num_tests
    b=rand(G.N,1);
    gLf_exact=G.U*diag(f(G.e))*G.U'*b;
    gLf_cheb=gsp_cheby_op(G,c,b);
   % gLf_leg=gsp_cheby_op(G,plegc,b);
    gLf_warped=gsp_cheby_op(G,p_warped_c,b);
    gLf_Lanwarped=gsp_cheby_op(G,Lan_warped_c,b);
    gLf_ls=G.U*diag(polyval(lsc,G.e))*G.U'*b;
    
    nmse_cheb(i)=sum((gLf_cheb-gLf_exact).^2)/sum(gLf_exact.^2);
   % nmse_leg(i)=sum((gLf_leg-gLf_exact).^2)/sum(gLf_exact.^2);
    nmse_warped(i)=sum((gLf_warped-gLf_exact).^2)/sum(gLf_exact.^2);
    nmse_Lanwarped(i)=sum((gLf_Lanwarped-gLf_exact).^2)/sum(gLf_exact.^2);
    nmse_ls(i)=sum((gLf_ls-gLf_exact).^2)/sum(gLf_exact.^2);
end

avg_nmse_cheb=mean(nmse_cheb)
%avg_nmse_leg=mean(nmse_leg)
avg_nmse_warped=mean(nmse_warped)
avg_nmse_Lanwarped=mean(nmse_Lanwarped)
avg_nmse_ls=mean(nmse_ls)