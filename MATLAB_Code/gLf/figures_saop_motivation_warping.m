% Figures for motivation and warping example of saop paper

close all;
clear all;
randn('seed', 18); rand('seed', 18)

% fig_motivation1: approximations for f=exp(-\lambda)
G=gsp_erdos_renyi(500,.2);
G2=G; param.force_svd=1;
G2=gsp_compute_fourier_basis(G2,param);

G=gsp_estimate_lmax(G);
G.lmin=G2.e(1)*(G2.e(1));
param.cdf_method='kpm';
param.num_pts=10;
param.num_vec=10;
param.order=30;
G=spectral_cdf_approx2(G,param); 
gi=@(s) G.spectrum_inv_cdf_approx((s-G.lmin)/(G.lmax-G.lmin));

tau=1; f=@(x) exp(-tau*x);
K=5;
c=gsp_cheby_coeff(G,f,K,1000);
y=f(G2.e);
[lsc,s_ls,mu_ls]=polyfit(G2.e,y,K);

y_cheb=gsp_cheby_eval(G2.e,c,[G.lmin,G.lmax]);
errors_cheb=y-y_cheb;
y_ls=polyval(lsc,G2.e,s_ls,mu_ls);
errors_ls=y-y_ls;

xx=0:G.lmax/20000:G.lmax;
xx=xx';
yy=f(xx);
yy_cheb=gsp_cheby_eval(xx,c,[G.lmin,G.lmax]);
yy_ls=polyval(lsc,xx,s_ls,mu_ls);

figure;
plot1=plot(xx,[yy,yy_cheb,yy_ls],'LineWidth',4);
set(gca,'FontSize',24);
legend(plot1,'f','Chebyshev','Discrete LS','Location','Northeast');  
xlabel('\lambda');
grid on;
hold on;
xlabel('\lambda');
ylabel('Function approximations');
plot(G2.e,zeros(G.N,1),'xk','LineWidth',2,'MarkerSize',6,'DisplayName','Eigenvalues');

% fig_motivation2: absolute errors in fig_motivation1
figure;
plot2=semilogy(xx,[abs(yy-yy_cheb),abs(yy-yy_ls)],'LineWidth',4);hold on;
sz=100;
scatter(G2.e,abs(errors_cheb),sz,'Marker','o','MarkerEdgeColor',[ 153 204 255]/256,'MarkerFaceColor',[ 153 204 255]/256);
scatter(G2.e,abs(errors_ls),sz,'Marker','o','MarkerEdgeColor',[255 179 153]/256,'MarkerFaceColor',[255 179 153]/256);
set(gca,'FontSize',24);
set(plot2, {'MarkerFaceColor'}, get(plot2,'Color')); 
hold on;
plot(G2.e,ones(G.N,1),'xk','LineWidth',2,'MarkerSize',6);
legend(plot2,'Chebyshev','Discrete LS','Location','SouthEast');  
xlabel('\lambda');
ylabel('$$|f(\lambda)-p_5(\lambda)|$$','Interpreter','latex');
grid on;

max(abs(errors_cheb))
max(abs(errors_ls))

max(abs(yy-yy_cheb))
max(abs(yy-yy_ls))

% fig_inverse
pts_interp=cos((0:K)*pi/K)'; 
pts_tx_interp=(pts_interp+1)/2;
pts_tx_interp=sort(pts_tx_interp,'ascend');
pts_warped_interp=gi(pts_tx_interp*(G.lmax-G.lmin)+G.lmin);
pts_warped_interp(1)=G.lmin;

figure;
Fbar=@(x)G.spectrum_cdf_approx(x);
fbarparam.plot_eigenvalues=0;
gsp_plot_filter(G,Fbar,fbarparam);
hold on;
plot(pts_warped_interp,zeros(length(pts_warped_interp),1),'xr','LineWidth',...
            4,'MarkerSize',15);
plot(zeros(length(pts_tx_interp),1),pts_tx_interp,'xr','LineWidth',...
            4,'MarkerSize',15);
xlim([0,70]);
ylim([0,1]);
set(gca,'FontSize',24);
set(gca,'XTick',0:10:70);
set(gca,'YTick',0:.1:1);
xlabel('\lambda');


