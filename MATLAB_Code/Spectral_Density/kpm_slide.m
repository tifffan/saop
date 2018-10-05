close all;
clear all;
randn('seed', 18); rand('seed', 18);

N=500;
p=.2;
G=gsp_erdos_renyi(N,p);
G=gsp_estimate_lmax(G);
lm=G.lmax;

order=30;
num_vec=30;
param.num_pts=20;
[cdf_approx,pts,vals]=spectral_cdf_approx(G,num_vec,order,param.num_pts);

param.force_svd=1;
G=gsp_compute_fourier_basis(G,param);

h=@(x) x <= pts(10);
figure;
gsp_plot_filter(G,h);
set(gca,'FontSize',24);
xlim([0,70]);
set(gca,'XTick',0:10:70);

[~, jch] = gsp_jackson_cheby_coeff(0,pts(10),[0,G.lmax], order);
htilde=@(x) gsp_cheby_eval(x,jch,[0,G.lmax]);
hh{1}=h;
hh{2}=htilde;
param.show_sum=0;
figure;
gsp_plot_filter(G,hh,param);
xlim([0,70]);
set(gca,'XTick',0:10:70);
set(gca,'FontSize',24);


for i=2:4
figure;
scatter(pts(1:i),vals(1:i),'r','LineWidth',4);
hold on;
plot(G.e,zeros(G.N,1),'xk','LineWidth',2,'MarkerSize',6);
grid on;
box on;
xlim([0,70]);
ylim([0,1.02]);
xlabel('\lambda');
set(gca,'XTick',0:10:70);
set(gca,'FontSize',24);

figure;
h=@(x) x <= pts(i);
[~, jch] = gsp_jackson_cheby_coeff(0,pts(i),[0,G.lmax], order);
htilde=@(x) gsp_cheby_eval(x,jch,[0,G.lmax]);
hh{1}=h;
hh{2}=htilde;
param.show_sum=0;
gsp_plot_filter(G,hh,param);
xlim([0,70]);
set(gca,'XTick',0:10:70);
set(gca,'FontSize',24);
end


figure;
scatter(pts,vals,'r','LineWidth',4);
hold on;
plot(G.e,zeros(G.N,1),'xk','LineWidth',2,'MarkerSize',6);
grid on;
box on;
xlim([0,70]);
ylim([0,1.02]);
xlabel('\lambda');
set(gca,'XTick',0:10:70);
set(gca,'FontSize',24);

figure;


ee=0:.001:lm;
plot(ee,cdf_approx(ee),'b','LineWidth',4);
hold on;
scatter(pts,vals,'r','LineWidth',4);
plot(G.e,zeros(G.N,1),'xk','LineWidth',2,'MarkerSize',6);
grid on;
box on;
xlim([0,70]);
ylim([0,1.02]);
xlabel('\lambda');
set(gca,'XTick',0:10:70);
set(gca,'FontSize',24);