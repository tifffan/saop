close all;
clear all;
randn('seed', 18); rand('seed', 18)

% Graph
N=500;
T=100;
p=.2;
G=gsp_erdos_renyi(N,p);

% Add Time Dimension
G = gsp_jtv_graph(G, T);
G.lmin = 0;

% G2 = G with Exact Eigenvalues
G2=G;
G2=gsp_compute_fourier_basis(G2);

% Graph Spectral CDF and inverse CDF
G=gsp_estimate_lmax(G);
param.cdf_method='kpm'; % can change to 'kpm' or 'lanczos' or 'ldlt'
param.num_pts=10; % for kpm and ldlt
param.num_vec=10; % for lanczos only
param.order=30; % up to 100 works for kpm and ldlt
G=spectral_cdf_approx2(G,param);

% Filter
filter_type = "lowpass";

r = 100;
lcut = G.lmax/2;
wcut = max(G.jtv.omega)/2; % change to a fraction of max 

switch filter_type
    case "lowpass" % ideal lowpass filter
        hlp = @( lambda,omega ) double(and(abs(lambda)<lcut,abs(omega)<wcut));
    case "approx-lowpass" % a smooth approximation to lowpass filter
        hlp = @(l,w) (1-1./(1+exp(-r*(l-lcut)))).* (1-1./(1+exp(-r*(abs(w)-wcut))));
    case "wave"
        hlp = @(l,w) (exp(-abs(2*pi*abs(w)-acos(1-l/G.lmax/2)).^2*r))*sqrt(T)/2;
    case "joint-lowpass"
        t1 = 1;
        t2 = 1;
        hlp = @(l,w) 1 ./ (1 + t1*l + 2*t2*(1-cos(2*pi*w)));
end

[X, Y] = meshgrid(G.jtv.omega, G2.e);
Hlp = hlp(Y, X);

% Plot filter
y=min(G.jtv.omega):.001:max(G.jtv.omega);
x=0:.001:G.lmax;
y2 = 2-2*cos(2*pi*y);
%figure;mesh(y,x,hlp(x',y));colormap autumn;
figure; mesh(y2,x, hlp(x',y), 'EdgeColor','interp','FaceColor','interp');
colormap autumn; hold on;
scatter3( 2-2*cos(2*pi* vec(X)), vec(Y), vec(Hlp), 'filled');
set(gca,'fontsize',24)
xlabel('$\lambda_R$','Interpreter','latex');
ylabel('$\lambda_G$','Interpreter','latex');
zlabel('$h(\lambda_G, \lambda_R)$','Interpreter','latex');