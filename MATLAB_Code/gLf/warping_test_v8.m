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
graph='net25';
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
        load('/Users/lfan/Documents/MATLAB/git/spectral-warping/MATLAB_Code/Data/net25.mat');
        %load('/Users/lifan/Desktop/Research/git/spectral-warping/MATLAB_Code/Data/net25.mat');
        %load('/Users/davidshuman/Dropbox/Current Research Work/MCSFB/Shuni_Thesis/GitHub/mcsfb2018/net25 graph data/net25.mat');
        A=Problem.A;
        A = A - diag(diag(A)); 
        A(4228,6327) = 1;
        A(6327,4228) = 1;
        G=gsp_graph(A);
        G.L=G.L./20;
    otherwise
        error('graph type not recognized');
end

G2=G;
tic
G2=gsp_compute_fourier_basis(G2);
time_exact_spec=toc;

%cdf_approx
G=gsp_estimate_lmax(G);
param.num_pts=50; % for approximating spectral cdf 

param.cdf_method='kpm'; % can change to 'kpm' or 'lanczos' or 'ldlt'
param.num_vec=30;
param.order=100;
%param.pts=linspace(.1,G.lmax,param.num_pts);
G=spectral_cdf_approx2(G,param); 

% CDF and inverse CDF
g=chebfun(@(s) G.lmax*G.spectrum_cdf_approx(s),[0,G.lmax],'splitting','on'); 
gi=inv(g,'splitting','on'); % warping function is the inverse of the spectral CDF

% Filters - testing mostly for analytic filters. If we have a
% discontinuity, we need more interpolation points near the discontinuity
filter='mid';

switch filter
    case 'inverse'
        tau=.5; 
        f=@(x) tau./(tau+x);
    case 'low'
        f=@(x) (x<=(G.lmax/2));
        %f=@(x) (x<=gi(G.lmax/2));
    case 'mid'
        f=@(x) (x>=(G.lmax/3) & x<=(G.lmax*2/3)); % usually needs higher polynomial order
    case 'high'
        f=@(x) (x>=(G.lmax/2));
    case 'heat'
        tau=1;
        f=@(x) exp(-tau*x);
    otherwise
        error('filter type not recognized');
end


%--------------------------------------------------------------------------
% Poly approx order
K=100;
start_pts_interp='cheb'; % K pts for warped chebyshev interpolation
start_pts_ls='cheb'; % G.N/10 pts for warped LS fitting

% Chebyshev
h=chebfun(@(s) f(s),[0,G.lmax],'splitting','on');
c=chebcoeffs(h,K+1); 
c(1)=c(1)*2;

% Legendre
pleg = polyfit(h,K);
plegc=pleg.coeffs;
plegc(1)=plegc(1)*2;

% Regular Lanczos
lanc_param.method='lanczos';
lanc_param.order=K;           % set K = larger constant for testing

% % Warped
% g=chebfun(@(s) G.lmax*G.spectrum_cdf_approx(s),[0,G.lmax],'splitting','on'); 
% gi=inv(g,'splitting','on'); % warping function is the inverse of the spectral CDF

% Starting points for interpolation
%start_pts_interp='cheb';
switch start_pts_interp
    case 'cheb'
        pts_interp=cos((0:K)*pi/K)'; 
        pts_tx_interp=(pts_interp+1)*G.lmax/2;
        pts_tx_interp=sort(pts_tx_interp,'ascend');
    case 'even'
        pts_tx_interp=linspace(0,G.lmax,K+1);
end

pts_warped_interp=gi(pts_tx_interp);
pts_warped_interp(1)=0;

% Starting points for LS fitting
%npts_ls=ceil(G.N/10);
npts_ls=max(ceil(G.N/10),K+1);
%start_pts_ls='even';
switch start_pts_ls
    case 'cheb'
        pts_ls=cos((0:npts_ls)*pi/npts_ls)';
        pts_tx_ls=(pts_ls+1)*G.lmax/2;    
        pts_tx_ls=sort(pts_tx_ls,'ascend');
    case 'even'
        pts_tx_ls=linspace(0,G.lmax,npts_ls);
    case 'unif'
        pts_tx_ls=G.lmax*rand(npts_ls,1);
        pts_tx_ls=sort(pts_tx_ls,'ascend');
    otherwise
        error('grid type not recognized');
end

pts_warped_ls=gi(pts_tx_ls);
pts_warped_ls(1)=0;

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
dom=domain(0,G.lmax);
p_warped=chebfun.interp1(pts_warped_interp,f(pts_warped_interp),dom);
p_warped_c=p_warped.coeffs;
p_warped_c(1)=p_warped_c(1)*2;

% Legendre interpolation on warped points
%p_warped_c_leg=polyfit(pts_warped_interp,f(pts_warped_interp),K);

% Piecewise cubic hermite polynomial interpolation + Chebyshev approximation
p_warped_pchip=pchip(pts_warped_interp,f(pts_warped_interp));
p_warped_pchip_fun=@(x) pchip(pts_warped_interp,f(pts_warped_interp),x); 
p_warped_pchip_interp_c = gsp_cheby_coeff(G, p_warped_pchip_fun, K, 2000);

% Include damping factors against Gibbs oscillations
damping='jackson';
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
        p_warped_pchip_interp_c=p_warped_pchip_interp_c.*gamma;
    case  'sigma'
        sigma=ones(K+1,1);
        for k=1:K
            sigma(k+1,1)=sin(k*pi/(K+1))/(k*pi/(K+1));
            % Damping factors calculated from 3.1, Napoli et.al, Efficient 
            % estimation of eigenvalue counts in an interval
        end
        p_warped_c=p_warped_c.*sigma;
        p_warped_pchip_interp_c=p_warped_pchip_interp_c.*sigma;
    case 'none'
    otherwise
        error('damping type not recognized');
end

% Many sample points + LS fitting
p_warped_lsc=polyfit(pts_warped_ls,f(pts_warped_ls),K); 


% Figure 1: warped points on approximated CDF
figure;
Fbar=@(x)G.lmax*G.spectrum_cdf_approx(x);
fbarparam.plot_eigenvalues=0;
gsp_plot_filter(G,Fbar,fbarparam);
hold on;
plot(pts_warped_interp,zeros(length(pts_warped_interp),1),'xr','LineWidth',...
            4,'MarkerSize',15);
plot(zeros(length(pts_warped_interp),1),pts_tx_interp,'xr','LineWidth',...
            4,'MarkerSize',15);
%xlim([0,70]);
%ylim([0,70]);
set(gca,'FontSize',24);
%set(gca,'XTick',0:10:70);
%set(gca,'YTick',0:10:70);

% Matrix/Spectrum adapted orthogonal polynomials
mop_param=struct;
mop_param.grid_order=1000;
mop_param.init_poly_order=100;

[ab,absc,weights,Pi]= matrix_adapted_ortho_poly(G,K,mop_param);
c_spec_adapted_ortho = matrix_adapted_poly_coeff(G, f, absc', weights', Pi, K);

% Least squares (assumes full knowledge of eigenvalues; just for comparison
% to ideal; not scalable)
y=f(G2.e);
lsc=polyfit(G2.e,y,K);
% Sampling from points, then LS fitting, for higher stability
%J=G.N/10;
%I=randperm(G.N);
%e_sub=sort(G2.e(I(1:J)),'ascend');
%lsc=polyfit(e_sub,f(e_sub),K);

% Compute polynomial approximation values at the actual eigenvalues and 
% the corresponding superior and squared error
leaveout_end=0;    % Option to leave out the last 1/4 of spectrum

y_cheb=gsp_cheby_eval(G2.e,c,[0,G.lmax]);
errors_cheb=y-y_cheb;
if leaveout_end
    errors_cheb=errors_cheb(1:ceil(3*G.N/4));
end
sup_err_cheb=max(abs(errors_cheb));
se_cheb=sum(errors_cheb.^2);

y_leg=gsp_cheby_eval(G2.e,plegc,[0,G.lmax]);
errors_leg=y-y_leg;
if leaveout_end
    errors_leg=errors_leg(1:ceil(3*G.N/4));
end
sup_err_leg=max(abs(errors_leg));
se_leg=sum(errors_leg.^2);

y_lanczos=G2.U'*gsp_filter(G,f,sum(G2.U')',lanc_param);
errors_lanczos=y-y_lanczos;
if leaveout_end
    errors_lanczos=errors_lanczos(1:ceil(3*G.N/4));
end
sup_err_lanczos=max(abs(errors_lanczos));
se_lanczos=sum(errors_lanczos.^2);

y_cheb_warped=p_warped(G2.e); % original warped chebyshev
errors_cheb_warped=y-y_cheb_warped;
if leaveout_end
    errors_cheb_warped=errors_cheb_warped(1:ceil(3*G.N/4));
end
sup_err_cheb_warped=max(abs(errors_cheb_warped));
se_cheb_warped=sum(errors_cheb_warped.^2);

y_cheb_warped_pchip=ppval(p_warped_pchip,G2.e); % piecewise cubic hermite interpolation (ideal)
errors_cheb_warped_pchip=y-y_cheb_warped_pchip;
if leaveout_end
    errors_cheb_warped_pchip=errors_cheb_warped_pchip(1:ceil(3*G.N/4));
end
sup_err_cheb_warped_pchip=max(abs(errors_cheb_warped_pchip));
se_cheb_warped_pchip=sum(errors_cheb_warped_pchip.^2);

y_cheb_warped_pchip_interp=gsp_cheby_eval(G2.e,p_warped_pchip_interp_c,[0,G.lmax]);  % approx by chebyshev
errors_cheb_warped_pchip_interp=y-y_cheb_warped_pchip_interp;
if leaveout_end
    errors_cheb_warped_pchip_interp=errors_cheb_warped_pchip_interp(1:ceil(3*G.N/4));
end
sup_err_cheb_warped_pchip_interp=max(abs(errors_cheb_warped_pchip_interp));
se_cheb_warped_pchip_interp=sum(errors_cheb_warped_pchip_interp.^2);

y_cheb_warped_ls=polyval(p_warped_lsc,G2.e); % many sample points warped + ls fitting
errors_cheb_warped_ls=y-y_cheb_warped_ls;
if leaveout_end
    errors_cheb_warped_ls=errors_cheb_warped_ls(1:ceil(3*G.N/4));
end
sup_err_cheb_warped_ls=max(abs(errors_cheb_warped_ls));
se_cheb_warped_ls=sum(errors_cheb_warped_ls.^2);

y_spec_adapted_ortho=three_term_eval(G,G2.e,ab,c_spec_adapted_ortho);
errors_spec_adapted_ortho=y-y_spec_adapted_ortho;
if leaveout_end
    errors_spec_adapted_ortho=errors_spec_adapted_ortho(1:ceil(3*G.N/4));
end
sup_err_spec_adapted_ortho=max(abs(errors_spec_adapted_ortho));
se_spec_adapted_ortho=sum(errors_spec_adapted_ortho.^2);

y_ls=polyval(lsc,G2.e);
errors_ls=y-y_ls;
if leaveout_end
    errors_ls=errors_ls(1:ceil(3*G.N/4));
end
sup_err_ls=max(abs(errors_ls));
se_ls=sum(errors_ls.^2);

Method = ["Chebyshev"; "Legendre"; "Regular Lanczos";"Warped Cheby";"Warped PCHIP+Cheby";...
    "Warped LS";"Spectrum-Adapted Ortho. Poly.";"Discrete LS";"Exact"];
Sup_Err = [sup_err_cheb;sup_err_leg;sup_err_lanczos;sup_err_cheb_warped;...
    sup_err_cheb_warped_pchip_interp;...
    sup_err_cheb_warped_ls;sup_err_spec_adapted_ortho;sup_err_ls;0];
Sq_Err = [se_cheb;se_leg;se_lanczos;se_cheb_warped;...
    se_cheb_warped_pchip_interp;se_cheb_warped_ls;se_spec_adapted_ortho;se_ls;0];

% Plots
xmax=G.lmax;
delta=xmax/5000;
xx=0:delta:xmax;
xx=xx';
yy=f(xx);

yy_cheb=gsp_cheby_eval(xx,c,[0,G.lmax]);
yy_leg=gsp_cheby_eval(xx,plegc,[0,G.lmax]);
yy_cheb_warped=p_warped(xx);              % chebyshev on warped pts
%yy_leg_warped=polyval(p_warped_c_leg,xx); % legendre on warped pts (should be identical to above)
yy_cheb_warped_pchip=ppval(p_warped_pchip,xx);       % piecewise cubic interpolation, ideal
yy_cheb_warped_pchip_interp=gsp_cheby_eval(xx,p_warped_pchip_interp_c,[0,G.lmax]); % approx by chebyshev
yy_cheb_warped_ls=polyval(p_warped_lsc,xx);
yy_spec_adapted_ortho=three_term_eval(G,xx,ab,c_spec_adapted_ortho);
yy_ls=polyval(lsc,xx);

% Figure 2: approximations for function g
figure;
morecolors={[0 0.4470 0.7410]; [0.8500 0.3250 0.0980]; [0.9290 0.6940 0.1250]; [0.4940 0.1840 0.5560]; [0.4660 0.6740 0.1880]; [0.3010 0.7450 0.9330]; [0.6350 0.0780 0.1840]; [0.308 0.785 0.636]; [1 0.6 0.8]};
p1=plot(xx,[yy_cheb,pleg(xx),yy_cheb_warped,yy_cheb_warped_pchip,yy_cheb_warped_pchip_interp,yy_cheb_warped_ls,yy_spec_adapted_ortho,yy_ls],'LineWidth',4);
set(p1,{'Color'},morecolors(1:8));
set(gca,'FontSize',20);
legend(p1,'Chebyshev','Legendre','Warped Chebyshev','Warped PCHIP','Warped PCHIP+Chebyshev',...
     'Warped LS','Spectrum-Adapted Ortho. Poly.','Discrete LS','Location','NorthEast');
grid on;
hold on;
xlabel('\lambda');
ylabel('Filter approximations');
plot(G2.e, y_lanczos, 'LineWidth',4,'DisplayName','Regular Lanczos','Color',[1 0.6 0.8]);
plot(xx, yy, 'LineWidth',4,'DisplayName','g','Color',[0.9 0.5 0.25]);
plot(G2.e,zeros(G.N,1),'xk','LineWidth',2,'MarkerSize',6,'DisplayName','Eigenvalues');
title(['Graph=',graph,', Filter=',filter,', K=',num2str(K),', CDF=',param.cdf_method,', StartPtsInterp=',start_pts_interp,', StartPtsLS=',start_pts_ls]);
%xlim([0,40]);

% Figure 3: warped interpolation (w/ damping)
figure;
cc=lines(5);
p1=plot(xx,[yy_cheb_warped,yy_cheb_warped_pchip,yy_cheb_warped_pchip_interp,yy_cheb_warped_ls,yy],'LineWidth',4);
set(gca,'FontSize',20);
legend(p1,'Warped Chebyshev','Warped PCHIP','Warped PCHIP+Chebyshev','Warped LS','g','Location','NorthEast');
%legend(p1,'Warped Interpolation w/ Damping','Warped Interpolation w/out Damping','g','Location','NorthEast');
%set(p1, {'color'}, {cc(3,:);cc(5,:)});
grid on;
hold on;
xlabel('\lambda');
title(['Graph=',graph,', Filter=',filter,', K=',num2str(K),', CDF=',param.cdf_method,', StartPtsInterp=',start_pts_interp,', StartPtsLS=',start_pts_ls]);
plot(G2.e,zeros(G.N,1),'xk','LineWidth',2,'MarkerSize',6,'DisplayName','Eigenvalues');
%plot(chebpts_warped,f(chebpts_warped),'xr','LineWidth',4,'MarkerSize',15,'DisplayName','Warped Pts');
%xlim([0,70]);
%set(gca,'XTick',0:10:70);

% Figure 4: absolute error comparison across methods
figure;
p2=semilogy(xx,[abs(yy-yy_cheb),abs(yy-yy_leg),abs(yy-yy_cheb_warped),abs(yy-yy_cheb_warped_pchip),abs(yy-yy_cheb_warped_pchip_interp),abs(yy-yy_cheb_warped_ls),abs(yy-yy_spec_adapted_ortho),abs(yy-yy_ls)],'LineWidth',4);
set(gca,'FontSize',20);
set(p2,{'Color'},morecolors(1:8));
hold on;
plot(G2.e,ones(G.N,1),'xk','LineWidth',2,'MarkerSize',6);
legend(p2,'Chebyshev','Legendre','Warped Chebyshev','Warped PCHIP','Warped PCHIP+Chebyshev','Warped LS','Spectrum-Adapted Ortho. Poly.','Discrete LS','Location','SouthEast');  
plot(G2.e, abs(y-y_lanczos), 'LineWidth',4,'DisplayName','Regular Lanczos','Color',[1 0.6 0.8]);
xlabel('\lambda');
ylabel('$$|\tilde{g}(\lambda)-g(\lambda)|$$','Interpreter','latex');
title(['Graph=',graph,', Filter=',filter,', K=',num2str(K),', CDF=',param.cdf_method,', StartPtsInterp=',start_pts_interp,', StartPtsLS=',start_pts_ls]);
grid on;

% Figure 5: absolute error comparison focused on eigenvalues
figure;
plot3=semilogy(G2.e,[abs(errors_cheb),abs(errors_leg),abs(errors_cheb_warped),abs(errors_cheb_warped_pchip),abs(errors_cheb_warped_pchip_interp),abs(errors_cheb_warped_ls),abs(errors_spec_adapted_ortho),abs(errors_ls),abs(y-y_lanczos)],'-o','LineWidth',1,'MarkerSize',10);
set(gca,'FontSize',20);
set(plot3,{'Color'},morecolors);
set(plot3, {'MarkerFaceColor'}, get(plot3,'Color')); 
hold on;
plot(G2.e,ones(G.N,1),'xk','LineWidth',2,'MarkerSize',6);
legend(plot3,'Chebyshev','Legendre','Warped Chebyshev','Warped PCHIP','Warped PCHIP+Chebyshev','Warped LS','Spectrum-Adapted Ortho. Poly.','Discrete LS','Lanczos','Location','SouthEast');  
xlabel('\lambda');
ylabel('$$|\tilde{g}(\lambda)-g(\lambda)|$$','Interpreter','latex');
title(['Graph=',graph,', Filter=',filter,', K=',num2str(K),', CDF=',param.cdf_method,', StartPtsInterp=',start_pts_interp,', StartPtsLS=',start_pts_ls]);
grid on;

% Test f(L)b on random signal b
num_tests=5;

ntime_cheb=zeros(num_tests,1);
ntime_leg=zeros(num_tests,1);
ntime_lanczos=zeros(num_tests,1);
ntime_warped=zeros(num_tests,1);
ntime_warped_pchip_interp=zeros(num_tests,1);
ntime_warped_ls=zeros(num_tests,1);
ntime_spec=zeros(num_tests,1);
ntime_ls=zeros(num_tests,1);
ntime_exact=zeros(num_tests,1);

nmse_cheb=zeros(num_tests,1);
nmse_leg=zeros(num_tests,1);
nmse_lanczos=zeros(num_tests,1);
nmse_warped=zeros(num_tests,1);
nmse_warped_pchip_interp=zeros(num_tests,1);
nmse_warped_ls=zeros(num_tests,1);
nmse_spec=zeros(num_tests,1);
nmse_ls=zeros(num_tests,1);

lanc_param.method='lanczos';
lanc_param.order=K;

for i=1:num_tests
    b=rand(G.N,1);
    bhat=G2.U'*b;
    tic
    gLf_exact=G2.U*(f(G2.e).*bhat);
    ntime_exact(i)=toc;
    % add time to compute exact spectrum for every new vector b
    % can change to add this time only once for all b's
    ntime_exact(i)=ntime_exact(i)+time_exact_spec;
    tic
    gLf_ls=G2.U*(polyval(lsc,G2.e).*bhat);
    ntime_ls(i)=toc;
    % add time to compute exact spectrum for every new vector b
    % can change to add this time only once for all b's
    ntime_ls(i)=ntime_ls(i)+time_exact_spec;
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
    tic
    gLf_warped_pchip_interp=gsp_cheby_op(G,p_warped_pchip_interp_c,b);
    ntime_warped_pchip_interp(i)=toc;
    tic
    gLf_warped_ls=poly_op(G,p_warped_lsc,b);
    ntime_warped_ls(i)=toc;
    tic
    gLf_spec=three_term_recurr_op(G2,ab,c_spec_adapted_ortho,b);
    ntime_spec(i)=toc;

    nmse_cheb(i)=sum((gLf_cheb-gLf_exact).^2)/sum(gLf_exact.^2);
    nmse_leg(i)=sum((gLf_leg-gLf_exact).^2)/sum(gLf_exact.^2);
    nmse_lanczos(i)=sum((gLf_lanczos-gLf_exact).^2)/sum(gLf_exact.^2);
    nmse_warped(i)=sum((gLf_warped-gLf_exact).^2)/sum(gLf_exact.^2);
    nmse_warped_pchip_interp(i)=sum((gLf_warped_pchip_interp-gLf_exact).^2)/sum(gLf_exact.^2);
    nmse_warped_ls(i)=sum((gLf_warped_ls-gLf_exact).^2)/sum(gLf_exact.^2);
    nmse_spec(i)=sum((gLf_spec-gLf_exact).^2)/sum(gLf_exact.^2);
    nmse_ls(i)=sum((gLf_ls-gLf_exact).^2)/sum(gLf_exact.^2);
    
end

avg_ntime_cheb=mean(ntime_cheb);
avg_ntime_leg=mean(ntime_leg);
avg_ntime_lanczos=mean(ntime_lanczos);
avg_ntime_warped=mean(ntime_warped);
avg_ntime_warped_pchip_interp=mean(ntime_warped_pchip_interp);
avg_ntime_warped_ls=mean(ntime_warped_ls);
avg_ntime_spec=mean(ntime_spec);
avg_ntime_ls=mean(ntime_ls);
avg_ntime_exact=mean(ntime_exact);

avg_nmse_cheb=mean(nmse_cheb);
avg_nmse_leg=mean(nmse_leg);
avg_nmse_lanczos=mean(nmse_lanczos);
avg_nmse_warped=mean(nmse_warped);
avg_nmse_warped_pchip_interp=mean(nmse_warped_pchip_interp);
avg_nmse_warped_ls=mean(nmse_warped_ls);
avg_nmse_spec=mean(nmse_spec);
avg_nmse_ls=mean(nmse_ls);


Avg_Time = [avg_ntime_cheb; avg_ntime_leg; avg_ntime_lanczos; avg_ntime_warped; avg_ntime_warped_pchip_interp; avg_ntime_warped_ls; avg_ntime_spec; avg_ntime_ls; avg_ntime_exact];
Avg_NMSE = [avg_nmse_cheb; avg_nmse_leg; avg_nmse_lanczos; avg_nmse_warped; avg_nmse_warped_pchip_interp; avg_nmse_warped_ls; avg_nmse_spec; avg_nmse_ls;0];
summary = table(Method,Sup_Err,Sq_Err,Avg_Time,Avg_NMSE)
parameters = cell2table({graph,param.cdf_method,filter,K,start_pts_interp,start_pts_ls,damping},'VariableNames',{'Graph' 'CDFMethod' 'Filter' 'K' 'StartPtsInterp','StartPtsLS','Damping'})
