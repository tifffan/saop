%% Denoising of dynamic mesh using joint time-vertex Tikhonov: Comparison between Fast Fourier Chebyshev and Fast Fourier Spectrum-Adapted Weighted Least Squares  
%   In this demo, we perform denoising of a dynamic mesh representing a dog
%   walking using a joint Tikhonov approach. The regularization parameters
%   are chosen to be the best for the exact filter
%   Dataset can be found at http://research.microsoft.com/en-us/um/redmond/events/geometrycompression/data/default.html
%

% Author: Original code by Francesco Grassi / Modified by David Shuman
% Date: April 2017 / July 2020

clear
close all
gsp_start
init_unlocbox

load('0_data/dog.mat');

rand('seed',0);
randn('seed',0);

%% Signal
X1 = X(:,:,1) - repmat(mean(X(:,:,1),1),[N,1]);
X2 = X(:,:,2) - repmat(mean(X(:,:,2),1),[N,1]);
X3 = X(:,:,3) - repmat(mean(X(:,:,3),1),[N,1]);

X = cat(3,X1,X2,X3);

%% Parameters
param.k = 5; % number of nearest neighbors for graph construction
param.transform = 'dft';
tau1 = 7.2; 
tau2 = .45;
itermax = 1;
err = @(x,y) norm(vec(x)-vec(y),'fro')/norm(vec(x),'fro');
orders=3:25;
ns_tik = @(l,w) 1./(1+tau1*l+2*tau2*(1-cos(2*pi*w)));
ns_tik_ring = @(l,r) 1./(1+tau1*l+tau2*r);
tplot=5;
iterplot=1;
kplot=6;
num_grid_pts=50;

%% 
T=size(X,2);
GT=gsp_ring(T);
GT=gsp_compute_fourier_basis(GT);
errjoint_cheb=zeros(itermax,length(orders));
errjoint_cheby2D=zeros(itermax,length(orders));
errjoint_wls=zeros(itermax,length(orders));
errjoint_exact=zeros(itermax,1);
for iter=1:itermax
    iter
    % Add noise to vertices position
    noise = randn(size(X));
    noise = 0.2 * noise * norm(X(:)) / norm(noise(:));
    Xn = X + noise;

    % Graph
    x0 = squeeze(mean(Xn,2));
    G = gsp_nn_graph(x0,param);
    G2 = G;
    G2 = gsp_compute_fourier_basis(G2);
    G2 = gsp_jtv_graph(G2,size(X,2),[],param);
    G = gsp_estimate_lmax(G);
    G = gsp_jtv_graph(G,size(X,2),[],param);
    v = gsp_jtv_fa(G,0);
    G.lmin=0;
    
    % Exact (G2)
    Y_joint_exact=gsp_jtv_filter_analysis(G2,ns_tik,'js',Xn);
    errjoint_exact(iter)=err(X,Y_joint_exact);
    
    % Fast Fourier Chebyshev (G)
    for k=1:length(orders)
        filter_param.order=orders(k);
        xfilt=gsp_jtv_filter_analysis(G,ns_tik,'js',Xn(:,:,1),filter_param);
        yfilt=gsp_jtv_filter_analysis(G,ns_tik,'js',Xn(:,:,2),filter_param);
        zfilt=gsp_jtv_filter_analysis(G,ns_tik,'js',Xn(:,:,3),filter_param);
        Y_joint_cheb=cat(3,xfilt,yfilt,zfilt);
        if iter==iterplot && orders(k)==kplot
            X_den_cheb=Y_joint_cheb;
        end
        errjoint_cheb(iter,k) = err(X,Y_joint_cheb);
    end
    
%     % Cheby2D
      [XX,YY]=meshgrid(G2.jtv.omega,G2.e); 
      xxyy=[YY(:) XX(:)];
      ns_tik2 = @(xy) 1./(1+tau1*xy(:,1)+2*tau2*(1-cos(2*pi*xy(:,2))));
      lower = [ min( G2.e ), min( G2.jtv.omega ) ];
      upper = [ max( G2.e ), max( G2.jtv.omega ) ];

     for k=1:length(orders)
        funobj = chebint( 2, ns_tik2, 'hyperrectangle', [ lower; upper ], 'degree', orders(k) );
        H_cheb2D = reshape( feval( funobj, xxyy ), N, T );
        xfilt_2D=gsp_ijft(G2,gsp_jft(G2,Xn(:,:,1)).*H_cheb2D);
        yfilt_2D=gsp_ijft(G2,gsp_jft(G2,Xn(:,:,2)).*H_cheb2D);
        zfilt_2D=gsp_ijft(G2,gsp_jft(G2,Xn(:,:,3)).*H_cheb2D);
        Y_joint_cheby2D=real(cat(3,xfilt_2D,yfilt_2D,zfilt_2D));
        if iter==iterplot && orders(k)==kplot
            X_den_cheby2D=Y_joint_cheby2D;
        end
        errjoint_cheby2D(iter,k) = err(X,Y_joint_cheby2D);
    end
    
    % Approximate CDF
    cdf_param.cdf_method='kpm';
    cdf_param.num_pts=10; 
    cdf_param.num_vec=10; 
    cdf_param.order=30; 
    G=spectral_cdf_approx2(G,cdf_param); 

    if iter==iterplot
        Gplot=G;
        Gplot2=G2;
    end
    
    % Apply DFT
    XnTilde1=Xn(:,:,1)*conj(GT.U); %gsp_tft(G,Xn(:,:,1));
    XnTilde2=Xn(:,:,2)*conj(GT.U);
    XnTilde3=Xn(:,:,3)*conj(GT.U);
    
    % Weighted LS (G)
    for k=1:length(orders)
        filter_param.order=orders(k);
        mop_param.grid_order=num_grid_pts;
        mop_param.absc_method='linear';
        mop_param.weights_method='pdf';
        [absc_weighted_ls,weights_weighted_ls]=gen_absc_weights(G,filter_param.order,mop_param);
        xfilt_wls=zeros(G.N,T);
        yfilt_wls=zeros(G.N,T);
        zfilt_wls=zeros(G.N,T);
        for r=1:T
            fi=@(x) ns_tik_ring(x,GT.e(r));
            [p_weighted_lsc,s_weighted_ls,mu_weighted_ls]=weighted_polyfit(absc_weighted_ls,fi(absc_weighted_ls),weights_weighted_ls,filter_param.order);
            pp=polyval(p_weighted_lsc,G2.e,s_weighted_ls,mu_weighted_ls);
            xfilt_wls(:,r)=G2.U*(pp.*(G2.U'*XnTilde1(:,r)));
            yfilt_wls(:,r)=G2.U*(pp.*(G2.U'*XnTilde2(:,r)));
            zfilt_wls(:,r)=G2.U*(pp.*(G2.U'*XnTilde3(:,r)));
        end
        xfilt_wls=xfilt_wls*transpose(GT.U); %gsp_itft(G,c1);
        yfilt_wls=yfilt_wls*transpose(GT.U); %gsp_itft(G,c2);
        zfilt_wls=zfilt_wls*transpose(GT.U); %gsp_itft(G,c3);
        Y_joint_wls=real(cat(3,xfilt_wls,yfilt_wls,zfilt_wls));
        if iter==iterplot && orders(k)==kplot
            X_den_wls=Y_joint_wls;
        end
        errjoint_wls(iter,k) = err(X,Y_joint_wls);
    end
end
avg_error_cheb=mean(errjoint_cheb,1); % vector the size of orders
avg_error_cheby2D=mean(errjoint_cheby2D,1); % vector the size of orders
avg_error_wls=mean(errjoint_wls,1); % vector the size of orders
avg_error_exact=mean(errjoint_exact); % single number

figure;
plot0=plot(orders,[avg_error_cheb',avg_error_cheby2D',avg_error_wls'],'o-','LineWidth',2,'MarkerSize',6);
hold on;
plot(orders,avg_error_exact*ones(size(orders)),'k:','LineWidth',2);
set(gca,'FontSize',24);
xlim([min(orders),max(orders)]);
ylim([0,.65]);
legend('Fast Fourier Chebyshev','Cheby2D','Fast Fourier Weighted LS','Exact','Location','Northeast');  
xlabel('K');
ylabel({'Average Relative Error:';'$$||{\bf X}_{\hbox{denoised}}-{\bf X}||_F~/~||{\bf X}||_F$$'},'Interpreter','latex');
grid on;
newcolors = [0, 0.4470, 0.7410
            0.8500,    0.3250,    0.0980
             0.4940, 0.1840, 0.5560];
colororder(newcolors);         

param.view=[0,90];
Gplot.plotting.vertex_size=25;
figure;
gsp_plot_graph(Gplot);
view(param.view)

% plot cdf
xx=linspace(Gplot.lmin,Gplot.lmax,10000);
mu=@(t) sum(Gplot2.e<=t)/Gplot.N;
figure;
plot1=plot(xx,[Gplot.spectrum_cdf_approx(xx)',mu(xx)'],'linewidth',4);hold on;
plot(xx,xx/Gplot.lmax,'k:','LineWidth',1.5);
grid on;
legend(plot1,'Estimated Spectral CDF','Actual Spectral CDF','Location','Southeast');  
xlabel('$$\lambda_G$$','Interpreter','latex');
xlim([Gplot.lmin,Gplot2.lmax]);
set(gca,'FontSize',24);

psize=72;

figure;
scatter3(X(:,tplot,1),X(:,tplot,2),X(:,tplot,3),psize,'k.')
axis([min(vec(X(:,tplot,1))) max(vec(X(:,tplot,1))) min(vec(X(:,tplot,2))) max(vec(X(:,tplot,2))) min(vec(X(:,tplot,3))) max(vec(X(:,tplot,3))) ])
view(param.view)
axis off;
grid off;
box off;

figure;
scatter3(Xn(:,tplot,1),Xn(:,tplot,2),Xn(:,tplot,3),psize,'k.')
axis([min(vec(Xn(:,tplot,1))) max(vec(Xn(:,tplot,1))) min(vec(Xn(:,tplot,2))) max(vec(Xn(:,tplot,2))) min(vec(Xn(:,tplot,3))) max(vec(Xn(:,tplot,3))) ])
view(param.view)
axis off;
grid off;
box off;

figure;
scatter3(X_den_cheb(:,tplot,1),X_den_cheb(:,tplot,2),X_den_cheb(:,tplot,3),psize,'k.')
axis([min(vec(X_den_cheb(:,tplot,1))) max(vec(X_den_cheb(:,tplot,1))) min(vec(X_den_cheb(:,tplot,2))) max(vec(X_den_cheb(:,tplot,2))) min(vec(X_den_cheb(:,tplot,3))) max(vec(X_den_cheb(:,tplot,3))) ])
view(param.view)
axis off;
grid off;
box off;

figure;
scatter3(X_den_wls(:,tplot,1),X_den_wls(:,tplot,2),X_den_wls(:,tplot,3),psize,'k.')
axis([min(vec(X_den_wls(:,tplot,1))) max(vec(X_den_wls(:,tplot,1))) min(vec(X_den_wls(:,tplot,2))) max(vec(X_den_wls(:,tplot,2))) min(vec(X_den_wls(:,tplot,3))) max(vec(X_den_wls(:,tplot,3))) ])
view(param.view)
axis off;
grid off;
box off;

[XX,YY]=meshgrid(GT.e,G2.e);
ZZ=ns_tik_ring(YY,XX);
xx=min(GT.e):.05:max(GT.e);
yy=0:.05:G.lmax;

figure; 
surf(xx,yy,ns_tik_ring(yy',xx)); 
colormap autumn; 
hold on;
scatter3(vec(XX),vec(YY),vec(ZZ),10,'filled','MarkerFaceColor','blue');
xlabel('$$\lambda_R$$','Interpreter','latex');
ylabel('$$\lambda_G$$','Interpreter','latex');
zlabel('$$h_{\hbox{tik}}\left(\lambda_G,\lambda_R\right)$$','Interpreter','latex');
set(gca,'FontSize',24);
view([-25,32]);