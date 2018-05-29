close all;
clear all;
randn('seed', 18); rand('seed', 18);

% Graph
graph='net25';

switch graph
    case 'gnp'
        N=100;
        p=.3;
        G=gsp_erdos_renyi(N,p);
    case 'sensor'
        N=500;
        G=gsp_david_sensor_network(N);
    case 'minnesota'
        G=gsp_minnesota(1);
    case 'comet'
        G=gsp_comet(100,20);
    case 'community'
        N=5000;
        G=gsp_community(N);
    case 'net25'
        load('/Users/davidshuman/Dropbox/Current Research Work/MCSFB/Shuni_Thesis/GitHub/mcsfb2018/net25 graph data/net25.mat');
        A=Problem.A;
        A = A - diag(diag(A)); 
        A(4228,6327) = 1;
        A(6327,4228) = 1;
        G=gsp_graph(A);
    case 'temperature'
        MyData=csvread('/Users/davidshuman/Dropbox/Current Research Work/MCSFB/Shuni_Thesis/GitHub/mcsfb2018/Data/MyData2.csv');
        coords=MyData(:,2:3);
        inds=MyData(:,4);

        % create 8 neighbor graph and find subset where we have data
        param.eight_neighbor=1;
        G0=gsp_2dgrid(1385,596,param);
        W=G0.W(inds,inds);

        % remove disconnected nodes
        disconnected=(sum(W,2)==0);
        inds=inds(~disconnected);
        coords=coords(~disconnected,:);
        W1=G0.W(inds,inds);

        % remove small components
        G1=gsp_graph(W1,coords);
        [V,D]=eigs(G1.L+1e-11*speye(G1.N),4,'sm'); % make sure this works for other machines. compute total number in each and choose the large one
        large_component=(abs(V(:,4))>.000000001);
        W=G1.W(large_component,large_component);
        coords=coords(large_component,:);

        % create graph
        G=gsp_graph(W,coords);
    otherwise
        error('graph type not recognized');
end

G=gsp_estimate_lmax(G);
ee=linspace(0,G.lmax,2000);
param.order=80;
param.num_vec=30;
param.num_pts=80;

tic % Lanczos
G2=spectral_cdf_approx2(G,param);
lanc_time=toc
tic
yy_lanc=G2.spectrum_cdf_approx(ee)';
lanc_plot_time=toc

tic % KPM
param.cdf_method='kpm';
G3=spectral_cdf_approx2(G,param);
kpm_time=toc
tic
yy_kpm=G3.spectrum_cdf_approx(ee)';
kpm_plot_time=toc

if G.N <= 10000
    tic % LDLt
    param.cdf_method='ldlt';
    G4=spectral_cdf_approx2(G,param);
    ldlt_time=toc
    tic
    yy_ldl=G4.spectrum_cdf_approx(ee)';
    ldlt_plot_time=toc

    figure;
    hold on;
    plot(ee,[yy_lanc,yy_kpm,yy_ldl],'LineWidth',4);
    plot(ee,ee/G.lmax,'k:','LineWidth',1.5);
    xlim([0,G.lmax]);
    ylim([0,1]);
    xlabel('\lambda');
    grid on;
    legend('Lanczos','KPM','LDLT','Location','NorthWest');
    set(gca,'FontSize',24);
    
else % too big to do ldlt
    figure;
    hold on;
    plot(ee,[yy_lanc,yy_kpm],'LineWidth',4);
    plot(ee,ee/G.lmax,'k:','LineWidth',1.5);
    xlim([0,G.lmax]);
    ylim([0,1]);
    xlabel('\lambda');
    grid on;
    legend('Lanczos','KPM','Location','NorthWest');
    set(gca,'FontSize',24);
end

