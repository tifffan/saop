clear all;
close all;
rand('seed',0);
randn('seed',0);

%% Parameters 1
graph_type='bunny';
num_filters=4; 
filter_type='wav_itersine';

%% Graph
switch graph_type
    case 'bunny'
        G=gsp_bunny();
        hadamard_size=2560;
    case 'brain'
        G0 = spg_construct_cerebellum_graph;
        G = gsp_graph(G0.W,G0.coords);
    case 'minnesota'
        G=gsp_minnesota();
        hadamard_size=3072;
    case 'gnp'
        N=500;
        p=.2;
        G=gsp_erdos_renyi(N,p);
    case 'community'
        N=2500;
        G=gsp_community(N);
    case 'net25'
        load('0_data/net25.mat');
        A=Problem.A;
        A=A-diag(diag(A)); 
        A(4228,6327)=1;
        A(6327,4228)=1;
        G=gsp_graph(A);
        hadamard_size=10240;
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
        A=A/2000;
        G=gsp_graph(A);
    otherwise
        error('unknown graph type');
end

G2=G;
G2=gsp_compute_fourier_basis(G2);
G=gsp_estimate_lmax(G);

%% Spectral filters
filts=spm_parseval_filters(G,num_filters,filter_type);
figure;
plot_param.show_sum=0;
gsp_plot_filter(G,filts,plot_param);
set(gca,'FontSize',24);
xlabel('$$\lambda$$','Interpreter','latex');

%% Exact atom norms
Dict=zeros(G2.N,G2.N*num_filters);
for i=1:num_filters
    Dict(:,(i-1)*G2.N+1:i*G2.N)=G2.U*diag(filts{i}(G2.e))*G2.U';
end
atom_norms=sqrt(sum(Dict.^2))';
[sorted_atom_norms,idx]=sort(atom_norms,'descend');
atom_colors=kron(1:num_filters,ones(1,G2.N));
sorted_atom_colors=atom_colors(idx);

figure;
scatter(1:G2.N*num_filters,atom_norms,5,atom_colors);
box on;
xlabel('Atom Index');
ylabel('Atom Norm');
set(gca,'FontSize',24);
title('Exact');

figure;
scatter(1:G2.N*num_filters,sorted_atom_norms,5,sorted_atom_colors);
box on;
xlabel('Sorted Atom Index');
ylabel('Atom Norm');
set(gca,'FontSize',24);
title('Exact');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Parameters 2 
grid_pts=100;
cdf_param.num_pts=10;
cdf_param.cdf_method='kpm';
KK=3:2:25; 
JJ=[15,50];
numK=length(KK);
numJ=length(JJ);
plot_figures=length(KK)<3 && length(JJ)<3; 
mop_param.grid_order=100;
mop_param.absc_method='linear';
mop_param.weights_method='pdf';

randvec_method='normal'; %'normal' or 'hadamard'
H=hadamard(hadamard_size);

mean_rel_error=zeros(numK,numJ);
mean_rel_error_saop=zeros(numK,numJ);

for k=1:numK
    for j=1:numJ
        K=KK(k);
        cdf_param.order=K;
        cdf_param.num_vec=JJ(j);
        G=spectral_cdf_approx2(G,cdf_param);
        X=H(1:G.N,1:cdf_param.num_vec);

        %% Approximate method 1: stochastic estimation
        % (no additional complexity beyond density estimation, which has KJ matvecs)

        cc=gsp_cheby_coeff(G,filts,K,grid_pts);
        filtered_noise=zeros(G2.N*num_filters,cdf_param.num_vec);
        for i=1:num_filters
            switch randvec_method
                case 'normal'
                    filtered_noise((i-1)*G2.N+1:i*G2.N,:)=gsp_cheby_opX(G,cc(:,i));
                case 'hadamard'
                    filtered_noise((i-1)*G2.N+1:i*G2.N,:)=gsp_cheby_op(G,cc(:,i),X);
            end
        end
        var_noise_tx=var(filtered_noise')';
        stochastic_atom_norm_est=sqrt(var_noise_tx);
        stochastic_errors=stochastic_atom_norm_est-atom_norms;
        rel_error=abs(stochastic_atom_norm_est-atom_norms)./atom_norms;
        mean_rel_error(k,j)=mean(rel_error);

        if plot_figures
            ratio=stochastic_atom_norm_est./atom_norms;
            %mean_abs_ratio_error=mean(abs(ratio-1))

            figure;
            hold on;
            scatter(1:G2.N*num_filters,ratio,5,atom_colors);
            box on;
            xlabel('Atom Index');
            ylabel('Atom Norm Estimate / Atom Norm');
            set(gca,'FontSize',24);
            ylim([0,2]);
            title('Chebyshev');

            figure;
            hold on;
            scatter(1:G2.N*num_filters,rel_error,5,atom_colors);
            box on;
            xlabel('Atom Index');
            ylabel('Relative Error of Atom Norm Estimate');
            set(gca,'FontSize',24);
            title('Chebyshev');

            figure;
            hold on;
            scatter(1:G2.N*num_filters,stochastic_atom_norm_est(idx),5,sorted_atom_colors);
            scatter(1:G2.N*num_filters,sorted_atom_norms,5,'k');
            box on;
            xlabel('Sorted Atom Index');
            ylabel('Atom Norm Estimate');
            set(gca,'FontSize',24);
            title('Chebyshev');
        end

        %% Approximate method 1a: stochastic estimation with spectrum-adapted polynomials
        % (no additional complexity beyond density estimation, which has KJ matvecs)
        [absc_weighted_ls,weights_weighted_ls]=gen_absc_weights(G,K,mop_param);
        filtered_noise_saop=zeros(G2.N*num_filters,cdf_param.num_vec);
        GXhat=G2.U'*G.X;
        Xhat=G2.U'*X; % hadamard
        for i=1:num_filters
            [p_weighted_lsc,s_weighted_ls,mu_weighted_ls]=weighted_polyfit(absc_weighted_ls,filts{i}(absc_weighted_ls),weights_weighted_ls,K);
            switch randvec_method
                case 'normal'
                    filtered_noise_saop((i-1)*G2.N+1:i*G2.N,:)=G2.U*(polyval(p_weighted_lsc,G2.e,s_weighted_ls,mu_weighted_ls).*GXhat);
                case 'hadamard'
                    filtered_noise_saop((i-1)*G2.N+1:i*G2.N,:)=G2.U*(polyval(p_weighted_lsc,G2.e,s_weighted_ls,mu_weighted_ls).*Xhat);
            end
        end
        var_noise_tx_saop=var(filtered_noise_saop')';
        stochastic_atom_norm_est_saop=sqrt(var_noise_tx_saop);
        rel_error_saop=abs(stochastic_atom_norm_est_saop-atom_norms)./atom_norms;
        mean_rel_error_saop(k,j)=mean(rel_error_saop);

        if plot_figures
            ratio_saop=stochastic_atom_norm_est_saop./atom_norms;
            mean_abs_ratio_error_saop=mean(abs(ratio_saop-1))
            %mean_square_ratio_error_saop=mean((ratio_saop-1).^2)

            figure;
            hold on;
            scatter(1:G2.N*num_filters,ratio_saop,5,atom_colors);
            box on;
            xlabel('Atom Index');
            ylabel('Atom Norm Estimate / Atom Norm');
            set(gca,'FontSize',20);
            ylim([0,2]);
            title('Weighted Least Squares');

            figure;
            hold on;
            scatter(1:G2.N*num_filters,rel_error_saop,5,atom_colors);
            box on;
            xlabel('Atom Index');
            ylabel('Relative Error of Atom Norm Estimate');
            set(gca,'FontSize',20);
            title('Weighted Least Squares');

            figure;
            hold on;
            scatter(1:G2.N*num_filters,stochastic_atom_norm_est_saop(idx),5,sorted_atom_colors);
            scatter(1:G2.N*num_filters,sorted_atom_norms,5,'k');
            box on;
            xlabel('Sorted Atom Index');
            title('Weighted Least Squares');
            ylabel('Atom Norm Estimate');
            set(gca,'FontSize',20);
        end
    end
end

newcolors = [0, 0.4470, 0.7410
             0.4940, 0.1840, 0.5560];
         
all_mean_rel_errors = reshape([mean_rel_error;mean_rel_error_saop], size(mean_rel_error,1), []); % alternate columns so they come in pairs for differnt Js
if length(KK)>1
    figure;
    hold on;
    eplot1=plot(KK,all_mean_rel_errors(:,1:2),'--o','LineWidth',2,'MarkerSize',10);
    colororder(newcolors)
    cc=get(eplot1,'Color');
    set(eplot1,{'MarkerFaceColor'},cc);
    set(eplot1,{'MarkerEdgeColor'},cc);
    for j=2:numJ
        eplot2=plot(KK,all_mean_rel_errors(:,(2*j-1):(2*j)),'-o','LineWidth',2,'MarkerSize',10);
        set(eplot2,{'Color'},cc);
        set(eplot2,{'MarkerFaceColor'},cc);
        set(eplot2,{'MarkerEdgeColor'},cc);
    end
    xlabel('K');
    ylabel('Mean Relative Error');
    set(gca,'FontSize',20);
    legend('Chebyshev, J=15','Weighted LS, J=15','Chebyshev, J=50','Weighted LS, J=50','Location','Northeast');
    xlim([min(KK),max(KK)]);
    box on;
    grid on;
end
