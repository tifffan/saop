
close all;
clear all;

N=64;
T=4;
G=gsp_david_sensor_network(N);
G=gsp_compute_fourier_basis(G);
figure;
gsp_plot_graph(G);
G.lmin=0;

param.weighted=0;
param.directed=0;
param.periodic=1;
Gm=create_multilayer_graph(G,T,param);
figure;
gsp_plot_graph(Gm);

Gm=gsp_compute_fourier_basis(Gm);
% figure;
% gsp_plot_signal(Gm,Gm.U(:,1));
% figure;
% gsp_plot_signal(Gm,Gm.U(:,2));
% figure;
% gsp_plot_signal(Gm,Gm.U(:,30));

% create an unweighted joint time_vertex graph
Gjtv=G;
Gjtv.W=Gjtv.W>0;
Gjtv.d=sum(Gjtv.W);
Gjtv.d=Gjtv.d(:);
Gjtv.L=diag(Gjtv.d)-Gjtv.W;
Gjtv=gsp_jtv_graph(Gjtv,T);

% plot delta
delta=zeros(N,T);
%delta(1,1)=1;
delta(30,2)=1;
delta_col=delta(:);
gsp_plot_signal(Gm,delta_col);



% param.method = 'cheby';
% param.order=10;
% delta1=gsp_jtv_filter_analysis(Gjtv, hlp, 'js', delta, param);
% delta1_col=delta1(:);
% gsp_plot_signal(Gm,delta1_col);



% joint laplacian
Gjtv.LJ=kron(Gjtv.jtv.LT, eye(N))+kron(eye(T), Gjtv.L); % equivalent to cartesian product
% product of LT and LG
max(max(abs(Gjtv.LJ-Gm.L))) % should be 0; happens for undirected time graph with periodic boundary conditions
% now the two laplacian matrices match
% but the eigendecomposition rank eigenvalues from zero to largest, which
% is not be the order we want

2-2*cos(2*(0:T-1)*pi/T) % eigenvalues of ring graph; checked

[VT,DT]=eigs(Gjtv.jtv.LT);

[dT,ind] = sort(diag(DT)); % sorted time eigenvalues
DTs = DT(ind,ind);
VTs = VT(:,ind); % sorted time eigenvectors

% filter
lcut = 1;
wcut = 1;
hlp = @( lambda,omega ) double(and(abs(lambda)<1,abs(omega)<wcut));

lcut = 0.8;
wcut2 = 0.5;
hlp2 = @( lambda,omega ) double(and(abs(lambda)<1,abs(omega)<wcut2));
% signal in matrix form
delta=zeros(N,T);
delta(30,3)=1;
signal=delta;

% evaluate filter coefficients
[X,Y] = meshgrid(dT,G.e); 
Hlp = hlp( Y,X ); % N*T matrix of filter coefficients, h(lambda, omega)

% perform the transform and scaling in matrix form
signal_hat=conj(G.U')*signal*VTs;
scaled_signal_hat=signal_hat.*Hlp;
filtered_signal=G.U*scaled_signal_hat*conj(VTs');

figure;
% convert filtered singal to vector form for plotting
filtered_signal_col=filtered_signal(:);
gsp_plot_signal(Gm,filtered_signal_col);

coef = gsp_jtv_filter_analysis(Gjtv, hlp2, 'js', delta)
coef_col=coef(:);
figure;gsp_plot_signal(Gm,coef_col);

figure;gsp_plot_signal(Gm,delta(:));

% evaluate filter coefficients in vector form
xy =[ Y(:) X(:) ]; % 2 columns; each row is one graph eigenvalue and one time eigenvalue; loop through graph eigenvalues first
Hlp_col = hlp(xy(:,1), xy(:,2) );

