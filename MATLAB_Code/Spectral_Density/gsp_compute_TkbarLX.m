function [G] = gsp_compute_TkbarLX(G, param)
% Computes \bar{T}_k(L)X for k=0,1,...,K (order)
% Output: Two things stored in G: 
% 1. X (random N x num_vec matrix)
% 2. TkbarLX: N x (K+1)*num_vec matrix [\bar{T}_0(L)X | \bar{T}_1(L)X | ... |
% \bar{T}_K(L)X ]

if nargin < 2
    param = struct;
end

if ~isfield(param,'verbose')
    param.verbose = 1; 
end

if ~isfield(param,'num_vec')
    num_vec = 30; 
else
    num_vec = param.num_vec;
end

if ~isfield(param,'order')
    order = 30; 
else
    order = param.order;
end

if ~isfield(G,'lmax')
    G = gsp_estimate_lmax(G);
    if param.verbose
    warning(['GSP_CHEBY_OP: The variable lmax is not ',...
        'available. The function will compute it for you. ',...
        'However, if you apply many time this function, you ',...
        'should precompute it using the function: ',...
        'gsp_estimate_lmax']);
    end
end

signal = randn(G.N, num_vec);
G.X = signal;

arange = [0, G.lmax];

a1 = (arange(2) - arange(1))/2;
a2 = (arange(2) + arange(1))/2;


TkbarLX = zeros(G.N, num_vec*(order+1));

Twf_old=signal;                     % j = 0;
TkbarLX(1:G.N, 1:num_vec) = Twf_old;

Twf_cur=(G.L*signal-a2*signal)/a1;  % j = 1;
TkbarLX(1:G.N, (1+num_vec):(2*num_vec)) = Twf_cur;

for k=2:order
    Twf_new = (2/a1) * (G.L*Twf_cur-a2*Twf_cur) - Twf_old;
    TkbarLX(1:G.N, (1+k*num_vec):((k+1)*num_vec)) = Twf_new;
    Twf_old=Twf_cur;
    Twf_cur=Twf_new;
end

G.TkbarLX = TkbarLX;


end

