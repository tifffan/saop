function [r] = gsp_cheby_opX(G, coeffs)
%GSP_CHEBY_OPX: Computes 1/2 \alpha_0 \bar{T}_0(L)X + \sum_{k=1}^K \alpha_k
% \bar{T}_k(L)X, an approximation of h(L)X for some filter h
%
% Usage: r = gsp_cheby_opX(G, coeffs);
%
% Input parameters:
%     G         : Graph structure
%     coeffs    : Chebyshev coefficients {\alpha_k}
%
% Output parameters:
%     r         : Matrix of 1/2 \alpha_0 \bar{T}_0(L)X + \sum_{k=1}^K 
%                 \alpha_k \bar{T}_k(L)X, where each column is the result
%                 for one column vector in X
	
if ~isfield(G,'TkbarLX') || ~isfield(G,'X')
   G=gsp_compute_TkbarLX(G);
end

num_vec = size(G.X, 2);
spI = speye(num_vec);
coeffs(1) = 1/2*coeffs(1);
coeffs_matrix = kron(coeffs, spI);

if size(G.TkbarLX,2) < size(coeffs_matrix,1)
   param.num_vec=size(coeffs_matrix,1);
   G=gsp_compute_TkbarLX(G,param);
end

r = G.TkbarLX(:,1:size(coeffs_matrix,1))*coeffs_matrix;