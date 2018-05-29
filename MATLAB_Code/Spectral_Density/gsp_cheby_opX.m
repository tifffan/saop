function [r] = gsp_cheby_opX(G, coeffs)

% Computes 1/2 \alpha_0 \bar{T}_0(L)X + \sum_{k=1}^K \alpha_k
% \bar{T}_k(L)X, an approximation of h(L)X for some filter h
	
if ~isfield(G,'TkbarLX') || ~isfield(G,'X')
   G=gsp_compute_TkbarLX(G);
end

num_vec = size(G.X, 2);
spI = speye(num_vec);
coeffs(1) = 1/2*coeffs(1);
coeffs_matrix = kron(coeffs, spI);

r = G.TkbarLX*coeffs_matrix;

end

