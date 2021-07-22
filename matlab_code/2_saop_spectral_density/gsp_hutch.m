function [m] = gsp_hutch(G, r)
% GSP_HUTCH: Computes \frac{1}{num_vec} \sum_{i=1}^{num_vec} x_i^t \tilde{h}(L) x_i to
% approximate tr(h(L))
%
% Usage: m = gsp_hutch(G, r);
%
% Input parameters:
%     G          : Graph structure
%     r          : Matrix of \tilde{h}(L)X, computed from gsp_cheby_opX
%
% Output parameters:
%     m          : Approximation to tr(h(L))

	num_vec = size(G.X,2);
    zz =G.X.*r;
	m = sum(sum(zz))/num_vec;
end

