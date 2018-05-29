function [m] = gsp_hutch(G, r)
% Computes \frac{1}{num_vec} \sum_{i=1}^{num_vec} x_i^t \tilde{h}(L) x_i to
% approximate tr(h(L))
% r is \tilde{h}(L)X, computed from gsp_cheby_opX
	num_vec = size(G.X,2);
    zz =G.X.*r;
	m = sum(sum(zz))/num_vec;
end

