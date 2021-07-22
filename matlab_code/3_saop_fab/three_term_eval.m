function r = three_term_eval(G,x,ab,c)
%GSP_CHEBY_EVAL Evaluate chebyshev polynomial
%
%   Usage:  r = gsp_cheby_eval(x,c,arange)
%
%   Input parameters:
%       x       : Points to evaluate the polynomial
%       c       : Chebyshef coefficients
%       arrange : arange (range to evaluate the polynomial)
%   Output parameters
%       r       : Result
%
%   In this function, *arrange* has to be [0, lmax]!

% Author: David K Hammond, Nathanael Perraudin
% Testing: test_dual
% Date: 30 December 2014

% By setting the operator L to mean (pointwise) multiplication by x,
% and f to be vector of ones, p(L)f will contain p(x) at each
% point. This lets us use gsp_cheby_eval to evaluate the Chebyshev polynomials.

% if arange(1)
%     error('This function will not work!')
% end

[N1,N2] = size(x);

L=spdiags(x(:),0,speye(numel(x)));
f=ones(size(x(:)));
N = length(f);

%G.lmax = arange(2);
 G.L = L;
 G.N = N;
param=struct;
r = three_term_recurr_op(G, ab, c, f,param);

r = reshape(r, N1, N2);

end

