function c = matrix_adapted_poly_coeff(G, filter, absc, weights, Pi, m)
%MATRIX_ADAPTED_POLY_COEFF : Computes expansion coefficients for a filterbank
%
%   Usage: c = matrix_adapted_poly_coeff(G, filter, absc, weights,Pi,m);
%
%   Input parameters:
%       G       : Graph structure or range of application
%       filter  : Filter or cell array of filters
%       m       : Maximum order Chebyshev coefficient to compute
%       Pi      : Values of the orthogonal polynomials at the abscissae
%                 (i-1)th polynomial in ith column of Pi
%   Output parameters
%       c       : Matrix of expansion coefficients
% 
%   This function compute the expansion coefficients for all of the filters
%   contained in the cell array filter. The coefficients are returned in a
%   matrix. Every column corresponds to a filter. The coefficients are
%   ordered such that c(j+1) is j'th expansion coefficient
%
%   This function is inspired by the GSPBox and sgwt_toolbox
%
%   See also: matrix_adapted_ortho_poly three_term_recurr_op gsp_cheby_op 
%

% Author: David Shuman, David K Hammond, Nathanael Perraudin
% Testing: 
% Date: 1 March 2018


if iscell(filter)
   Nf = length(filter);
   c = zeros(m+1,Nf);
   for ii = 1: Nf
       c(:,ii) = matrix_adapted_poly_coeff(G, filter{ii}, absc, weights, Pi, m);
   end
   return;
end

c = zeros(m+1,1);
weights=weights/sum(weights);

for ii=1:m+1
    c(ii) = sum(filter(absc).*weights.*Pi(:,ii))/sum(Pi(:,ii).^2.*weights);  
end

end
  
