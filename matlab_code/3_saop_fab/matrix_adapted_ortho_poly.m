function [recurr_coeffs,Pi]=matrix_adapted_ortho_poly(absc,weights,poly_order)
%MATRIX_ADAPTED_ORTHO_POLY: Computes the recurrence coefficients from the 
%discrete measure, using the Lanczos method described in Gautschi, 
%Section 2.2.3.2 (pp. 96-98)
%   Usage: [ab,Pi]=matrix_adapted_ortho_poly(absc,weights,poly_order);
%
%   Input parameters:
%       absc       : Abscissae of a discrete set of N points {x_i}, the support
%                    points of a discrete measure with respect to which the 
%                    polynomials {\pi_k} are orthogonal
%       weights    : Weights of the N support points
%       poly_order : Maximum order Chebyshev coefficient to compute

%   Output parameters
%       recurr_coeffs : Matrix of expansion coefficients
%       Pi            : Values of the orthogonal polynomials at the abscissae,
%                       (i-1)th polynomial in ith column of Pi

xw=[absc,weights/sum(weights)];
%tic
recurr_coeffs=opq_lanczos(poly_order+2,xw); % from https://www.cs.purdue.edu/archives/2002/wxg/codes/OPQ.html
%lanczos_time=toc

% Evaluate Pi at a discrete set of pts
Pi=eval_pi(recurr_coeffs,absc);

end