function [recurr_coeffs,Pi]=matrix_adapted_ortho_poly(absc,weights,poly_order)

% Compute the recurrence coefficients from the discrete measure, using the
% Lanczos method described in Gautschi, Section 2.2.3.2 (pp. 96-98)
xw=[absc,weights/sum(weights)]; %TODO: place weights at absc or in between?
%tic
recurr_coeffs=opq_lanczos(poly_order+2,xw); % from https://www.cs.purdue.edu/archives/2002/wxg/codes/OPQ.html
%lanczos_time=toc

% Evaluate Pi at a discrete set of pts
Pi=eval_pi(recurr_coeffs,absc);

end