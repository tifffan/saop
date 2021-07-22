function [p,S,mu] = weighted_polyfit(x,y,w,n)
%WEIGHTED_POLYFIT: Fits a degree n polynomial with least weighted least 
%squares error to data, adapted from the polyfit function built in Matlab
%
%   Usage: [p,S,mu] = weighted_polyfit(x,y,w,n);
%
%   Input parameters:
%       x      : x coordinates of the data
%       y      : y coordinates of the data
%       w      : Weights of the points (x,y)
%       n      : Order of polynomial fitting
%   Output parameters:
%       p      : (n+1) polynomial coefficients in descending powers
%       S      : a structure that can be used as an input to polyval to obtain error estimates
%       mu     : a two-element vector with mean and standard deviation of x

w=max(w, 1e-5);

if ~isequal(size(x),size(y)) 
    error('X and Y vectors must be the same size.') 
end 
 
x = x(:);
y = y(:);
 
if nargout > 2 
   mu = [mean(x); std(x)]; 
   x = (x - mu(1))/mu(2); 
end 
 
% Construct Vandermonde matrix. 
V(:,n+1) = sqrt(w); 
for j = n:-1:1 
   V(:,j) = x.*V(:,j+1); 
end 
 
% Solve least squares problem, and save the Cholesky factor. 
[Q,R] = qr(V,0); 
ws = warning('off','all');  
p = R\(Q'*(sqrt(w).*y));    % Same as p = V\y; 
warning(ws); 
if size(R,2) > size(R,1) 
   warning('MATLAB:polyfit:PolyNotUnique', ... 
       'Polynomial is not unique; degree >= number of data points.') 
elseif condest(R) > 1.0e10 
    if nargout > 2 
        warning('MATLAB:polyfit:RepeatedPoints', ... 
            'Polynomial is badly conditioned. Remove repeated data points.') 
    else 
        warning('MATLAB:polyfit:RepeatedPointsOrRescale', ... 
            ['Polynomial is badly conditioned. Remove repeated data points\n' ... 
            '         or try centering and scaling as described in HELP POLYFIT.']) 
    end 
end 
r = y - V*p; 
p = p.';          % Polynomial coefficients are row vectors by convention. 
 
% S is a structure containing three elements: the Cholesky factor of the 
% Vandermonde matrix, the degrees of freedom and the norm of the residuals. 
 
S.R = R; 
S.df = length(y) - (n+1); 
S.normr = norm(r); 