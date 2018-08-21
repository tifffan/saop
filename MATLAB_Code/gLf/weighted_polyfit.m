function [p,S,mu] = weighted_polyfit(x,y,w,n)

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
V(:,n+1) = w; 
for j = n:-1:1 
   V(:,j) = x.*V(:,j+1); 
end 
 
% Solve least squares problem, and save the Cholesky factor. 
[Q,R] = qr(V,0); 
ws = warning('off','all');  
p = R\(Q'*(w.*y));    % Same as p = V\y; 
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