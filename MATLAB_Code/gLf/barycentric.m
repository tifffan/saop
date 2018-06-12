function [p,y0]=barycentric(x,y,x0)
%BARYCENTRIC Barycentric interpolation
%   Usage: p=barycentric(x,y,x0);
%
%   Input parameters:
%       x       : vector of x-coordinates of interpolating points
%       y       : vector of y-coordinates of interpolating points
%       x0      : vector of x-coordinates of points to interpolate
%   Output parameters
%       p       : interpolated function through data (x, y)
%       y0      : interpolated values at x0
%
% Notations are from Thm 5.1, p.36, Approximation Theory and Approximation
% Practice by L. Trefethen

n=length(x);
if length(y) ~= n
    error('x and y vectors have different dimensions');
end

lambda=zeros(n,1);
for i=1:n
    diffs=x(i)*ones(n,1)-x;
    diffs(i)=1;
    lambda(i)=1e25/prod(diffs);
end

p=@(s) sum(lambda.*y./(s-x))/sum(lambda./(s-x));

n_interpolate=length(x0);
y0=zeros(size(x0));

for j=1:n_interpolate
    y0(j)=p(x0(j));
end


