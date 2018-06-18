function r = poly_op(G, c, signal)
%POLY_OP Evaluate filtered signal with polyfit coefficients
%   Usage:  r = poly_op(G, c, signal)
%
%   Input parameters:
%       G       : Graph structure
%       c       : Polyfit(LS fitting) coefficients
%       signal  : Signal to filter
%   Output parameters
%       r       : Result
%

n=size(c,2);
r=c(1)*signal;

for i=2:n
    r=c(i)*signal+G.L*r;
end
