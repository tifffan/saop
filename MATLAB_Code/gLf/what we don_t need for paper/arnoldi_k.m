function [H,V,f] = arnoldi_k(A,v0,k)
%
% k-step basic arnoldi iteration
%
% usage:    [H,V,f] = arnoldi_k(A,v0,k)
%

% INITIALIZATION
beta   = norm(v0);
v      = v0 / beta;
V(:,1) = v;
Av     = A*v;
h      = V(:,1)'*Av;
H(1,1) = h;
f      = Av - V(:,1)*h;
% ONE STEP OF REORTHOGONALIZATION
s    = V(:,1)'*f;
h    = h + s;
f    = f - V(:,1)*s;
%
%
% MAIN LOOP
for j = 2:k
	beta = norm(f);
	H(j,j-1) = beta;
	v    = f/beta;
	V(:,j) = v;
	Av   = A*v;
	h    = V'*Av;
	f    = Av - V*h;
	% ONE STEP OF REORTHOGONALIZATION
	s    = V'*f;
	h    = h + s;
	f    = f - V*s;
	%
	H(1:j,j) = h;
end;
