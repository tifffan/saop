% Notations are from W. Gautschi, "Orthogonal Polynomials"
% Weights computed using Lanczos CDF approach

function [recurr_coeffs,absc,weights,Pi]=matrix_adapted_ortho_poly_Lan(G,Mdeg,nvec,poly_order)

if ~isfield(G,'lmax')
    G=gsp_estimate_lmax(G);
end
 grid_order=poly_order+1;

a=0;
b=G.lmax;

% Create the grid points
absc=(b-a)/(2*(grid_order-1)):(b-a)/(grid_order-1):b-a; %linspace(a,b,grid_order);
absc=[a,absc];
%weights=zeros(size(absc));
%weights(1)=1; % know a priori that there is an eigenvalue at 0


% Obtain weights using the Lanczos CDF idea 
weights= LanczosCDOS(G.L, Mdeg, nvec,absc)';
weights(2:end)=weights(2:end)-weights(1:end-1);

% Compute the recurrence coefficients from the discrete measure, using the
% Lanczos method described in Gautschi, Section 2.2.3.2 (pp. 97-98)
xw=[absc',weights'];
%tic
recurr_coeffs=lanczos(poly_order+1,xw); % from https://www.cs.purdue.edu/archives/2002/wxg/codes/OPQ.html
Pi=eval_pi(recurr_coeffs,absc');
%lanczos_time=toc
%tic
%[recurr_coeffs,Pi]=stieltjes(poly_order+1,xw);
%stieltjes_time=toc
%diff=abs(rSecurr_coeffs-recurr_coeffs2)