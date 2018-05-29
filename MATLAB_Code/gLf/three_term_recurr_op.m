function r = three_term_recurr_op(G, ab, c, signal,param)
%THREE_TERM_RECURR_OP : Orthogonal polynomial of graph Laplacian applied to vector
%   Usage: r = three_term_recurr_op(G, recurr_coeffs,c, signal)
%
%   Input parameters:
%       G       : Graph structure
%       ab      : Recurrence coefficients alpha (first column) and beta
%                 (second column)
%       c       : Expansion coefficients in the chosen polynomials
%       signal  : Signal(s) to filter (one per column)
%   Output parameters
%       r       : Result of the filtering
% 
%   Compute (possibly multiple) polynomials of graph laplacian applied to input.
%
%   Coefficients for multiple polynomials may be passed as a matrix.
%   This is equivalent to setting::
%
%       r(1) = three_term_recurr_op(G, ab, c(:,1), signal);
%       r(2) = three_term_recurr_op(G, ab, c(:,2), signal);
%       ...
% 
%   but is more efficient as the polynomials of G.L applied
%   to *signal* can be computed once and shared.
%
%   The output *r* is a matrix with each column corresponding to a filter.
%
%   *param* contain only one field param.verbose to controle the verbosity.
%
%   Example:::
%
%         Nf = 4;
%         G = gsp_sensor(100);
%         G = gsp_estimate_lmax(G);
%         poly_order=30;
%         num_vec=20;
%         num_its=1;
%         [ab,absc,weights,Pi]= matrix_adapted_ortho_poly(G,poly_order,num_vec,num_its);
%         g = gsp_design_meyer(G, Nf);  
%         c = matrix_adapted_poly_coeff(G, g , absc, weights);
%         f = rand(G.N,1);
%         r = gsp_three_term_op(G, ab, c, f);
%
%   This function is inspired by the gsp_cheby_op function from the GSPBox, which was adapted from sgwt_toolbox
%
%   See also: matrix_adapted_poly_coeff matrix_adapted_ortho_poly gsp_cheby_op gsp_cheby_coeff gsp_filter_analysis
%

% Author: David Shuman, David K Hammond, Nathanael Perraudin
% Testing: 
% Date: 1 March 2018

if nargin < 5
    param = struct;
end

if ~isfield(param,'verbose'), param.verbose = 1; end;


Nscales=size(c,2);

M = size(c,1);
% To handle different order of approximation


assert(all(M>=2));

maxM=max(M);



if ~isfield(G,'lmax');
    G = gsp_estimate_lmax(G);
    if param.verbose
    warning(['THREE_TERM_RECURR_OP: The variable lmax is not ',...
        'available. The function will compute it for you. ',...
        'However, if you apply many time this function, you ',...
        'should precompute it using the function: ',...
        'gsp_estimate_lmax']);
    end
end

if isa(signal,'single')
    signal = double(signal); 
    bsingle = 1;
else
    bsingle = 0;
end

arange = [0, G.lmax];

a1 = (arange(2) - arange(1))/2;
a2 = (arange(2) + arange(1))/2;


%Twf_new = T_j(L) f
%Twf_cur T_{j-1}(L) f
%TWf_old T_{j-2}(L) f

Twf_old=signal;                     % j = 0;
Twf_cur=G.L*signal-ab(1,1)*signal;  % j = 1;

Nv = size(signal,2);
r = zeros(G.N*Nscales,Nv);

for ii=1:Nscales
    r((1:G.N)+G.N * (ii-1),:) = c(1,ii) * Twf_old + c(2,ii) * Twf_cur;
end

for k=2:maxM
    Twf_new = G.L*Twf_cur-ab(k,1)*Twf_cur - ab(k,2)*Twf_old; % this is the three-term recurrence relation in (1.3.2) in Gautschi. Should we be using the on in (1.3.13) instead?
    for ii=1:Nscales
        if 1+k <= M
            r((1:G.N)+G.N * (ii-1),:) =...
                r((1:G.N)+G.N * (ii-1),:) + c(k+1,ii)*Twf_new;
        end
    end
    Twf_old=Twf_cur;
    Twf_cur=Twf_new;
end


if bsingle
    r = single(r);
end

end
