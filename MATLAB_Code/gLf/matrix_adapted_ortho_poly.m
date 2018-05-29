% Notations are from W. Gautschi, "Orthogonal Polynomials"

function [recurr_coeffs,absc,weights,Pi]=matrix_adapted_ortho_poly(G,poly_order,num_vec,num_its,param)

if ~isfield(G,'lmax')
    G=gsp_estimate_lmax(G);
end

if ~isfield(param,'jackson')
    param.jackson=1;
end

if ~isfield(param,'grid_order')
    grid_order=poly_order+1;
else
    if param.grid_order< (poly_order+1)s
        error('Grid order must be greater than poly order');
    else
        grid_order=param.grid_order;
    end
end

if ~isfield(param,'init_poly_order')
    param.init_poly_order=30;
end

a=0;
b=G.lmax;

% Initialize discrete measure
% TODO: Experiment with initial placement of points and whether to force
% one at 0
absc=(b-a)/(2*(grid_order-1)):(b-a)/(grid_order-1):b-a; %linspace(a,b,grid_order);
absc=[a,absc];
weights=zeros(size(absc));
weights(1)=1; % know a priori that there is an eigenvalue at 0

ch=zeros(param.init_poly_order+1,grid_order-1);
jch=zeros(param.init_poly_order+1,grid_order-1);

Sig = randn(G.N,num_vec);
%Sig = Sig*diag(1./sqrt(sum(Sig.^2,1)));
    
for j=1:grid_order-1
    [ch(:,j), jch(:,j)] = gsp_jackson_cheby_coeff(0,absc(j+1),[a,b], param.init_poly_order);
end

% Filtering
if param.jackson
    X = gsp_cheby_op(G, jch, Sig);
else
    X = gsp_cheby_op(G, ch, Sig);
end

St=Sig';
quad_sums=zeros(grid_order-1,1);
for j=1:grid_order-1
    for i=1:num_vec
        quad_sums(j)=quad_sums(j)+St(i,:)*X((j-1)*G.N+1:j*G.N,i)/num_vec;
    end
end
weights(2)=quad_sums(1)-weights(1);
weights(3:end)=quad_sums(2:end)-quad_sums(1:grid_order-2);
weights=weights/sum(weights);

% TODO: add some iteration here to update the weights and maybe also update
% the support pts?


% Compute the recurrence coefficients from the discrete measure, using the
% Lanczos method described in Gautschi, Section 2.2.3.2 (pp. 97-98)
xw=[absc',weights']; %TODO: place weights at absc or in between?
%tic
recurr_coeffs=lanczos(poly_order+1,xw); % from https://www.cs.purdue.edu/archives/2002/wxg/codes/OPQ.html
Pi=eval_pi(recurr_coeffs,absc');
%lanczos_time=toc
%tic
%[recurr_coeffs,Pi]=stieltjes(poly_order+1,xw);
%stieltjes_time=toc
%diff=abs(recurr_coeffs-recurr_coeffs2)