function [xx, yy,xx2,yy2] = landos(A, Mdeg, Npts, nvec, lmin, lMax, sigma)
%% function [xx, yy] = landos(A, Mdeg, Npts, nvec)
%% Lanczos approximation metohd  for computing the DOS.
%% A = input symmetric matrix
%% Mdeg = degree of polynomial to be used
%% Npts = Number of points in DOS (discretization)
%% nvec = number of starting vectors
%% lmin = Smallest eigenvalue of A
%% lMax  = Largest eigenvlaue of A
%% sigma = scaling for Gaussian blurring
%% Output:
% xx = discretized points in the interval [lmin,lMax]
% yy = the corresponding values of DOS
%% Discretization of the interval
x = linspace( lmin, lMax, Npts )';
y = zeros(Npts,1);
h = x(2)-x(1);

%% Parameters for Gaussian blurring
%%----------- if gaussian smaller than tol ignore point.
tol    = 1.e-02;
width  = sigma*sqrt(-2 *log(tol));
sigma2 = 2*sigma^2;
x2=[];y2=[];
%% Mdeg Lanczos steps for different starting vectors
for m = 1 : nvec
	%% Mdeg steps of Lanczos (full reorthogonalization) 
    %% by Arnoldi 
    w = randn(size(A,1),1); % starting vector
	v0 = w /norm(w);
	[H,V,f] = arnoldi_k(A,v0,Mdeg); 
	H = (H+H')/2;
	[eigvec,D]=eig(H);
	theta  = real(diag(D));  % eigenvalues of Tridiagonal matrix
	gamma2 = eigvec(1,:).^2; % square of top entry of eigenvectors
    eta2=cumsum(gamma2);
	%% Accumulate DOS
	for i=1:size(theta)
		t = theta(i);        
		ind = find( abs(x - t) < width );
		y(ind) = y(ind) + eta2(i) * exp(-(x(ind)-t).^2/sigma2);
    end 
    x2=[x2,theta'];
    y2=[y2,eta2];
end

%% -------------To plot the DOS
[xx2,idx2]=sort(x2,'ascend');
yy2=y2(idx2);
%y2=sort(y2,'ascend');
xx  = x; 
yy  = y';
yy2= gsp_mono_cubic_warp_fn(xx2',y2',xx);

%yy  = yy / (sum(yy)*(xx(2)-xx(1)));
end