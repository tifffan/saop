function yy = InvLanczosCDOS(A, Mdeg, nvec,lmax,x,tx)
%% function y = InvLanczosCDOS(A, Mdeg, nvec,x)
%% Lanczos approximation for computing the Cummulative DOS.
% A = input symmetric matrix
% Mdeg = degree of polynomial to be used
% nvec = number of starting vectors
% x = discretized points in the interval [lmin,lMax]
%% Output:
% y = the corresponding values of CDOS

%acc_theta=[];acc_eta=[];
acc_ym=[];
order=length(tx);
%% Mdeg Lanczos steps for different starting vectors
for m = 1 : nvec
	%% Mdeg steps of Lanczos (full reorthogonalization) 
    %% by Arnoldi 
    w = randn(size(A,1),1); % starting vector
	v0 = w /norm(w);
	[H,V,f] = arnoldi_k(A,v0,Mdeg); % Arnoldi works well for symmetrix matrices
	H = (H+H')/2;
	[eigvec,D]=eig(H);
	theta  = real(diag(D));  % eigenvalues of Tridiagonal matrix
	gamma2 = eigvec(1,:).^2; % square of top entry of eigenvectors
    eta2=cumsum(gamma2)';
	%% Accumulate CDOS 
    y_m=gsp_mono_cubic_warp_fn(theta,eta2,x);
   acc_ym=[acc_ym,y_m];
%     acc_theta=[acc_theta;theta];
%     acc_eta=[acc_eta;eta2'];
end

%% -------------To plot the DOS
%[xx,idx2]=sort(acc_theta,'ascend');
%yy=acc_eta(idx2);
yhat=mean(acc_ym,2);
xfull=0:0.01:lmax;
yfull=gsp_mono_cubic_warp_fn(x,yhat,xfull');

%% --Inverse
yy=zeros(order,1);
for i=1:length(tx)
    [~,idx]=min(abs(yfull-tx(i,1)));
    yy(i,1)=xfull(idx(1));
end
