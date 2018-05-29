function [ G, vals ] = spectral_cdf_approx2( G , param)

if nargin<2
   param = struct;
end

if ~isfield(G,'lmax')
    G=gsp_estimate_lmax(G);
end

if ~isfield(param, 'order')
    param.order = 30;
end

if ~isfield(param, 'num_vec')
    param.num_vec = 30;
end

if ~isfield(param, 'num_pts')
    param.num_pts = 50;
end

if ~isfield(param, 'pts')
    param.pts=linspace(0,G.lmax,param.num_pts);
end

if ~isfield(param,'cdf_method')
    cdf_method='lanczos';
else
    cdf_method=param.cdf_method; % 'lanczos' or 'kpm' or 'ldl'
end

switch cdf_method

    case 'kpm'
    G=gsp_compute_TkbarLX(G, param);

    vals=zeros(param.num_pts,1);
    vals(1)=0; % force to be 0 or 1?
    vals(end)=G.N;
    %TODO: force last one to be N?

    for j=2:param.num_pts-1
        [~, jch] = gsp_jackson_cheby_coeff(0,param.pts(j),[0,G.lmax], param.order);
        r=gsp_cheby_opX(G,jch);
        vals(j)=gsp_hutch(G,r);
    end
    vals=vals/G.N;
    if(vals(1)>vals(2))
        vals(1)=vals(2);
    end
    vals=min(vals,1);
    G.spectrum_cdf_approx = @(s) gsp_mono_cubic_warp_fn(param.pts',vals,s);
    
    case 'ldlt'
        if ~isfield(param, 'use_speedup')
            if ( exist('ldlsymbol_extra','file')==3 && exist('ldlnumeric','file')==3 )
                use_speedup=1;
            else
                use_speedup=0;
            end
        else
            use_speedup = param.use_speedup;
        end

        if ~isfield(param, 'use_ldl_package')
            if exist('ldlsparse','file')==3
                use_ldl_package=1;
            else
                use_ldl_package=0;
            end
        else
            use_ldl_package = param.use_ldl_package;
        end

        if ~isfield(param, 'use_permutation')
            use_permutation=1;
        else
            use_permutation = param.use_permutation;
        end
        
        if ~isfield(param, 'ldl_thresh')
            ldl_thresh=.001;
        else
            ldl_thresh = param.ldl_thresh;
        end
        
        counts=zeros(param.num_pts,1);
        counts(param.num_pts)=G.N-1;

        interp_x=(0:param.num_pts-1)*G.lmax/(param.num_pts-1);
        interp_x=interp_x';

        identity=speye(G.N);

        if use_speedup
            if use_permutation
                P=symamd(G.L);
                [Parent, Lp, PO, PIn, flopcount] = ldlsymbol_extra(G.L,P);
                for i=1:param.num_pts-2 
                    mat=G.L-interp_x(i+1)*identity;
                    [~, HD]=ldlnumeric(mat,Lp,Parent,PO,PIn);
                    counts(i+1)=sum(diag(HD)<0);
                end
            else
                [Parent, Lp, flopcount] = ldlsymbol_extra(G.L);
                for i=1:param.num_pts-2 
                    mat=G.L-interp_x(i+1)*identity;
                    [~, HD]=ldlnumeric(mat,Lp,Parent);            
                    counts(i+1)=sum(diag(HD)<0);
                end
            end
        else
            if use_permutation
                P=symamd(G.L);
                for i=2:param.num_pts-1 
                    mat=G.L-interp_x(i)*identity;
                    if use_ldl_package
                        [~,HD]=ldlsparse(mat,P);
                    else
                        [~,HD,~]=ldl(mat(P,P),ldl_thresh);
                    end
                    counts(i)=sum(diag(HD)<0);
                end
            else
                for i=2:param.num_pts-1 
                    mat=G.L-interp_x(i)*identity;
                    if use_ldl_package
                        [~,HD]=ldlsparse(mat);
                    else
                        [~,HD,~]=ldl(mat,ldl_thresh);
                    end
                    counts(i)=sum(diag(HD)<0);
                end
            end
        end
        interp_y=counts/(G.N-1);

        G.spectrum_cdf_approx = @(s) gsp_mono_cubic_warp_fn(interp_x,interp_y,s);


    case 'lanczos' % cdf_method = 'lanczos', adapted from LanczosCDOS; see also Appendix C of "Approximating Spectral Densities of Large Matrices" by Lin, Saad, and Yang
        G.spectrum_cdf_approx=@(x) 0;
        for m =1:param.num_vec
            %% Mdeg steps of Lanczos (full reorthogonalization) 
            %% by Arnoldi 
            w = randn(G.N,1); % starting vector
            v0 = w /norm(w);
            [H,~,~] = arnoldi_k(G.L,v0,param.order); % Arnoldi works well for symmetrix matrices
            H = (H+H')/2;
            [eigvec,D]=eig(H);
            theta  = real(diag(D));  % eigenvalues of Tridiagonal matrix
            gamma2 = eigvec(1,:).^2; % square of top entry of eigenvectors
            eta2=cumsum(gamma2)';
            %% Accumulate CDOS 
            G.spectrum_cdf_approx=@(x)G.spectrum_cdf_approx(x)+(1/param.num_vec)*gsp_mono_cubic_warp_fn(theta,eta2,x);
        end
        
    otherwise
        error('Unknown spectral cdf approximation method');
end


end

