function [absc,weights]=gen_absc_weights(G,poly_order,param)

if ~isfield(G,'lmax')
    G=gsp_estimate_lmax(G);
end

if ~isfield(param,'jackson')
    param.jackson=1;
end

if ~isfield(param,'grid_order')
    param.grid_order=poly_order+1;
else
    if param.grid_order< (poly_order+1)
        error('Grid order must be greater than polynomial order');
    else
        param.grid_order=param.grid_order;
    end
end

if ~isfield(param,'init_poly_order')
    param.init_poly_order=30;
end

if ~isfield(param,'num_vec')
    param.num_vec=30;
end

if (~isfield(G,'TkbarLX') || (poly_order > (size(G.TkbarLX,2)/size(G.X,2)-1)))
    Tk_param.order=poly_order;
    G=gsp_compute_TkbarLX(G,Tk_param);
end

a=G.lmin;
b=G.lmax;
grid_order=param.grid_order;

switch param.absc_method
    case 'linear'
        absc=linspace(G.lmin,G.lmax,grid_order)';
        
    case 'even'
        absc=(b-a)/(2*(grid_order-2)):(b-a)/(grid_order-2):b-a;
        absc=[a,a+absc,b];
        absc=absc';
        
    case 'warp'
        gi=@(s) G.spectrum_inv_cdf_approx((s-G.lmin)/(G.lmax-G.lmin));
        absc=gi((b-a)/2*sort((1+cos((0:(grid_order-1))*pi/(grid_order-1))),'ascend'));
        absc=absc';
        
    case 'spline'
        nodes=G.cdf_nodes;
        n=length(nodes);
        m=ceil(grid_order/(n-1))+1; % number of chebyshev extrema in each interval, one of them removed for repetition
        absc=zeros((n-1)*(m-1)+1,1);
        cheb_pts=sort(cos((0:m)*pi/m),'ascend');
        cheb_pts=cheb_pts(1:end-1);
        for i=1:n-1
            aa=nodes(i);
            bb=nodes(i+1);
            shift_cheb_pts=(cheb_pts+1)/2*(bb-aa)+aa;
            absc((m-1)*(i-1)+1:(m-1)*i)=shift_cheb_pts(1:end-1);
        end
        absc(end)=nodes(end);
        
    otherwise
        error('absc method not recognized');
end



switch param.weights_method
    case 'count'
        weights=zeros(1,grid_order);
        weights(1)=1; % know a priori that there is an eigenvalue at G.lmin

        ch=zeros(param.init_poly_order+1,grid_order-2);
        jch=zeros(param.init_poly_order+1,grid_order-2);

        for j=1:grid_order-3
            [ch(:,j), jch(:,j)] = gsp_jackson_cheby_coeff(a,(absc(j+1)+absc(j+2))/2,[a,b], param.init_poly_order);
        end
        [ch(:,grid_order-2), jch(:,grid_order-2)] = gsp_jackson_cheby_coeff(a,b,[a,b], param.init_poly_order);

        if param.jackson
            r = gsp_cheby_opX(G,jch);
        else
            r = gsp_cheby_opX(G,ch);
        end

        num_vec=size(G.X,2);
        quad_sums=zeros(grid_order-2,1);
        for j=1:grid_order-2
            for i=1:num_vec
                quad_sums(j)=quad_sums(j)+gsp_hutch(G,r(:,((j-1)*num_vec+1):j*num_vec));
            end
        end
        weights(2)=quad_sums(1)-weights(1);
        weights(3:grid_order-1)=max(0,quad_sums(2:end)-quad_sums(1:grid_order-3));
        weights(grid_order)=1;
        weights(grid_order-1)=weights(grid_order-1)-weights(grid_order);
        weights=real(weights');
        
    case 'constant'
        weights=ones(grid_order,1)/grid_order;
        
    case 'pdf'
        weights=G.spectrum_pdf_approx(absc);
 
    case 'pdf2'
        absc(1)=[];
        absc(end)=[];
        weights=pdf_est(G,absc,param.poly_order);
        
    otherwise
        error('weights method not recognized');
end

end
