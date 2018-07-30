% Adapted from Lin et.al, 'Approximating spectral densities of large
% matrices'
% The approximated pdf is not well defined at G.lmin and G.lmax

function phi=pdf_est(G,absc,poly_order)

num_vec=size(G.X,2);
mu=zeros(poly_order+1,1);

for k=1:poly_order+1
    r=G.TkbarLX(:,(k-1)*num_vec+1:k*num_vec);
    for i=1:size(G.X,2)
        r(:,i)=r(:,i);
    end
    s=gsp_hutch(G,r);
    mu(k)=s*2/G.N/pi;
end
mu(1)=mu(1)/2;

% Include damping factors against Gibbs oscillations
damping='sigma';
K=poly_order;
switch damping
    case 'jackson'
        gamma=ones(K+1,1);
        for k=1:K
            gamma(k+1,1)=((1-k/(K+2))*sin(pi/(K+2))*cos(k*pi/(K+2))...
                    +1/(K+2)*cos(pi/(K+2))*sin(k*pi/(K+2)))/sin(pi/(K+2));
        end
        mu=mu.*gamma;
    case  'sigma'
        sigma=ones(K+1,1);
        for k=1:K
            sigma(k+1,1)=sin(k*pi/(K+1))/(k*pi/(K+1));
        end
        mu=mu.*sigma;
    case 'none'
    otherwise
        error('damping type not recognized');
end


a=G.lmin; %lambda_min
b=G.lmax; %lambda_max
c=(a+b)/2;
d=(b-a)/2;

t=(absc-c)/d; %shift absc to [-1,1]
grid_order=length(absc);
Tkt_old=ones(grid_order,1);
Tkt_cur=t;
zeta=mu(1)*Tkt_old+mu(2)*Tkt_cur;
for k=3:poly_order+1
    Tkt_new=2*t.*Tkt_cur-Tkt_old;
    zeta=zeta+mu(k)*Tkt_new;
    Tkt_old=Tkt_cur;
    Tkt_cur=Tkt_new; 
end

phi=zeta./sqrt(1-t.^2);
% phi=phi/sum(phi);
% figure;stem(t,phi);
% figure;plot(t,cumsum(phi));

%Compare to derivative of CDF est from KPM
% xx = 0:0.001:G.lmax;
% delta=.1;
% G.spectrum_pdf_approx = @(x) (G.spectrum_cdf_approx(x+delta) - G.spectrum_cdf_approx(x-delta)) / (2*delta);
% figure;plot(xx,G.spectrum_pdf_approx(xx));

end