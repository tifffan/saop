function [V,H,orth] = lanczos(A,order,x)

[N,M] = size(x);

% normalization
norm2vec = @(x) (sum(x.^2,1)).^0.5;
q = x./repmat(norm2vec(x),N,1);

% Initialization
hiv =0:order:(order*M-1); % helping indice vector

V = zeros(N,M*order);
V(:,1+hiv) = q;


H = zeros(order+1,order*M);

r = A*q;
H(1,1+hiv) = sum(q .* r, 1 );
r = r - repmat(H(1,1+hiv),N,1).*q; 
H(2,1+hiv) = norm2vec(r);

if (nargout > 2)
    orth = zeros(M,1);
    orth(1) = norm(V'*V - M);
end

for k = 2:order
    
    if (sum(abs(H(k,k-1+hiv))) <= eps)
        H = H(1:k-1,sum_ind(1:k-1,hiv));
        V = V(:,sum_ind(1:k-1,hiv));
        if (nargout > 2)
            orth = orth(1:k-1);
        end
        return;
    end
    
    H(k-1,hiv+k) = H(k,hiv+k-1);
    v = q;
    q = r./repmat(H(k-1,k+hiv),N,1);
    V(:,k+hiv) = q;
  
    r = A*q;
    r = r - repmat(H(k-1,k+hiv),N,1).*v;
    H(k,k+hiv) = sum(q .* r, 1 );
    
    r = r - repmat(H(k,k+hiv),N,1).*q;
    % The next line has to be checked
    r = r - V*(V'*r); % full reorthogonalization
    H(k+1,k+hiv) = norm2vec(r);
    
    if (nargout > 2)
        orth(k) = [orth, norm(V'*V - M)];
    end
end
   
H = H(1:order,1:order);