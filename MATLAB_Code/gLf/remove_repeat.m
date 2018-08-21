function [x2, weights2, ind]=remove_repeat(x, delta, weights)

n=size(x,1);
x2=x;

if nargin < 3
    weights2=ones(n,1);
else
    weights2=weights;
end

ind=1:n;
i=2;
while i<= size(x2,1)
    if abs(x2(i)-x2(i-1))< delta
        weights2(i-1)=weights2(i-1)+weights2(i);
        x2(i)=[];
        weights2(i)=[];
        ind(i)=[];
        i=i-1;
    end
    i=i+1;
end
