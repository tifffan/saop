function r=roots_cubic(coef)
%ROOTS_CUBIC: Solves the cubic equation ax^3+bx^2+cx+d=0 for x
%
%   Usage: r = three_term_recurr_op(G, recurr_coeffs, c, signal);
%
%   Input parameters:
%       coef    : Coefficients a, b, c, d
%   Output parameters:
%       r       : Three cubic roots

a=coef(1);
b=coef(2);
c=coef(3);
d=coef(4);
r=zeros(3,1);

delta=18*a*b*c*d-4*b^3*d+b^2*c^2-4*a*c^3-27*a^2*d^2;
delta0=b^2-3*a*c;
delta1=2*b^3-9*a*b*c+27*a^2*d;
X=(delta1+sqrt(delta1^2-4*delta0^3))/2;
C=X.^(1/3);

if delta==0
    if delta0==0
        r=-b/3/a*ones(3,1);
    else 
        r(1)=(9*a*d-b*c)/(2*delta0);
        r(2)=(9*a*d-b*c)/(2*delta0);
        r(3)=(4*a*b*c-9*a^2*d-b^3)/a/delta0;
    end
else
    r(1)=-(b+C+delta0/C)/3/a;
    r(2)=-(b+(-1/2+sqrt(3)*1i/2)*C+delta0/C/(-1/2+sqrt(3)*1i/2))/3/a;
    r(3)=-(b+(-1/2+sqrt(3)*1i/2)^2*C+delta0/C/(-1/2+sqrt(3)*1i/2)^2)/3/a;
end




