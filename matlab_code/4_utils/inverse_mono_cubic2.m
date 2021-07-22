function [x0_values]=inverse_mono_cubic2(x,y,y0)
%INVERSE_MONO_CUBIC2: Interpolates a monotonic cubic function f
%through points (x,y) and computes the x0 values given y0=f(x0)
%
% Usage: x0 = inverse_mono_cubic2(x,y,y0);
%
% Input parameters:
%     x         :a vector of x coordinates of the known points
%     y         :a vector of y coordinates of the known points
%     y0        :a vector of y coordinates of the points to be interpolated
%
% Output parameters:
%     x0        :a vector of x coordinates such that y0=f(x0)

cut=1e-4;

[x,x_ind]=sort(x,'ascend');
y=y(x_ind);
if ( isequal(sort(y,'ascend'),y)==0 && isequal(sort(y,'descend'),y)==0 )
    error('Data points are not monotonic');
end

% Monotonic cubic interpolation using the Fritsch-Carlson method
num_pts=length(x);
if length(y) ~= num_pts
    error('x and y vectors have different dimensions');
end

% 1. Compute slopes of secant lines
Delta=(y(2:end)-y(1:num_pts-1))./(x(2:end)-x(1:num_pts-1));

% 2. Initialize tangents m at every data point
m = (Delta(1:num_pts-2)+Delta(2:num_pts-1))/2;
m = [Delta(1);m;Delta(end)];

% 3. Check for equal y's to set slopes equal to zero
for k=1:num_pts-1
    if Delta(k)==0
        m(k)=0;
        m(k+1)=0;
    end
end

% 4. Initialize alpha and beta
alpha = m(1:num_pts-1)./Delta;
beta = m(2:num_pts)./Delta;

% 5. Make monotonic
for k=1:num_pts-1
    if alpha(k)^2+beta(k)^2 > 9
        tau=3/sqrt(alpha(k)^2+beta(k)^2);
        m(k)=tau*alpha(k)*Delta(k);
        m(k+1)=tau*beta(k)*Delta(k);
    end
end

% 6. Cubic interpolation
coef_t=zeros(4,num_pts-1);
for lower_ind=1:num_pts-1
    h=x(lower_ind+1)-x(lower_ind);
    coef_t(1,lower_ind)=y(lower_ind)*2+h*m(lower_ind)+y(lower_ind+1)*(-2)+h*m(lower_ind+1);
    coef_t(2,lower_ind)=y(lower_ind)*(-3)+h*m(lower_ind)*(-2)+y(lower_ind+1)*3+h*m(lower_ind+1)*(-1);
    coef_t(3,lower_ind)=h*m(lower_ind);
    coef_t(4,lower_ind)=y(lower_ind);
end

num_pts_to_interpolate=length(y0);
x0_values=zeros(size(y0));
%derivative_values=zeros(size(x0));


%tic
for ii=1:num_pts_to_interpolate
    if y0(ii) < y(1)
        x0_values(ii)=x(1);
    elseif y0(ii) >= y(end)-cut
        x0_values(ii)=x(end);
    else
    [~,closest_ind]=min(abs(y-y0(ii)));
    %if sign(x(closest_ind)-x0(i))<0 || ( sign(x(closest_ind)-x0(i))==0 && closest_ind < num_pts)
    if (y(closest_ind)-y0(ii))<-cut || ( abs(y(closest_ind)-y0(ii))<cut && closest_ind < num_pts)
        lower_ind=closest_ind;
    else
        lower_ind=closest_ind-1;
    end
    
    % knowing that y is between y(lower_ind) and y(lower_ind+1), solve for
    % t and then convert to x
    adj_coef_t=coef_t(:,lower_ind)';
    adj_coef_t(1,4)=adj_coef_t(1,4)-y0(ii);
    
% Symbolic Math Toolbox needed for this part
%     syms x
%     eqn= adj_coef_t(1,1)*x^3+adj_coef_t(1,2)*x^2+adj_coef_t(1,3)*x+adj_coef_t(1,4) ==0;
%     r=solve(eqn,x,'Real',true);
%     t=max(r.*(r>=0).*(r<1-cut));
    r=roots_cubic(adj_coef_t);
%    r=roots(adj_coef_t);
%     [~,i]=min(abs(imag(r)));
%     t=r(i);
    t=real(max((abs(imag(r))<1e-6).*r.*(r>=0).*(r<1-cut)));
    h=x(lower_ind+1)-x(lower_ind);
    x0_values(ii) = t*h+x(lower_ind);
    end
end

%time_newmethod=toc



