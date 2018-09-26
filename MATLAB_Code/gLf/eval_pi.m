function [ Pi ] = eval_pi(ab,absc)
%EVAL_PI: Evaluates the first K orthogonal polynomials {\pi_k(x)} for a 
%discrete set of N points via the three term recursion
%
%   Usage: Pi = eval_pi(ab,absc);
%
%   Input parameters:
%       ab      : Recurrence coefficients alpha (first column) and beta
%                 (second column)
%       absc    : Abscissae of a discrete set of N points {x_i}, the support
%                 points of a discrete measure with respect to which the 
%                 polynomials {\pi_k} are orthogonal
%
%   Output parameters:
%       Pi      : values of the orthogonal polynomials at the abscissae,
%                 (i-1)th polynomial in ith column of Pi
% 

Ncap=length(absc);
N=size(ab,1);
Pi=zeros(Ncap,N);
Pi(:,1)=ones(Ncap,1);
Pi(:,2)=((absc-ab(1,1)).*Pi(:,1))/sqrt(ab(2,2));
for i=3:N
    Pi(:,i)=((absc-ab(i-1,1)).*Pi(:,i-1)-sqrt(ab(i-1,2)).*Pi(:,i-2))/sqrt(ab(i,2));
end
end

