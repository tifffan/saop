function [ Pi ] = eval_pi(ab,absc)

Ncap=length(absc);
N=size(ab,1);
Pi=zeros(Ncap,N);
Pi(:,1)=ones(Ncap,1);
Pi(:,2)=(absc-ab(1,1)).*Pi(:,1);
for i=3:N
    Pi(:,i)=(absc-ab(i-1,1)).*Pi(:,i-1)-ab(i-1,2).*Pi(:,i-2);
   % Pi(:,i)=Pi(:,i)/sum(Pi(:,ii).^2.*weights)
end
end

