%
% The function kappa(x)
%

function y=kappa(x)
global par 
n=size(x,2);
y=zeros(1,n);
for i=1:n
    y(i)=par.kappa0 * 1.0/(1.0 + power(par.b1 * x(i), 6));
end
end