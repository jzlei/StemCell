%
% This function defines the differentiation rate kappa = kappa(t, x)
%

function y=kappa(t,x)
global par 
n=size(x,2);
y=zeros(1,n);
k0=fkappa0(t);
for i=1:n
    y(i) = k0 * 1.0/(1.0 + power(par.b1 * x(i), 6));
end
end