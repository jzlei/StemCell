function b0=beta0(x)
global par
n=size(x,2);
b0=zeros(1,n);
for i=1:n
    b0(i) = par.betabar * (par.a1*x(i) + power(par.a2*x(i), 6))/(1.0 + (par.a3*x(i))^6);
end
end