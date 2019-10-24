function p0=p(x,y)
n=size(x,2);
p0=zeros(n,n);
for i=1:n
    for j=1:n
        a=eta(y(i))*phi(y(i));
        b=eta(y(i))*(1-phi(y(i)));
        p0(i,j)=betapdf(x(j),a,b);
    end
end
end