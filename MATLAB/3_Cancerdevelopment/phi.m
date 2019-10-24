function z=phi(x)
global t T0
%if (t<T0)
    z=phidev(x);
%else
%    z=phiTissueGrowth(x);
%end
end

function z=phidev(x)
a=1.65;
s=1.8;
z=0.08 + 1.0 *power(a*x,s)/(1 + power(a*x,s));
end

function z=phiaging(x)
a=1.65;
s=1.8;
z=0.04 + 0.9 *power(a*x,s)/(1 + power(a*x,s));
end

function z=phiabnormalgrowth(x)
a=1.65;
s=1.8;
z=0.08 + 1.2 *power(a*x,s)/(1 + power(a*x,s));
end

function z=phiTissueGrowth(x)
a=1.65;
s=2.0;
z=0.08 + 1.05 *power(a*x,s)/(1 + power(a*x,s));
end
