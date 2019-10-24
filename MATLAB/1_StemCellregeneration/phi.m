function z=phi(x)
global t T0
global Flag

switch Flag
    case 1
        z=phidev(x);
    case 2
        z=phiaging(x);
    case 3
        z=phiabnormalgrowth(x);
    case 4
        z=phidev(x);
end
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


