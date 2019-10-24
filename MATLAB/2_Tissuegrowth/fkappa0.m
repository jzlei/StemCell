%
% This function defines the time dependent differentiation rate 
%

function k0=fkappa0(t)
global T0 par
if T0 < 0
    k0 = par.kappa0;
else
    k0 = par.kappa0 * 1.0/(1.0 + exp(-(t-T0)/1000)); 
end
end