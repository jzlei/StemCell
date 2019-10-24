%
% The clean rate of terminatlly differentiated cells
%
function nu=fnu(x)
global par
n=size(x,2);

nu=par.nu0 * ones(1,n);
end