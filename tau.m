%
% The function for the delay tau(x)
%

function z=tau(x)
global par
z=par.tau * ones(1,size(x,2));
end