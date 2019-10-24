%
% The function for the apoptosis rate mu(x)
%

function z=mu(x)
global par
z=par.mu * ones(1,size(x,2));
end