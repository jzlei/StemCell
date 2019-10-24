%
% The function to set up the model parameters
% Set the parameter Flag for different applications 
%

function parameter()
global par
global Flag

par=struct;
par.theta=1e3;
par.a1=5.8;
par.a2=2.2;
par.a3=3.75;
par.b1=4.0;
par.mu=2.0e-4;
par.tau=20;
par.nu0=2.0e-5;

switch Flag
    case 1  %%  For normal development (Figure 2 B-C)
        par.betabar = 0.12;
        par.kappa0 = 0.02;
    case 2  %% For the aging process (Figure 2 D-E)
        par.betabar = 0.12;
        par.kappa0 = 0.02;
    case 3  %% For the abnormal growth (Figure 2 F-G)
        par.betabar = 0.12;
        par.kappa0 = 0.0002;
    case 4  %% For the situation of subculture (Figure 3)
        par.betabar = 4.8;
        par.kappa0 = 0.4;
end
        
end