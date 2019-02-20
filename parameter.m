%
% The function to set up the model parameters
%

function parameter()
global par
global t T0
beta_dev=0.12;  % Rates for normal development
kappa_dev=0.02;
%kappa_dev=0.00005;  % The rate for the simulation of tissue growth.

beta_aging=0.12; % Rates for aging process
kappa_aging=0.02;

beta_ag=0.12;   % Rates for abnormal growth
kappa_ag=0.0002; 

beta_HSC=0.12;
kappa_HSC=0.002;

par=struct;
par.theta=1e3;
par.a1=5.8;
par.a2=2.2;
par.a3=3.75;
par.b1=4.0;
par.mu=2.0e-4;
par.tau=20;
par.nu0=2.0e-5;


%beta_dev=4.8;  % Rates for subculture
%kappa_dev=0.4;


if (t< T0)
    par.betabar=beta_dev;
    par.kappa0=kappa_dev;
else
    par.betabar=beta_HSC;
    par.kappa0=kappa_HSC;
end

% Mutation rates for cancer development
p0=1e-6;
par.p01=0.5*p0;
par.p02=0.5*p0;
par.p13=p0;
par.p23=p0;
end