%
% This function set up the parameters used in model simulations.
% Different sets of parameters are used for different applications.
%

function parameter()
global par

par=struct;
par.theta=1e3;
par.a1=5.8;
par.a2=2.2;
par.a3=3.75;
par.b1=4.0;
par.mu=2.0e-4;
par.tau=20;
par.nu0=2.0e-5;
par.kappa0=0.02;
par.betabar=0.12;

end