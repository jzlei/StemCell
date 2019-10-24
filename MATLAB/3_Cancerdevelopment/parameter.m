%
% The function to set up the model parameters
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
par.betabar=0.12;
par.kappa0=0.02;


% Mutation rates for cancer development
p0=1e-6;
par.p12=0.5*p0;
par.p13=0.5*p0;
par.p24=p0;
par.p34=p0;

end