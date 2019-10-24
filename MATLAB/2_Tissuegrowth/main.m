%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This code contains the program for tissue growth including the terminally 
% differentiation cells analogous to a cell lineages dynamics ?Figure 6 in the paper).
% Copyright: Jinzhi Lei
% Contact: jzlei@tsinghua.edu.cn
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function main()
global par              % The structure variable for model parameters
global T dt Td fname dx ntpr T0 fmd     % Control variables. 
global gi g0 f0
global t                % The current time

t = 0;
Control();      % Set the simulation controls;
parameter();    % Set the model parameters

x=0:dx:1;
x0=0.7;         % The boundary between stem cell and progenitor cells.
J0=find(x<x0);
J1=find(x>=x0);
n=size(x,2);
Td=round(par.tau/dt)+1;
Qdelay=zeros(Td,n);     % The matrix to store the delay states. 

% Pre-define the valuables used for simulation
Q=Initialization(x);    % Initialize the system.
Q(J0)=zeros(1,size(J0,2));
P=zeros(1,n);           % Initialize the terminally differentiated cell with 0 cells.
Qsum=sum(Q);            % Total number
Q0=sum(Q(J0));          % The total number of progenctor cells
Q1=sum(Q(J1));          % The total numbe of stem cells.
Psum=sum(P);            % The total number of terminallary differentiated cells
for i=1:Td
    Qdelay(i,1:n)=Q(1:n);   % Initialized the delays
end
tau_index=floor(tau(x)/dt); % The index for the delays 

% Here, we predefine the rate constant that are inddpendent to the time t

b0=beta0(x);
gamma = mu(x).*tau(x);  % gamma = mu * tau
nu = fnu(x);
p0 = p(x,x);

fp=fopen(char(fname),'w');
fmdip=fopen(char(fmd),'w');

OutPut(fp,Q,P,t,n);
fprintf(fmdip,'%f %f %f %f %f\n',t,Q0, Q1, Qsum, Psum);


% Perform the iteration
step=0;
for t=dt:dt:T
    % Calculate the differentiation rate.
    fkappa = kappa(t,x);
    
    c=Integrate(Q,x);   % Calculate the integral of Q(t, x) over all x from 0 to 1.
    fbeta = beta(c,b0);
    f0 = -1.0 * Q.*(fbeta + fkappa); % The coefficient -1.0*(beta(c,x) + kappa(x))* Q
    
    Qtau=FindDelay(Qdelay,tau_index,n);  % Find Qtau = Q(t-tau)
    ctau=Integrate(Qtau,x);             % The integral c_\tau
    gi = 2*(beta(ctau,b0).*Qtau).*exp(-gamma);  % The coefficient 2*exp(-mu*tau) * beta(c_tau, y) * Q(t-tau,y)
    g0 = gi * p0 * dx;
    
    Q = Q + (f0+g0)*dt; 
    P = P + (Q.*fkappa - nu.*P)*dt;  % Here P for T(t,x) in the equaiton
    
    % Update the Qdelay
    for i = 1:Td-1
        Qdelay(i,1:n) = Qdelay(i+1,1:n);
    end
    Qdelay(Td,1:n) = Q(1:n);
    
    step=step + 1;
    % Output the results
    if(mod(step,ntpr)==0)
        OutPut(fp,Q,P,t,n);
    end

    if(mod(step,ntpr)==4)   
        Qsum=sum(Q);            % Total number
        Q0=sum(Q(J0));          % The total number of progenctor cells
        Q1=sum(Q(J1));          % The total numbe of stem cells.
        Psum=sum(P);            % The total number of terminallary differentiated cells 
        fprintf(fmdip,'%f %f %f %f %f\n',t,Q0, Q1, Qsum, Psum);
    end

end
fclose(fp);
fclose(fmdip);
end

% The function to calculate the integration Integrate[Q(t,x)*xi(x),{x,0,1}]

function Q1=Integrate(Q,x)
Q1=sum(Q.*xi(x));
end

% The function to find the delay of Qdelay with delay tau
function Qt=FindDelay(Qdelay, tau,n)
global Td
Qt=zeros(1,n);
for i=1:n
    Qt(i)=Qdelay(Td-tau(i),i);
end
end

% This function output the simulation data

function OutPut(fp,Q,P,t,n)

fprintf(fp,'%f ',t);
for i=1:n
    fprintf(fp,'%f ',Q(i));
end
fprintf(fp,'\n');
fprintf(fp,'%f ',t);
for i=1:n
    fprintf(fp,'%f ',P(i));
end
fprintf(fp,'\n');

end
