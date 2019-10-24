%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This code contains the basic program for the model of heterogeneous stem
% cell regeneration. Figures 2, 3 in the paper.
% Copyright: Jinzhi Lei
% Contact: jzlei@tsinghua.edu.cn
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function main()
global par              % The structure variable for model parameters
global T dt Td fname dx ntpr T0 fmd    % Control variables. 
global gi g0 f0
global t                % The current time
global Flag

Flag = 2;   % Set Flag = 1 for Figure 2 (B-C)
            %     Flag = 2 for Figure 2 (D-E)
            %     Flag = 3 for Figure 2 (F-G)
            %     Flag = 4 for Figure 3
            % See the function parameter(Flag) for details. 
t = 0;
Control();      % Set the simulation controls;
parameter();    % Set the model parameters

x=0:dx:1;
n=size(x,2);
Td=round(par.tau/dt)+1;
Qdelay=zeros(Td,n);     % The matrix to store the delay states. 

% Pre-define the valuables used for simulation
Q=Initialization(Flag,x);    % Initialize the system.
Qsum=sum(Q);            % Total number
for i=1:Td
    Qdelay(i,1:n)=Q(1:n);   % Initialized the delays
end
tau_index=floor(tau(x)/dt); % The index for the delays 
gamma = mu(x).*tau(x);
b0 = beta0(x);
p0 = p(x,x);
fp=fopen(char(fname),'w');
fmdip=fopen(char(fmd),'w');
OutPut(fp,Q,t,n);
fprintf(fmdip,'%f %f\n',t, Qsum);

% Perform the iteration
step=1;
for t=dt:dt:T
    % Calculate the differentiation rate.
    fkappa = kappa(t,x);

    c=Integrate(Q,x);   % Calculate the integral of Q(t,x) over all x from 0 to 1.
    fbeta = beta(c,b0);
    f0 = -1.0 * Q.*(fbeta + fkappa); % The coefficient -1.0*(beta(c,x) + kappa(x))*Q
    
    Qtau=FindDelay(Qdelay,tau_index,n);  % Find Qtau = Q(t-tau)
    ctau=Integrate(Qtau,x);             % The integral c_\tau
    gi = 2*(beta(ctau,b0).*Qtau).*exp(-gamma); % The coefficient 2*exp(-mu*tau) * beta(c_tau,y) * Q(t-tau, y)
    g0 = gi * p0 * dx;
    
    Q = Q + (f0+g0)*dt; 
    
    % Update the Qdelay
    for i = 1:Td-1
        Qdelay(i,1:n) = Qdelay(i+1,1:n);
    end
    Qdelay(Td,1:n) = Q(1:n);
    
    step=step + 1;
    % Output the results
    if(mod(step,ntpr)==0)
        OutPut(fp,Q,t,n);
    end
    
    if(mod(step,ntpr)==4)   
        Qsum=sum(Q);            % Total number
        fprintf(fmdip,'%f %f\n',t, Qsum);
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

function OutPut(fp,Q,t,n)

fprintf(fp,'%f ',t);
for i=1:n
    fprintf(fp,'%f ',Q(i));
end
fprintf(fp,'\n');
end