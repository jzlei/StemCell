%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This code contains the program for the dynamics of cancer development due 
% to mutations to the driver genes related to the proliferation rate and
% the differentiation rate.
% Copyright: Jinzhi Lei
% Contact: jzlei@tsinghua.edu.cn
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function main()
global par              % The structure variable for model parameters
global T dt Td fname dx ntpr fmd     % Control variables. 
global gi g0 f0
global t                % The current time

t = 0;
Control();      % Set the simulation controls;
parameter();    % Set the model parameters

x=0:dx:1;
n=size(x,2);
Td=round(par.tau/dt)+1;
Qdelay1=zeros(Td,n);     % The matrix to store the delay states. 
Qdelay2=zeros(Td,n);
Qdelay3=zeros(Td,n);
Qdelay4=zeros(Td,n);

% Pre-define the valuables used for simulation
Q=Initialization(x);    % Initialize the system.
Qsum=zeros(4,1);
for i=1:4
    Qsum(i)=sum(Q(i,:));
end
for i=1:Td
    Qdelay1(i,1:n)=Q(1,1:n);   % Initialized the delays
    Qdelay2(i,1:n)=Q(2,1:n);
    Qdelay3(i,1:n)=Q(3,1:n);
    Qdelay4(i,1:n)=Q(4,1:n);
end
tau_index=floor(tau(x)/dt); % The index for the delays 
fkappa=zeros(4,n);
fkappa(1,:) = kappa(0,x);     % wile type
fkappa(2,:) = fkappa(1,:);  % mutation to the proliferation rate
fkappa(3,:) = 0.1 * fkappa(1,:);    % mutation to the differentiation rate
fkappa(4,:) = 0.1 * fkappa(1,:);
b0=zeros(4,n);
b0(1,:)=beta0(x);           % wilde type
b0(3,:) = b0(1,:);        % mutation to the differentiation rate
b0(2,:) = 2 * b0(1,:);    % mutation to the proliferation rate
b0(4,:) = 2 * b0(1,:);

gamma = mu(x).*tau(x);
p0 = p(x,x);
fp=fopen(char(fname),'w');
fmdip=fopen(char(fmd),'w');
OutPut(fp,Q,t,n);
fprintf(fmdip,'%f %f %f %f %f\n',t, Qsum(1), Qsum(2), Qsum(3), Qsum(4));

% Perform the iteration
step=0;
for t=dt:dt:T
    for i=1:4
        Qsum(i)=Integrate(Q(i,:),x);
    end
    c=sum(Qsum);
    
    f0=zeros(4,n);
    gi=zeros(4,n);
    g0=zeros(4,n);
    for i=1:4
        fbeta = beta(c, b0(i,:));
        f0(i,:) = -1.0 * Q(i,:).*(fbeta + fkappa(i,:));
    end
    
    Qtau=zeros(4,n);
    Qtau1=zeros(4,1);
    
    Qtau(1,:)=FindDelay(Qdelay1,tau_index,n);  % Find Qtau = Q(t-tau)
    Qtau(2,:)=FindDelay(Qdelay2,tau_index,n);
    Qtau(3,:)=FindDelay(Qdelay3,tau_index,n);
    Qtau(4,:)=FindDelay(Qdelay4,tau_index,n);
    for i=1:4
        ctau1(i)=Integrate(Qtau(i,:),x);   
    end
    ctaus1=sum(ctau1);
    for i=1:4
        fbeta = beta(ctaus1, b0(i,:));
        gi(i,:) = 2*(fbeta.*Qtau(i,:)).*exp(-gamma);
        g0(i,:) = gi(i,:) * p0 * dx;
    end
    
    Q(1,:)=Q(1,:) + dt * (f0(1,:) + (1-par.p12-par.p13)*g0(1,:));
    Q(2,:)=Q(2,:) + dt * (f0(2,:) + (1-par.p24) * g0(2,:) + par.p12 * g0(1,:));
    Q(3,:)=Q(3,:) + dt * (f0(3,:) + (1-par.p34) * g0(3,:) + par.p13 * g0(1,:));
    Q(4,:)=Q(4,:) + dt * (f0(4,:) + g0(4,:) + (par.p24 * g0(2,:) + par.p34 * g0(3,:)));
    
    % Update the Qdelay
    for i = 1:Td-1
        Qdelay1(i,1:n) = Qdelay1(i+1,1:n);
        Qdelay2(i,1:n) = Qdelay2(i+1,1:n);
        Qdelay3(i,1:n) = Qdelay3(i+1,1:n);
        Qdelay4(i,1:n) = Qdelay4(i+1,1:n); 
    end
    Qdelay1(Td,1:n) = Q(1,1:n);
    Qdelay2(Td,1:n) = Q(2,1:n);
    Qdelay3(Td,1:n) = Q(3,1:n);
    Qdelay4(Td,1:n) = Q(4,1:n);
    
    
    step=step + 1;
    % Output the results
    if(mod(step,ntpr)==0)
        OutPut(fp,Q,t,n);
    end
    if(mod(step,40)==0)
        fprintf(fmdip,'%f %f %f %f %f\n',t, sum(Q(1,:)), sum(Q(2,:)), sum(Q(3,:)), sum(Q(4,:)));
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

for j=1:4
    fprintf(fp,'%f ',t);
    for i=1:n
        fprintf(fp,'%f ',Q(j,i));
    end
    fprintf(fp,'\n');
end
end