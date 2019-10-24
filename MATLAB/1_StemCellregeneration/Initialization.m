%
% The function to initialize the state.
%

function Q=Initialization(Flag, x)
switch Flag
    case 1      % Initialize the system from high stemness cells
        Q=IniEmbryo(x);        
    case 2
        Q=IniNormal(x); % Initialize the system from the stady state at normal growth
    case 3
        Q=IniNormal(x); 
    case 4
        Q=IniSubcolture(x); % Inititlize the system from the subcolture of steady state at normal growth
end
end

function Q=IniEmbryo(x) % Initialize the system as an early embryo development, so that most cells are naive stem cells.
n=size(x,2);
Q=zeros(1,n);

for i=floor(n/2):n
    Q(i)=randi(4);
end
end

function Q=IniNormal(x) % Initialize the system as the normal steady state, which calculated from the developmental progress.
A=load('output/md-dev.dat');
[m,n]=size(A);
Q=A(m,2:end);
end

function Q=IniSubcolture(x) % Initialize the system to mimi the sub-colture from normal cells
Q=IniNormal(x);
x0=0.75;
x1=1.0;
J=find(x<x0 | x>x1);
Q(J)=0;
end


