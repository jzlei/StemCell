%
% The function to initialize the state.
%

function Q=Initialization(x)
Q=IniEmbryo(x);
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

function Q=IniEven(x)   % Initialize the system with even distribution.
n=size(x,2);
Q=zeros(1,n);

for i=1:n
    Q(i)=randi(1000);
end
end

% Initialize the system for the program of cancer development. 
% The state Q include 4 mutant types:
% Q(1,:) for the wilde type cells
% Q(2,:)--Q(4,:) for the mutant type cells.
function Q=IniCancerDev(x)
A=load('output/md-dev.dat');
[m,n]=size(A);
Q=zeros(4,n-1);
Q(1,:)=A(m,2:end);
end
