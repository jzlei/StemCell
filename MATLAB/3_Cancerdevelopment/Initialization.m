%
% The function to initialize the state.
%

function Q=Initialization(x)
Q=IniCancerDev(x);
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
