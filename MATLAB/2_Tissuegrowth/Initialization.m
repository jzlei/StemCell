%
% The function to initialize the state.
% Here, we initialize the system as an early embryo development, so that most cells are naive stem cells.
%

function Q=Initialization(x)
n=size(x,2);
Q=zeros(1,n);

for i=floor(n/2):n
    Q(i)=randi(4);
end
end

