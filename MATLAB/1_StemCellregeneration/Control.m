%
% The function to define the control parameters
%

function Control()
global T dt dx fname ntpr T0 fmd
T0=-1;        % The timing for system trainning. 
T=240000;    % The timing for total simulation
dt=0.25;
fname='output/md.dat';   % The fine name to output the function Q(t, x)
fmd='output/md-ntpx.dat';    % The file name to output the dynamics 
dx=0.01;
ntpr=100;                       % Steps to write the outpu files.
end