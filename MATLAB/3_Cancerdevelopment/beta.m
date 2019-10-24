function y=beta(Q,b0)
global par
y=b0 * par.theta/(par.theta + Q);
end

