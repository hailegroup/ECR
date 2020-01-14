function [d] =fun1(t,S,x,t_offset)
count=size(t);       
S_est = 1 - (2 * (tan(x(1)))^2 / (x(1)^2 + x(1)^2 * (tan(x(1)))^2 + x(1)*tan(x(1))) * exp( - (xdata-t_offset).*60*x(2)*x(1)^2));
d = sqrt( sum ((S_est -S).^2)/count(1) );


