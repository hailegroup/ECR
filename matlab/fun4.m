function [d] =fun4(t,S,x,t_offset)
count=size(t);
S_est =  1 - (2 * (tan(x(1)))^2 / (x(1)^2 + x(1)^2 * (tan(x(1)))^2 + x(1)*tan(x(1))) * exp( - (t-t_offset).*60*x(5)*x(1)^2) + 2*(tan(x(2)))^2 / (x(2)^2 + x(2)^2 * (tan(x(2)))^2 + x(2)*tan(x(2))) * exp( -(t-t_offset).*60*x(5)*x(2)^2)+ 2*(tan(x(3)))^2 / (x(3)^2 + x(3)^2 * (tan(x(3)))^2 + x(3)*tan(x(3))) * exp( -(t-t_offset).*60*x(5)*x(3)^2)+2*(tan(x(4)))^2 / (x(4)^2 + x(4)^2 * (tan(x(4)))^2 + x(4)*tan(x(4))) * exp( -(t-t_offset).*60*x(5)*x(4)^2)  ) ;
d = sqrt( sum ((S_est -S).^2)/count(1) );


