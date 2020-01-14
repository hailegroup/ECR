function [d] =fun3(t,S,x,t_offset)
%t=data(:,1);
%S=data(:,2);
count=size(t);
S_est =  1 - (2 * (tan(x(1)))^2 / (x(1)^2 + x(1)^2 * (tan(x(1)))^2 + x(1)*tan(x(1))) * exp( - (t-t_offset).*60*x(4)*x(1)^2) + 2*(tan(x(2)))^2 / (x(2)^2 + x(2)^2 * (tan(x(2)))^2 + x(2)*tan(x(2))) * exp( -(t-t_offset).*60*x(4)*x(2)^2)+ 2*(tan(x(3)))^2 / (x(3)^2 + x(3)^2 * (tan(x(3)))^2 + x(3)*tan(x(3))) * exp( -(t-t_offset).*60*x(4)*x(3)^2)  ) ;
d = sqrt( sum ((S_est -S).^2)/count(1) );


