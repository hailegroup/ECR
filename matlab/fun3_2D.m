function [d] =fun3_2D(t,S,x,t_offset,thickness,thickness2)
%t=data(:,1);
%S=data(:,2);
count=size(t);
S_est= 1 - ((2 * (tan(x(1)))^2 / (x(1)^2 + x(1)^2 * (tan(x(1)))^2 + x(1)*tan(x(1))) * exp( - (t-t_offset).*60*x(4)*4/thickness^2*x(1)^2))+( 2 * (tan(x(2)))^2 / (x(2)^2 + x(2)^2 * (tan(x(2)))^2 + x(2)*tan(x(2))) * exp( - (t-t_offset).*60*x(4)*4/thickness^2*x(2)^2))+( 2 * (tan(x(3)))^2 / (x(3)^2 + x(3)^2 * (tan(x(3)))^2 + x(3)*tan(x(3))) * exp( - (t-t_offset).*60*x(4)*4/thickness^2*x(3)^2))).*((2*(tan(x(5)))^2 / (x(5)^2 + x(5)^2 * (tan(x(5)))^2 + x(5)*tan(x(5))) * exp( - (t-t_offset).*60*x(4)*4/thickness2^2*x(5)^2))+( 2 * (tan(x(6)))^2 / (x(6)^2 + x(6)^2 * (tan(x(6)))^2 + x(6)*tan(x(6))) * exp( - (t-t_offset).*60*x(4)*4/thickness2^2*x(6)^2))+( 2 * (tan(x(7)))^2 / (x(7)^2 + x(7)^2 * (tan(x(7)))^2 + x(7)*tan(x(7))) * exp( - (t-t_offset).*60*x(4)*4/thickness2^2*x(7)^2))); 
d = sqrt( sum ((S_est -S).^2)/count(1) );


