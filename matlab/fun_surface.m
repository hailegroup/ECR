function [d] =fun_surface(t,S,x)
count=size(t,1);
S_est =  1 - exp( - (t-x(2)).*60*x(1)) ;
d = sqrt(sum ((S_est -S).^2)/count);


