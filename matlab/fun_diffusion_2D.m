function [d] =fun_diffusion_2D(t,S,thickness,thickness2,x)
x
count=size(t);       
S_est= 1 - ((8/(pi^2))*exp(-(t-x(2)).*60*x(1)*pi^2/(thickness^2))+ (8/(9*pi^2))*exp(-(t-x(2)).*60*x(1)*9*pi^2/(thickness^2))+ (8/(25*pi^2))*exp(-(t-x(2)).*60*x(1)*25*pi^2/(thickness^2))).*((8/(pi^2))*exp(-(t-x(2)).*60*x(1)*pi^2/(thickness2^2))+ (8/(9*pi^2))*exp(-(t-x(2)).*60*x(1)*9*pi^2/(thickness2^2))+ (8/(25*pi^2))*exp(-(t-x(2)).*60*x(1)*25*pi^2/(thickness2^2))) ;
%S_est= 1 - ((8/(pi^2))*exp(-(t-x(2)).*60*x(1)*pi^2/(thickness^2)));
d = sqrt( sum ((S_est -S).^2)/count(1) );

