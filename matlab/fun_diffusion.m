function [d] =fun_diffusion(t,S,x)
count=size(t);
S_est= 1 - (8/(pi^2))*exp(-(t-x(2)).*60*x(1)*pi^2/4)- (8/(9*pi^2))*exp(-(t-x(2)).*60*x(1)*9*pi^2/4)- (8/(25*pi^2))*exp(-(t-x(2)).*60*x(1)*25*pi^2/4);
d = sqrt( sum ((S_est -S).^2)/count(1) );

