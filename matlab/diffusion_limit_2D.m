function [param]=diffusion_limit_2D(datafile,D,t_offset,thickness,algo)
data=load (datafile);
xdata=data(1:end,1);
ydata=data(1:end,2);
thickness2=.55;

x0=[D,t_offset] ;
lb=[0,0]; 
ub=[1,80];
for j=1:5

c = @(x)[];
ceq = @(x)[];
nonlinfcn = @(x)deal(c(x),ceq(x));
if algo=='s'
	opts = optimset('Algorithm','sqp','Display','on');
end
if algo=='i'
	opts = optimset('TolFun',1e-12,'Algorithm','interior-point','Display','on');
end
if algo=='a'
	opts = optimset('Algorithm','active-set','Display','on');
	end
[z,fval,exitflag,output,lambda] = fmincon(@(x) fun_diffusion_2D(xdata,ydata,thickness,thickness2,x),x0,[],[],[],[],lb,ub,nonlinfcn,opts)
 yest= 1 - ((8/(pi^2))*exp(-(xdata-z(2)).*60*z(1)*pi^2/thickness^2)+ (8/(9*pi^2))*exp(-(xdata-z(2)).*60*z(1)*9*pi^2/thickness^2)+ (8/(25*pi^2))*exp(-(xdata-z(2)).*60*z(1)*25*pi^2/thickness^2)).*((8/(pi^2))*exp(-(xdata-z(2)).*60*z(1)*pi^2/thickness2^2)+ (8/(9*pi^2))*exp(-(xdata-z(2)).*60*z(1)*9*pi^2/thickness2^2)+ (8/(25*pi^2))*exp(-(xdata-z(2)).*60*z(1)*25*pi^2/thickness2^2)) ;
Do(j)=z(1);
t0(j)=z(2);                                                 
D=Do(j);
t_offset=t0(j);
residual(j) = sqrt( sum ((yest-ydata).^2)/size(xdata,1) )

end
Do'
residual'
t0'
%save('D_only','Do','residual','-ascii');
close all 
plot(xdata,ydata,'ko');
hold on
plot(xdata,yest,'k-','Linewidth',2);
%title('  T=750 C. \Delta pO_2  2.6x10^{-21} atm to 3.5x10^{-21} atm  ','FontWeight','bold','FontSize',13)
ylabel('  Normalized conductivity  ','FontSize',24)
xlabel('  Time(min) ','FontSize',24)            
legend1=sprintf('%s','actual');   
legend2=sprintf('%s','fit');
leg1=legend(legend1,legend2) ;
set(leg1,'Location','SouthEast','FontSize',24);
set(gca,'FontSize',24);
