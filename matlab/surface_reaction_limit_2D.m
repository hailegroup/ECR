function [param]=surface_reaction_limit_2D(datafile,K,t_offset,thickness,algo)
data=load(datafile);
xdata=data(1:1:end,1)-data(1,1);
ydata=data(1:1:end,2);
thickness2=0.55;
%thickness2=0.493;
for j=1:4
x0=[K,t_offset] ;
lb=[10^-10,0];    
ub=[1,20];
c = @(x)[];
ceq = @(x)[];
nonlinfcn = @(x)deal(c(x),ceq(x));

if algo=='a'
opts = optimset('TolFun',1e-8,'Algorithm','active-set');
end
if algo=='i'
opts = optimset('TolFun',1e-8,'Algorithm','interior-point');
end
if algo=='s'
opts = optimset('TolFun',1e-8,'Algorithm','sqp');
end

[z,fval,exitflag,output,lambda] = fmincon(@(x) fun_surface_2D(xdata,ydata,x,thickness,thickness2),x0,[],[],[],[],lb,ub,nonlinfcn, opts)

yest= 1 -   exp(-(xdata-z(2)).*60*z(1)*2/thickness).*exp(-(xdata-z(2)).*60*z(1)*2/thickness2);

K1(j)=z(1);
t0(j)=z(2);
residual(j) = sqrt( sum ((yest-ydata).^2)/size(xdata,1) );
K=K1(j);
t_offset =t0(j);
end

K1'
residual'
t0'
%save('K_only','K1','t0','residual','-ascii');
close all 
plot(xdata-t_offset,ydata,'bx');
hold on
plot(xdata-t_offset,yest,'k-','LineWidth',2);
%title('  T=850 C.  pO_2  7.76x10^{-19} atm   ','FontWeight','bold','FontSize',18)
ylabel('  Normalized conductivity  ','FontSize',24)
xlabel('  Time(min) ','FontSize',24)
set(gca,'FontSize',24);

legend1=sprintf('%s','actual');
legend2=sprintf('%s','fit');
leg1=legend(legend1,legend2) ;
set(leg1,'Location','SouthEast');



