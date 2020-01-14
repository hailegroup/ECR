function [param]=fit_multiple_fmincon(t_offset,datafile,thickness,nexp,ntrials,algo)
data=load (datafile);
xdata=data(1:1:end,1)-data(1,1);
ydata=data(1:1:end,2);
%%%%%%%%%%%%%%% Constants %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%ntrials=8;   %(number of starting K for each L)
%nl=5;  %(number of L)
ndD=600;%(no of bins for D histogram) 
ndK=600;%(no of bins for D histogram) 
%n=input('Number of exponentials to fit to(max=5):');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

f=fopen('data_after_fit.txt','w');
fprintf(f,'%s\n','%%%%%%%%%%L(initial) K(initial) D(initial) L(fit)  K(fit) D(fit) residual exitflag%%%%%%%%%%%');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if algo=='s'
	opts = optimset('Algorithm','sqp','Display','off');
end
if algo=='i'
	opts = optimset('TolFun',1e-8,'Algorithm','interior-point','Display','off');
end			
if algo=='a'
	opts = optimset('Algorithm','active-set','Display','off');
end			

if nexp==2
	for itrials=1:ntrials
		K_ini =10^-7*10.^(4*rand(1,1));
		D_ini =10^-7*10.^(4*rand(1,1));
	   	Ktemp=K_ini;
		Dtemp=D_ini;
		Ltemp= thickness/2*K_ini/D_ini;
		alpha=get_eigenval(Ltemp,nexp);
		for j=1:3
			x0=[alpha(1),alpha(2),Dtemp*4/thickness^2] ;
			lb=[0,pi,10^(-8)]; 
			ub=[pi*0.5,1.5*pi,10^(2)];

			c = @(x)[];
			ceq = @(x) [x(1)*tan(x(1)) - x(2)*tan(x(2)) ];
			nonlinfcn = @(x)deal(c(x),ceq(x));


			[z,fval,exitflag,output,lambda] = fmincon(@(x) fun2(xdata,ydata,x,t_offset),x0,[],[],[],[],lb,ub,nonlinfcn, opts)
	        yest= 1 - (2 * (tan(z(1)))^2 / (z(1)^2 + z(1)^2 * (tan(z(1)))^2 + z(1)*tan(z(1))) * exp( - (xdata-t_offset).*60*z(3)*z(1)^2))-( 2 * (tan(z(2)))^2 / (z(2)^2 + z(2)^2 * (tan(z(2)))^2 + z(2)*tan(z(2))) * exp( - (xdata-t_offset).*60*z(3)*z(2)^2));
           %[gf]= gfit(ydata,yest,'7')
           alpha(1)=z(1);
           alpha(2)=z(2);
           Ltemp=z(1)*tan(z(1));
           Dtemp=z(3)* thickness^2/4;
           Ktemp =2*Dtemp*Ltemp/thickness;
           residual = sqrt( sum ((yest-ydata).^2)/size(xdata,1) )
		end
       fprintf(f,'%f\t%e\t%e\t%f\t%e\t%e\t%f\t%d\n',thickness/2*K_ini/D_ini,K_ini,D_ini,Ltemp,Ktemp,Dtemp,residual,exitflag);
	end   
end
	
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nexp==3
	for itrials=1:ntrials
		K_ini =10^-6*10.^(6*rand(1,1));
		D_ini =10^-7*10.^(6*rand(1,1));
	   	Ktemp=K_ini;
		Dtemp=D_ini;
		Ltemp= thickness/2*K_ini/D_ini;
		alpha=get_eigenval(Ltemp,nexp);
		for j=1:3
			x0=[alpha(1),alpha(2),alpha(3),Dtemp*4/thickness^2] ;
			lb=[0,pi,2*pi,1e-10]; 
			ub=[pi*0.5-0.000001,1.5*pi-0.000001,2.5*pi-0.000001,1e7];

			c = @(x)[];
			ceq = @(x) [x(1)*tan(x(1)) - x(2)*tan(x(2)), x(1)*tan(x(1)) - x(3)*tan(x(3)) ];
			nonlinfcn = @(x)deal(c(x),ceq(x));

			[z,fval,exitflag,output,lambda] = fmincon(@(x) fun3(xdata,ydata,x,t_offset),x0,[],[],[],[],lb,ub,nonlinfcn, opts)
			yest= 1 - (2 * (tan(z(1)))^2 / (z(1)^2 + z(1)^2 * (tan(z(1)))^2 + z(1)*tan(z(1))) * exp( - (xdata-t_offset).*60*z(4)*z(1)^2))-( 2 * (tan(z(2)))^2 / (z(2)^2 + z(2)^2 * (tan(z(2)))^2 + z(2)*tan(z(2))) * exp( - (xdata-t_offset).*60*z(4)*z(2)^2))-( 2 * (tan(z(3)))^2 / (z(3)^2 + z(3)^2 * (tan(z(3)))^2 + z(3)*tan(z(3))) * exp( - (xdata-t_offset).*60*z(4)*z(3)^2)) ;
			%[gf]= gfit(ydata,yest,'7')
			alpha(1)=z(1)
			alpha(2)=z(2)
			alpha(3)=z(3)
			Ltemp=z(1)*tan(z(1))
	      	Dtemp=z(4)* thickness^2/4
	       	Ktemp =2*Dtemp*Ltemp/thickness
			residual = sqrt( sum ((yest-ydata).^2)/size(xdata,1) )
		end
       fprintf(f,'%f\t%e\t%e\t%f\t%e\t%e\t%f\t%d\n',thickness/2*K_ini/D_ini,K_ini,D_ini,Ltemp,Ktemp,Dtemp,residual,exitflag);
	parameters(itrials,:)=[Ktemp Dtemp Ltemp residual exitflag];
	end   
	save('statistics_from_fit','parameters','-ascii' );

	
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nexp==4
for itrials=1:ntrials
K_ini =10^-9*10.^(6*rand(1,1))
D_ini =10^-10*10.^(6*rand(1,1))
Ktemp=K_ini;
Dtemp=D_ini;
Ltemp= thickness/2*K_ini/D_ini;
alpha=get_eigenval(Ltemp,nexp);
for j=1:3
    
x0=[alpha(1),alpha(2),alpha(3),alpha(4),Dtemp*4/thickness^2] ;
lb=[0,pi,2*pi,3*pi,10^(-14)];
ub=[pi*0.5-0.000001,1.5*pi-0.000001,2.5*pi-0.000001,3.5*pi-.000001,10^(0)];

c = @(x)[];
ceq = @(x) [x(1)*tan(x(1)) - x(2)*tan(x(2)), x(1)*tan(x(1)) - x(3)*tan(x(3)),x(1)*tan(x(1)) - x(4)*tan(x(4)) ];
nonlinfcn = @(x)deal(c(x),ceq(x));

[z,fval,exitflag,output,lambda] = fmincon(@(x) fun4(xdata,ydata,x,t_offset),x0,[],[],[],[],lb,ub,nonlinfcn, opts)
yest= 1 - (2 * (tan(z(1)))^2 / (z(1)^2 + z(1)^2 * (tan(z(1)))^2 + z(1)*tan(z(1))) * exp( - (xdata-t_offset).*60*z(5)*z(1)^2))-( 2 * (tan(z(2)))^2 / (z(2)^2 + z(2)^2 * (tan(z(2)))^2 + z(2)*tan(z(2))) * exp( - (xdata-t_offset).*60*z(5)*z(2)^2))-( 2 * (tan(z(3)))^2 / (z(3)^2 + z(3)^2 * (tan(z(3)))^2 + z(3)*tan(z(3))) * exp( - (xdata-t_offset).*60*z(5)*z(3)^2))-( 2 * (tan(z(4)))^2 / (z(4)^2 + z(4)^2 * (tan(z(4)))^2 + z(4)*tan(z(4))) * exp( - (xdata-t_offset).*60*z(5)*z(4)^2)) ;
%[gf]= gfit(ydata,yest,'7')
alpha(1)=z(1)
alpha(2)=z(2)
alpha(3)=z(3)
alpha(4)=z(4)
Ltemp=z(1)*tan(z(1))
Dtemp=z(5)* thickness^2/4
Ktemp =2*Dtemp*Ltemp/thickness
residual = sqrt( sum ((yest-ydata).^2)/size(xdata,1) )
end
fprintf(f,'%f\t%e\t%e\t%f\t%e\t%e\t%f\t%d\n',thickness/2*K_ini/D_ini,K_ini,D_ini,Ltemp,Ktemp,Dtemp,residual,exitflag);
parameters(itrials,:)=[Ktemp Dtemp Ltemp residual exitflag];
end
save('statistics_from_fit','parameters','-ascii' );

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fclose(f);

close all
%stat_and_plot(datafile,'statistics_from_fit',thickness);

toplot=load('data_after_fit.txt');
minD=min(toplot(:,6))
maxD=max(toplot(:,6))
minK=min(toplot(:,5))
maxK=max(toplot(:,5))
%dD=linspace(minD,maxD,ndD);
%dK=linspace(minK,maxK,ndK);
dD=logspace(log10(minD),log10(maxD),ndD);
dK=logspace(log10(minK),log10(maxK),ndK);
%loglog(5.45e-5,2.14e-6,'kv','MarkerSize',8,'MarkerFaceColor','k');
%hold on
loglog(toplot(:,2),toplot(:,3),'ko','MarkerSize',8);
hold on
loglog(toplot(:,5),toplot(:,6),'ko','MarkerSize',8,'MarkerFaceColor','k');
ylabel('D_{Chem} (cm^2/s)','FontSize',24)
xlabel('k_S (cm/s)','FontSize',24)
set(gca,'FontSize',24);
legend1=sprintf('%s','Initial guess');
legend2=sprintf('%s','Optimized');
%legend3=sprintf('%s','Input');
%leg1=legend(legend3,legend1,legend2) ;
leg1=legend(legend1,legend2) ;
set(leg1,'Location','Best');
hold off
pause
hist(toplot(:,5),dK,'FaceColor','k');
xlabel('k_S (cm/s)','FontSize',24)
ylabel('Count','FontSize',24)
pause
hist(toplot(:,6),dD,'FaceColor','k');
xlabel('D_{Chem} (cm^2/s)','FontSize',24)
ylabel('Count','FontSize',24)
pause
loglog(toplot(:,1),toplot(:,4),'bv','MarkerSize',6,'MarkerFaceColor','b' )
ylabel('L(final)','FontSize',24)
xlabel('L(initial)','FontSize',24)
pause
plot(toplot(:,4),toplot(:,7),'bv','MarkerSize',6,'MarkerFaceColor','b' );
xlabel('L(final)','FontSize',24)
ylabel('residual','FontSize',24)
 
nelem=histc(toplot(:,5),dK);
index=min(find(nelem>max(nelem)-1));
%Kmode=1/2*(dK(index)+dK(index+1));
Kmode=dK(index) ;

nelem=histc(toplot(:,6),dD);
index=min(find(nelem>max(nelem)-1));
%Dmode=1/2*(dK(index)+dK(index+1));
Dmode=dD(index)  ;            

t=xdata;
[r,c]=find(toplot == min(toplot(:,7)));
%%k_minerror= toplot(r,5)
%%D_minerror= toplot(r,6)
%%L_minerror=thickness/2*k_minerror/D_minerror
Lmode= thickness/2*Kmode/Dmode   ;
%Lmode= thickness/2*Kmode/Dmode   
z=get_eigenval(Lmode,nexp);
%z(nexp+1) = D_minerror*4/thickness^2;
z(nexp+1) = Dmode*4/thickness^2;
yest=1;
for j=1:nexp
	yest= yest - 1*((2 * (tan(z(j)))^2 / (z(j)^2 + z(j)^2 * (tan(z(j)))^2 + z(j)*tan(z(j))) * exp( - (t-t_offset).*60*z(nexp+1)*z(j)^2)));
end

plot(t,ydata,'bo');
hold on
plot(t,yest,'k-','Linewidth',2);
ylabel('  Normalized conductivity  ','FontWeight','bold','FontSize',24)
xlabel('  Time(min) ','FontWeight','bold','FontSize',24)
fit = [t,ydata,yest]
dlmwrite('fit.txt',fit);
legend1=sprintf('%s','actual');
legend2=sprintf('%s','fit');
leg1=legend(legend1,legend2,'FontSize',24) ;

Lmode
Kmode
Dmode 
