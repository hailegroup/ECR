function generate_data_2D(K,D,thickness,totaltime_hours,t0_minutes,noise_amp)
format shortEng;
token=1; i=1;
L = thickness/2*K/D
%thickness2=.55;
thickness2=.493;
L2 = thickness2/2 *K/D

n=4;
npts=1200;
weight=totaltime_hours*3600.00/npts;
alpha=get_eigenval(L,n);
alpha2=get_eigenval(L2,n);
filename=sprintf('%2.2ecm_k%2.2e_D%2.2e_1D',thickness,K,D);
filename2=sprintf('%2.2ecm_k%2.2e_D%2.2e',thickness,K,D);
f=fopen(filename,'w');
f2=fopen(filename2,'w');
i=1;     
while i< npts
		c(i)=0;
		c2(i)=0;
        
	 t(i)=(i-1)*weight;
	for j=1:n    
		c(i) = c(i) + (2*(tan(alpha(j)))^2 / (alpha(j)^2 + alpha(j)^2 * (tan(alpha(j)))^2 + alpha(j)*tan(alpha(j))) * exp( -(t(i)-t0_minutes*60)*D*4/thickness^2*alpha(j)^2));
		c2(i) = c2(i) + (2*(tan(alpha2(j)))^2 / (alpha2(j)^2 + alpha2(j)^2 * (tan(alpha2(j)))^2 + alpha2(j)*tan(alpha2(j))) * exp( -(t(i)-t0_minutes*60)*D*4/thickness2^2*alpha2(j)^2));
	 end
	    c2(i) = 1 + (-0.5+ rand(1,1))*noise_amp - c(i)*c2(i); 
       c(i) = 1 + (-0.5+ rand(1,1))*noise_amp - c(i); 
    if (c(i)<0)  
		c(i)=0;
	end
	fprintf(f,'%e\t%e\n',t(i)/60-t0_minutes,c(i));
	fprintf(f2,'%e\t%e\n',t(i)/60-t0_minutes,c2(i));
	i=i+1;
end
    plot(t(1:5:end)/60,c(1:5:end),'ro','MarkerSize',5,'MarkerFaceColor','r');
    hold on
    plot(t(1:5:end)/60,c2(1:5:end),'bo','MarkerSize',5,'MarkerFaceColor','b');
%    hold on
%	plot(t(1:5:end)/60,c(1:5:end).*c2(1:5:end),'ro','MarkerSize',5);
  	xlabel('Time in min','FontSize',16)
	ylabel('normalized conductivity','FontSize',16)
	name=sprintf('%s%3.2e%s%3.2e%s%s%f','K=',K,'cm/s  D=',D,'cm ^{  2}/s','  L=',L);
    title(name,'FontSize',18);
    set(gca,'FontSize',18)
	fclose(f);
	fclose(f2);
