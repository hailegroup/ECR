%% Script to generate normalized conductivity vs time data from raw file.
%% Input  :  starting and ending line numbers, initial and final conductivity
%% Output :  matrix called 'writeme' with col1:time in minutes and col2:normalized conductivity

%display('Pass in the matrix of raw data');
%data=input('Pass in the matrix of raw data : ');
data=a;
%cond_col=input('Which column is conductivity ? ');
cond_col=5;
line1=input('Enter starting line number : ');
cond1=input('Enter initial equilibrated conductivity : ');
line2=input('Enter ending line number : ');
t=data(line1:line2,1);
t=t-t(1);
cond2=input('Enter final equilibrated conductivity : ');
s=(data(line1:line2,cond_col)-cond1)/(cond2-cond1);
writeme=t;
writeme(:,2)=s;
