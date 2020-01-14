function eigenval = get_eigenval(L,n)
format long; 
for i=1:n 
   eigenval(i)=fsolve(@(v) v*tan(v)-L, (i-1)*pi+0.1);
end
