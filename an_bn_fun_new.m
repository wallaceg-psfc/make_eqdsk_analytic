function [a,b] = an_bn_fun_new(nseries, lam, epsbar, k, a_in, b_in)

nterm = nseries+1;
a = NaN(nterm,1);
b = a;

a(1:3) = a_in;
b(1:3) = b_in;

p=lam^2-epsbar*k^2;

for nn=4:nterm
    n=nn-1;
    
    T1=epsbar*(n-1)*(n-2)*a(nn-1)+p*a(nn-3);
    T2=2*(n-1)*k*b(nn-1)+2*epsbar*(n-2)*k*b(nn-2);
    a(nn)=-(T1+T2)/n/(n-1);
    T3=epsbar*(n-1)*(n-2)*b(nn-1)+p*b(nn-3);
    T4=-2*(n-1)*k*a(nn-1)-2*epsbar*(n-2)*k*a(nn-2);
    b(nn)=-(T3+T4)/n/(n-1);

end


end