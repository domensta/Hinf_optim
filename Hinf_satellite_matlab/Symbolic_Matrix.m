kd = sym('kd');
ki = sym('ki');
kp = sym('kp');
T = sym('T');
p = sym('p');
but = p^4 + 0.008568*p^2 + 1.835e-05/(p^4 + 0.007071*p^3 + 0.008593*p^2 + 3.029e-05*p + 1.835e-05);
pidf=kp+ki/p+kd*p/(1+T*p);
pidtot = but*pidf;
[n,d] = numden(pidtot);
b = coeffs(n,p);
a = coeffs(d,p);
nbstate = length(a)-1;

for i=0:nbstate-1
    Actrl(nbstate,i+1) = -a(nbstate-i);
end
Bctrl = zeros(nbstate,1);
Bctrl(end) = 1;
for i=0:nbstate-1
    Cctrl(1,i+1) = b(nbstate-i)-a(nbstate-i)*b(1);
end
Dctrl = b(1);




load('Ass');
Aev =subs(Ass,[kp,ki,kd,T],[-9.95,9.45,9.57,1.006]);
cpol_ev = charpoly(Aev);
zeta = roots(cpol_ev);
center = 1; % where to locate eigenvalue;
nbeigen = 1; % how many eigenvalue
[delta,index] = sort(abs(zeta-center)) 
c = mean(zeta(index(1:nbeigen))) % best center of disc
c = 1;

cvec = [c^10,c^9,c^8,c^7,c^6,c^5,c^4,c^3,c^2,c,1];
pol_evatc = cpol_ev*cvec'

prod=1;
for i=nbeigen+1:length(index)
    prod = prod*(c-zeta(index(i)));
end