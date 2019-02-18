clc;clear
m=4;
n=2;
j=2;
N=100;

syms x u1 u2 u3 u4 u5 u6 udd1 udd2 udd3 udd4 udd5 udd6 V real
u=[u1 u2 u3 u4 u5 u6];
udd=[udd1 udd2 udd3 udd4 udd5 udd6];
b=30e-6;
L=250e-6;
d0=3e-6;
ep=8.854187817e-12;
h=2e-6;
ro=2400;
E=120e9;
Ib=1/12*b*h^3;
dx=L/N;


M=vpa(ro*b*h*int(mode_function(x,j,L)^2,0,L));
K1=vpa(E*Ib*int(diff(mode_function(x,j,L),4)*mode_function(x,j,L),0,L));
z=0;
for i=1:n
    z=z+mode_function(x,i,L)*u(i);
end

C=0;

x_int=0:dx:L;
fun=sym(zeros(size(x_int)));

phi_j=mode_function(x,j,L);
for k=0:m
    mm=0;
    for x_val=0:dx:L
        mm=mm+1;
%         x_vals=sym(x_val);
%         
%         x_vals
        
        fun(mm)=subs(phi_j*z^k,x,x_val);
        if rem(mm,10)<0.1
            mm;
        end
    end
    
    C=C+simplify(vpa((k+1)/(d0^(k+2))*trapz(x_int,fun)));
    k;

end
    
udd(j)=-1/M*(K1*u(j))+V^2*vpa(1/M*b*ep/2*C);
udd=udd(j)

% pretty(DiffEq)
