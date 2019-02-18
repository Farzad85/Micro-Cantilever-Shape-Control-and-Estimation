clc;clear
%% Important Parameters
m=4;
n=2;
j=1;
N=100;
L=250e-6;
x0=L; % Imposing position of electrostatic excitation;

%% Parmaters

b=30e-6;
d0=3e-6;
ep=8.854187817e-12;
h=2e-6;
ro=2400;
E=120e9;
Ib=1/12*b*h^3;
dx=L/N;

%% Finding DE
syms x u1 u2 u3 u4 u5 u6 udd1 udd2 udd3 udd4 udd5 udd6 V real
u=[u1 u2 u3 u4 u5 u6];
udd=[udd1 udd2 udd3 udd4 udd5 udd6];

M=vpa(ro*b*h*mode_function(x0,j,L)^2);
K1=vpa(E*Ib*subs(diff(mode_function(x,j,L),4)*mode_function(x,j,L),x,x0));
z=0;
for i=1:n
    z=z+mode_function(x0,i,L)*u(i);
end

C=0;
phi_j=mode_function(x0,j,L);
for k=0:m    
    C=C+simplify(vpa((k+1)/(d0^(k+2))*phi_j*z^k));
    k;
end
    
% DiffEq=udd(j)+1/M*(K1*u(j))-V^2*vpa(1/M*b*ep/2*C);
% pretty(DiffEq)
udd(j)=-1/M*(K1*u(j))+V^2*vpa(1/M*b*ep/2*C);
udd=udd(j)
