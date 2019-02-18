function [phi,d_phi,dd_phi,ddd_phi,dddd_phi]=mode_function(x,n,L)


switch n
    case 1
        kn=1.875104069/L;
    case 2
        kn=4.694091133/L;
    case 3
        kn=7.854757438/L;
    case 4
        kn=10.99554073/L;
    case 5
        kn=14.13716839/L;
    case 6
        kn=17.27875953/L;
end
        
a=0.5;        
phi=a*((cos(kn*x)-cosh(kn*x))+(-cos(kn*L)-cosh(kn*L))/(sin(kn*L)+sinh(kn*L))*...
    (sin(kn*x)-sinh(kn*x)));

d_phi=a*kn*((-sin(kn*x)-sinh(kn*x))+(-cos(kn*L)-cosh(kn*L))/(sin(kn*L)+sinh(kn*L))*...
    (cos(kn*x)-cosh(kn*x)));
dd_phi=a*kn^2*((-cos(kn*x)-cosh(kn*x))+(-cos(kn*L)-cosh(kn*L))/(sin(kn*L)+sinh(kn*L))*...
    (-sin(kn*x)-sinh(kn*x)));
ddd_phi=a*kn^3*((sin(kn*x)-sinh(kn*x))+(-cos(kn*L)-cosh(kn*L))/(sin(kn*L)+sinh(kn*L))*...
    (-cos(kn*x)-cosh(kn*x)));
dddd_phi=a*kn^4*((cos(kn*x)-cosh(kn*x))+(-cos(kn*L)-cosh(kn*L))/(sin(kn*L)+sinh(kn*L))*...
    (sin(kn*x)-sinh(kn*x)));

