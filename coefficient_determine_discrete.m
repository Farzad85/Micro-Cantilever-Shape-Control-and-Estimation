function [M,L,D,S,E,G,J,P,K,B]=coefficient_determine_discrete(xa,n)
b=30e-6;
l=250e-6;
d0=3e-6;
ep=8.854187817e-12;
h=2e-6;
ro=2400;

E_modul=120e9;
Ib=1/12*b*h^3;


[space_part,d_space_pare,dd_space_part,ddd_space_part...
    ,dddd_space_part]=mode_function(l,1,l);


M=ro*b*h*space_part.^2;
L=-0.5*ep*b*space_part/d0^2;
B=-ep*b*space_part.^2/d0^3;
D=-1.5*ep*b*space_part.^3/d0^4;
E=-2.5*ep*b*space_part.^4/d0^5;
G=-2.5*ep*b*space_part.^5/d0^6;
J=-3*ep*b*space_part.^6/d0^7;
P=-3.5*ep*b*space_part.^7/d0^8;
S=0;
K=E_modul*Ib*space_part.*dddd_space_part;
