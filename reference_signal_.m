function dy=reference_signal_(t,y)
dy=zeros(2,1);
K1=600;
K2=90000;
a=1e-6;
b=.5e-6;
t1=0.03;
u=(a*heaviside(t)-b*heaviside(t-t1));
dy(1)=y(2);
dy(2)=K2*u-K1*y(2)-K2*y(1);