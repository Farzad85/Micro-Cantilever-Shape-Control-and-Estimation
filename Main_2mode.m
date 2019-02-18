% Model of a micro cantilever subjected to the electrosatatic actuation. 
% The solver is fixed step Runge Kutta 4.
% States are estimated by an Extended Kalman Filter
% The controller is a state feedbal linearization controller. 

%%
clc
clear 
close all
tf=0.1;
N=100000;
t=linspace(0,tf,N);
y=zeros(4,N);
y(1,1)=0e-6;
yr=zeros(2,N);


% [T,Y]=ode45(@(t,x)beam_time_ode(t,x,4),[0 .001],[1e-6 0],odeset('Reltol',1e5*eps,'Abstol',1e5*eps));
h=t(2)-t(1);
% n=1;
[M,L,D,S,E,G,J,P,K,B]=coefficient_determine(1);
%% ODE solver loop
for i=2:1:N
        [a,yrr(i)] = feval(@beam_time_ode_2mode, t(i-1)    , y(:,i-1)       ,M,L,D,S,E,G,J,P,K,B,yr(:,i-1));
        
        yr(:,i)=yr(:,i-1)+h*feval(@reference_signal_,t(i-1),yr(:,i-1));
        k1 = feval(@beam_time_ode_2mode, t(i-1)    , y(:,i-1)       ,M,L,D,S,E,G,J,P,K,B,yr(:,i-1));
        k2 = feval(@beam_time_ode_2mode, t(i-1)+h/2, y(:,i-1)+k1*h/2,M,L,D,S,E,G,J,P,K,B,yr(:,i-1));
        k3 = feval(@beam_time_ode_2mode, t(i-1)+h/2, y(:,i-1)+k2*h/2,M,L,D,S,E,G,J,P,K,B,yr(:,i-1));
        k4 = feval(@beam_time_ode_2mode, t(i-1)+h  , y(:,i-1)+k3*h  ,M,L,D,S,E,G,J,P,K,B,yr(:,i-1));
        y(:,i) = y(:,i-1) + (k1 + 2*k2 + 2*k3 + k4)*h/6;
end
% plot(T,Y(:,1))


% plot(t,-z)
% hold on
% plot(t,yrr);

%% Post Processing
x=linspace(0,250e-6,100);
z=0;
for i=1:2
    z=z+mode_function(x(end),i,x(end)).*y(i,:);
end

u_cont=zeros(size(t));

for i=1:length(t)
[ydot_tmp,xr_tmp,u_cont(i)]=beam_time_ode_2mode(t(i),y(:,i),M,L,D,S,E,G,J,P,K,B);
end
%% Plots
close all
figure;
subplot(2,1,1)
plot(t,y(1,:))
ylabel('u_1','FontSize',12,'FontWeight','Bold')
subplot(2,1,2)
plot(t,y(2,:))
ylabel('u_2','FontSize',12,'FontWeight','Bold')
xlabel('Time(s)','FontSize',12,'FontWeight','Bold')

figure;
subplot(2,1,1)
plot(t,-z)
hold on 
plot(t,yrr,'r')
ylabel('z','FontSize',12,'FontWeight','Bold')
legend('z','z_d')
subplot(2,1,2)
plot(t,u_cont)
ylabel('Control Signal','FontSize',12,'FontWeight','Bold')
xlabel('Time(s)','FontSize',12,'FontWeight','Bold')

% 3D plots
figure;
t2=linspace(0,tf,N/100);
[T,X]=meshgrid(t2,x);
Z=zeros(size(T));
for i=1:length(t2)
    for j=1:length(x)
        for nn=1:2
            Z(j,i)=Z(j,i)+mode_function(x(j),nn,x(end)).*y(nn,100*i);
        end
    end
    plot(x,Z(:,i),'linewidth',2);
    ylim([-2 2]*1e-6)    
%     mov(i)=getframe;
end
surf(T,X,Z)
ylabel('x','FontSize',12,'FontWeight','Bold')
xlabel('Time(s)','FontSize',12,'FontWeight','Bold')
zlabel('z(x,t)','FontSize',12,'FontWeight','Bold')

% shading interp
% movie2avi(mov, 'microbeam.avi', 'compression', 'None');
% movie(mov,10)
