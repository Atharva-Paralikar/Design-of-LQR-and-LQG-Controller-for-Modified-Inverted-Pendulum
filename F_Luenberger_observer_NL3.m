clc;

M=1000;
m1=100;
m2=100;
l1=20;
l2=10;
g=9.81;


x0 = [0;0;pi/6;0;pi/3;0];

time_span = 0:0.01:400;

% initial external force
force = 0;

A=[0 1 0 0 0 0; 
    0 0 -(m1*g)/M 0 -(m2*g)/M 0;
    0 0 0 1 0 0;
    0 0 -((M+m1)*g)/(M*l1) 0 -(m2*g)/(M*l1) 0;
    0 0 0 0 0 1;
    0 0 -(m1*g)/(M*l2) 0 -(g*(M+m2))/(M*l2) 0];

B=[0; 1/M; 0; 1/(M*l1); 0; 1/(M*l2)];

C = [1 0 0 0 0 0;
     0 0 1 0 0 0;
     0 0 0 0 1 0];

D = 0;
nlterm = 0;

%% Kalman Estimator Design
Bd = 0.1*eye(6);                % Process Disturbance
Vn = 0.01;                      % Measurement Noise
[Lue4,P,E] = lqe(A,Bd,C,Bd,Vn*eye(3));
% Running the state variables and updating them for 400 steps
[t1,x1,F1] = ode45(@NLObs4,time_span,x0,force);
x = x1(:, 1);
x_dot = x1(:, 2);
theta1 = x1(:, 3);
theta1dot = x1(:, 4);
theta2 = x1(:, 5);
theta2dot = x1(:, 6);

figure(5)
k = tiledlayout(3,2);
title(k,"Observer 4")
nexttile
plot(t1,x);
nexttile
plot(t1,x_dot);
nexttile
plot(t1,theta1);
nexttile
plot(t1,theta1dot);
nexttile
plot(t1,theta2);
nexttile
plot(t1,theta2dot);
grid on;

function dxdt = NLObs4(t,x,force)
M=1000;
m1=100;
m2=100;
l1=20;
l2=10;
g=10;

A=[0 1 0 0 0 0; 
    0 0 -(m1*g)/M 0 -(m2*g)/M 0;
    0 0 0 1 0 0;
    0 0 -((M+m1)*g)/(M*l1) 0 -(m2*g)/(M*l1) 0;
    0 0 0 0 0 1;
    0 0 -(m1*g)/(M*l2) 0 -(g*(M+m2))/(M*l2) 0];

B=[0; 1/M; 0; 1/(M*l1); 0; 1/(M*l2)];

% LQR Parameters
R = [0.0001];
Q=[0.05 0 0 0 0 0;
   0 1 0 0 0 0;
   0 0 1000 0 0 0;
   0 0 0 10 0 0;
   0 0 0 0 500 0;
   0 0 0 0 0 10];
[K, P, Poles] = lqr(A,B,Q,R);
Lue4 = [1.42853787432444 -0.0730849372050412 0.0167575382445034;
0.973171340756836 -0.316168695665875 -0.172164093106471;
-0.0730849372050378 0.463573037997319 -0.0450754880501592;
0.177128534238979 0.0611365846136433 -0.0382715173564502;
0.0167575382444999 -0.0450754880501546 0.386886721863738;
0.205880479420738 -0.00128809501661267 0.0259969751327212];
C = [1 0 0 0 0 0;
     0 0 1 0 0 0;
     0 0 0 0 1 0];

nlterm = Lue4*([x(1); x(3); x(5)] - C*x);
F =- K*x ;
dxdt=zeros(6,1);

dxdt(1) = x(2) + nlterm(1); 
dxdt(2)=(F-(g/2)*(m1*sin(2*x(3))+m2*sin(2*x(5)))-(m1*l1*(x(4)^2)*sin(x(3)))-(m2*l2*(x(6)^2)*sin(x(5))))/(M+m1*((sin(x(3)))^2)+m2*((sin(x(5)))^2))+ nlterm(2);
dxdt(3)= x(4)+ nlterm(3); 
dxdt(4)= (dxdt(2)*cos(x(3))-g*(sin(x(3))))/l1+ nlterm(4); 
dxdt(5)= x(6)+ nlterm(5); 
dxdt(6)= (dxdt(2)*cos(x(5))-g*(sin(x(5))))/l2+ nlterm(6); 
end