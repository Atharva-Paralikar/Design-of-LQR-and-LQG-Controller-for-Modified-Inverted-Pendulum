clc
clear


M=1000; %Mass of the cart
m1=100; %Mass of Pendulum 1
m2=100; %Mass of Pendulum 2
l1=20;  %Length of the string of Pendulum 1
l2=10;  %length of the string of Pendulum 2
g=9.81; %declaring the value of the accelertaion due to gravity in m/s^2

A=[0 1 0 0 0 0; 
    0 0 -(m1*g)/M 0 -(m2*g)/M 0;
    0 0 0 1 0 0;
    0 0 -((M+m1)*g)/(M*l1) 0 -(m2*g)/(M*l1) 0;
    0 0 0 0 0 1;
    0 0 -(m1*g)/(M*l2) 0 -(g*(M+m2))/(M*l2) 0];
B=[0; 1/M; 0; 1/(M*l1); 0; 1/(M*l2)];

if (rank(ctrb(A,B))==6)
    disp("System is controllable")
else
    disp("System is uncontrollable")
end

x_initial = [0;0;pi/6;0;pi/3;0]; 

Q=[0.05 0 0 0 0 0;
   0 1 0 0 0 0;
   0 0 1000 0 0 0;
   0 0 0 10 0 0;
   0 0 0 0 500 0;
   0 0 0 0 0 10];
R=0.00001;

C = eye(6);
D = 0; 

disp("Now, seeing the results using an LQR controller")
[K_val, P_mat, Poles] = lqr(A,B,Q,R);
K_val %computes the K matrix and displays
P_mat %positive definite matrix calculated for the same
Poles %To see the poles of the given equation
sys2 = ss(A-(B*K_val),B,C,D); %Using the K matrix to define ss
figure
initial(sys2,x_initial)
grid on
