clc
clear
M=1000;
m1=100;
m2=100;
l1=20;
l2=10;
g=9.81; 

A=[0 1 0 0 0 0; 
    0 0 -(m1*g)/M 0 -(m2*g)/M 0;
    0 0 0 1 0 0;
    0 0 -((M+m1)*g)/(M*l1) 0 -(m2*g)/(M*l1) 0;
    0 0 0 0 0 1;
    0 0 -(m1*g)/(M*l2) 0 -(g*(M+m2))/(M*l2) 0];

B=[0; 1/M; 0; 1/(M*l1); 0; 1/(M*l2)];

C1 = [1 0 0 0 0 0];  
C2 = [1 0 0 0 0 0; 0 0 0 0 1 0]; 
C3 = [1 0 0 0 0 0; 0 0 1 0 0 0; 0 0 0 0 1 0]; 
D = 0;

Q=[0.01 0 0 0 0 0;
   0 1 0 0 0 0;
   0 0 10000 0 0 0;
   0 0 0 1 0 0;
   0 0 0 0 5000 0;
   0 0 0 0 0 1];
R=0.001; 


x0=[0,0,30,0,60,0,0,0,0,0,0,0];

poles=[-1;-2;-3;-4;-5;-6];

K=lqr(A,B,Q,R);

L1 = place(A',C1',poles)' 
L2 = place(A',C2',poles)' 
L3 = place(A',C3',poles)' 

Al1 = [(A-B*K) B*K; 
        zeros(size(A)) (A-L1*C1)];
Bl = [B;zeros(size(B))];
Cl1 = [C1 zeros(size(C1))];

Al2 = [(A-B*K) B*K;
        zeros(size(A)) (A-L2*C2)];
Cl2 = [C2 zeros(size(C2))];

Al3 = [(A-B*K) B*K;
        zeros(size(A)) (A-L3*C3)];
Cl3 = [C3 zeros(size(C3))];

sys1 = ss(Al1, Bl, Cl1, D);
figure 
initial(sys1,x0)
figure
step(sys1)

sys2 = ss(Al2, Bl, Cl2, D);
figure 
initial(sys2,x0)
figure 
step(sys2)

sys3 = ss(Al3, Bl, Cl3, D);
figure 
initial(sys3,x0)
figure 
step(sys3)

grid on