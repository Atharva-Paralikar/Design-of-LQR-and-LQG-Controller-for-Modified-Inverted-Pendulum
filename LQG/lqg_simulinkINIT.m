clear all, close all, clc

M=1000;
m1=100;
m2=100;
l1=20;
l2=10;
g=9.81;
s = 1; 

A=[0 1 0 0 0 0; 
    0 0 -(m1*g)/M 0 -(m2*g)/M 0;
    0 0 0 1 0 0;
    0 0 -s*((M+m1)*g)/(M*l1) 0 -s*(m2*g)/(M*l1) 0;
    0 0 0 0 0 1;
    0 0 -s*(m1*g)/(M*l2) 0 -s*(g*(M+m2))/(M*l2) 0];


B=[0; 1/M; 0; 1/(M*l1); 0; s*1/(M*l2)];
 
C = [1 0 0 0 0 0];

D = zeros(size(C,1),size(B,2));

Q=[0.05 0 0 0 0 0;
   0 1 0 0 0 0;
   0 0 1000 0 0 0;
   0 0 0 10 0 0;
   0 0 0 0 500 0;
   0 0 0 0 0 10];
R=0.00001;

K = lqr(A,B,Q,R); 

Vd = .001*eye(6);
Vn = .001;


BF = [B Vd 0*B];

sysC = ss(A,BF,C,[0 0 0 0 0 0 0 Vn]);


sysFullOutput = ss(A,BF,eye(6),zeros(6,size(BF,2)));


 [L,P,E] = lqe(A,Vd,C,Vd,Vn);  
 Kf = (lqr(A',C',Vd,Vn))';  

 sysKF = ss(A-Kf*C,[B Kf],eye(6),0*[B Kf]);  
