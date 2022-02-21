clc
clear

syms M m1 m2 l1 l2 g;

A=[0 1 0 0 0 0; 
    0 0 -(m1*g)/M 0 -(m2*g)/M 0;
    0 0 0 1 0 0;
    0 0 -((M+m1)*g)/(M*l1) 0 -(m2*g)/(M*l1) 0;
    0 0 0 0 0 1;
    0 0 -(m1*g)/(M*l2) 0 -(g*(M+m2))/(M*l2) 0];
B=[0; 1/M; 0; 1/(M*l1); 0; 1/(M*l2)];

C1 = [1 0 0 0 0 0];  

O1 = [C1' A'*C1' A'*A'*C1' A'*A'*A'*C1' A'*A'*A'*A'*C1' A'*A'*A'*A'*A'*C1'];
if rank(O1)==6
    disp('x(t):System is observable')
else
    disp('x(t):System is not observable')
end

C2 = [0 0 1 0 0 0; 0 0 0 0 1 0]; 

O2 = [C2' A'*C2' A'*A'*C2' A'*A'*A'*C2' A'*A'*A'*A'*C2' A'*A'*A'*A'*A'*C2'];
if rank(O2)==6 
    disp('Theta_1(t)Theta_2(t):System is observable')
else
    disp('Theta_1(t)Theta_2(t):System is not observable')
end

C3 = [1 0 0 0 0 0; 0 0 0 0 1 0]; 
 
O3 = [C3' A'*C3' A'*A'*C3' A'*A'*A'*C3' A'*A'*A'*A'*C3' A'*A'*A'*A'*A'*C3'];
if rank(O3)==6
    disp('x(t)Theta_2(t):System is observable')
else
    disp('x(t)Theta_2(t):System is not observable')
end

C4 = [1 0 0 0 0 0; 0 0 1 0 0 0; 0 0 0 0 1 0];

O4 = [C4' A'*C4' A'*A'*C4' A'*A'*A'*C4' A'*A'*A'*A'*C4' A'*A'*A'*A'*A'*C4'];
if rank(O4)
    disp('x(t)Theta_1(t)Theta_2(t):System is observable')
else
    disp('x(t)Theta_1(t)Theta_2(t):System is not observable')
end