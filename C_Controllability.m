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

disp("Controllability:");
C= [B A*B A*A*B A*A*A*B A*A*A*A*B A*A*A*A*A*B]
disp("Determinant:");
disp(simplify(det(C)));
disp("Rank:");
rank(C)

disp("Conditional substitution: (l1=l2)")
Ct1 = subs(C,l1,l2) 
disp("Rank :")

rank(Ct1)

if (rank(Ct1)==rank(C))
    disp("System is controllable")
else
    disp("System is not controllable")
end