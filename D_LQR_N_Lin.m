clear 
clc

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
Q=[0.05 0 0 0 0 0;
   0 1 0 0 0 0;
   0 0 1000 0 0 0;
   0 0 0 10 0 0;
   0 0 0 0 500 0;
   0 0 0 0 0 10];
R=0.00001;

[K, P_mat, Poles] = lqr(A,B,Q,R);
y0 = [0; 0; pi/6; 0; pi/6; 0]

tspan = 0:0.3:60;
[t,y] = ode45(@(t,y)lqrc(y,K,M,m1,m2,l1,l2,g),tspan,y0);

figure
plot(t,y)
grid on

figure
for k = 1:length(t)
    draw_cart_pendulum(y(k, :), 1, 5, 2 ,1);  
end

figure
plot(t,y)
grid on


function dydt = lqrc(y,K,M,m1,m2,l1,l2,g)
F=-K*y;
dydt=zeros(6,1);
% y(1)=x; y(2)=xdot; y(3)=theta1;   y(4)=theta1dot;  y(5)=theta2;    y(6)=theta2dot;
dydt(1) = y(2); 
dydt(2)=(F-(g/2)*(m1*sin(2*y(3))+m2*sin(2*y(5)))-(m1*l1*(y(4)^2)*sin(y(3)))-(m2*l2*(y(6)^2)*sin(y(5))))/(M+m1*((sin(y(3)))^2)+m2*((sin(y(5)))^2));
dydt(3)= y(4); 
dydt(4)= (dydt(2)*cos(y(3))-g*(sin(y(3))))/l1'; 
dydt(5)= y(6); 
dydt(6)= (dydt(2)*cos(y(5))-g*(sin(y(5))))/l2; 
end

function draw_cart_pendulum(y,m, M, L ,l)
  x = y(1);
  theta = y(3);
  theta2 = y(5);
  
  W = 1*sqrt(M/5);   
  H = 0.5*sqrt(M/5); 
  wr = 0.2;          
  mr = 0.3*sqrt(m);    
  
  y = wr/2 + H/2;
  w1x = x - 0.9*W/2;
  w1y = 0;
  w2x = x + 0.9*W/2 - wr;
  w2y = 0;
  
  px = x + L*sin(theta);
  py = y - L*cos(theta);

  px1 = x + l*sin(theta2);
  py1 = y - l*cos(theta2);
   
  hold on;
  clf;
  axis equal;
  rectangle('Position',[x-W/2,y-H/2,W,H],'Curvature',.1,'FaceColor',[1 0.1 0.1])
  rectangle('Position',[w1x,w1y,wr,wr],'Curvature',1,'FaceColor',[0 0 0])
  rectangle('Position',[w2x,w2y,wr,wr],'Curvature',1,'FaceColor',[0 0 0])
  
  line ([-10 10], [0 0], "linestyle", "-", "color", "k");
  line ([x px], [y py], "linestyle", "-", "color", "k");
  rectangle('Position',[px-mr/2,py-mr/2,mr,mr],'Curvature',1,'FaceColor',[.1 0.1 1])

  line ([-10 30], [0 0], "linestyle", "-", "color", "k");
  line ([x px1], [y py1], "linestyle", "-", "color", "k");
  rectangle('Position',[px1-mr/2,py1-mr/2,mr,mr],'Curvature',1,'FaceColor',[.1 0.1 1])
  
  xlim([-5 10]);
  ylim([-2 3]);
  set(gcf, 'Position', [300 300 1000 400]);
  drawnow
  hold off
end