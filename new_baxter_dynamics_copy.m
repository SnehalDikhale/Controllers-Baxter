
function [M,C,G] = new_baxter_dynamics_copy(q,dq)
%%
% Performing forward kinematics for the robot

% Initializing symbolic thetas for differentiation
syms q1 q2 q3 q4 q5 q6 q7;

% Intitializing symbolic masses for links
m1=5.70044; m2=3.22698; m3=4.31272; m4=2.07206; m5=2.24665; m6=1.60979; m7=0.54218;  
% Initializing gravity
syms g;
% Intializing lengths of links
d1=270.35; d3=364.35; d5=374.29; d7= 229.525; a1=69; a3=69; a5=10;

% Getting individual homogenous transformation matrices
T_0_1 = calcdh_deg(q1, d1, a1, -pi/2);
T_1_2 = calcdh_deg((q2-pi/2), 0, 0, -pi/2);
T_2_3 = calcdh_deg(q3, d3, -a3, pi/2);
T_3_4 = calcdh_deg(q4, 0, 0, -pi/2);
T_4_5 = calcdh_deg(q5, d5, -a5, pi/2);
T_5_6 = calcdh_deg(q6, 0, 0, -pi/2);
T_6_7 = calcdh_deg(q7, d7, 0, 0);

%Getting all rotation matrices
R_0_2 = T_0_1(1:3,1:3);
R_0_3 = T_0_1(1:3,1:3)*T_1_2(1:3,1:3);
R_0_4 = T_0_1(1:3,1:3)*T_1_2(1:3,1:3)*T_2_3(1:3,1:3);
R_0_5 = T_0_1(1:3,1:3)*T_1_2(1:3,1:3)*T_2_3(1:3,1:3)*T_3_4(1:3,1:3);
R_0_6 = T_0_1(1:3,1:3)*T_1_2(1:3,1:3)*T_2_3(1:3,1:3)*T_3_4(1:3,1:3)*...
        T_4_5(1:3,1:3);
R_0_7 = T_0_1(1:3,1:3)*T_1_2(1:3,1:3)*T_2_3(1:3,1:3)*T_3_4(1:3,1:3)*...
        T_4_5(1:3,1:3)*T_5_6(1:3,1:3);

% Final Homogenous transformation matrix
T_base_tip = simplify(T_0_1*T_1_2*T_2_3*T_3_4*T_4_5*T_5_6*T_6_7);

% Getting position of center of masses
cm1 = T_0_1(1:3,4)/2;
cm2 = T_0_1(1:3,4) + R_0_2*T_1_2(1:3,4)/2;
cm3 = T_0_1(1:3,4) + R_0_2*T_1_2(1:3,4) + R_0_3*T_2_3(1:3,4)/2;
cm4 = T_0_1(1:3,4) + R_0_2*T_1_2(1:3,4) + R_0_3*T_2_3(1:3,4) + ...
    R_0_4*T_3_4(1:3,4)/2;
cm5 = T_0_1(1:3,4) + R_0_2*T_1_2(1:3,4) + R_0_3*T_2_3(1:3,4) + ...
    R_0_4*T_3_4(1:3,4) + R_0_5*T_4_5(1:3,4)/2;
cm6 = T_0_1(1:3,4) + R_0_2*T_1_2(1:3,4) + R_0_3*T_2_3(1:3,4) + ...
    R_0_4*T_3_4(1:3,4) + R_0_5*T_4_5(1:3,4) + R_0_6*T_5_6(1:3,4)/2;
cm7 = T_0_1(1:3,4) + R_0_2*T_1_2(1:3,4) + R_0_3*T_2_3(1:3,4) + ...
    R_0_4*T_3_4(1:3,4) + R_0_5*T_4_5(1:3,4) + R_0_6*T_5_6(1:3,4)+...
    R_0_7*T_6_7(1:3,4)/2;


 syms t1(t) t2(t) t3(t) t4(t) t5(t) t6(t) t7(t); 

 %syms cm1(q1,q2,q3,q4,q5,q6,q7);

% Substituting symbolic theta variables (as a function of time) into 
% center of mass expressions.
[k1,p1]= dynamics(cm1,m1);
[k2,p2]= dynamics(cm2,m2);
[k3,p3]= dynamics(cm3,m3);
[k4,p4]= dynamics(cm4,m4);
[k5,p5]= dynamics(cm5,m5);
[k6,p6]= dynamics(cm6,m6);
[k7,p7]= dynamics(cm7,m7);

%Total kinetic energy 
k = k1 + k2 + k3 + k4 + k5 + k6 + k7;
%Total potential energy 
p = p1 + p2 + p3 + p4 + p5 + p6 + p7;

%%
%Finding lagrangian
L = k - p;

% Finding partial derivatives
% % First substituting so that diff(t),t is differentiable
syms th1_d th2_d th3_d th4_d th5_d th6_d th7_d;
syms th1 th2 th3 th4 th5 th6 th7;

L = subs(L,[diff(t1(t),t), diff(t2(t),t),diff(t3(t),t),...
    diff(t4(t),t),diff(t5(t),t),diff(t6(t),t),diff(t7(t),t),t1(t), t2(t), t3(t), t4(t), t5(t), t6(t) t7(t)],[th1_d,...
    th2_d, th3_d, th4_d, th5_d, th6_d, th7_d,th1, th2,...
    th3, th4, th5, th6, th7]);

%Reversing the substitutions to get time derivatives
partial_wrt_th1_d = subs(diff(L, th1_d),[th1_d, th2_d, th3_d, th4_d,...
    th5_d, th6_d, th7_d,th1, th2, th3, th4, th5, th6,...
    th7],[diff(t1(t),t), diff(t2(t),t) diff(t3(t),t),...
    diff(t4(t),t),diff(t5(t),t),diff(t6(t),t),diff(t7(t),t),t1(t), t2(t), t3(t), t4(t), t5(t), t6(t), t7(t)]);

partial_wrt_th2_d = subs(diff(L, th2_d),[th1_d, th2_d, th3_d, th4_d,...
    th5_d, th6_d, th7_d,th1, th2, th3, th4, th5, th6,...
    th7],[diff(t1(t),t), diff(t2(t),t) diff(t3(t),t),...
    diff(t4(t),t),diff(t5(t),t),diff(t6(t),t),diff(t7(t),t),t1(t), t2(t), t3(t), t4(t), t5(t), t6(t), t7(t)]);

partial_wrt_th3_d = subs(diff(L, th3_d),[th1_d, th2_d, th3_d, th4_d,...
    th5_d, th6_d, th7_d,th1, th2, th3, th4, th5, th6,...
    th7],[diff(t1(t),t), diff(t2(t),t) diff(t3(t),t),...
    diff(t4(t),t),diff(t5(t),t),diff(t6(t),t),diff(t7(t),t),t1(t), t2(t), t3(t), t4(t), t5(t), t6(t), t7(t)]);

partial_wrt_th4_d = subs(diff(L, th4_d),[th1_d, th2_d, th3_d, th4_d,...
    th5_d, th6_d, th7_d,th1, th2, th3, th4, th5, th6,...
    th7],[diff(t1(t),t), diff(t2(t),t) diff(t3(t),t),...
    diff(t4(t),t),diff(t5(t),t),diff(t6(t),t),diff(t7(t),t),t1(t), t2(t), t3(t), t4(t), t5(t), t6(t), t7(t)]);

partial_wrt_th5_d = subs(diff(L, th5_d),[th1_d, th2_d, th3_d, th4_d,...
    th5_d, th6_d, th7_d,th1, th2, th3, th4, th5, th6,...
    th7],[diff(t1(t),t), diff(t2(t),t) diff(t3(t),t),...
    diff(t4(t),t),diff(t5(t),t),diff(t6(t),t),diff(t7(t),t),t1(t), t2(t), t3(t), t4(t), t5(t), t6(t), t7(t)]);

partial_wrt_th6_d = subs(diff(L, th6_d),[th1_d, th2_d, th3_d, th4_d,...
    th5_d, th6_d, th7_d,th1, th2, th3, th4, th5, th6,...
    th7],[diff(t1(t),t), diff(t2(t),t) diff(t3(t),t),...
    diff(t4(t),t),diff(t5(t),t),diff(t6(t),t),diff(t7(t),t),t1(t), t2(t), t3(t), t4(t), t5(t), t6(t), t7(t)]);

partial_wrt_th7_d = subs(diff(L, th7_d),[th1_d, th2_d, th3_d, th4_d,...
    th5_d, th6_d, th7_d,th1, th2, th3, th4, th5, th6,...
    th7],[diff(t1(t),t), diff(t2(t),t) diff(t3(t),t),...
    diff(t4(t),t),diff(t5(t),t),diff(t6(t),t),diff(t7(t),t),t1(t), t2(t), t3(t), t4(t), t5(t), t6(t), t7(t)]);

%The Torque vector
Tau =  [diff(partial_wrt_th1_d, t) - diff(L, th1); 
         diff(partial_wrt_th2_d, t) - diff(L, th2); 
         diff(partial_wrt_th3_d, t) - diff(L, th3);
         diff(partial_wrt_th4_d, t) - diff(L, th4) ; 
         diff(partial_wrt_th5_d, t) - diff(L, th5); 
         diff(partial_wrt_th6_d, t) - diff(L, th6); 
         diff(partial_wrt_th7_d, t) - diff(L, th7)];
% 
% %%
%Performing substitutions to obtain M, C, G matrices
syms th1_dd th2_dd th3_dd th4_dd th5_dd th6_dd th7_dd
Tau = subs(Tau,[t1(t), t2(t),t3(t), t4(t), t5(t), t6(t), t7(t),diff(t1(t),t,t),diff(t2(t),t,t),diff(t3(t),t,t),...
    diff(t4(t),t,t),diff(t5(t),t,t),diff(t6(t),t,t),diff(t7(t),t,t),diff(t1(t),t),diff(t2(t),t),diff(t3(t),t),diff(t4(t),t),...
    diff(t5(t),t),diff(t6(t),t),diff(t7(t),t)],[th1,th2,th3,th4,th5,th6,th7,th1_dd, th2_dd, th3_dd, th4_dd, th5_dd, th6_dd, th7_dd,th1_d, th2_d, th3_d,th4_d, th5_d, th6_d, th7_d]);                                                                                                                                                                            
% %disp('The updated Tau');
% disp(Tau);
% 
% %%
%The Inertia matrix(M)
% syms IM(g,th1, th2,th3, th4, th5, th6, th7,th1_d, th2_d, th3_d,th4_d, th5_d, th6_d, th7_d)
 IM = equationsToMatrix(Tau,[th1_dd th2_dd th3_dd th4_dd th5_dd th6_dd,th7_dd]);

%M = subs(IM,[g,th1, th2,th3, th4, th5, th6, th7,th1_d, th2_d, th3_d,th4_d, th5_d, th6_d, th7_d],[9.8,q(1),q(2),q(3),q(4),q(5),q(6),q(7),dq(1),dq(2),dq(3),dq(4),dq(5),dq(6),dq(7)]);
%M = IM(9.8,q(1),q(2),q(3),q(4),q(5),q(6),q(7),dq(1),dq(2),dq(3),dq(4),dq(5),dq(6),dq(7))
M =vpa(M,3);


%%
%The Gravity matrix(G)
GM = equationsToMatrix(Tau,g);
G = subs(GM,[g,th1, th2,th3, th4, th5, th6, th7,th1_d, th2_d, th3_d,th4_d, th5_d, th6_d, th7_d],[9.8,q(1),q(2),q(3),q(4),q(5),q(6),q(7),dq(1),dq(2),dq(3),dq(4),dq(5),dq(6),dq(7)]);
G = vpa(G,3);

%%
%The Coriolis matrix(C)
CM = (subs(Tau,[th1_dd, th2_dd, th3_dd,th4_dd, th5_dd, th6_dd, th7_dd],[0,0,0,0,0,0,0]))- G;
C = subs(CM,[g,th1, th2,th3, th4, th5, th6, th7,th1_d, th2_d, th3_d,th4_d, th5_d, th6_d, th7_d],[9.8,q(1),q(2),q(3),q(4),q(5),q(6),q(7),dq(1),dq(2),dq(3),dq(4),dq(5),dq(6),dq(7)]);
C = vpa(C,3);


