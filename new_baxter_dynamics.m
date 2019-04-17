
function [M,C,G] = new_baxter_dynamics(q,dq)
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

% Intializing time dependent theta values (to be used for differentiation) 
syms t1(t) t2(t) t3(t) t4(t) t5(t) t6(t) t7(t); 

% Substituting symbolic theta variables (as a function of time) into 
% center of mass expressions.
cm1 = subs(cm1,[q1,q2,q3,q4,q5,q6,q7],[t1(t),t2(t),t3(t),t4(t),t5(t),...
        t6(t),t7(t)]);
cm2 = subs(cm2,[q1,q2,q3,q4,q5,q6,q7],[t1(t),t2(t),t3(t),t4(t),t5(t),...
        t6(t),t7(t)]);
cm3 = subs(cm3,[q1,q2,q3,q4,q5,q6,q7],[t1(t),t2(t),t3(t),t4(t),t5(t),...
        t6(t),t7(t)]);
cm4 = subs(cm4,[q1,q2,q3,q4,q5,q6,q7],[t1(t),t2(t),t3(t),t4(t),t5(t),...
        t6(t),t7(t)]);
cm5 = subs(cm5,[q1,q2,q3,q4,q5,q6,q7],[t1(t),t2(t),t3(t),t4(t),t5(t),...
        t6(t),t7(t)]);
cm6 = subs(cm6,[q1,q2,q3,q4,q5,q6,q7],[t1(t),t2(t),t3(t),t4(t),t5(t),...
        t6(t),t7(t)]);
cm7 = subs(cm7,[q1,q2,q3,q4,q5,q6,q7],[t1(t),t2(t),t3(t),t4(t),t5(t),...
        t6(t),t7(t)]) ;  
    
%Finding individual velocities
v1 = diff(cm1, t);
v2 = diff(cm2, t);
v3 = diff(cm3, t);
v4 = diff(cm4, t);
v5 = diff(cm5, t);
v6 = diff(cm6, t);
v7 = diff(cm7, t);

%Finding the Kinetic energies
k1 = (1/2)*m1*(transpose(v1)*v1);
k2 = (1/2)*m2*(transpose(v2)*v2);
k3 = (1/2)*m3*(transpose(v3)*v3);
k4 = (1/2)*m4*(transpose(v4)*v4);
k5 = (1/2)*m5*(transpose(v5)*v5);
k6 = (1/2)*m6*(transpose(v6)*v6);
k7 = (1/2)*m7*(transpose(v7)*v7);

%Total kinetic energy 
k = k1 + k2 + k3 + k4 + k5 + k6 + k7;

%%
%Finding the potential energies
p1 = m1 * g * (transpose(cm1) * [0; 0; 1]);
p2 = m2 * g * (transpose(cm2) * [0; 0; 1]);
p3 = m3 * g * (transpose(cm3) * [0; 0; 1]);
p4 = m4 * g * (transpose(cm4) * [0; 0; 1]);
p5 = m5 * g * (transpose(cm5) * [0; 0; 1]);
p6 = m6 * g * (transpose(cm6) * [0; 0; 1]);
p7 = m7 * g * (transpose(cm7) * [0; 0; 1]);

%Total potential energy 
p = p1 + p2 + p3 + p4 + p5 + p6 + p7;

%%
%Finding lagrangian
%disp('The Lagrangian');
L = k - p;

% Finding partial derivatives
 
% First substituting so that diff(t),t is differentiable
syms th1_d th2_d th3_d th4_d th5_d th6_d th7_d;
syms th1 th2 th3 th4 th5 th6 th7;

L = subs(L,[diff(t1(t),t), diff(t2(t),t),diff(t3(t),t),...
    diff(t4(t),t),diff(t5(t),t),diff(t6(t),t),diff(t7(t),t),t1(t), t2(t), t3(t), t4(t), t5(t), t6(t) t7(t)],[th1_d,...
    th2_d, th3_d, th4_d, th5_d, th6_d, th7_d,th1, th2,...
    th3, th4, th5, th6, th7]);

partial_wrt_th1 = diff(L, th1);
partial_wrt_th2 = diff(L, th2);
partial_wrt_th3 = diff(L, th3);
partial_wrt_th4 = diff(L, th4);
partial_wrt_th5 = diff(L, th5);
partial_wrt_th6 = diff(L, th6);
partial_wrt_th7 = diff(L, th7);

partial_wrt_th1_d = diff(L, th1_d);
partial_wrt_th2_d = diff(L, th2_d);
partial_wrt_th3_d = diff(L, th3_d);
partial_wrt_th4_d = diff(L, th4_d);
partial_wrt_th5_d = diff(L, th5_d);
partial_wrt_th6_d = diff(L, th6_d);
partial_wrt_th7_d = diff(L, th7_d);

%Reversing the substitutions to get time derivatives
partial_wrt_th1_d = subs(partial_wrt_th1_d,[th1_d, th2_d, th3_d, th4_d,...
    th5_d, th6_d, th7_d,th1, th2, th3, th4, th5, th6,...
    th7],[diff(t1(t),t), diff(t2(t),t) diff(t3(t),t),...
    diff(t4(t),t),diff(t5(t),t),diff(t6(t),t),diff(t7(t),t),t1(t), t2(t), t3(t), t4(t), t5(t), t6(t), t7(t)]);

partial_wrt_th2_d = subs(partial_wrt_th2_d,[th1_d, th2_d, th3_d, th4_d,...
    th5_d, th6_d, th7_d,th1, th2, th3, th4, th5, th6,...
    th7],[diff(t1(t),t), diff(t2(t),t) diff(t3(t),t),...
    diff(t4(t),t),diff(t5(t),t),diff(t6(t),t),diff(t7(t),t),t1(t), t2(t), t3(t), t4(t), t5(t), t6(t), t7(t)]);

partial_wrt_th3_d = subs(partial_wrt_th3_d,[th1_d, th2_d, th3_d, th4_d,...
    th5_d, th6_d, th7_d,th1, th2, th3, th4, th5, th6,...
    th7],[diff(t1(t),t), diff(t2(t),t) diff(t3(t),t),...
    diff(t4(t),t),diff(t5(t),t),diff(t6(t),t),diff(t7(t),t),t1(t), t2(t), t3(t), t4(t), t5(t), t6(t), t7(t)]);

partial_wrt_th4_d = subs(partial_wrt_th4_d,[th1_d, th2_d, th3_d, th4_d,...
    th5_d, th6_d, th7_d,th1, th2, th3, th4, th5, th6,...
    th7],[diff(t1(t),t), diff(t2(t),t) diff(t3(t),t),...
    diff(t4(t),t),diff(t5(t),t),diff(t6(t),t),diff(t7(t),t),t1(t), t2(t), t3(t), t4(t), t5(t), t6(t), t7(t)]);

partial_wrt_th5_d = subs(partial_wrt_th5_d,[th1_d, th2_d, th3_d, th4_d,...
    th5_d, th6_d, th7_d,th1, th2, th3, th4, th5, th6,...
    th7],[diff(t1(t),t), diff(t2(t),t) diff(t3(t),t),...
    diff(t4(t),t),diff(t5(t),t),diff(t6(t),t),diff(t7(t),t),t1(t), t2(t), t3(t), t4(t), t5(t), t6(t), t7(t)]);

partial_wrt_th6_d = subs(partial_wrt_th6_d,[th1_d, th2_d, th3_d, th4_d,...
    th5_d, th6_d, th7_d,th1, th2, th3, th4, th5, th6,...
    th7],[diff(t1(t),t), diff(t2(t),t) diff(t3(t),t),...
    diff(t4(t),t),diff(t5(t),t),diff(t6(t),t),diff(t7(t),t),t1(t), t2(t), t3(t), t4(t), t5(t), t6(t), t7(t)]);

partial_wrt_th7_d = subs(partial_wrt_th7_d,[th1_d, th2_d, th3_d, th4_d,...
    th5_d, th6_d, th7_d,th1, th2, th3, th4, th5, th6,...
    th7],[diff(t1(t),t), diff(t2(t),t) diff(t3(t),t),...
    diff(t4(t),t),diff(t5(t),t),diff(t6(t),t),diff(t7(t),t),t1(t), t2(t), t3(t), t4(t), t5(t), t6(t), t7(t)]);

%Finding time derivatives
time_derivative_tau1 = diff(partial_wrt_th1_d, t);
time_derivative_tau2 = diff(partial_wrt_th2_d, t);
time_derivative_tau3 = diff(partial_wrt_th3_d, t);
time_derivative_tau4 = diff(partial_wrt_th4_d, t);
time_derivative_tau5 = diff(partial_wrt_th5_d, t);
time_derivative_tau6 = diff(partial_wrt_th6_d, t);
time_derivative_tau7 = diff(partial_wrt_th7_d, t);

%%
%Finding individual torques
tau1 = time_derivative_tau1 - partial_wrt_th1;
tau2 = time_derivative_tau2 - partial_wrt_th2;
tau3 = time_derivative_tau3 - partial_wrt_th3;
tau4 = time_derivative_tau4 - partial_wrt_th4;
tau5 = time_derivative_tau5 - partial_wrt_th5;
tau6 = time_derivative_tau6 - partial_wrt_th6;
tau7 = time_derivative_tau7 - partial_wrt_th7;

%The Torque vector
Tau =  [tau1; tau2; tau3; tau4; tau5; tau6; tau7];

%%
%Performing substitutions to obtain M, C, G matrices
syms th1_dd th2_dd th3_dd th4_dd th5_dd th6_dd th7_dd
Tau = subs(Tau,[t1(t), t2(t),t3(t), t4(t), t5(t), t6(t), t7(t),diff(t1(t),t,t),diff(t2(t),t,t),diff(t3(t),t,t),...
    diff(t4(t),t,t),diff(t5(t),t,t),diff(t6(t),t,t),diff(t7(t),t,t),diff(t1(t),t),diff(t2(t),t),diff(t3(t),t),diff(t4(t),t),...
    diff(t5(t),t),diff(t6(t),t),diff(t7(t),t)],[th1,th2,th3,th4,th5,th6,th7,th1_dd, th2_dd, th3_dd, th4_dd, th5_dd, th6_dd, th7_dd,th1_d, th2_d, th3_d,th4_d, th5_d, th6_d, th7_d]);
                                                                                                                                                                                
%disp('The updated Tau');
disp(Tau);

% %%
% %The Inertia matrix(M)
% IM = equationsToMatrix(Tau,[th1_dd th2_dd th3_dd th4_dd th5_dd th6_dd,th7_dd]);
% M = subs(IM,[g,th1, th2,th3, th4, th5, th6, th7,th1_d, th2_d, th3_d,th4_d, th5_d, th6_d, th7_d],[9.8,q(1),q(2),q(3),q(4),q(5),q(6),q(7),dq(1),dq(2),dq(3),dq(4),dq(5),dq(6),dq(7)]);
% M =vpa(M)
% 
% 
% %%
% %The Gravity matrix(G)
% GM = equationsToMatrix(Tau,g);
% G = subs(GM,[g,th1, th2,th3, th4, th5, th6, th7,th1_d, th2_d, th3_d,th4_d, th5_d, th6_d, th7_d],[9.8,q(1),q(2),q(3),q(4),q(5),q(6),q(7),dq(1),dq(2),dq(3),dq(4),dq(5),dq(6),dq(7)]);
% G = vpa(G)
% 
% %%
% %The Coriolis matrix(C)
% CM = Tau - IM*[th1_dd;th2_dd;th3_dd;th4_dd;th5_dd;th6_dd;th7_dd] - GM;
% C = subs(CM,[g,th1, th2,th3, th4, th5, th6, th7,th1_d, th2_d, th3_d,th4_d, th5_d, th6_d, th7_d],[9.8,q(1),q(2),q(3),q(4),q(5),q(6),q(7),dq(1),dq(2),dq(3),dq(4),dq(5),dq(6),dq(7)]);
% C = vpa(C)


