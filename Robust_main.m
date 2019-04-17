%% Robust Controller for 1 link Arm
clc
clear all;
close all;

%% Initial values
tf=5;
q=[10,10,10,10,10,10,10]; %Joint Angle
dq=[0.5,0.5,0.5,0.5,0.5,0.5,0.5]; %Joint Angle
[M,C,G]= new_baxter_dynamics(q,dq)
% q0=-10; q0_f=80;
% dq0 =0.5; dq0_f = 0.4;
% 
% %% Desired state Trajectory
% %Using cubic function for trajectory generation
% a(:,1) = cubic(q0,dq0,q0_f,dq0_f,tf); 
% T = linspace(0,tf,200);
% X(1,:) = a(1,1) + a(2,1).*T +a(3,1).*T.^2 + a(4,1).*T.^3 ;
% 
% %% Ode Function: Robust
% I_est = 8;
% mgd_est = 5;
% fv_est= 2.5;
% kp=90; kd=40;
% z0=[1,0.5]; %Small Deviation
% z_l=[50,0.5]; %Large Deviation
% param = {kp,kd,a,I_est,mgd_est,fv_est};
% options = odeset('RelTol',1e-4,'AbsTol',[1e-4, 1e-4]);
% [T2,Z] = ode45(@(t,z) robust_link(t,z,param),[0 tf],z0, options);
% [T22,Z2] = ode45(@(t,z1) robust_link(t,z1,param),[0 tf],z_l, options);
% %% Plotting the result
% figure(1);
% hold on
% title('Robust Controller: 1 Link Robot Arm');
% xlabel('Time');
% ylabel('Angle');
% plot(T2, Z(:,1));
% hold on
% plot(T22, Z2(:,1));
% hold on
% plot(T, X(1,:));
% legend('Actual(Small Deviation)','Actual(Large Deviation)','Desired');
% grid on
% 
