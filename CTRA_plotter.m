close all

GNDT=importdata('groundtruth_logfile.mat');
GPS=importdata('GPS_logfile.mat');
IMU=importdata('IMU_logfile.mat');
EKF=importdata('EKF_Estimator.mat');

set(0,'defaultTextInterpreter','latex');


%X Position plots

figure('Name','Position x');
f=figure(1);
p=get(f,'position');
set(f,'Position',[50,50, p(3)*1.35,p(4)*1.85]) %maybe 1.5 or 1.8 to getfull page

ax1=subplot(6,1,1:3);hold on; grid on; box on;
hold all;
plot(GPS.x_gps,'color',[0.6,0.6,0.6],'Displayname','GPS','LineWidth',4);
%plot(GPS.x_gps,'color',[0.8,0.8,0.8],'Displayname','GPS','LineWidth',2);
plot(GNDT.x,'color',[0.4,0.4,0.4],'Displayname','Ground Truth','LineWidth',1,'LineStyle','--');
plot(EKF.Time,EKF.Data(:,1),'color',[0,0,0],'Displayname','Estimator','LineWidth',1,'LineStyle',':');
leg1=legend('location','northwest');
set(leg1,'Interpreter','latex')
ylabel('x displacement [m]') 
title('Position x')


% x absolute error 
subplot(6,1,4:6); hold on; grid on; box on;
plot(EKF.Time, abs(GNDT.x.Data - EKF.Data(:,1)))
title('Position x absolute error')
xlabel('Time [s]') 
ylabel('Absolute error [m]')
xticklabels(ax1,{})
%PNG File output
saveas(gcf,'f1_x.png');

%Y Position plots
figure('Name','Position y');
f=figure(2);
p=get(f,'position');
set(f,'Position',[50,50, p(3)*1.35,p(4)*1.85]) %maybe 1.5 or 1.8 to getfull page
ax1=subplot(6,1,1:3);hold on; grid on; box on;
hold all;
plot(GPS.y_gps,'color',[0.6,0.6,0.6],'Displayname','GPS','LineWidth',4);
%plot(GPS.y_gps, 'color',[0.8,0.8,0.8],'Displayname','GPS','LineWidth',2);
plot(GNDT.y,'color',[0.4,0.4,0.4],'Displayname','Ground Truth','LineWidth',1,'LineStyle','--');
plot(EKF.Time,EKF.Data(:,2),'color',[0,0,0],'Displayname','Estimator','LineWidth',1,'LineStyle',':');
leg1=legend('location','northwest');
set(leg1,'Interpreter','latex')
ylabel('y displacement [m]') 
title('Position y')

% y absolute error 
subplot(6,1,4:6); hold on; grid on; box on;
plot(EKF.Time, abs(GNDT.y.Data - EKF.Data(:,2)))
title('Position y absolute error')
xlabel('Time [s]') 
ylabel('Absolute error [m]')
xticklabels(ax1,{})
%PNG File output
saveas(gcf,'f2_y.png');

%x y plot figure 3
figure('Name','x y position');hold on; grid on; box on;
f=figure(3);
p=get(f,'position');
set(f,'Position',[50,50, p(3)*1.35,p(4)*1.85]) %maybe 1.5 or 1.8 to getfull page
hold all;
plot(GPS.x_gps.Data,GPS.y_gps.Data,'.','color',[0.6,0.6,0.6],'Displayname','GPS','LineWidth',2)
plot(GNDT.x.Data,GNDT.y.Data,'color',[0.4,0.4,0.4],'Displayname','Ground Truth','LineWidth',1,'LineStyle','--');
plot(EKF.Data(:,1),EKF.Data(:,2),'color',[0,0,0],'Displayname','Estimation','LineWidth',0.05,'LineStyle','-');
leg1=legend('location','northwest');
set(leg1,'Interpreter','latex')
daspect([2 2 1])
title('x y plot')
xlabel('x displacement [m]')
ylabel('y displacement [m]')
%PNG File output
saveas(gcf,'f3_xy.png');

%Theta orientation plots
figure('Name','Orientation theta');
f=figure(4);
p=get(f,'position');
set(f,'Position',[50,50, p(3)*1.35,p(4)*1.85]) %maybe 1.5 or 1.8 to getfull page
ax1=subplot(6,1,1:3);hold on; grid on; box on;
hold all;
plot(GPS.theta_gps*180/pi, 'color',[0.8,0.8,0.8],'Displayname','GPS','LineWidth',2);
plot(GNDT.theta*180/pi,'color',[0.4,0.4,0.4],'Displayname','Ground Truth','LineWidth',1,'LineStyle','--');
plot(EKF.Time,EKF.Data(:,3)*180/pi,'color',[0,0,0],'Displayname','Estimator','LineWidth',1,'LineStyle',':');
leg1=legend('location','northwest');
set(leg1,'Interpreter','latex')
ylabel('$\theta$ [$^{\circ}$]') 
title('Yaw angle $\theta$ ')

%Theta absolute error
subplot(6,1,4:6);hold on; grid on; box on;
plot(EKF.Time, abs((GNDT.theta.Data - EKF.Data(:,3)))*180/pi)
title('Yaw angle $\theta$ absolute error')
xlabel('Time [s]') 
ylabel('Absolute error [$^{\circ}$]')
xticklabels(ax1,{})
%PNG File output
saveas(gcf,'f4_theta.png');

%Speed plots
figure('Name','Speed v');
f=figure(5);
p=get(f,'position');
set(f,'Position',[50,50, p(3)*1.35,p(4)*1.85]) %maybe 1.5 or 1.8 to getfull page
ax1=subplot(6,1,1:3);hold on; grid on; box on;
hold all;
plot(GNDT.vf*3.6,'color',[0.4,0.4,0.4],'Displayname','Ground Truth','LineWidth',1,'LineStyle','--');
plot(EKF.Time,EKF.Data(:,4)*3.6,'color',[0,0,0],'Displayname','Estimator','LineWidth',1,'LineStyle',':');
leg1=legend('location','northwest');
set(leg1,'Interpreter','latex')
xlabel('')
ylabel('Speed [kph]') 
title('Speed v')

% v absolute error 
subplot(6,1,4:6); hold on; grid on; box on;
plot(EKF.Time, abs(GNDT.vf.Data - EKF.Data(:,4))*3.6)
title('Speed absolute error')
xlabel('Time [s]') 
ylabel('Absolute error [kph]')
xticklabels(ax1,{})
%PNG File output
saveas(gcf,'f5_v.png');


%Yaw rate plots
figure('Name','Yaw rate');
f=figure(6);
p=get(f,'position');
set(f,'Position',[50,50, p(3)*1.35,p(4)*1.85]) %maybe 1.5 or 1.8 to getfull page
ax1=subplot(6,1,1:3);hold on; grid on; box on;
hold all;
plot(IMU.theta_dot_imu*180/pi, 'color',[0.8,0.8,0.8],'Displayname','IMU','LineWidth',2);
plot(GNDT.theta_dot*180/pi,'color',[0.4,0.4,0.4],'Displayname','Ground Truth','LineWidth',1,'LineStyle','--');
plot(EKF.Time,EKF.Data(:,5)*180/pi,'color',[0,0,0],'Displayname','Estimator','LineWidth',1,'LineStyle',':');
leg1=legend('location','northwest');
set(leg1,'Interpreter','latex')
ylabel('$\dot{\theta}$ [$^\circ$/s]') 
title('Yaw rate $\dot{\theta}$')

%Theta dot absolute error
subplot(6,1,4:6);hold on; grid on; box on;
plot(EKF.Time, abs((GNDT.theta_dot.Data - EKF.Data(:,5)))*180/pi)
title('Yaw rate $\dot{\theta}$ Absolute error ')
xlabel('Time [s]') 
ylabel('Absolute error [$^\circ$/s]')
xticklabels(ax1,{})
%PNG File output
saveas(gcf,'f6_theta_dot.png');

%acceleration plots
figure('Name','Acceleration a');
figure(7);
ax1=subplot(6,1,1:3);hold on; grid on; box on;
hold all;
plot(IMU.x_ddot_imu, 'b', 'color',[0.8,0.8,0.8],'Displayname','IMU','LineWidth',2);
plot(GNDT.X_ddot,'k', 'color',[0.4,0.4,0.4],'Displayname','Ground Truth','LineWidth',1,'LineStyle','--');
plot(EKF.Time,EKF.Data(:,6),'r', 'color',[0,0,0],'Displayname','Estimator','LineWidth',1,'LineStyle',':');
leg1=legend('location','northwest');
set(leg1,'Interpreter','latex')
xlabel('')
ylabel('Acceleration [m/s^2]') 
title('Acceleration')

% a absolute error 
subplot(6,1,4:6); hold on; grid on; box on;
plot(EKF.Time, abs(GNDT.X_ddot.Data - EKF.Data(:,6)))
title('Acceleration absolute error')
xlabel('Time [s]')
ylabel('Absolute error [m/s^2]')
xticklabels(ax1,{})

%PNG File output
saveas(gcf,'f7_a.png');

% RMSE Values calculation
%x RMSE
x_RMSE = sqrt(mean((GNDT.x.Data - EKF.Data(:,1)).^2));

%y RMSE
y_RMSE = sqrt(mean((GNDT.y.Data - EKF.Data(:,2)).^2));

%Theta RMSE
theta_RMSE = sqrt(mean((GNDT.theta.Data - EKF.Data(:,3)).^2));

%v RMSE
v_RMSE = sqrt(mean((GNDT.vf.Data - EKF.Data(:,4)).^2));

%theta dot RMSE
theta_dot_RMSE = sqrt(mean((GNDT.theta_dot.Data - EKF.Data(:,5)).^2));