%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% init()                                                                  %
%                                                                         %              
% Set initial parameters for part1.slx and part2.slx                      %
%                                                                         %
% Created:      2018.07.12	Jon Bjørnø                                    %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all;
close all

load('supply.mat');
load('supplyABC.mat');
load('thrusters_sup.mat')

%% INITIAL CONDITION
% Initial position x, y, z, phi, theta, psi
eta0 = [0,0,0,0,0,0]';
% Initial velocity u, v, w, p, q, r
nu0 = [0,0,0,0,0,0]';




%% CURRENT INPUT
Current_Speed = 0;            % m/s
Current_input_type = 1;         % 0 = constant heading, 1 = linear varied heading

CurrHead.Slope = 0.01;          %[rad/s] slop of the linear varied heading
CurrHead.TimeStart = 300;       %[s] start at which the linear varied heading starts
% current heading given aaccording N-W coodinates from whee the current is
% flowing. An angle of pi() is added in the code to uniform it with the
% reference sstem
CurrHead.HeadStart = 0;         %[rad] initial heading value (alway to be defined, it 
                                % is used also for he case of "consat
                                % heading"
CurrHead.HeadEnd = pi/2;          %[rad] final heading value


%% CONTROLLER INPUT
Kp_surge = (2*pi/30)^2*(6.4794^5+6.362*10^6);   % based on natural period of abot 30 seconds from RAOs
Ki_surge = 18000;
Kd_surge = 2*0.7*sqrt(6.4794^5+6.362*10^6)*sqrt(Kp_surge);
% Kp_surge = 0;   
% Ki_surge = 0;
% Kd_surge = 0;


Kp_sway = (2*pi/20)^2*(2.157*10^6+6.362*10^6);        % based on natural period of 25 seconds from RAOs
Ki_sway = 80000;
Kd_sway = 2*0.7*sqrt(2.157*10^6+6.362*10^6)*sqrt(Kp_sway);
% Kp_sway = 0;       
% Ki_sway = 0;
% Kd_sway = 0;

Kp_yaw = (2*pi/10)^2*(1.0711*10^9+2.724*10^9);    % based on natural period of 10 seconds from RAOs
Ki_yaw = 1.5*10^8;
Kd_yaw = 2*0.7*sqrt(1.0711*10^9+2.724*10^9)*sqrt(Kp_yaw);
% Kp_yaw = 0;  
% Ki_yaw = 0;
% Kd_yaw = 0;



%% SET-POINTS INPUT
typeOfSimulation = 1; % 0=constant set-point , 1 sequence of set-points

%if constant set-point
ConstSPSuge = 0;  %m
ConstSPSway = 0;  %m
ConstSPYaw = 0;   %rad

% if sequence of set-points
% Used to tune the controller
% n0 = [0 0 0];
% n1 = [50 0 0];
% n2 = [50 0 0];
% n3 = n0;
% n4 = n0;
% n5 = n0;

n0 = [0 0 0];
n1 = [50 0 0];
n2 = [50 -50 0];
n3 = [50 -50 -pi/4];
n4 = [0 -50 -pi/4];
n5 = [0 0 0];

nEnd = [0 0 0];
TimeSteady = 200;       % seconds , time for which the system has to stay  steady
TimeTansition = 200;    % seconds , transidition time between one set-point to the next


% do not change anything from now on
setPoints = [n0;n0;n1;n1;n2;n2;n3;n3;n4;n4;n5;n5;nEnd];

timeVector = [0;TimeSteady;TimeTansition;TimeSteady;TimeTansition;TimeSteady;TimeTansition;...
    TimeSteady;TimeTansition;TimeSteady;TimeTansition;TimeSteady;TimeTansition];

for idx =1:length(timeVector)
    timeSetPoints(idx)=sum(timeVector(1:idx));
end   
% RateOfChangeSuge = 10;  % m/s
% RateOfChangeSway = 10;  % m/s
% RateOfChangeYaw = 1;  % rad/s


SurgeSP = [timeSetPoints' setPoints(:,1)];
SwaySP = [timeSetPoints' setPoints(:,2)];
YawSP = [timeSetPoints' setPoints(:,3)];



%% SIMULATION

%sim('part1a','StartTime','0.0','StopTime','2400')
sim('part1')


%% Plot results

figure()
linkx(1)=subplot(3,1,1)
plot(SetPointPos.Time,SetPointPos.Data(:,1))
hold on
plot(Eta.Time,Eta.Data(:,1))
hold off
grid on
xlabel('Time [s]')
ylabel('X_E [m]')
legend('Set-point','Vessel')
title('Vessel position - Earth reference frame')
linkx(2)=subplot(3,1,2)
plot(SetPointPos.Time,SetPointPos.Data(:,2))
hold on
plot(Eta.Time,Eta.Data(:,2))
hold off
grid on
xlabel('Time [s]')
ylabel('Y_E [m]')
legend('Set-point','Vessel')
linkx(3)=subplot(3,1,3)
plot(SetPointPos.Time,SetPointPos.Data(:,3))
hold on
plot(Eta.Time,Eta.Data(:,3))
hold off
grid on
legend('Set-point','Vessel')
xlabel('Time [s]')
ylabel('Heading [rad]')

figure()
linkx(4)=subplot(3,1,1)
plot(SetPointSpeed.Time,SetPointSpeed.Data(:,1))
hold on
plot(Nu.Time,Nu.Data(:,1))
hold off
grid on
xlabel('Time [s]')
ylabel('Surge velocity [m/s]')
legend('Set-point','Vessel')
title('Vessel velocities - Earth reference frame')
linkx(5)=subplot(3,1,2)
plot(SetPointSpeed.Time,SetPointSpeed.Data(:,2))
hold on
plot(Nu.Time,Nu.Data(:,2))
hold off
grid on
xlabel('Time [s]')
ylabel('Sway velocity [m/s]')
legend('Set-point','Vessel')
linkx(6)=subplot(3,1,3)
plot(SetPointSpeed.Time,SetPointSpeed.Data(:,3))
hold on
plot(Nu.Time,Nu.Data(:,3))
hold off
grid on
legend('Set-point','Vessel')
xlabel('Time [s]')
ylabel('Yaw rate [rad/s]')


figure()
plot(SetPointPos.Data(:,2),SetPointPos.Data(:,1));
hold on
plot(Eta.Data(:,2),Eta.Data(:,1))
%plot(Eta.Data(:,1),Eta.Data(:,2))
hold off
grid on
legend('Set-point','Vessel')
ylabel('X_E [m] Earth reference frame')
xlabel('Y_E [m] Earth reference frame')
% xlim([min(SetPointPos.Data(:,1))-abs(max(SetPointPos.Data(:,1))*0.3) max(SetPointPos.Data(:,1))+abs(max(SetPointPos.Data(:,1))*0.3)])
% ylim([min(SetPointPos.Data(:,2))-abs(max(SetPointPos.Data(:,2))*0.3) max(SetPointPos.Data(:,2))+abs(max(SetPointPos.Data(:,2))*0.3)])
axis equal

figure()
linkx(7)=subplot(3,1,1)
plot(Tau_Surge.Time,Tau_Surge.Data(:,1));
hold on
plot(Tau_Surge.Time,Tau_Surge.Data(:,2));
plot(Tau_Surge.Time,Tau_Surge.Data(:,3));
hold off
grid on
title('Surge')
legend('\tau_P','\tau_D','\tau_I')
xlabel('Time [s]')
ylabel('\tau ')
linkx(8)=subplot(3,1,2)
plot(Tau_Sway.Time,Tau_Sway.Data(:,1));
hold on
plot(Tau_Sway.Time,Tau_Sway.Data(:,2));
plot(Tau_Sway.Time,Tau_Sway.Data(:,3));
hold off
grid on
title('Sway')
legend('\tau_P','\tau_D','\tau_I')
xlabel('Time [s]')
ylabel('\tau ')
linkx(9)=subplot(3,1,3)
plot(Tau_Yaw.Time,Tau_Yaw.Data(:,1));
hold on
plot(Tau_Yaw.Time,Tau_Yaw.Data(:,2));
plot(Tau_Yaw.Time,Tau_Yaw.Data(:,3));
hold off
grid on
title('Yaw')
legend('\tau_P','\tau_D','\tau_I')
xlabel('Time [s]')
ylabel('\tau ')
% xlim([min(SetPointPos.Data(:,1))-abs(max(SetPointPos.Data(:,1))*0.3) max(SetPointPos.Data(:,1))+abs(max(SetPointPos.Data(:,1))*0.3)])
% ylim([min(SetPointPos.Data(:,2))-abs(max(SetPointPos.Data(:,2))*0.3) max(SetPointPos.Data(:,2))+abs(max(SetPointPos.Data(:,2))*0.3)])

linkaxes(linkx,'x');
