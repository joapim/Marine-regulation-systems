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
Current_input_type = 0;         % 0 = constant heading, 1 = linear varied heading

CurrHead.Slope = 0.01;          %[rad/s] slop of the linear varied heading
CurrHead.TimeStart = 100;       %[s] start at which the linear varied heading starts
% current heading given aaccording N-W coodinates from whee the current is
% flowing. An angle of pi() is added in the code to uniform it with the
% reference sstem
CurrHead.HeadStart = 0;         %[rad] initial heading value (alway to be defined, it 
                                % is used also for he case of "consat
                                % heading"
CurrHead.HeadEnd = pi;          %[rad] final heading value


%% CONTROLLER INPUT
Kp_surge = 200000;
Ki_surge = 0;
Kd_surge = 0;

Kp_sway = 200000;
Ki_sway = 0;
Kd_sway = 0;

Kp_yaw = 2000000;
Ki_yaw = 0;
Kd_yaw = 0;



%% SET-POINTS INPUT
typeOfSimulation = 1; % 0=constant set-point , 1 sequence of set-points

%if constant set-point
ConstSPSuge = 1;  %m
ConstSPSway = 0;  %m
ConstSPYaw = 0;   %rad

%if sequence of set-points
n0 = [0 0 0];
n1 = [50 0 0];
n2 = [50 50 0];
n3 = [50 50 pi/4];
n4 = [0 50 pi/4];
n5 = [0 0 0];
nEnd = [0 0 0];
TimeSteady = 200;      %seconds , time for which the system has to stay  steady
TimeTansition = 300;    %seconds , transidition time between one set-point to the next


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
subplot(3,1,1)
plot(SetPointPos.Time,SetPointPos.Data(:,1))
hold on
plot(Eta.Time,Eta.Data(:,1))
hold off
grid on
ylabel('X-pos [m]')
legend('Set-point','Vessel')
title('Vessel position - Earth-Ref. Frame')
subplot(3,1,2)
plot(SetPointPos.Time,SetPointPos.Data(:,2))
hold on
plot(Eta.Time,Eta.Data(:,2))
hold off
grid on
ylabel('Y-pos [m]')
legend('Set-point','Vessel')
subplot(3,1,3)
plot(SetPointPos.Time,SetPointPos.Data(:,3))
hold on
plot(Eta.Time,Eta.Data(:,3))
hold off
grid on
legend('Set-point','Vessel')
ylabel('Time [s]')
ylabel('Yaw [rad]')

figure()
plot(SetPointPos.Data(:,1),SetPointPos.Data(:,2));
hold on
plot(Eta.Data(:,1),Eta.Data(:,2))
hold off
grid on
legend('Set-point','Vessel')
xlabel('X-pos [m] Earth-Ref. Frame')
ylabel('Y-pos [m] Earth-Ref. Frame')
% xlim([min(SetPointPos.Data(:,1))-abs(max(SetPointPos.Data(:,1))*0.3) max(SetPointPos.Data(:,1))+abs(max(SetPointPos.Data(:,1))*0.3)])
% ylim([min(SetPointPos.Data(:,2))-abs(max(SetPointPos.Data(:,2))*0.3) max(SetPointPos.Data(:,2))+abs(max(SetPointPos.Data(:,2))*0.3)])
%axis equal

