%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% init()                                                                  %
%                                                                         %              
% Set initial parameters for part1.slx and part2.slx                      %
%                                                                         %
% Created:      2018.07.12	Jon Bjørnø                                    %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all;

load('supply.mat');
load('supplyABC.mat');
load('thrusters_sup.mat')

% Initial position x, y, z, phi, theta, psi
eta0 = [0,0,0,0,0,0]';
% Initial velocity u, v, w, p, q, r
nu0 = [0,0,0,0,0,0]';




% CURRENT INPUT
Current_Speed = 0.5;            % m/s
Current_input_type = 1;         % 0 = constant heading, 1 = linear varied heading

CurrHead.Slope = 0.01;          %[rad/s] slop of the linear varied heading
CurrHead.TimeStart = 100;       %[s] start at which the linear varied heading starts
% current heading given aaccording N-W coodinates from whee the current is
% flowing. An angle of pi() is added in the code to uniform it with the
% reference sstem
CurrHead.HeadStart = 0;         %[rad] initial heading value (alway to be defined, it 
                                % is used also for he case of "consat
                                % heading"
CurrHead.HeadEnd = pi;          %[rad] final heading value





% SIMULATION
sim('part1')