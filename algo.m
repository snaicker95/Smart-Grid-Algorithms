% algo.m finds the global minima of fuction ga_func.m & ga_func2.m
clc; clear all; close all;
% Task parameters
global p; global r; global a; global d;
p = xlsread ('ProjectData.xlsx','O4:O16'); % Power rating
r = xlsread ('ProjectData.xlsx','P4:P16'); % Duration
a = xlsread ('ProjectData.xlsx','R4:R16')-7; % Arrival time
d = xlsread ('ProjectData.xlsx','Q4:Q16'); % Deadline

% Input variables
global S; global C; global G;
S = xlsread ('ProjectData.xlsx','I4:I27'); % Solar power generation
C = xlsread ('ProjectData.xlsx','J4:J27'); % Tariffs for buying from grid
G = xlsread ('ProjectData.xlsx','K4:K27'); % Tariffs for selling to grid
global Bga;
Bga = zeros(1,25); % Battery charge level
Bga(1) = 6; % Initialise battery capacity as half of 12kWh battery
global Bga2; % (not selling)
Bga2 = zeros(1,25); % Battery charge level
Bga2(1) = 6; % Initialise battery capacity as half of 12kWh battery

starttime = tic; % Start stopwatch
% Bounds
s_ub = d-r; % Upper bound for s
b_max = 5; % Max discharge/charge
% Total
T = 24; % Time slots (24 hours)
N = 13; % Number of tasks

% Genetic algorithm input
nvars = 2*T+N; % Number of variables
% Inequalities
Aueq = [-diag(S),eye(T),zeros(T,N);diag(S),-eye(T),zeros(T,N)];
bueq = [Bga(1:end-1)';Bga(2:end)'];
% Inequalities (not selling)
Aueq2 = [-diag(S),eye(T),zeros(T,N);diag(S),-eye(T),zeros(T,N)];
bueq2 = [Bga2(1:end-1)';Bga2(2:end)'];
% Equalities
Aeq = [];
beq = [];
% Linear contraints
lb = [0*ones(1,T),-b_max*ones(1,T),0*ones(1,N)];
ub = [1*ones(1,T),b_max*ones(1,T),s_ub'.*ones(1,N)];
% Non-linear contraints
nonlcon = [];
% Integer
IntCon = T*2+1:T*2+N;
% Setting the genetic algorithm
options = gaoptimset;
rng default

% GA FUNCTION
% Output:
% x:            Solution
% fval:         Objective function value using x
% exitflag:     Reason algorithm stops
% output:       Information about the process
% population:   Final population
% scores:       Final scores
[x,fval,exitflag,output,population,scores] = ga(@ga_func,2*T+N,Aueq,bueq,Aeq,beq,lb,ub,nonlcon,IntCon,options);

% Decision variables
% x is split into the different decision variables
k = x(1:T); % Proportion of solar power siphoned to battery
b = x(T+1:T*2); % Battery power bought or used for tasks
s = x(T*2+1:end); % Delay time for tasks

% Not selling
[x2,fval2,exitflag2,output2,population2,scores2] = ga(@ga_func2,2*T+N,Aueq2,bueq2,Aeq,beq,lb,ub,nonlcon,IntCon,options);
% Decision variables
k2 = x2(1:T);
b2 = x2(T+1:T*2);
s2 = x2(T*2+1:end);

finishtime = toc(starttime); % Stop stopwatch

