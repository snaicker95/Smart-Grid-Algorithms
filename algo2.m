% algo2.m reduces the average cost according to the tarrifs 
% algo2.m plots the results
% Input variables
global S; global C; global G; global p;

starttime2 = tic; % Start stopwatch
C_sw = [C(:)',C(1)]; % Needs to loop back to first tariff

T = 24; % Total time slots (24 hours)
% n_sw{t) is a set of tasks in time t
n_sw = {[9 10 11 12], [9 10 11 12], [7 9 10 11 12], [5 7 9 11 12],...
    [7 9 11 12], [9 11 12], [11 12], [11 12], [11 12], [3 6 13 11 12],...
    [6 8 13 11 12], [1 2 4 6 8 13 11 12], [6 8 13 11 12], [6 8 13 11 12],...
    [8 13 11 12], [13 11 12], [11 12], [11 12], [11 12], [11 12], [11 12],...
    [11 12], [11 12], [11 12],};

% Getting the array power consumption in time t (W(t))
for t = 1:T
    w = [];
    for i = n_sw{t}      
        w(i) = (p(i));
    end
    W(t) = sum(w);
end

% The decision making
b_max = 5; % Max discharge/charge
% Initialising the variables
k_sw = zeros(1,T+1); % Proportion of solar power siphoned to battery
b_sw = zeros(1,T); % Battery power bought or used for tasks
B = zeros(1,T+1); % Battery charge level
k_sw(1) = 0.5; B(1) = 6; b_sw(1) = 0.5*b_max; 
% Loop to go through all time slots
for t = 1:T
    % When next time slot has higher buying tariff
    if (C_sw(t+1) > C_sw(t))
        % -b(t) increases proportionally to C(t+1)
        % Buy power to store
        b_sw(t+1) = ((C_sw(t+1))./(C_sw(t)))*(b_max/2);
        % Solar siphoned to battery increases proportionally to C(t+1)
        k_sw(t+1) = k_sw(t).*(C_sw(t))./C_sw(t+1);
        % Battery level at end of 1 time slot
        B(t+1) = B(t) + k_sw(t).*S(t) - b_sw(t);
    % When next time slot has lower buying tarrif    
    elseif (C_sw(t+1) < C_sw(t))
        % b(t) increases proportionally to C(t+1)
        % Use power from store
        b_sw(t+1) = -((C_sw(t))./(C_sw(t+1)))*(b_max/2);
        % Solar siphoned to battery decreases proportionally to C(t+1)
        k_sw(t+1) = k_sw(t).*(C_sw(t))./C_sw(t+1);
        % Battery level at end of 1 time slot
        B(t+1) = B(t) + k_sw(t).*S(t) - b_sw(t);
    % When next time slot has no change in tariff    
    else
        % Maintain all variables
        b_sw(t+1) = b_sw(t);
        k_sw(t+1) = k_sw(t);
        B(t+1) = B(t);
    end
end

% Getting the total costs over all time slots
for t = 1:T
    cost_sw(t) = W(t).*C_sw(t) - b_sw(t).*C_sw(t) - (1-k_sw(t)).*S(t).*G(t);
end
% Average cost
avg_cost = sum(cost_sw)/T;

% Getting the total costs over all time slots (not selling)
for t = 1:T
    cost_sw2(t) = W(t).*C_sw(t) - b_sw(t).*C_sw(t);
end
% Average cost
avg_cost2 = sum(cost_sw2)/T;

% Getting the total costs over all time slots (controlled)
for t = 1:T
    cost_sw3(t) = W(t).*C_sw(t);
end
% Average cost
avg_cost3 = sum(cost_sw3)/T;

finishtime2 = toc(starttime2); % Stop stopwatch

% ANALYSING DATA
% Data from genetic algorithm
% % n_ga{t) is a set of tasks in time t with delay
nga = {[10 11 12], [9 10 11 12], [9 10 11 12], [5 7 9 11 12], [7 9 11 12],...
    [7 9 11 12], [9 11 12], [11 12], [11 12], [3 13 11 12], [8 13 11 12],...
    [8 13 11 12], [8 13 11 12], [6 8 13 11 12], [2 6 8 13 11 12],...
    [4 6 13 11 12], [1 6 11 12], [6 11 12], [11 12], [11 12], [11 12],...
    [11 12], [11 12], [11 12]};

% Array power consumption in time t (W_ga(t))
for t = 1:T
    wga = [];
    for i = nga{t}      
        wga(i) = (p(i));
    end
    Wga(t) = sum(wga);
end

% Net Power Consumption
% Removing Negative Values
for t = 1:T
    if b_sw(t) < 0
        b_sw_new(t) = 0;
    else
        b_sw_new(t) = b_sw(t);
    end
    if b (t) < 0
        b_new(t) = 0;
    else
        b_new(t) = b(t);
    end
end

% Effect of Battery
W_net = W - b_sw_new;
Wga_net = Wga - b_new;
% minW_net = min(W_net);
% minWga_net = min(Wga_net);
% W_netOff = W_net - minW_net*ones(1,T);
% Wga_netOff = Wga_net - minWga_net*ones(1,T);
% Offset
for t = 1:T
    if W_net(t) < 0
        W_netOff(t) = 0;
    else
        W_netOff(t) = W_net(t);
    end
    if Wga_net(t) < 0
        Wga_netOff(t) = 0;
    else
        Wga_netOff(t) = Wga_net(t);
    end 
end

% Peak-to-average Ratios (PAR)
PAR1 = max(W_netOff)/mean(W_netOff); % Tariff-proportional Algorithm
PAR2 = max(Wga_netOff)/mean(Wga_netOff); % Genetic Algorithm

% Direct Effect of Solar
% Getting power profiles
S_de = W-S';
S_de2 = Wga-S';
% minS_de = min(S_de);
% minS_de2 = min(S_de2);
% S_denew = S_de - minS_de*ones(1,T);
% S_denew2 = S_de2 - minS_de2*ones(1,T);
% Offset
for t = 1:T
    if S_de(t) < 0
        S_denew(t) = 0;
    else
        S_denew(t) = S_de(t);
    end
    if S_de2(t) < 0
        S_denew2(t) = 0;
    else
        S_denew2(t) = S_de2(t);
    end
end

% PAR with direct solar
PARS_de = max(S_denew)/mean(S_denew); % Tariff-Proportional Algorithm
PARS_de2 = max(S_denew2)/mean(S_denew2); % Genetic Algorithm

% PAR with task scheduling
PARTs = max(W)/mean(W); % Tariff-Proportional Algorithm
PARTs2 = max(Wga)/mean(Wga); % Genetic Algorithm

T = 24; N = 13;
t = 1:T; % Time slots
i = 1:N; % Tasks

% Utility Function
Ureal = s.^2; % Dissatisfaction for each task
% Plotting function
sU = 0:0.5:max(s);
U = sU.^2;
figure (1)
plot(sU,U)
title ('Utility Function')
xlabel ('Task Delay, s')
ylabel ('Dissatisfaction, U')

% % Power Demand Profiles (Task Scheduling)
% figure (2)
% plot (t, W, t, Wga);
% title ('Power Demand Profiles (Task Scheduling)');
% xlabel ('Time Slots (hours)');
% ylabel ('Power Demand (kW)');
% legend ('Tariff-Proportional', 'Genetic');
% 
% % Power Demand Profiles (Battery Offset)
% figure (3)
% title ('Battery Offset')
% subplot(2,1,1)
% plot (t, W,'*-', t, W_netOff,'o-')
% title ('Tariff-Proportional (Battery Offset)')
% ylabel ('Power Demand (kW)')
% xlabel ('Time Slots (hours)')
% legend ('No Offset', 'Battery Offset')
% 
% subplot(2,1,2)
% plot (t, Wga,'*-', t, Wga_netOff,'o-')
% title ('Genetic (Battery Offset)')
% ylabel ('Power Demand (kW)')
% xlabel ('Time Slots (hours)')
% legend ('No Offset', 'Battery Offset')
% 
% % Power Demand Profiles (Solar Offset)
% figure (4)
% subplot(2,1,1)
% plot (t, W,'*-', t, S_denew,'o-')
% title ('Tariff-Proportional (Solar Offset)')
% ylabel ('Power Demand (kW)')
% xlabel ('Time Slots (hours)')
% legend ('No Offset', 'Solar Offset')
% 
% subplot(2,1,2)
% plot (t, Wga,'*-', t, S_denew2,'o-')
% title ('Genetic (Solar Offset)')
% ylabel ('Power Demand (kW)')
% xlabel ('Time Slots (hours)')
% legend ('No Offset', 'Solar Offset')

% % DECISION VALUES
% % Task Delays
% figure (5)
% bar (i, s);
% title ('Task Delays');
% xlabel ('Tasks');
% ylabel ('Delay');
% legend ('Genetic');
% 
% % Battery charged/discharged
% figure (6)
% plot (t, b_sw(1:end-1), t, b);
% title ('Battery Power Charged/Discharged');
% xlabel ('Time Slots (hours)');
% ylabel ('Battery Power Charged/Discharged (kW)');
% legend ('Tariff-Proportional', 'Genetic');
% 
% % Proportion of solar power siphoned to battery
% figure (7)
% plot (t, k_sw(1:end-1), t, k);
% title ('Proportion of solar power siphoned to battery');
% xlabel ('Time Slots (hours)');
% ylabel ('Solar Power (kW)');
% legend ('Tariff-Proportional', 'Genetic');
% 
% % Battery power levels
% figure (8)
% plot (t, B(1:end-1), t, Bga(1:end-1));
% title ('Battery Power Levels');
% xlabel ('Time Slots (hours)');
% ylabel ('Battery Power Levels (kW)');
% legend ('Tariff-Proportional', 'Genetic');

% % CONSTANT VALUES
% % Tariffs
% figure (9)
% plot (t, C, t, G);
% title ('Rates');
% xlabel ('Time Slots (hours)');
% ylabel ('Price (cents/kWh)');
% legend ('Buying','Selling');
% 
% % Solar Power Harnessed
% figure (10)
% plot (t, S, t, mean(S)*ones(1,T));
% title ('Solar Power Harnessed');
% xlabel ('Time Slots (hours)');
% ylabel ('Power Generated (kW)');
% legend ('Solar Power Harnessed', 'Mean Solar Power Harnessed');
% 
% % Relative Power & Tariff
% WRel = W./max(W);
% W_gaRel = W_ga./max(W_ga);
% CRel = C./max(C);
% GRel = G./max(G);
% % Buying 
% figure (11)
% plot (t, WRel, t, CRel);
% title ('Comparing Tariff to Power Demand');
% xlabel ('Time Slots (hours)');
% ylabel ('Relative Difference Between Tarrif & Power Demand');
% legend ('Relative Power Demand (No Task Scheduling', 'Relative Tariff for Buying');
% % Selling
% figure (12)
% plot (t, WRel, t, GRel);
% title ('Comparing Tariff to Power Demand');
% xlabel ('Time Slots (hours)');
% ylabel ('Relative Difference Between Tarrif & Power Demand');
% legend ('Relative Power Demand (No Task Scheduling', 'Relative Tariff for Selling');
% 
% % Solar vs Power
% figure (13)
% plot (t, S, t, W);
% title ('Comparing Solar Power Harnessed to Power Demand');
% xlabel ('Time Slots (hours)');
% ylabel ('Power (kW)');
% legend ('Solar Power Harnessed', 'Power Demand');

