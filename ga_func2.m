function [avg_cost] = ga_func2(x)
% Function finds the average cost in a day but without ability to sell
% Total time (T) and tasks (N)
T = 24; N = 13;
% Input: x is an array of k, b and s
k = x(1:T); % Proportion of solar power siphoned to battery
b = x(T+1:T*2); % Battery power bought or used for tasks
s = x(T*2+1:T*2+N); % Delay time for tasks
% Output: avg_cost is the average cost in a day

% n{t) is a set of tasks arriving in time t
n = {[9 10 11 12], [], [7], [5], [], [], [], [], [], [3 6 13],...
    8, [1 2 4], [], [], [], [], [], [], [], [], [], [], [], []};

% Using variables from algo.m
global p; global r; global a; 
global S; global C;
global Bga2;

% Loop to go through all time slots
for t = 1:T
    % Loop to go through tasks in time slot
    for i = n{t}
        % Only for time slots which have tasks arriving
        if i ~= 0
            % If task takes longer than 1 time slot to complete
            if r(i) > 1
                % Running through all the time slots in task i 
                for j = 1:r(i)-1
                    t1 = a(i) + j + s(i);
                    cost_j(j) = p(i).*C(t1) - b(t1).*C(t1);
                end
                t2 = a(i) + s(i);
                % Adding the first time slot for task i
                cost_i(i) = sum(cost_j) + p(i).*C(t2) - b(t2).*C(t2);
            else
                % If task takes only 1 time slot
                t3 = a(i) + s(i);
                cost_i(i) = p(i).*C(t3) - b(t3).*C(t3);
            end
        end
    end
    % Battery level at end of 1 time slot
    Bga2(t+1) = Bga2(t) + k(t).*S(t) - b(t);
end
% Average cost taken as sum of all costs for each time slot over total time
avg_cost = sum(cost_i)/T;

end

