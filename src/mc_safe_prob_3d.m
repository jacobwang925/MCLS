function [prob] = mc_safe_prob_3d(x_t, h, sigma)
% compute safe probability of state x in time horizon h using Monte Carlo
% safe set is h(x) = x-1 >= 0
Nt = h;
loop_num = 10000;
dt = 0.05;
K = 2.5;
bnd = 1; % boundary of safe set
A = [2, 1.9, 1.8];

% store the first exit time each loop
first_exit_time = (Nt+2)*ones(loop_num,1);

for loop_count = 1:loop_num
    x = zeros(Nt,3);
    x(1, :) = x_t;
    for i = 1:(Nt-1)
        x(i+1, 1) = exp((A(1)-K)*dt) * x(i, 1) + randn*sigma;
        x(i+1, 2) = exp((A(2)-K)*dt) * x(i, 2) + randn*sigma;
        x(i+1, 3) = exp((A(3)-K)*dt) * x(i, 3) + randn*sigma;
        if x(i+1, 1) <= bnd || x(i+1, 2) <= bnd || x(i+1, 3) <= bnd % unsafe condition
            first_exit_time(loop_count) = i+1;
            break
        end
    end
end

% first_exit_time
prob = sum(first_exit_time > h) / size(first_exit_time, 1); % safe probability for fixed time horizon
end