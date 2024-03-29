function [prob] = mc_safe_prob_nonlinear(x_t, h, sigma)
% compute safe probability of state x in time horizon h using Monte Carlo
% safe set is h(x) = x-1 >= 0
Nt = 100;
x = zeros(Nt,1);
loop_num = 10000;
dt = 0.1;
A = 2;
K = 2.5;
K_trap = 5;
bnd = 1.5; % system dynamics switching boundary
safe_bnd = 1; % boundary of safe set

% store the first exit time each loop
first_exit_time = (Nt+2)*ones(loop_num,1);

for loop_count = 1:loop_num
    x = zeros(Nt,1);
    x(1) = x_t;
    for i = 1:(Nt-1)
        if x(i) > bnd
            x(i+1) = exp((A-K)*dt) * x(i) + randn*sigma;
        else
            x(i+1) = exp((A-K_trap)*dt) * x(i) + randn*sigma;
        end
        if x(i+1) <= safe_bnd % unsafe condition
            first_exit_time(loop_count) = i+1;
            break
        end
    end
end

% first_exit_time
prob = sum(first_exit_time > h) / size(first_exit_time, 1); % safe probability for fixed time horizon

end

