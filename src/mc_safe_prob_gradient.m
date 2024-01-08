function [prob_grad, prob, safe_num, all_num] = mc_safe_prob_gradient(x_t, h, dt, sigma, x_next, df_x)
% compute safe probability of state x in time horizon h using Monte Carlo
% safe set is h(x) = x-1 >= 0
Nt = h+2;
x = zeros(Nt,1);
loop_num = 100000;
% dt = 0.1;
A = 2;
K = 2.5;
bnd = 1; % boundary of safe set

% store the first exit time each loop
first_exit_time = (Nt+2)*ones(loop_num,1);

sigma_cont = sigma/sqrt(dt);

x_1 = zeros(loop_num, 1);

for loop_count = 1:loop_num
    x = zeros(Nt,1);
    x(1) = x_t;
    for i = 1:(Nt-1)
%         x(i+1) = exp((A-K)*dt) * x(i) + (exp(2*(A-K)*dt)-1)/(2*(A-K))*randn*sigma;
        x(i+1) = exp((A-K)*dt) * x(i) + randn*sigma;
        if i == 1
            x_1(loop_count, 1) = x(i+1);
        end
        if x(i+1) <= bnd % unsafe condition
            first_exit_time(loop_count) = i+1;
            break
        end
    end
end

% first_exit_time
prob = sum(first_exit_time > h) / size(first_exit_time, 1); % safe probability for fixed time horizon
safe_num = sum(first_exit_time > h);
all_num = size(first_exit_time, 1);

prob_grad = (first_exit_time > h).*(x_1-x_next)/(sigma^2)*df_x;
% prob_grad = (first_exit_time > h).*(x_1-x_next)/sigma_cont*df_x;

prob_grad = sum(prob_grad) / size(first_exit_time, 1);

end

