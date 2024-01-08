function [prob, value_func] = mc_safe_prob_K(x_t, h, dt, sigma, K)
% compute safe probability of state x in time horizon h using Monte Carlo
% safe set is h(x) = x-1 >= 0
Nt = h+2;
x = zeros(Nt,1);
loop_num = 100000;
% dt = 0.1;
% dt = 0.05;
A = 2;
bnd = 1; % boundary of safe set

% store the first exit time each loop
first_exit_time = (Nt+2)*ones(loop_num,1);

% test continuous function
value_func = zeros(loop_num, 1);

x_1 = zeros(loop_num, 1);

for loop_count = 1:loop_num
    x = zeros(Nt,1);
    exit = false;
    x(1) = x_t;
    for i = 1:(Nt-1)
%         x(i+1) = exp((A-K)*dt) * x(i) + (exp(2*(A-K)*dt)-1)/(2*(A-K))*randn*sigma;
        x(i+1) = exp((A-K)*dt) * x(i) + randn*sigma;
%         x(i+1) = x(i) + (A-K)*x(i)*dt + randn*sigma;
        if i == 1
            x_1(loop_count, 1) = x(i+1);
        end
        if x(i+1) <= bnd && exit == false % unsafe condition
            first_exit_time(loop_count) = i+1;
%             break
            exit = true;
        end
    end
    value_func(loop_count) = x(3);
end

% first_exit_time
prob = sum(first_exit_time > h) / size(first_exit_time, 1); % safe probability for fixed time horizon

value_func = sum(value_func) / loop_num;

end

