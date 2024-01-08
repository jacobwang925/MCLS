function [prob, weight, value_func] = mc_safe_prob_IS(x_t, h, dt, sigma, K1, K2)
% compute safe probability of state x in time horizon h using Monte Carlo
% safe set is h(x) = x-1 >= 0
% using importance sampling
% K1 is behavioural controller
% K2 is target controller
Nt = h+1;
x = zeros(Nt,1);
loop_num = 100000;
% dt = 0.1;
% dt = 0.05;
A = 2;
bnd = 1; % boundary of safe set

sigma_cont = sigma/sqrt(dt);

% one from reformulation, one from discretization
% this is wrong
% sigma_cont = sigma_cont^2;


% sigma_cont = 3.3;

% store the first exit time each loop
first_exit_time = (Nt+2)*ones(loop_num,1);

% test continuous function
value_func = zeros(loop_num, 1);

x_1 = zeros(loop_num, 1);
IS_weight = zeros(loop_num, 1);

for loop_count = 1:loop_num
    x = zeros(Nt,1);
    exit = false; 
    noise = zeros(Nt, 1);
    x(1) = x_t;
    for i = 1:(Nt-1)
%         x(i+1) = exp((A-K)*dt) * x(i) + (exp(2*(A-K)*dt)-1)/(2*(A-K))*randn*sigma;
        noise(i) = randn*sigma;
        x(i+1) = exp((A-K1)*dt) * x(i) + noise(i);
%         x(i+1) = x(i) + (A-K1)*x(i)*dt + noise(i);

        if x(i+1) <= bnd && exit == false % unsafe condition
            first_exit_time(loop_count) = i+1;
            exit = true;
        end
    end
    % continuous function just equals to the value of x(3)
    value_func(loop_count) = x(3);
    % Girsanov theorem
    % TODO: check this
%     x(1:end-1)
%     noise(1:end-1)
%     x(1:end-1)'*noise(1:end-1)

    u = -(K1-K2)*x(1:end-1);
%     disp(-sum(u.^2)*dt/(sigma_cont)^2)
%     disp(sum(u'*noise(1:end-1))/sigma_cont)
%     IS_weight(loop_count) = sum(u.^2)*dt/(sigma_cont)^2 + sum(u'*noise(1:end-1))/sigma_cont; % importance sampling weight
%     IS_weight(loop_count) = -sum(u.^2)*dt/(sigma_cont)^2 - sum(u'*noise(1:end-1))*sqrt(dt)/sigma_cont/sigma; % IS according to wiki
    IS_weight(loop_count) = -sum(u.^2)*dt/2/(sigma_cont)^2 - sum(u'*noise(1:end-1))/sigma_cont;
    % in cross term, one sigma from control, another from noise
    IS_weight(loop_count) = -sum(u.^2)*dt/2/(sigma_cont)^2 - sum(u'*noise(1:end-1))/sigma_cont^2;
%     IS_weight(loop_count) = -sum(u.^2)*dt/2/(sigma_cont)^2 - sum(u'*noise(1:end-1))/sigma_cont/sigma;
%     IS_weight(loop_count) = -sum(u.^2)*dt/2 - sum(u'*noise(1:end-1));
%     IS_weight(loop_count) = -sum(u.^2)*dt/2/(sigma)^2 - sum(u'*noise(1:end-1))/sigma;
%     u
%     x
%     noise
%     sum(u'*noise(1:end-1))
    %     disp(IS_weight(loop_count))
    IS_weight(loop_count) = exp(IS_weight(loop_count));
end

disp(mean(IS_weight))

weight = mean(IS_weight);


% first_exit_time
weighted_prob = (first_exit_time > h).*IS_weight;

prob = sum(weighted_prob) / size(first_exit_time, 1);
prob = prob / weight; % regularize the weight or not

% prob = sum(weighted_prob) / sum(IS_weight); % equivalent expression

value_func = sum(value_func.*IS_weight) / sum(IS_weight);

end

