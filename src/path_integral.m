clear; clc;
Nt = 50; % time points
dt = 0.1; % time step size
% dt = 0.05;
traj_num = 50; % number of trajectories
dx = 0.1; % state step size (for derivative calculations)
eps = 0.1; % epsilon
K = 2.5; % proportional controller gain
A = 2; % system dynamics f(x) = A
sigma = 1; % magnitude of noise
h = 10; % safe prob time horizon
bnd = 1; % safe when x > 1
Nx = 40; 

% W \sim N(0, sigma^2 dt)
sigma = sigma*sqrt(dt); % equivalent magnitude with discretized system

% x_0 = 1.5; % initial state
% x_next = exp((A-K)*dt) * x_0; % expected value for x_1
% df_x = exp((A-K)*dt); % gradient of discretized system of the deterministic dynamics f_cl(x) = (A-K)x
% 
% prob_grad = mc_safe_prob_gradient(x_0, h, sigma, x_next, df_x)

safe_prob = zeros(Nx,1);
prob_grad = zeros(Nx,1);
x_init = zeros(Nx, 1);

for i = 1:Nx
    x_0 = dx*i; % initial state
    x_init(i) = x_0;
    x_next = exp((A-K)*dt) * x_0; % expected value for x_1
    df_x = exp((A-K)*dt); % gradient of discretized system of the deterministic dynamics f_cl(x) = (A-K)x
    
    [prob_grad(i), safe_prob(i)] = mc_safe_prob_gradient(x_0, h, dt, sigma, x_next, df_x);
end

prob_grad_FD = zeros(Nx,1); % gradient from finite difference
for i = 2:Nx-1
    prob_grad_FD(i) = (safe_prob(i+1) - safe_prob(i-1)) / (2*dx);
end

figure

plot(x_init, prob_grad_FD)
hold on
plot(x_init, prob_grad)
legend('Finite difference', 'Proposed computation')
title('estimated probability gradient')
xlabel('$x_0$', 'Interpreter','latex')
set(gca, 'FontSize', 19)

figure

plot(x_init, safe_prob)
title('safety probability')
xlabel('$x_0$', 'Interpreter','latex')
set(gca, 'FontSize', 19)

%% importance sampling
clear; clc;
Nt = 50; % time points
dt = 0.1; % time step size
dt = 0.05;
traj_num = 50; % number of trajectories
dx = 0.1; % state step size (for derivative calculations)
eps = 0.1; % epsilon
K1 = 3; % proportional controller gain
A = 2; % system dynamics f(x) = A
sigma = 0.3; % magnitude of noise
h = 20; % safe prob time horizon
bnd = 1; % safe when x > 1
Nx = 40; 


K2 = 2.5; % controller for importance sampling

sigma = sigma*sqrt(dt); % equivalent magnitude with discretized system

safe_prob_K1 = zeros(Nx,1);
safe_prob_K2 = zeros(Nx,1);
safe_prob_IS = zeros(Nx,1); % estimated safe prob with importance sampling
IS_weight = zeros(Nx,1); % averaged IS weight for each initial state
prob_grad = zeros(Nx,1);
x_init = zeros(Nx, 1);

vaule_func_K1 = zeros(Nx,1);
vaule_func_K2 = zeros(Nx,1);
vaule_func_IS = zeros(Nx,1);

for i = 1:Nx
    x_0 = dx*i; % initial state
    x_init(i) = x_0;
    
    [safe_prob_K1(i), vaule_func_K1(i)] = mc_safe_prob_K(x_0, h, dt, sigma, K1);
    [safe_prob_K2(i), vaule_func_K2(i)] = mc_safe_prob_K(x_0, h, dt, sigma, K2);
    [safe_prob_IS(i), IS_weight(i), vaule_func_IS(i)] = mc_safe_prob_IS(x_0, h, dt, sigma, K1, K2);
end

figure
plot(x_init, safe_prob_K1)
hold on
plot(x_init, safe_prob_K2)
plot(x_init, safe_prob_IS)
legend('Data sample on controller K1', 'Direct estimation for K2', 'Importance sampling')
title('safety probability')
xlabel('$x_0$', 'Interpreter','latex')
set(gca, 'FontSize', 19)

figure
plot(x_init, IS_weight)
title('averaged IS weight')
xlabel('$x_0$', 'Interpreter','latex')
set(gca, 'FontSize', 19)

figure
plot(x_init, vaule_func_K1)
hold on
plot(x_init, vaule_func_K2)
plot(x_init, vaule_func_IS)
% legend('Data on K1', 'Naive estimation', 'Importance sampling')
legend('Data sample on controller K1', 'Direct estimation for K2', 'Importance sampling')
title('continuous function')
xlabel('$x_0$', 'Interpreter','latex')
set(gca, 'FontSize', 19)



