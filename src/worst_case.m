% worst case control

clear; clc;
Nt = 100; % time points
dt = 0.1; % time step size
traj_num = 100; % number of trajectories
dx = 0.1; % state step size (for derivative calculations)
eps = 0.1; % epsilon
K = 2.5; % proportional controller gain
A = 2; % A matrix in system dynamics
sigma = 1; % magnitude of noise
h = 10; % safe prob time horizon
x_0 = 3; % initial state

sigma = sigma * sqrt(dt); % discretization

%% Proposed controller
x = zeros(Nt, traj_num); % initialization
x(1, :) = x_0;
F = zeros(Nt, 1); % to store expected value of safe probability

alpha = 1;

for j = 1:traj_num
    for i = 1:Nt
        P = mc_safe_prob(x(i, j), h, sigma);
        P_dx1 = mc_safe_prob(x(i, j)+dx, h, sigma);
        P_dx2 = mc_safe_prob(x(i, j)-dx, h, sigma);
        dP_x = (P_dx1 - P_dx2) / (2*dx); % gradient
        if abs(dP_x) < 0.05
            dP_x = 0.05;
        end
        u = (-alpha * (P - (1-eps))) / dP_x - A*x(i, j);
        x(i+1, j) = exp((A)*dt) * x(i, j) + u*(exp(A*dt)-1)/A + randn*sigma;
    end
end
x_ave = mean(x, 2);


%% Proportional controller
x_r = zeros(Nt, traj_num); % initialization
x_r(1, :) = x_0;
F_r = zeros(Nt, 1); % to store expected value of safe probability

for j = 1:traj_num
    for i = 1:Nt
        x_r(i+1, j) = exp((A-K)*dt) * x_r(i, j) + randn*sigma; % nominal controller
    end
end

x_ave_r = mean(x_r, 2);


%% StoCBF (Clark, 2019)
x_clark = zeros(Nt, traj_num); % initialization
x_clark(1, :) = x_0;
F_clark = zeros(Nt, 1); % to store expected value of safe probability

alpha = 1;

for j = 1:traj_num
    for i = 1:Nt
        u = ((-alpha*dt + 1 - exp(A*dt)) * x_clark(i, j) + alpha*dt) * A / (exp(A*dt)-1);
        x_clark(i+1, j) = exp((A)*dt) * x_clark(i, j) + randn*sigma + u*(exp(A*dt)-1)/A; % zero-hold control
    end
end

x_ave_clark = mean(x_clark, 2);


%% CVaR (Ahmadi et al, 2020)
x_cvar = zeros(Nt, traj_num); % initialization
x_cvar(1, :) = x_0;
u_cvar = zeros(Nt, traj_num);
F_cvar = zeros(Nt, 1); % to store expected value of safe probability

beta_risk = 0.1;
sigma_risk = sqrt(dt);
alpha_risk = 0.65;

epsilon = norminv(beta_risk, 0, sigma_risk);
l_risk = quadgk(@(x) x.*normpdf(x, 0, sigma_risk), -inf, epsilon)/normcdf(epsilon, 0, sigma_risk); % expectation of h(x) given h(x) < epsilon

for j = 1:traj_num
    for i = 1:Nt
        u_cvar(i, j) = ((alpha_risk - exp((A)*dt)) * x_cvar(i, j) - l_risk) * A / (exp(A*dt)-1);
        x_cvar(i+1, j) = exp((A)*dt) * x_cvar(i, j) + randn*sigma + u_cvar(i, j)*(exp(A*dt)-1)/A; % zero-hold control
    end
end

x_ave_cvar = mean(x_cvar, 2);


%% PrSBC (Luo et al, 2019)
x_prsbc = zeros(Nt, traj_num); % initialization
x_prsbc(1, :) = x_0;
u_prsbc = zeros(Nt, traj_num);
F_prsbc = zeros(Nt, 1); % to store expected value of safe probability

alpha_risk = 1;
sigma_risk = sqrt(dt);
epsilon = 0.9;
l_risk = norminv(1-epsilon, 0, sigma_risk);

for j = 1:traj_num
    for i = 1:Nt
        u_prsbc(i, j) = ((- alpha_risk * x_prsbc(i, j) - l_risk + alpha_risk) * dt + x_prsbc(i, j) - exp((A)*dt) * x_prsbc(i, j)) * A / (exp(A*dt)-1);
        x_prsbc(i+1, j) = exp((A)*dt) * x_prsbc(i, j) + randn*sigma + u_prsbc(i, j)*(exp(A*dt)-1)/A; % zero-hold control
    end
end

x_ave_prsbc = mean(x_prsbc, 2);


%% Plots

% safe probability
for i = 1:Nt
    F(i) = mc_safe_prob(x_ave(i), h, sigma);
    F_r(i) = mc_safe_prob(x_ave_r(i), h, sigma);
    F_clark(i) = mc_safe_prob(x_ave_clark(i), h, sigma);
    F_cvar(i) = mc_safe_prob(x_ave_cvar(i), h, sigma);
    F_prsbc(i) = mc_safe_prob(x_ave_prsbc(i), h, sigma);
end

figure
plot(F, 'linewidth', 1.5)
hold on
plot(0)
plot(F_clark, 'linewidth', 1.5)
plot(F_prsbc, 'linewidth', 1.5)
plot(F_cvar, 'linewidth', 1.5)
legend('Proposed controller', '', 'StoCBF', 'PrSBC', 'CVaR')
set(gca, 'FontSize', 16)
set(gca, 'FontSize', 19)
legend('Location', 'best')
xlabel('Time')
ylabel('Expected value of safe probability')
xt = get(gca, 'XTick'); 
set(gca, 'XTick', xt, 'XTickLabel', xt/10)
set(gcf, 'position', [200 200 600 469])

% averaged state
figure
plot(x_ave, 'linewidth', 1.5)
hold on
plot(0)
plot(x_ave_clark, 'linewidth', 1.5)
plot(x_ave_prsbc, 'linewidth', 1.5)
plot(x_ave_cvar, 'linewidth', 1.5)
yline(1, 'LineStyle', '--', 'color', 'red', 'linewidth', 1.5)
legend('Proposed controller', '', 'StoCBF', 'PrSBC', 'CVaR')
set(gca, 'FontSize', 16)
set(gca, 'FontSize', 19)
legend('Location', 'best')
xlabel('Time')
ylabel('Averaged state over ' + string(traj_num) + ' trajectories')
xlim([0,100])
xt = get(gca, 'XTick'); 
set(gca, 'XTick', xt, 'XTickLabel', xt/10)
set(gcf, 'position', [200 200 600 469])

save('worst_case.mat')

