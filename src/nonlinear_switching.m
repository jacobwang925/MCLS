clear; clc;
Nt = 100; % time points
dt = 0.1; % time step size
traj_num = 50;
dx = 0.1; % state step size (for derivative calculations)
% eps = 0.03; % epsilon
eps = 0.15;
eps = 0.1;
K = 2.5; % proportional controller gain
A = 2; % A
K_trap = 5; % trapping nonlinear dynamics equivalent gain
sigma = 2; % magnitude of noise
h = 10; % safe prob time horizon
bnd = 1.5; % system dynamics switching boundary

%% our controller
x = zeros(Nt, traj_num); % initialization
x(1, :) = 3;
u = zeros(Nt, traj_num);
F = zeros(Nt, 1); % to store expected value of safe probability

alpha = 2;
alpha = 1;

for j = 1:traj_num
    for i = 1:Nt
        P = mc_safe_prob_nonlinear(x(i, j), h, sigma);
        if x(i, j) > bnd
            if P > 1 - eps
                u(i, j) = -K * x(i, j);
                x(i+1, j) = exp((A-K)*dt) * x(i, j) + (exp(2*(A-K)*dt)-1)/(2*(A-K))*randn*sigma; % nominal controller
            else
                P_dx1 = mc_safe_prob_nonlinear(x(i, j)+dx, h, sigma);
                P_dx2 = mc_safe_prob_nonlinear(x(i, j)-dx, h, sigma);
                dP_x = (P_dx1 - P_dx2) / (2*dx); % gradient
               
                if dP_x*(A-K)*x(i, j) >= -alpha * (P - (1-eps))
                    u(i, j) = -K * x(i, j);
                    x(i+1, j) = exp((A-K)*dt) * x(i, j) + (exp(2*(A-K)*dt)-1)/(2*(A-K))*randn*sigma; % nominal controller
                else
    %                 u(i, j) = 1/dP_x-A*x(i, j);
                    u(i, j) = (-alpha * (P - (1-eps))) / dP_x - A*x(i, j);
                    x(i+1, j) = exp((A)*dt) * x(i, j) + (exp(2*(A)*dt)-1)/(2*A)*randn*sigma + u(i, j)*(exp(A*dt)-1)/A; % zero-hold control
                end
            end   
        else
            x(i+1, j) = exp((A-K_trap)*dt) * x(i, j) + (exp(2*(A-K_trap)*dt)-1)/(2*(A-K_trap))*randn*sigma; % uncontrollable dynamics
        end
    end
end
x_ave = mean(x, 2);


%% Proportional controller
x_r = zeros(Nt, traj_num); % initialization
x_r(1, :) = 3;
u_r = zeros(Nt, traj_num);
F_r = zeros(Nt, 1); % to store expected value of safe probability


for j = 1:traj_num
    for i = 1:Nt
        if x_r(i, j) > bnd
            u_r(i, j) = -K*x_r(i, j);
            x_r(i+1, j) = exp((A-K)*dt) * x_r(i, j) + (exp(2*(A-K)*dt)-1)/(2*(A-K))*randn*sigma; % nominal controller
        else
            x_r(i+1, j) = exp((A-K_trap)*dt) * x_r(i, j) + (exp(2*(A-K_trap)*dt)-1)/(2*(A-K_trap))*randn*sigma; % uncontrollable dynamics
        end 
    end
end

x_ave_r = mean(x_r, 2);


%% Clarks
x_clark = zeros(Nt, traj_num); % initialization
x_clark(1, :) = 3;
u_clark = zeros(Nt, traj_num);
F_clark = zeros(Nt, 1); % to store expected value of safe probability

alpha = 9;
alpha = 1;

for j = 1:traj_num
    for i = 1:Nt
        P = mc_safe_prob_nonlinear(x_clark(i, j), h, sigma);
        if x_clark(i, j) > bnd
            if P > 1 - eps
                u_clark = -K*x_clark(i, j);
                x_clark(i+1, j) = exp((A-K)*dt) * x_clark(i, j) + (exp(2*(A-K)*dt)-1)/(2*(A-K))*randn*sigma; % nominal controller
            else
                if exp((A-K)*dt) * x_clark(i, j) - x_clark(i, j) >= -alpha * dt * (x_clark(i, j)-1)
                    u_clark = -K*x_clark(i, j);
                    x_clark(i+1, j) = exp((A-K)*dt) * x_clark(i, j) + (exp(2*(A-K)*dt)-1)/(2*(A-K))*randn*sigma; % nominal controller
                else
                    u_clark(i, j) = ((-alpha*dt + 1 - exp(A*dt)) * x_clark(i, j) + alpha*dt) * A / (exp(A*dt)-1);
                    x_clark(i+1, j) = exp((A)*dt) * x_clark(i, j) + (exp(2*(A)*dt)-1)/(2*A)*randn*sigma + u_clark(i, j)*(exp(A*dt)-1)/A; % zero-hold control
                end
            end
        else
            x_clark(i+1, j) = exp((A-K_trap)*dt) * x_clark(i, j) + (exp(2*(A-K_trap)*dt)-1)/(2*(A-K_trap))*randn*sigma; % uncontrollable dynamics
        end     
    end
end

x_ave_clark = mean(x_clark, 2);


%% CVaR
x_cvar = zeros(Nt, traj_num); % initialization
x_cvar(1, :) = 3;
u_cvar = zeros(Nt, traj_num);
F_cvar = zeros(Nt, 1); % to store expected value of safe probability

beta_risk = 0.03;
beta_risk = 0.1;
sigma_risk = sqrt(dt);
% alpha_risk = 0.5;
alpha_risk = 0.65;
% alpha_risk = 0.7;

epsilon = norminv(beta_risk, 0, sigma_risk);
l_risk = quadgk(@(x) x.*normpdf(x, 0, sigma_risk), -inf, epsilon)/normcdf(epsilon, 0, sigma_risk); % expectation of h(x) given h(x) < epsilon


for j = 1:traj_num
    for i = 1:Nt
        P = mc_safe_prob_nonlinear(x_cvar(i, j), h, sigma);
        if x_cvar(i, j) > bnd
            if P > 1 - eps
                u_cvar(i, j) = -K*x_cvar(i, j);
                x_cvar(i+1, j) = exp((A-K)*dt) * x_cvar(i, j) + (exp(2*(A-K)*dt)-1)/(2*(A-K))*randn*sigma; % nominal controller
            else
                if exp((A-K)*dt) * x_cvar(i, j) + l_risk >= alpha_risk * x_cvar(i, j)
                    u_cvar(i, j) = -K*x_cvar(i, j);
                    x_cvar(i+1, j) = exp((A-K)*dt) * x_cvar(i, j) + (exp(2*(A-K)*dt)-1)/(2*(A-K))*randn*sigma; % nominal controller
                else
                    u_cvar(i, j) = ((alpha_risk - exp((A)*dt)) * x_cvar(i, j) - l_risk) * A / (exp(A*dt)-1);
                    x_cvar(i+1, j) = exp((A)*dt) * x_cvar(i, j) + (exp(2*(A)*dt)-1)/(2*A)*randn*sigma + u_cvar(i, j)*(exp(A*dt)-1)/A; % zero-hold control
                end
            end
        else
            x_cvar(i+1, j) = exp((A-K_trap)*dt) * x_cvar(i, j) + (exp(2*(A-K_trap)*dt)-1)/(2*(A-K_trap))*randn*sigma; % uncontrollable dynamics
        end  
    end
end

x_ave_cvar = mean(x_cvar, 2);


%% PrSBC
x_prsbc = zeros(Nt, traj_num); % initialization
x_prsbc(1, :) = 3;
u_prsbc = zeros(Nt, traj_num);
F_prsbc = zeros(Nt, 1); % to store expected value of safe probability

alpha_risk = 0.1;
alpha_risk = 1;
sigma_risk = sqrt(dt);
epsilon = 0.9;
% epsilon = 0.1;
l_risk = norminv(1-epsilon, 0, sigma_risk);

for j = 1:traj_num
    for i = 1:Nt
        P = mc_safe_prob_nonlinear(x_prsbc(i, j), h, sigma);
        if x_prsbc(i, j) > bnd
            if P > 1 - eps
                u_prsbc(i, j) = -K*x_prsbc(i, j);
                x_prsbc(i+1, j) = exp((A-K)*dt) * x_prsbc(i, j) + (exp(2*(A-K)*dt)-1)/(2*(A-K))*randn*sigma; % nominal controller
            else
                if exp((A-K)*dt) * x_prsbc(i, j) - x_prsbc(i, j) >= (- alpha_risk * x_prsbc(i, j) - l_risk + alpha_risk) * dt
                    u_prsbc(i, j) = -K*x_prsbc(i, j);
                    x_prsbc(i+1, j) = exp((A-K)*dt) * x_prsbc(i, j) + (exp(2*(A-K)*dt)-1)/(2*(A-K))*randn*sigma; % nominal controller
                else
                    u_prsbc(i, j) = ((- alpha_risk * x_prsbc(i, j) - l_risk + alpha_risk) * dt + x_prsbc(i, j) - exp((A)*dt) * x_prsbc(i, j)) * A / (exp(A*dt)-1);
                    x_prsbc(i+1, j) = exp((A)*dt) * x_prsbc(i, j) + (exp(2*(A)*dt)-1)/(2*A)*randn*sigma + u_prsbc(i, j)*(exp(A*dt)-1)/A; % zero-hold control
                end
            end
        else
            x_prsbc(i+1, j) = exp((A-K_trap)*dt) * x_prsbc(i, j) + (exp(2*(A-K_trap)*dt)-1)/(2*(A-K_trap))*randn*sigma; % uncontrollable dynamics
        end  
    end
end

x_ave_prsbc = mean(x_prsbc, 2);


for i = 1:Nt
    F(i) = mc_safe_prob_nonlinear(x_ave(i), h, sigma);
    F_r(i) = mc_safe_prob_nonlinear(x_ave_r(i), h, sigma);
    F_clark(i) = mc_safe_prob_nonlinear(x_ave_clark(i), h, sigma);
    F_cvar(i) = mc_safe_prob_nonlinear(x_ave_cvar(i), h, sigma);
    F_prsbc(i) = mc_safe_prob_nonlinear(x_ave_prsbc(i), h, sigma);
end

figure
plot(F)
hold on
plot(0)
plot(F_r)
plot(F_clark)
plot(F_cvar)
plot(F_prsbc)
legend('proposed controller', '', 'StoCBF', 'PrSBC', 'CVaR', 'proportional controller')
set(gca, 'FontSize', 16)
legend('Location', 'best')
xlabel('timestep')
ylabel('expected value of safe probability')


x_all = zeros(5, Nt+1, traj_num);
x_all(1,:,:) = x;
x_all(2,:,:) = x_clark;
x_all(3,:,:) = x_prsbc;
x_all(4,:,:) = x_cvar;
x_all(5,:,:) = x_r;
colors = [0, 0.4470, 0.7410; 0.9290, 0.6940, 0.1250; 0.4940, 0.1840, 0.5560; 0.4660, 0.6740, 0.1880; 0.5, 0.5, 0.5];

figure
for i = 1:5
    stdshade(squeeze(x_all(i,:,:))', 0.3, colors(i,:,:));
    hold on
end
yline(1, 'LineStyle', '--', 'color', 'red', 'LineWidth', 1.5)
legend('', 'Proposed controller', '', 'StoCBF', '', 'PrSBC', '', 'CVaR', '', 'Nominal controller')
set(gca, 'FontSize', 16)
set(gca, 'FontSize', 19)
legend('Location', 'best')
xlabel('Time')
ylabel('Averaged state over ' + string(traj_num) + ' trajectories')
xlim([0,100])
xt = get(gca, 'XTick'); 
set(gca, 'XTick', xt, 'XTickLabel', xt/10)  
set(gcf, 'position', [200 200 600 469])

safe = zeros(5, traj_num, Nt);

for k = 1:5
    for j = 1:traj_num
        for i = 1:Nt
            if x_all(k, i, j) >= 1
                safe(k, j, i) = 1;
            else
                break
            end
        end
    end
end

safe = sum(safe, 2)/traj_num;

figure
for k = 1:5
    plot(squeeze(safe(k,:,:)), LineWidth=1.5, Color=colors(k,:,:))
    hold on
end
legend('Proposed controller', 'StoCBF', 'PrSBC', 'CVaR', 'Nominal controller')
set(gca, 'FontSize', 16)
set(gca, 'FontSize', 19)
legend('Location', 'best')
xlabel('Time')
ylabel('Empirical safe probability')
xt = get(gca, 'XTick'); 
set(gca, 'XTick', xt, 'XTickLabel', xt/10)  
set(gcf, 'position', [200 200 600 469])

save('nonlinear_switching.mat')






