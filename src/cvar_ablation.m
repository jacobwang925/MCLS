% switching control with multiple gamma values

clear; clc;
Nt = 100; % time points
dt = 0.1; % time step size
traj_num = 20; % number of trajectories
dx = 0.1; % state step size (for derivative calculations)
eps = 0.1; % epsilon
K = 2.5; % proportional controller gain
A = 2; % system dynamics f(x) = A
sigma = 1; % magnitude of noise
h = 10; % safe prob time horizon

sigma = sigma * sqrt(dt); % discretization

%% CVaR parameters (renaming alpha_risk to gamma for clarity)
gamma_values = [0.25, 0.45, 0.65, 0.75, 0.85];

beta_risk = 0.1;
sigma_risk = sqrt(dt);
epsilon = norminv(beta_risk, 0, sigma_risk);
l_risk = quadgk(@(x) x.*normpdf(x, 0, sigma_risk), -inf, epsilon)/normcdf(epsilon, 0, sigma_risk); 

% Arrays to store results:
x_all = zeros(length(gamma_values), Nt+1, traj_num);
u_all = zeros(length(gamma_values), Nt, traj_num);

for a_idx = 1:length(gamma_values)
    gamma_val = gamma_values(a_idx);

    x_cvar = zeros(Nt+1, traj_num);
    u_cvar = zeros(Nt, traj_num);
    x_cvar(1, :) = 3;

    for j = 1:traj_num
        for i = 1:Nt
            P = mc_safe_prob(x_cvar(i, j), h, sigma);
            if P > 1 - eps
                u_cvar(i, j) = -K*x_cvar(i, j);
                x_cvar(i+1, j) = exp((A-K)*dt)*x_cvar(i, j) + randn*sigma; 
            else
                if exp((A-K)*dt)*x_cvar(i, j) + l_risk >= gamma_val * x_cvar(i, j)
                    u_cvar(i, j) = -K*x_cvar(i, j);
                    x_cvar(i+1, j) = exp((A-K)*dt)*x_cvar(i, j) + randn*sigma; 
                else
                    u_cvar(i, j) = ((gamma_val - exp(A*dt))*x_cvar(i, j) - l_risk)*A/(exp(A*dt)-1);
                    x_cvar(i+1, j) = exp(A*dt)*x_cvar(i, j) + randn*sigma + u_cvar(i, j)*(exp(A*dt)-1)/A; 
                end
            end
        end
    end

    x_all(a_idx,:,:) = x_cvar;
    u_all(a_idx,:,:) = u_cvar;
end

%% Plotting the averaged states with std
colors = lines(length(gamma_values)); 
figure; hold on

% We'll collect line handles for the legend and insert empty entries as needed
hlines = gobjects(length(gamma_values),1);
for a_idx = 1:length(gamma_values)
    [hl, hp] = stdshade(squeeze(x_all(a_idx,:,:))', 0.3, colors(a_idx,:), 1:Nt+1);
    set(get(get(hp,'Annotation'),'LegendInformation'),'IconDisplayStyle','off'); % hide patches from legend
    hlines(a_idx) = hl;
    % For the requested legend style: alternate empty legend entries was previously done by user as placeholders.
    % However, since we are using a single legend now with \gamma, we don't need the empty placeholders.
    % If you still need them, you can reintroduce them by adding extra plot calls or strings.
end

yline(1, 'LineStyle', '--', 'color', 'red', 'LineWidth', 1.5)

legend_entries = arrayfun(@(x) ['\gamma = ', num2str(x)], gamma_values, 'UniformOutput', false);
legend(hlines, legend_entries{:}, 'Location', 'best')

set(gca, 'FontSize', 25)
xlabel('Time')
% ylabel(['Averaged state over ' num2str(traj_num) ' trajectories'])
ylabel(['Averaged state'])
xlim([0, Nt])
xt = get(gca, 'XTick'); 
set(gca, 'XTick', xt, 'XTickLabel', xt*dt)  
set(gcf, 'position', [200 200 600 469])

%% Compute empirical safe probability as per the user's original method
safe = zeros(length(gamma_values), traj_num, Nt);

for k = 1:length(gamma_values)
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
safe = squeeze(safe); % Now safe is (length(gamma_values), Nt)

figure; hold on
for k = 1:length(gamma_values)
    plot(safe(k,:), 'LineWidth', 1.5, 'Color', colors(k,:))
end

legend(legend_entries{:}, 'Location', 'best')
set(gca, 'FontSize', 25)
xlabel('Time')
ylabel('Empirical safe probability')
xt = get(gca, 'XTick'); 
set(gca, 'XTick', xt, 'XTickLabel', xt*dt)  
set(gcf, 'position', [200 200 600 469])
