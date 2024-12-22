% switching control with multiple alpha values (proposed method)

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

% Multiple alpha values
alpha_values = [1, 1.2, 1.5, 2];

% Arrays to store results: (#alpha_values, Nt+1, traj_num)
x_all = zeros(length(alpha_values), Nt+1, traj_num);
u_all = zeros(length(alpha_values), Nt, traj_num);

for a_idx = 1:length(alpha_values)
    alpha = alpha_values(a_idx);

    x = zeros(Nt+1, traj_num); % initialization: note now we store Nt+1 steps
    x(1, :) = 3; % initial state
    u = zeros(Nt, traj_num);

    for j = 1:traj_num
        for i = 1:Nt
            P = mc_safe_prob(x(i, j), h, sigma);
            P_dx1 = mc_safe_prob(x(i, j)+dx, h, sigma);
            P_dx2 = mc_safe_prob(x(i, j)-dx, h, sigma);
            dP_x = (P_dx1 - P_dx2) / (2*dx); % gradient

            if dP_x*(A-K)*x(i, j) >= -alpha * (P - (1-eps))
                u(i, j) = -K * x(i, j);
                x(i+1, j) = exp((A-K)*dt)*x(i, j) + randn*sigma; % nominal controller
            else
                u(i, j) = (-alpha * (P - (1-eps)))/dP_x - A*x(i, j);
                x(i+1, j) = exp(A*dt)*x(i, j) + randn*sigma + u(i, j)*(exp(A*dt)-1)/A; % zero-hold control
            end
        end
    end

    x_all(a_idx,:,:) = x;
    u_all(a_idx,:,:) = u;
end

% Plot the averaged states with std
colors = lines(length(alpha_values)); 
figure; hold on

% Plot each alpha scenario with stdshade
hlines = gobjects(length(alpha_values),1);
for a_idx = 1:length(alpha_values)
    [hl, hp] = stdshade(squeeze(x_all(a_idx,:,:))', 0.3, colors(a_idx,:), 1:Nt+1);
    set(get(get(hp,'Annotation'),'LegendInformation'),'IconDisplayStyle','off'); % hide patch from legend
    hlines(a_idx) = hl;
end
yline(1, 'LineStyle', '--', 'color', 'red', 'LineWidth', 1.5)

ylim([0, 5])

% Legend entries using \alpha
legend_entries = arrayfun(@(x) ['\alpha = ' num2str(x)], alpha_values, 'UniformOutput', false);
legend(hlines, legend_entries{:}, 'Location', 'best')

set(gca, 'FontSize', 25)
xlabel('Time')
% ylabel(['Averaged state over ' num2str(traj_num) ' trajectories'])
ylabel(['Averaged state'])
xlim([0, Nt])
xt = get(gca, 'XTick'); 
set(gca, 'XTick', xt, 'XTickLabel', xt*dt)  
set(gcf, 'position', [200 200 600 469])

% Compute empirical safe probability
safe = zeros(length(alpha_values), traj_num, Nt);

for k = 1:length(alpha_values)
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

safe = sum(safe, 2)/traj_num; % Now safe is (length(alpha_values), 1, Nt)
safe = squeeze(safe); % shape: (length(alpha_values), Nt)

figure; hold on
for k = 1:length(alpha_values)
    plot(safe(k,:), 'LineWidth', 1.5, 'Color', colors(k,:))
end

ylim([0, 1])

legend(legend_entries{:}, 'Location', 'best')
set(gca, 'FontSize', 25)
xlabel('Time')
ylabel('Empirical safe probability')
xt = get(gca, 'XTick'); 
set(gca, 'XTick', xt, 'XTickLabel', xt*dt)  
set(gcf, 'position', [200 200 600 469])
