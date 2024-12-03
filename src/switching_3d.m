% worst case control

clear; clc;
Nt = 200; % time points
dt = 0.05; % time step size
traj_num = 10; % number of trajectories
dx = 0.1; % state step size (for derivative calculations)
eps = 0.1; % epsilon
K = 2.5; % proportional controller gain
A = [2, 1.9, 1.8]; % A matrix in system dynamics
sigma = 1; % magnitude of noise
h = 10; % safe prob time horizon
x_0 = [5, 5, 5]; % initial state

sigma = sigma * sqrt(dt); % discretization

%% Proposed controller
x = zeros(Nt, 3, traj_num); % initialization
u_all = zeros(Nt, traj_num); 
x(1, :, :) = repmat(x_0, [1, 1, traj_num]);
F = zeros(Nt, 1); % to store expected value of safe probability

alpha = 1;

for j = 1:traj_num
    for i = 1:Nt
        P = mc_safe_prob_3d(x(i, :, j), h, sigma);
        P_dx1_p = mc_safe_prob_3d(x(i, :, j)+[dx,0,0], h, sigma);
        P_dx1_m = mc_safe_prob_3d(x(i, :, j)-[dx,0,0], h, sigma);
        P_dx2_p = mc_safe_prob_3d(x(i, :, j)+[0,dx,0], h, sigma);
        P_dx2_m = mc_safe_prob_3d(x(i, :, j)-[0,dx,0], h, sigma);
        P_dx3_p = mc_safe_prob_3d(x(i, :, j)+[0,0,dx], h, sigma);
        P_dx3_m = mc_safe_prob_3d(x(i, :, j)-[0,0,dx], h, sigma);
        dP_x1 = (P_dx1_p - P_dx1_m) / (2*dx); % gradient
        dP_x2 = (P_dx2_p - P_dx3_m) / (2*dx);
        dP_x3 = (P_dx2_p - P_dx3_m) / (2*dx);
        if abs(dP_x1) < 0.005
            dP_x1 = 0.005;
        end
        if abs(dP_x2) < 0.005
            dP_x2 = 0.005;
        end
        if abs(dP_x3) < 0.005
            dP_x3 = 0.005;
        end
        if P > 1-eps || dP_x1*(A(1)-K)*x(i, 1, j) + dP_x2*(A(2)-K)*x(i, 2, j) + dP_x3*(A(3)-K)*x(i, 3, j)>= -alpha * (P - (1-eps))
            % nominal controller
            u_all(i, j) = -K * x(i, 1, j);
            x(i+1, 1, j) = exp((A(1)-K)*dt) * x(i, 1, j) + randn*sigma;
            x(i+1, 2, j) = exp((A(2)-K)*dt) * x(i, 2, j) + randn*sigma;
            x(i+1, 3, j) = exp((A(3)-K)*dt) * x(i, 3, j) + randn*sigma;
        else
            u = (-alpha * (P - (1-eps)) - ((A(1)-K)*x(i,1,j)*dP_x1 + (A(2)-K)*x(i,2,j)*dP_x2 + (A(3)-K)*x(i,3,j)*dP_x3)) / (dP_x1 + dP_x2 + dP_x3);
            u_all(i, j) = u;
            x(i+1, 1, j) = exp(A(1)*dt) * x(i, 1, j) + u*(exp(A(1)*dt)-1)/A(1) + randn*sigma;
            x(i+1, 2, j) = exp(A(2)*dt) * x(i, 2, j) + u*(exp(A(2)*dt)-1)/A(2) + randn*sigma;
            x(i+1, 3, j) = exp(A(3)*dt) * x(i, 3, j) + u*(exp(A(3)*dt)-1)/A(3) + randn*sigma;
        end
    end
end
x_ave = mean(x, 3);


%%
% Define colors for each plot (assuming you have 5 different controllers or categories)
colors = lines(5);

% Plot for averaged state trajectories with shaded standard deviation
figure
hold on

% Loop through the first 3 state variables (assuming x is [Nt, 3, traj_num])
for i = 1:3
    stdshade(squeeze(x(:,i,:))', 0.3, colors(i,:)); % Plot mean and shaded std for each state
end

% Add a red dashed line at y = 1 for reference
yline(1, 'LineStyle', '--', 'color', 'red', 'LineWidth', 1.5)

% Set up the legend for each controller or state
legend('', '$x_1$', '', '$x_2$', '', '$x_3$', 'Location', 'best', 'interpreter', 'latex')

% Adjust the font size and labels
set(gca, 'FontSize', 16)
xlabel('Time (s)')
% ylabel('Averaged state over ' + string(traj_num) + ' trajectories')
ylabel('Averaged state')

% Set x and y limits
xlim([0, Nt])
ylim([-2, 20]) % Adjust based on your data

% Adjust X-tick labels to represent actual time
xt = get(gca, 'XTick'); 
set(gca, 'XTick', xt, 'XTickLabel', xt/20)  

% Set figure size
set(gcf, 'position', [200 200 600 469])

hold off


% Initialize a vector to store the safety probabilities of the average trajectory
safe_prob_ave = zeros(Nt, 1);

% Loop through each time step and compute safety probability for the average trajectory
for i = 1:Nt
    safe_prob_ave(i) = mc_safe_prob_3d(x_ave(i, :), h, sigma);
end

% Plot the safety probability of the mean trajectory over time
figure
plot((1:Nt) * dt, safe_prob_ave, 'LineWidth', 2, 'Color', [0, 0.4470, 0.7410]) % Blue line with increased width

% Add labels and title
xlabel('Time (s)', 'FontSize', 14)
ylabel('Safety Probability', 'FontSize', 14)
% title('Safety Probability of the Mean Trajectory', 'FontSize', 16)

% Adjust y-axis limits to match the probability range
ylim([0, 1])

% Add a grid for better readability
grid on
set(gca, 'YTick', 0:0.1:1) % Set y-axis ticks with a 0.1 interval
set(gca, 'XTick', 0:2.5:Nt*dt) % Set x-axis ticks with a 0.1 interval

% Set the font size for axes
set(gca, 'FontSize', 16)

% Set the figure size for better aspect ratio
set(gcf, 'Position', [200 200 600 469])



% Initialize a matrix to store the safety status for each trajectory at each time step
safe = zeros(traj_num, Nt);

% Define safety threshold
safety_threshold = 1;  % Change this value based on your safety criterion

% Loop through each trajectory and each time step to check safety
for j = 1:traj_num
    for i = 1:Nt
        % Check if all state components (x_1, x_2, x_3) are above the safety threshold
        if all(x(i, :, j) >= safety_threshold)
            safe(j, i) = 1;  % Mark the trajectory as safe at this time step
        else
            break  % If a trajectory becomes unsafe, stop checking further time steps for this trajectory
        end
    end
end

% Initialize a vector to store the empirical safety probabilities
empirical_safe_prob = zeros(Nt, 1);

% Compute the empirical safety probability at each time step
for i = 1:Nt
    % The empirical safety probability is the fraction of safe trajectories at each time step
    empirical_safe_prob(i) = sum(safe(:, i)) / traj_num;
end

% Plot the empirical safety probability over time
figure
plot((1:Nt) * dt, empirical_safe_prob, 'LineWidth', 2, 'Color', [0.8500, 0.3250, 0.0980]) % Orange line

% Add labels and title
xlabel('Time (s)', 'FontSize', 14)
ylabel('Empirical Safety Probability', 'FontSize', 14)
% title('Empirical Safety Probability Over Time', 'FontSize', 16)

% Adjust y-axis limits to match the probability range
ylim([0, 1])

% Add a grid for better readability with a 0.1 interval
grid on
set(gca, 'XTick', 0:2.5:Nt*dt) % Set x-axis ticks to match time scale
set(gca, 'YTick', 0:0.1:1) % Set y-axis ticks with a 0.1 interval

% Set the font size for axes
set(gca, 'FontSize', 16)

% Set the figure size for better aspect ratio
set(gcf, 'Position', [200 200 600 469])
