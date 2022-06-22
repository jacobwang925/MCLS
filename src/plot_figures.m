%% worst case
clear; clc;
load("worst_case.mat")
% load('nonlinear_worst_case.mat')
fontsz = 20;

figure
plot(F, 'linewidth', 1.5)
hold on
plot(0)
plot(F_clark, 'linewidth', 1.5)
plot(F_prsbc, 'linewidth', 1.5)
plot(F_cvar, 'linewidth', 1.5)
legend('Proposed controller', '', 'StoCBF', 'PrSBC', 'CVaR')
% set(gca, 'FontSize', 16)
set(gca, 'FontSize', 25)
legend('Location', 'best')
xlabel('Time')
% ylabel('Expected value of safe probability')
ylabel('Expected safe probability')
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
% set(gca, 'FontSize', 16)
set(gca, 'FontSize', 25)
legend('Location', 'best')
xlabel('Time')
% ylabel('Averaged state over ' + string(traj_num) + ' trajectories')
ylabel('Averaged state')
xlim([0,100])
xt = get(gca, 'XTick'); 
set(gca, 'XTick', xt, 'XTickLabel', xt/10)
set(gcf, 'position', [200 200 600 469])


%% switching
clear; clc;
load('switching.mat')
% load('nonlinear_switching.mat')
fontsz = 20;

figure
for i = 1:5
    stdshade(squeeze(x_all(i,:,:))', 0.3, colors(i,:,:));
    hold on
end
yline(1, 'LineStyle', '--', 'color', 'red', 'LineWidth', 1.5)
yline(1.5, 'LineStyle', '--', 'color', 'black', 'linewidth', 1.5)
% legend('', 'Proposed controller', '', 'StoCBF', '', 'PrSBC', '', 'CVaR', '', 'Nominal controller')
legend('', 'Proposed controller', '', 'Clark, (2019).', '', 'Luo et al., (2019)', '', 'Ahmadi et al., (2020)', '', 'Nominal controller')
set(gca, 'FontSize', 25)
legend('Location', 'best')
xlabel('Time')
% ylabel('Averaged state over ' + string(traj_num) + ' trajectories')
ylabel('Averaged state')
xlim([0,100])
ylim([-1,3.5])
xt = get(gca, 'XTick'); 
set(gca, 'XTick', xt, 'XTickLabel', xt/10)  
set(gcf, 'position', [200 200 600 469])

figure
for k = 1:5
    plot(squeeze(safe(k,:,:)), LineWidth=1.5, Color=colors(k,:,:))
    hold on
end
legend('Proposed controller', 'Clark, (2019).', 'Luo et al., (2019)', 'Ahmadi et al., (2020)', 'Nominal controller')
% legend('Proposed controller', 'StoCBF', 'PrSBC', 'CVaR', 'Nominal controller')
set(gca, 'FontSize', 25)
legend('Location', 'best')
xlabel('Time')
ylabel('Empirical safe probability')
xt = get(gca, 'XTick'); 
set(gca, 'XTick', xt, 'XTickLabel', xt/10)  
set(gcf, 'position', [200 200 600 469])