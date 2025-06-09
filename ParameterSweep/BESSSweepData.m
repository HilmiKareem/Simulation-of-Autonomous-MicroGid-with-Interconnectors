% Data
bessRatio = [1.94, 1.95, 1.96, 1.97, 1.98, 1.99, 2.0, 2.1, 2.2];
economicProfit = [-109.96, -3.62, 11.96, 30.15, 144.09, 28.55, 21.08, 118.05, 270.13];
exportRevenue = [6510, 6510, 6510, 6510, 6510, 6615, 6615, 6720, 6825];

%% Economic Profit vs BESS Ratio
figure('Name','BESS Analysis Dashboard','NumberTitle','off');

subplot(2,2,1)
scatter(bessRatio, economicProfit, 60, 'filled')
grid on
xlabel('BESS Ratio')
ylabel('Economic Profit (£)')
title('Economic Profit vs BESS Ratio')

%% Export Revenue vs BESS Ratio
subplot(2,2,2)
scatter(bessRatio, exportRevenue, 60, 'filled')
grid on
xlabel('BESS Ratio')
ylabel('Export Revenue (£)')
title('Export Revenue vs BESS Ratio')

%% Line Plot: Economic Profit Trend
subplot(2,2,[3 4])
plot(bessRatio, economicProfit, '-o', 'LineWidth', 2)
hold on
yline(0, 'r', 'LineWidth', 1.5)  % <-- Adds the red horizontal line at y = 0
hold off
grid on
xlabel('BESS Ratio')
ylabel('Economic Profit (£)')
title('Combined Analysis: Economic Profit Trend')
