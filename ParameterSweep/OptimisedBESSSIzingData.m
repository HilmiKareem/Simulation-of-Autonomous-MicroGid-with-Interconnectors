% Data
bessRatio = [1.4, 1.425, 1.44, 1.45, 1.55, 1.65, 1.75, 1.8, 1.85, 1.95];
bessSize = [924, 940.5, 950.4, 957, 1023, 1089, 1155, 1188, 1221, 1287];
economicProfit = [-312.06, -202.58, -12.32, 40.73, 626.85, 992.89, 1320.85, 1759.1, 2742.75, 2958.12];

% Create figure
figure('Name', 'Economic Profit vs BESS Sizing', 'Color', 'w');
hold on;
grid on;
box on;

% Plot economic profit vs BESS size
plot(bessSize, economicProfit, '-o', 'LineWidth', 2, 'Color', [0 0.45 0.74], ...
    'MarkerFaceColor', [0 0.45 0.74]);

% Add red dashed line for break-even point
yline(0, 'r--', 'LineWidth', 2);

% Labels and title
xlabel('BESS Size (kWh)');
ylabel('Economic Profit (£)');
title('Economic Profit vs BESS Sizing');

% Improve readability
set(gca, 'FontSize', 12);
ytickformat('£%,.0f');


