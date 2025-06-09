
%% 1. Define the data
threshold = [  0,   5,   10,   15,   20,   25,   30,   35,   40,   45,   50,   60,   70,   80,   90,  100 ]';
economicProfit      = [ -830.62, -747.39, -675.46, -849.38, -388.52, -471.88,  304.39,  304.39,  238.41,  590.94, 1040.15, 1350.14, 2296.38, 1988.38, 2958.12, 2834.07 ]';
unservedEnergyCost  = [ 7419.83, 7288.58, 7114.30, 7110.22, 6526.85, 6526.85, 5383.09, 5383.09, 5383.09, 4785.55, 4165.75, 3592.30, 2248.53, 2248.53,  613.66,   49.52 ]';
exportRevenue       = [ 7678.38, 7630.35, 7528.00, 7350.00, 7227.50, 7105.00, 6737.50, 6737.50, 6492.50, 6247.50, 6002.50, 5635.00, 5145.00, 4655.00, 3920.00, 3185.00 ]';
finalSoC            = [   6.85,   7.30,  10.89,  19.65,  19.65,  26.84,  26.84,  26.84,  24.75,  24.75,  26.18,  26.18,  27.38,  23.73,  37.33,  48.94 ]';

%% 2. Create a 2×2 tiled layout of line charts
figure('Units', 'normalized', 'Position', [0.1 0.1 0.8 0.8])
tiledlayout(2,2, 'Padding', 'compact', 'TileSpacing', 'compact');

% 2.1 Economic Profit vs Export Threshold
nexttile
plot(threshold, economicProfit, '-o', 'LineWidth', 2, 'MarkerSize', 6)
grid on
xlabel('Min Export Threshold (%)')
ylabel('Economic Profit (£)')
title('Economic Profit vs Export Threshold')
% Optional: annotate the point of maximum economic profit
[ep_max, idx_ep_max] = max(economicProfit);
hold on


% 2.2 Unserved Energy Cost vs Export Threshold
nexttile
plot(threshold, unservedEnergyCost, '-o', 'LineWidth', 2, 'MarkerSize', 6)
grid on
xlabel('Min Export Threshold (%)')
ylabel('Unserved Energy Cost (£)')
title('Unserved Energy Cost vs Export Threshold')
% Annotate the point of minimum unserved energy cost
[uec_min, idx_uec_min] = min(unservedEnergyCost);
hold on


% 2.3 Export Revenue vs Export Threshold
nexttile
plot(threshold, exportRevenue, '-o', 'LineWidth', 2, 'MarkerSize', 6)
grid on
xlabel('Min Export Threshold (%)')
ylabel('Export Revenue (£)')
title('Export Revenue vs Export Threshold')
% Annotate the point of maximum export revenue
[er_max, idx_er_max] = max(exportRevenue);
hold on



% 2.4 Final State of Charge vs Export Threshold
nexttile
plot(threshold, finalSoC, '-o', 'LineWidth', 2, 'MarkerSize', 6)
grid on
xlabel('Min Export Threshold (%)')
ylabel('Final SoC (%)')
title('Final State of Charge vs Export Threshold')
% Annotate the point of maximum final SoC
[fsc_max, idx_fsc_max] = max(finalSoC);
hold on

