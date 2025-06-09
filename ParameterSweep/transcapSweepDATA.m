% Transmission Capacity Parameter Analysis in MATLAB
% This script reproduces the Seller and Buyer Capacity analyses
% with dual y‐axes for Cost/Profit (or Revenue/Profit) and Final SoC.
% Each line uses a distinct color matching the original Recharts version.

%% Seller Capacity Analysis (Buyer fixed at 300 kW)
% Data
sellerCapacity   = [0, 50, 100, 150, 200, 250, 300];
importCost       = [0, 566.69, 840.53, 908.50, 909.61, 909.61, 909.66];
economicProfit   = [-1282.87, 23.74, 65.46, -2.51, -3.62, -3.62, -5.36];
finalSoC_seller  = [6.85, 6.85, 15.86, 27.00, 27.18, 27.18, 27.18];

% Define RGB colors (normalized 0–1) corresponding to hex codes:
% Import Cost:   #2563eb → [0.1451, 0.3882, 0.9216]
% Economic Profit: #dc2626 → [0.8627, 0.1490, 0.1490]
% Final SoC:     #16a34a → [0.0863, 0.6392, 0.2902]

colorImportCost    = [37, 99, 235] ./ 255;
colorEconomicProfit = [220, 38, 38] ./ 255;
colorFinalSoC      = [22, 163, 74] ./ 255;

% Create Figure
figure('Name', 'Seller Capacity Analysis', 'NumberTitle', 'off');
yyaxis left
    p1 = plot(sellerCapacity, importCost, '-o', ...
        'Color', colorImportCost, 'LineWidth', 2, 'MarkerSize', 6);
    hold on
    p2 = plot(sellerCapacity, economicProfit, '-o', ...
        'Color', colorEconomicProfit, 'LineWidth', 2, 'MarkerSize', 6);
    ylabel('Cost / Profit (£)', 'FontSize', 12)
yyaxis right
    p3 = plot(sellerCapacity, finalSoC_seller, '-o', ...
        'Color', colorFinalSoC, 'LineWidth', 2, 'MarkerSize', 6);
    ylabel('Final SoC (%)', 'FontSize', 12)
hold off

% Labels, Grid, Legend, Title
xlabel('Seller Capacity (kW)', 'FontSize', 12)
title('Seller Capacity Analysis (Buyer Fixed at 300 kW)', 'FontSize', 14, 'FontWeight', 'bold')
grid on
legend([p1, p2, p3], {'Import Cost', 'Economic Profit', 'Final SoC'}, ...
       'Location', 'best', 'FontSize', 10)

%% Buyer Capacity Analysis (Seller fixed at 300 kW)
% Data
buyerCapacity    = [0, 50, 100, 150, 200, 250, 300, 350, 400, 450, 500];
exportRevenue    = [0, 1557.5, 2835.0, 4042.5, 5180.0, 5950.0, 6510.0, 6737.5, 7000.0, 7087.5, 7350.0];
economicProfitB  = [-350.93, -922.20, -627.50, -275.99, -204.17, 97.70, -5.36, 303.29, 19.23, -442.65, -453.84];
finalSoC_buyer   = [48.94, 25.83, 25.83, 25.32, 23.70, 25.42, 27.18, 27.02, 25.65, 27.06, 17.88];

% Define RGB colors:
% Export Revenue:   #7c3aed → [0.4863, 0.2275, 0.9294]
% Economic Profit:  #dc2626 → [0.8627, 0.1490, 0.1490]
% Final SoC:        #16a34a → [0.0863, 0.6392, 0.2902]

colorExportRevenue  = [124, 58, 237] ./ 255;
% Reuse colorEconomicProfit and colorFinalSoC from above

% Create Figure
figure('Name', 'Buyer Capacity Analysis', 'NumberTitle', 'off');
yyaxis left
    q1 = plot(buyerCapacity, exportRevenue, '-s', ...
        'Color', colorExportRevenue, 'LineWidth', 2, 'MarkerSize', 6);
    hold on
    q2 = plot(buyerCapacity, economicProfitB, '-s', ...
        'Color', colorEconomicProfit, 'LineWidth', 2, 'MarkerSize', 6);
    ylabel('Revenue / Profit (£)', 'FontSize', 12)
yyaxis right
    q3 = plot(buyerCapacity, finalSoC_buyer, '-s', ...
        'Color', colorFinalSoC, 'LineWidth', 2, 'MarkerSize', 6);
    ylabel('Final SoC (%)', 'FontSize', 12)
hold off

% Labels, Grid, Legend, Title
xlabel('Buyer Capacity (kW)', 'FontSize', 12)
title('Buyer Capacity Analysis (Seller Fixed at 300 kW)', 'FontSize', 14, 'FontWeight', 'bold')
grid on
legend([q1, q2, q3], {'Export Revenue', 'Economic Profit', 'Final SoC'}, ...
       'Location', 'best', 'FontSize', 10)

