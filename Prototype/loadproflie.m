% Filepath for the Excel file
filename = 'C:\Users\ASUS\OneDrive - Imperial College London\Documents\MATLAB Drive\FYP\Clover Load 7days.xlsx'; 

% Read data from each sheet
PDomestic = readmatrix(filename, 'Sheet', 'domestic');
PCommerical = readmatrix(filename, 'Sheet', 'commercial');
PPublic = readmatrix(filename, 'Sheet', 'public');

% Reshape data if necessary
PDomestic = reshape(PDomestic, [], 1);
PCommerical = reshape(PCommerical, [], 1);
PPublic = reshape(PPublic, [], 1);

% Ensure all data vectors have 168 points (7 days x 24 hours)
if numel(PDomestic) ~= 168 || numel(PCommerical) ~= 168 || numel(PPublic) ~= 168
    error('The load data does not contain the expected 168 data points (7 days x 24 hours).');
end

% Create a time vector for the 7 days (assuming hourly data)
time = 0:167; % 168 points for 7 days (7 * 24 hours)

%% Telcommunication Load Calculation
% Input parameters
Total_Area = 56000; % m² (total area to be covered)
Range_per_tower = 40000; % m (coverage radius of one tower)
Consumption_per_tower = 1e3; % W (power consumption per tower)
Overlap_factor = 1.2; % Accounts for overlapping coverage (e.g., 20% overlap)

% Calculate coverage area per tower (assuming circular coverage)
Coverage_per_tower = pi * Range_per_tower^2; % m²

% Adjust for overlapping coverage
Effective_Coverage_per_tower = Coverage_per_tower / Overlap_factor;

% Calculate number of towers needed
Number_Tower_Needed = ceil(Total_Area / Effective_Coverage_per_tower);

% Calculate total power consumption
Total_Consumption = (Number_Tower_Needed+2) * Consumption_per_tower; % W, the plus 2 is added to cover some areas not in the centre town of the island

% Create a constant telecommunication load profile for 168 hours (7 days)
PTelco = Total_Consumption * ones(168, 1); % Constant load for 24/7 operation

%% Healthcare 
% Provided data and parameters
daily_profile = [250, 250, 250, 250, 250, 700, 750, 810, 910, 1200, 1200, 1600, ...
                 1700, 1600, 950, 900, 700, 700, 500, 500, 400, 400, 300, 300];  % in watts
daily_energy_target = 11.5 * 1000;  % Convert kWh to Wh
day_to_day_variability = 0.10;  % 10%
hour_to_hour_variability = 0.15;  % 15%
days = 7;  % Week-long profile

% Normalize and scale the daily profile to match the energy target
scaling_factor = daily_energy_target / sum(daily_profile);
scaled_daily_profile = daily_profile * scaling_factor;

% Initialize weekly profile
PHealthcare = zeros(1, days * 24);

% Generate the weekly load profile
for day = 1:days
    % Apply day-to-day variability
    day_variation_factor = 1 + (rand * 2 - 1) * day_to_day_variability;
    daily_load = scaled_daily_profile * day_variation_factor;

    % Apply hour-to-hour variability
    hourly_loads = daily_load .* (1 + (rand(1, 24) * 2 - 1) * hour_to_hour_variability);

    % Store the result in the weekly profile array
    PHealthcare((day-1)*24 + 1 : day*24) = hourly_loads;
end

% Optional: Ensure the total energy over the week matches the expected value
total_weekly_energy = sum(PHealthcare) / 1000;  % Convert back to kWh
if abs(total_weekly_energy - 11.5 * days) > 0.1
    warning('Total weekly energy deviates significantly from the target.');
end

% Plot the weekly load profile
figure;
plot(1:168, weekly_profile);
title('Clinic Weekly Load Profile (1x168 Array)');
xlabel('Hour');
ylabel('Load (Watts)');
grid on;

% Display total weekly energy consumption in kWh
total_energy_weekly = sum(weekly_profile) / 1000;  % Convert Wh to kWh
disp('Total weekly energy consumption (kWh):');
disp(total_energy_weekly);

%% School 

% Load profile data for the secondary school
total_daily_energy_target = 13.56 * 1000; % Convert kWh to Wh
peak_load = 3050; % Peak load in watts during school hours

% Define active consumption hours (6 AM to 5 PM)
active_hours = 6:16; % Hours from 6 AM (hour 6) to 5 PM (hour 16)

% Base load profile (non-active hours set to zero)
hourly_load = zeros(1, 24);
for h = active_hours
    hourly_load(h) = 1; % Uniform load within active period (will be scaled)
end

% Scale the load profile to match the daily energy target
scaling_factor = total_daily_energy_target / sum(hourly_load);
scaled_hourly_load = hourly_load * scaling_factor;

% Initialize weekly profile (1x168 array)
num_days = 7;
weekly_profile = zeros(1, num_days * 24);

% Generate the weekly profile with variability and weekend adjustment
index = 1;
day_to_day_variability = 0.10; % 10%
hour_to_hour_variability = 0.15; % 15%
weekend_load_reduction = 0.90; % Load reduced to 90% on weekends

for day = 1:num_days
    % Weekend load adjustment (90% lower consumption on weekends)
    if day >= 6 % Saturday (day 6) and Sunday (day 7)
        load_reduction_factor = 1 - weekend_load_reduction;
    else
        load_reduction_factor = 1;
    end

    % Apply day-to-day variability
    day_variation_factor = 1 + (rand * 2 - 1) * day_to_day_variability;
    daily_load = scaled_hourly_load * day_variation_factor * load_reduction_factor;

    % Apply hour-to-hour variability
    hourly_loads = daily_load .* (1 + (rand(1, 24) * 2 - 1) * hour_to_hour_variability);

    % Store in weekly profile array
    weekly_profile(index:index + 23) = hourly_loads;
    index = index + 24;
end

% Plot the weekly profile
figure;
plot(1:168, weekly_profile);
title('School Weekly Load Profile (1x168 Array with Weekend Load Reduction)');
xlabel('Hour');
ylabel('Load (Watts)');
grid on;

% Plot daily load profile for visual validation
figure;
bar(0:23, scaled_hourly_load / 1000, 'b'); % Convert to kW
title('Daily Load Profile (Scaled)');
xlabel('Hour');
ylabel('Load (kW)');
grid on;

% Display total weekly energy consumption in kWh
total_energy_weekly = sum(weekly_profile) / 1000; % Convert Wh to kWh
disp('Total weekly energy consumption (kWh):');
disp(total_energy_weekly);

%% Plotting the Load Profiles
figure;

% Plot domestic, commercial, and public loads
plot(time, PDomestic, '-o', 'DisplayName', 'Domestic', 'LineWidth', 1.5);
hold on;
plot(time, PCommerical, '-s', 'DisplayName', 'Commercial', 'LineWidth', 1.5);
plot(time, PPublic, '-^', 'DisplayName', 'Public', 'LineWidth', 1.5);

% Plot telecommunication load
plot(time, PTelco, '-d', 'DisplayName', 'Telecommunication', 'LineWidth', 1.5);

% Enhance the plot
xlabel('Time (Hour)');
ylabel('Load (W)'); % Updated to Watt (W)
title('Load Data for 7 Days');
legend('show');
grid on;

% Set the x-axis labels to show hours: 0000, 0100, ..., 2300
xticks(0:24:167); % Every 24th data point corresponds to a new day (24, 48, ..., 168)
xticklabels(repmat(0:23, 1, 7)); % Repeat 0000-2300 for 7 days

