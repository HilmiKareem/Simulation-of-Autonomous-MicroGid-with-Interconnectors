clear all

%% Define Value of Lost Load (VoLL) Parameters (in $/kWh)
VoLL_Domestic             = 2;    % 2 to 4
VoLL_HighImpactCommercial = 8;    % 8 to 12
VoLL_LowImpactCommercial  = 4;    % 4 to 6
VoLL_Education            = 1.5;  % 1.5 to 3
VoLL_Healthcare           = 10;   % 10 to 15
VoLL_Public               = 1;    % 1 to 2
VoLL_Telecommunications   = 12;   % 12 to 20

%% System Parameters
PV_Module_Price   = 1481.43;       % $/kW
Battery_Cost      = 2160.43;       % $/kWh
Rs                = 1.2;           % Sizing ratio
P_max_inverter    = 230e3;         % Maximum inverter power (W)
PV_peak           = Rs * P_max_inverter; % Peak PV generation (W)
n_inverter        = 0.9;           % Inverter efficiency (AC output)
n_charger         = 0.9;           % Charger efficiency (DC charging)
BESS_size         = 500e3;         % Battery capacity (Wh)
SoC_max           = 100;           % Maximum State of Charge (%)
SoC_min           = 0;             % Minimum State of Charge (%)
SoC               = 50;            % Initial State of Charge (%)
P_max_charger     = 500e3;         % Maximum charger power (W)
time_steps        = 168;           % Simulation duration (hours; 7 days)
delta_t           = 1;             % Time step (hours)
min_Trading_SoC   = 30;            % Minimum SoC for trading (%)
DOD               = 5;             % Depth of Discharge (%)
%% Time Settings
hours_per_day = 24;
days          = 7;
total_hours   = hours_per_day * days;

%% External Grid Parameters
% External Grid 1 (Seller)
ext1_total_energy = 100;           % kWh over 7 days (positive for selling)
ext1_energy = ext1_total_energy * ones(1, time_steps);
ext1_price = 0.1 * ones(1, time_steps);  % $0.5/kWh
ext1_trans_cap = 100;              % kWh per hour

% External Grid 2 (Buyer)
ext2_total_energy = -100;          % kWh over 7 days (negative for buying)
ext2_energy = ext2_total_energy * ones(1, time_steps);
ext2_price = 0.8 * ones(1, time_steps);  % $0.8/kWh
ext2_trans_cap = 100;              % kWh per hour

%% Initialize Variables
energy_sold_ext1 = zeros(1, time_steps);
energy_sold_ext2 = zeros(1, time_steps);
energy_bought_ext1 = zeros(1, time_steps);
energy_bought_ext2 = zeros(1, time_steps);
trading_cost = zeros(1, time_steps);
trading_profit = zeros(1, time_steps);
Cost_of_Unserved_array = zeros(1, time_steps);
SoC_array = zeros(1, time_steps);
SoC_array(1) = SoC;
E_storage_array = zeros(1, time_steps);
E_demand_effective_array = zeros(1, time_steps);

%% Load Profiles from Excel (Replace with actual data)
filename = 'C:\Users\ASUS\OneDrive - Imperial College London\Documents\MATLAB Drive\FYP\Clover Load 7days.xlsx'; 

% Read and reshape data from Excel sheets
PDomestic   = reshape(readmatrix(filename, 'Sheet', 'domestic'), [], 1);
PCommerical = reshape(readmatrix(filename, 'Sheet', 'commercial'), [], 1);
PPublic     = reshape(readmatrix(filename, 'Sheet', 'public'), [], 1);

% Split Commercial Load into High and Low Impact
PCommerical_HighImpact = 0.2 * PCommerical;  % 20%
PCommerical_LowImpact  = 0.8 * PCommerical;  % 80%

%% Telecommunication Load Calculation
Total_Area            = 56000;        % m²
Range_per_tower       = 40000;        % m (coverage radius)
Consumption_per_tower = 5e3;          % W per tower
Overlap_factor        = 1.2;          % accounts for overlap

Coverage_per_tower    = pi * Range_per_tower^2;                % m² per tower
Effective_Coverage    = Coverage_per_tower / Overlap_factor;   % m² per tower
Number_Tower_Needed   = ceil(Total_Area / Effective_Coverage);
Total_Consumption     = (Number_Tower_Needed + 2) * Consumption_per_tower;

% Create constant telecommunication load for 168 hours
PTelco = Total_Consumption * ones(time_steps, 1);

%% Healthcare Load Profile
daily_profile        = [250, 250, 250, 250, 250, 700, 750, 810, 910, 1200, 1200, 1600, ...
                        1700, 1600, 950, 900, 700, 700, 500, 500, 400, 400, 300, 300];  % in W
daily_energy_target  = 11.5 * 1000;  % Wh (11.5 kWh)
day_variability      = 0.10;          % 10%
hour_variability     = 0.15;          % 15%
scaling_factor       = daily_energy_target / sum(daily_profile);
scaled_daily_profile = daily_profile * scaling_factor;

PHealthcare = zeros(1, days * 24);
for day = 1:days
    day_factor   = 1 + (rand * 2 - 1) * day_variability;
    daily_load   = scaled_daily_profile * day_factor;
    hourly_loads = daily_load .* (1 + (rand(1, 24) * 2 - 1) * hour_variability);
    PHealthcare((day-1)*24 + 1 : day*24) = hourly_loads;
end

%% Education Facility Load Profile
total_daily_energy_target = 13.56 * 1000;  % Wh
peak_load                 = 3050;           % W (during school hours)
active_hours              = 6:16;           % 6 AM to 5 PM

hourly_load = zeros(1, 24);
hourly_load(active_hours) = 1;
scaling_factor = total_daily_energy_target / sum(hourly_load);
scaled_hourly_load = hourly_load * scaling_factor;

PEducation = zeros(1, days * 24);
index = 1;
for day = 1:days
    if day >= 6  % Saturday and Sunday
        load_reduction_factor = 1 - 0.90;
    else
        load_reduction_factor = 1;
    end
    day_factor   = 1 + (rand * 2 - 1) * day_variability;
    daily_load   = scaled_hourly_load * day_factor * load_reduction_factor;
    hourly_loads = daily_load .* (1 + (rand(1, 24) * 2 - 1) * hour_variability);
    PEducation(index:index+23) = hourly_loads;
    index = index + 24;
end

%% Reshape Loads
PTelco      = reshape(PTelco, [], 1);
PHealthcare = reshape(PHealthcare, [], 1);
PEducation  = reshape(PEducation, [], 1);

P1_effective = PTelco;                   % Telecommunication (Highest Priority)
P2_effective = PHealthcare;              % Healthcare
P3_effective = PCommerical_HighImpact;     % High Impact Businesses
P4_effective = PDomestic;                % Household
P5_effective = PCommerical_LowImpact;      % Low Impact Businesses
P6_effective = PEducation;               % Education
P7_effective = PPublic;                  % Public (Lowest Priority)

%% Create Time Vectors
t_hourly = 1:168;                  % 168 hourly time points
% To obtain 336 half-hourly data points, extend the range to 168.5
t_halfHourly = 1:0.5:168.5;          % This yields 336 points

%% Perform Linear Interpolation for each load profile
PTelco_30mins               = interp1(t_hourly, PTelco, t_halfHourly, 'linear');
PHealthcare_30mins          = interp1(t_hourly, PHealthcare, t_halfHourly, 'linear');
PCommerical_HighImpact_30mins = interp1(t_hourly, PCommerical_HighImpact, t_halfHourly, 'linear');
PDomestic_30mins            = interp1(t_hourly, PDomestic, t_halfHourly, 'linear');
PCommerical_LowImpact_30mins = interp1(t_hourly, PCommerical_LowImpact, t_halfHourly, 'linear');
PEducation_30mins           = interp1(t_hourly, PEducation, t_halfHourly, 'linear');
PPublic_30mins              = interp1(t_hourly, PPublic, t_halfHourly, 'linear');

%% Compute Total Energy Demand
% Sum the original hourly profiles
total_hourly = PTelco + PHealthcare + PCommerical_HighImpact + ...
               PDomestic + PCommerical_LowImpact + PEducation + PPublic;

%% Compute Total Energy Demand at Half-Hourly Resolution
total_halfHourly = PTelco_30mins + PHealthcare_30mins + ...
                   PCommerical_HighImpact_30mins + PDomestic_30mins + ...
                   PCommerical_LowImpact_30mins + PEducation_30mins + PPublic_30mins;

%% Calculate Maximum Deviation from Hourly Average
% Reshape the half-hourly data into a 2-row matrix (each column represents one hour)
half_hourly_matrix = reshape(total_halfHourly, 2, []);

% Compute the hourly averages from the half-hourly data
hourly_averages = mean(half_hourly_matrix, 1);

% Expand these hourly averages back to the half-hourly resolution:
expanded_hourly_averages = repelem(hourly_averages, 2);

% Calculate the absolute deviation for each half-hourly point
deviations = abs(total_halfHourly - expanded_hourly_averages);

% The maximum deviation (in Wh) across the week
max_deviation_Wh = max(deviations);

% Convert to kWh for reporting:
max_deviation_kWh = max_deviation_Wh / 1000;

%% Calculate the Minimum SoC Reserve Required
% BESS_size is in Wh. The reserve required (in percentage) is the maximum deviation 
% divided by the total battery capacity, times 100.
slack_SoC_required = (max_deviation_Wh / BESS_size) * 100;

%% Display the Results
fprintf('Maximum deviation: %.2f kWh\n', max_deviation_kWh);
fprintf('Minimum SoC reserve required: %.2f%% of BESS capacity\n', slack_SoC_required);

%% Plot Hourly Total Energy Demand
figure;
plot(t_hourly, total_hourly, 'LineWidth',1.5);
xlabel('Time (hours)');
ylabel('Total Energy Demand (W)');
title('Hourly Total Energy Demand');
grid on;

%% Plot 30-Minute Total Energy Demand
figure;
plot(t_halfHourly, total_halfHourly, 'LineWidth',1.5);
xlabel('Time (hours)');
ylabel('Total Energy Demand (W)');
title('30-Minute Total Energy Demand');
grid on;
