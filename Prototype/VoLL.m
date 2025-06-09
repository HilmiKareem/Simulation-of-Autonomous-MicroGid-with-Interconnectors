% Defining VoLL

VoLL_Domestic = 2; %2 to 4
VoLL_HighImpactCommercial = 8; %8 to 12
VoLL_LowImpactCommercial = 4; % 4 to 6
VoLL_Education = 1.5; %1.5 to 3
VoLL_Healthcare = 10; % 10 to 15
VoLL_Public = 1; % 1 to 2
VoLL_Telecommunications = 12; % 12 to 20 

% Initialization of parameters
PV_Module_Price = 490; % $/kW
Battery_Cost = 90; % $/kWh
Rs = 1.2; % Sizing ratio
P_max_inverter = 200e3; % Maximum inverter power in Watts
PV_peak = Rs * P_max_inverter ; % Peak PV generation in watts
n_inverter = 0.9; % Efficiency of the inverter (90%)
n_charger = 0.9; % Efficiency of the charger (90%)
BESS_size = 1e3; % Battery Energy Storage System size in Wh
SoC_max = 100; % Maximum State of Charge (%)
SoC_min = 0; % Minimum State of Charge (%)
SoC = 50; % Initial State of Charge (%)
P_max_charger = 500e3; % Maximum charger power in Watts
time_steps = 168; % Simulation time for 7 days (hours)
delta_t = 1; % Time interval in hours

% Define SoC thresholds for load prioritization
Threshold_1 = 45; % Threshold 1 (%)
Threshold_2 = 35; % Threshold 2 (%)
Threshold_3 = 25; % New threshold for additional load levels
Threshold_4 = 20;
Threshold_5 = 15;
Threshold_6 = 10; % New threshold for Healthcare
DOD = 5;    % Depth of Discharge set to minimum SoC (%)

% Time settings
hours_per_day = 24;
days = 7;
total_hours = hours_per_day * days;

%% Crest & Clover Load Profiles
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

% Split Commercial Load into High Impact and Low Impact Businesses
PCommerical_HighImpact = 0.2 * PCommerical; % 20% High Impact Businesses
PCommerical_LowImpact = 0.8 * PCommerical; % 80% Low Impact Businesses

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

%% Healthcare Load Profile
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

%% Education Facility
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
PEducation = zeros(1, num_days * 24);

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
    PEducation(index:index + 23) = hourly_loads;
    index = index + 24;
end

% Reshape data if necessary
PTelco = reshape(PTelco, [], 1);
PHealthcare = reshape(PHealthcare, [], 1);
PEducation = reshape(PEducation, [], 1);

% Initialize effective demand variables
P1_effective = PTelco; % Telecommunication (Highest Priority)
P2_effective = PHealthcare; % Healthcare
P3_effective = PCommerical_HighImpact; % High Impact Businesses
P4_effective = PDomestic; % Household
P5_effective = PCommerical_LowImpact; % Low Impact Businesses
P6_effective = PEducation; % Education
P7_effective = PPublic; % Public (Lowest Priority)

% Total energy demand
E_demand = (P1_effective + P2_effective + P3_effective + P4_effective + P5_effective + P6_effective + P7_effective) * delta_t;

% Inverter capacity limit
E_max_inverter = P_max_inverter * delta_t;
E_demand(E_demand > E_max_inverter) = E_max_inverter;

%% Read and Process Irradiance Data
filename = 'irradiance_data.xlsx'; 
sheet = 1; 
irradiance_data = readmatrix(filename, 'Sheet', sheet);
irradiance_data = irradiance_data(:)'; % Ensure correct shape

% Check data length
num_hours = length(irradiance_data);
if num_hours > time_steps
    irradiance_data = irradiance_data(1:time_steps);
elseif num_hours < time_steps
    error('The irradiance data is shorter than the required simulation time.');
end

% Remove negative values or NaNs
irradiance_data(irradiance_data < 0) = 0;
irradiance_data(isnan(irradiance_data)) = 0;

%% Compute PV Generation
EPV = irradiance_data * PV_peak;
EPV = EPV * n_inverter; 
EPV = min(EPV, P_max_inverter);

%% Cost Calculation
Solar_Panel_Capital = (PV_peak / 1000) * PV_Module_Price; 
Battery_Capital = Battery_Cost*(BESS_size/1000);

Total_Capital = Solar_Panel_Capital+Battery_Capital;% Total capital cost

i = 0.03; % Interest Rate (3%)
OnM_Cost = 25; % Fixed O&M cost ($/kW-yr)
Variable_OnM_Cost = 0.4; % Variable O&M cost ($/kWh)
Operational_Years = 20;
Capacity_Factor = 0.151; % Assumed capacity factor

% Capital Recovery Factor (CRF) Calculation
CRF = (i * (1 + i)^Operational_Years) / ((1 + i)^Operational_Years - 1);

% Annual Energy Output (in kWh)
Annual_Energy = 8760 * Capacity_Factor * (PV_peak / 1000); % Wh/year

% sLCOE formula
sLCOE = ((Total_Capital * CRF + OnM_Cost * ((PV_peak/1000) + (BESS_size/1000))) / Annual_Energy) + Variable_OnM_Cost;

%% Dynamic Electricity Pricing Based on PV Output
Price_min = sLCOE * 0.8; % Minimum price (low demand)
Price_max = sLCOE * 1.5; % Maximum price (high demand)

%% Compute Electricity Price with Storage Consideration
b = 0.000001; % Nonlinear coefficient for marginal cost variation
c = 0.05;  % Sensitivity factor for storage impact on pricing
electricity_price = zeros(1, time_steps);

%% Storage Energy Flow and SoC Update
E_storage = zeros(1, time_steps);
SoC_array = zeros(1, time_steps);
SoC_array(1) = SoC;
delta_E = zeros(1, time_steps);
Cost_of_Unserved = zeros(1, time_steps);

for t = 1:time_steps
    % Observe current SoC
    if t == 1
        current_SoC = SoC;
    else
        current_SoC = SoC_array(t-1);
    end

    % Normalize PV output (0 to 1) for price scaling
    PV_utilization = EPV / max(EPV);
    Electricity_Price = Price_max - PV_utilization * (Price_max - Price_min);

    % Compute electricity price considering storage level
    storage_factor = (1 - (SoC / SoC_max)); 
    electricity_price(t) = sLCOE + 2 * b * (EPV(t) / 1000) + c * storage_factor;

    %temporary
    actual_P1_Demand = P1_effective(t);
    actual_P2_Demand = P2_effective(t);
    actual_P3_Demand = P3_effective(t);
    actual_P4_Demand = P4_effective(t);
    actual_P5_Demand = P5_effective(t);
    actual_P6_Demand = P6_effective(t);
    actual_P7_Demand = P7_effective(t);
    
    % Apply load prioritization logic
    if current_SoC <= DOD
        % Shed all loads (critical SoC level)
        P1_effective(t) = 0; % Telecommunication
        P2_effective(t) = 0; % Healthcare
        P3_effective(t) = 0; % High Impact Businesses
        P4_effective(t) = 0; % Household
        P5_effective(t) = 0; % Low Impact Businesses
        P6_effective(t) = 0; % Education
        P7_effective(t) = 0; % Public
    elseif current_SoC <= Threshold_6
        % Shed Healthcare and all lower-priority loads
        P2_effective(t) = 0; % Healthcare
        P3_effective(t) = 0; % High Impact Businesses
        P4_effective(t) = 0; % Household
        P5_effective(t) = 0; % Low Impact Businesses
        P6_effective(t) = 0; % Education
        P7_effective(t) = 0; % Public
    elseif current_SoC <= Threshold_5
        % Shed High Impact Businesses and all lower-priority loads
        P3_effective(t) = 0; % High Impact Businesses
        P4_effective(t) = 0; % Household
        P5_effective(t) = 0; % Low Impact Businesses
        P6_effective(t) = 0; % Education
        P7_effective(t) = 0; % Public
    elseif current_SoC <= Threshold_4
        % Shed Household and all lower-priority loads
        P4_effective(t) = 0; % Household
        P5_effective(t) = 0; % Low Impact Businesses
        P6_effective(t) = 0; % Education
        P7_effective(t) = 0; % Public
    elseif current_SoC <= Threshold_3
        % Shed Low Impact Businesses and all lower-priority loads
        P5_effective(t) = 0; % Low Impact Businesses
        P6_effective(t) = 0; % Education
        P7_effective(t) = 0; % Public
    elseif current_SoC <= Threshold_2
        % Shed Education and Public
        P6_effective(t) = 0; % Education
        P7_effective(t) = 0; % Public
    elseif current_SoC <= Threshold_1
        % Shed Public
        P7_effective(t) = 0; % Public
    end

    % Recalculate effective total demand
    E_demand_effective = (P1_effective + P2_effective + P3_effective + P4_effective + P5_effective + P6_effective + P7_effective) * delta_t;
    delta_E(t) = EPV(t) - E_demand_effective(t);

    % Calculate unserved energy
    unserved_P1 = (actual_P1_Demand - P1_effective(t))* delta_t;
    unserved_P2 = (actual_P2_Demand - P2_effective(t))* delta_t;
    unserved_P3 = (actual_P3_Demand - P3_effective(t))* delta_t;
    unserved_P4 = (actual_P4_Demand - P4_effective(t))* delta_t;
    unserved_P5 = (actual_P5_Demand - P5_effective(t))* delta_t;
    unserved_P6 = (actual_P6_Demand - P6_effective(t))* delta_t;
    unserved_P7 = (actual_P7_Demand - P7_effective(t))* delta_t;

    Cost_of_Unserved(t) = electricity_price(t) * (unserved_P1*VoLL_Telecommunications + unserved_P2*VoLL_Healthcare + unserved_P3*VoLL_HighImpactCommercial + unserved_P4*VoLL_Domestic + unserved_P5*VoLL_LowImpactCommercial + unserved_P6*VoLL_Education + unserved_P7*VoLL_Public) / 1000;

    % Apply maximum charger power limit during charging
    if delta_E(t) > 0 
        delta_E(t) = min(delta_E(t), P_max_charger);
    end

    % Ensure charging does not exceed SoC_max
    DeltaAE_t = SoC_max - current_SoC;
    if delta_E(t) > DeltaAE_t
        delta_E(t) = DeltaAE_t;
    end

    % Ensure discharging does not fall below SoC_min
    if delta_E(t) < 0
        AE_t = current_SoC - SoC_min;
        if abs(delta_E(t)) > AE_t
            delta_E(t) = -AE_t;
        end
    end

    % Update storage energy and SoC
    E_storage(t) = delta_E(t);
    SoC = current_SoC + (n_charger*delta_E(t) / BESS_size) * 100;
    SoC_array(t) = SoC;

end

total_unserved_energy_cost = sum(Cost_of_Unserved);

%% Plot Results
figure;

% Subplot 1: Original Power Demand (in W)
subplot(7, 1, 1);
plot(1:time_steps, PTelco/ 1000, '-r', 'LineWidth', 1.5); % Telecommunication (P1)
hold on;
plot(1:time_steps, PHealthcare/ 1000, '-g', 'LineWidth', 1.5); % Healthcare (P2)
plot(1:time_steps, PCommerical_HighImpact/ 1000, '-b', 'LineWidth', 1.5); % High Impact Businesses (P3)
plot(1:time_steps, PDomestic/ 1000, '-m', 'LineWidth', 1.5); % Household (P4)
plot(1:time_steps, PCommerical_LowImpact/ 1000, '-c', 'LineWidth', 1.5); % Low Impact Businesses (P5)
plot(1:time_steps, PEducation/ 1000, '-k', 'LineWidth', 1.5); % Education (P6)
plot(1:time_steps, PPublic/ 1000, '-y', 'LineWidth', 1.5); % Public (P7)
xlabel('Time (hours)');
ylabel('Power (kW)');
title('Original Power Demand (No Load Shedding)');
legend('Telecommunication (P1)', 'Healthcare (P2)', 'High Impact Businesses (P3)', ...
       'Household (P4)', 'Low Impact Businesses (P5)', 'Education (P6)', 'Public (P7)');
grid on;

% Subplot 2: Effective Power Demand with Load Shedding (in W)
subplot(7, 1, 2);
plot(1:time_steps, P1_effective/ 1000, '-r', 'LineWidth', 1.5); % Telecommunication (P1)
hold on;
plot(1:time_steps, P2_effective/ 1000, '-g', 'LineWidth', 1.5); % Healthcare (P2)
plot(1:time_steps, P3_effective/ 1000, '-b', 'LineWidth', 1.5); % High Impact Businesses (P3)
plot(1:time_steps, P4_effective/ 1000, '-m', 'LineWidth', 1.5); % Household (P4)
plot(1:time_steps, P5_effective/ 1000, '-c', 'LineWidth', 1.5); % Low Impact Businesses (P5)
plot(1:time_steps, P6_effective/ 1000, '-k', 'LineWidth', 1.5); % Education (P6)
plot(1:time_steps, P7_effective/ 1000, '-y', 'LineWidth', 1.5); % Public (P7)
xlabel('Time (hours)');
ylabel('Power (kW)');
title('Effective Power Demand with Load Shedding');
legend('Telecommunication (P1)', 'Healthcare (P2)', 'High Impact Businesses (P3)', ...
       'Household (P4)', 'Low Impact Businesses (P5)', 'Education (P6)', 'Public (P7)');
grid on;

% Subplot 3: Energy Demand vs PV Generation
subplot(7, 1, 3);
plot(1:time_steps, EPV / 1000, '-g', 'LineWidth', 1.5); % Convert to kW
hold on;
plot(1:time_steps, E_demand_effective / 1000, '-r', 'LineWidth', 1.5); % Convert to kW
xlabel('Time (hours)'); ylabel('Energy (kWh)');
title('Energy Demand vs PV Generation');
legend('PV Generation', 'Effective Demand');
grid on;

% Subplot 4: Storage Energy Flow
subplot(7, 1, 4);
plot(1:time_steps, E_storage / 1000, '-b', 'LineWidth', 1.5); % Convert to kW
xlabel('Time (hours)'); ylabel('Energy (kWh)');
title('Storage Energy Flow');
grid on;

% Subplot 5: Battery SoC
subplot(7, 1, 5);
plot(1:time_steps, SoC_array, '-m', 'LineWidth', 1.5);
xlabel('Time (hours)'); ylabel('SoC (%)');
title('Battery SoC');
grid on;

% Subplot 6: Electricity Price
subplot(7, 1, 6);
plot(1:time_steps, electricity_price, '-r', 'LineWidth', 1.5);
xlabel('Time (hours)'); ylabel('Price ($/kWh)');
title('Electricity Price');
grid on;

% Compute cumulative cost
cumulative_Cost = cumsum(Cost_of_Unserved);

% Subplot 7: Cumulative Unserved Energy Cost 
subplot(7, 1, 7);
plot(1:time_steps, cumulative_Cost, '-b', 'LineWidth', 1.5);
xlabel('Time (hours)');
ylabel('Cumulative Cost ($)');
title('Cumulative Unserved Energy Cost');
grid on;

disp(['Total Unserved Energy Cost: $', num2str(total_unserved_energy_cost)]);
