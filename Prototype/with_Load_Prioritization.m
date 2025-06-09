clear;

% Initialization of parameters

PV_Module_Price = 490; % $/kW
Rs = 1.4; % Sizing ratio
P_max_inverter = 228000; % Maximum inverter power in Watts
PV_peak = Rs * P_max_inverter ; % Peak PV generation in watts
n_inverter = 0.9; % Efficiency of the inverter (90%)
n_charger = 0.9; % Efficiency of the charger (90%)
BESS_size = 500000; % Battery Energy Storage System size in Wh
SoC_max = 100; % Maximum State of Charge (%)
SoC_min = 0; % Minimum State of Charge (%)
SoC = 50; % Initial State of Charge (%)
P_max_charger = 1000; % Maximum charger power in Watts
time_steps = 168; % Simulation time for 7 days (hours)
delta_t = 1; % Time interval in hours

% Define SoC thresholds for load prioritization
Threshold_1 = 50; % Threshold 1 (%)
Threshold_2 = 30; % Threshold 2 (%)
DOD = 10;    % Depth of Discharge set to minimum SoC (%)

% Time settings
hours_per_day = 24;
days = 7;
total_hours = hours_per_day * days;

%%

filename = 'load_data_CREST.xlsx'; 
sheet = 1; 
P1_data = readmatrix(filename, 'Sheet', sheet);
P1_data_minute_daily = P1_data(:)'; % Ensure correct shape

% Number of minutes in a day (1440 minutes = 24 hours * 60 minutes)
minutes_per_day = 1440;

% Step 1: Convert minute-level data to hourly data
% Reshape the data to have 60 minutes per row (1 hour) and sum for each hour
P1_data_hourly = reshape(P1_data_minute_daily, 60, []); % Reshape data into 60-minute blocks
P1_data_hourly = sum(P1_data_hourly, 1) / 60; % Sum the minutes per hour and divide by 60 to get average per hour

% Step 2: Replicate hourly data for 7 days
P1_daily = P1_data_hourly; % This is the hourly data for 1 day
P1_demand = repmat(P1_daily, 1, days); % Repeat it for 7 days

% Daily demand profiles (24 hours each)
%P1_daily = [20 20 20 20 30 30 40 40 50 50 50 60 40 40 30 30 50 60 70 80 80 50 40 30];
%P2_daily = [0 0 0 0 0 0 50 50 100 150 200 200 100 100 50 50 50 150 200 300 200 100 50 0];
%P3_daily = [0 0 0 0 0 0 0 0 0 200 300 400 300 200 0 0 0 200 400 500 400 300 100 0];

%P2_daily = P1_daily;
%P3_daily = P1_daily;

% Replicate daily demand for 7 days
P2_demand = repmat(P2_daily, 1, days);
P3_demand = repmat(P3_daily, 1, days);

% Initialize effective demand variables
P1_effective = P1_demand;
P2_effective = P2_demand;
P3_effective = P3_demand;

% Total energy demand
E_demand = (P1_demand + P2_demand + P3_demand) * delta_t;

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

Solar_Panel_Capital = (PV_peak / 1000) * PV_Module_Price; % Total capital cost

i = 0.03; % Interest Rate (3%)
OnM_Cost = 25; % Fixed O&M cost ($/kW-yr)
Variable_OnM_Cost = 0.002; % Variable O&M cost ($/kWh)
Operational_Years = 20;
Capacity_Factor = 0.3; % Assumed capacity factor

% Capital Recovery Factor (CRF) Calculation
CRF = (i * (1 + i)^Operational_Years) / ((1 + i)^Operational_Years - 1);

% Annual Energy Output (in kWh)
Annual_Energy = 8760 * Capacity_Factor * (PV_peak / 1000); % kWh/year

% Corrected sLCOE formula
sLCOE = ((Solar_Panel_Capital * CRF + OnM_Cost) / Annual_Energy) + Variable_OnM_Cost;

%% Dynamic Electricity Pricing Based on PV Output
Price_min = sLCOE * 0.8; % Minimum price (low demand)
Price_max = sLCOE * 1.5; % Maximum price (high demand)

%% Compute Electricity Price with Storage Consideration

b = 0.001; % Nonlinear coefficient for marginal cost variation
c = 0.05;  % Sensitivity factor for storage impact on pricing
electricity_price_marginal = zeros(1, time_steps);

%% Storage Energy Flow and SoC Update
E_storage = zeros(1, time_steps);
SoC_array = zeros(1, time_steps);
SoC_array(1) = SoC;
delta_E = zeros(1, time_steps);

for t = 1:time_steps
    % Observe current SoC
    if t == 1
        current_SoC = SoC;
    else
        current_SoC = SoC_array(t-1);
    end
    
    % Apply load prioritization logic
    if current_SoC <= DOD
        P1_effective(t) = 0;
        P2_effective(t) = 0;
        P3_effective(t) = 0;
    elseif current_SoC <= Threshold_2
        P2_effective(t) = 0;
        P3_effective(t) = 0;
    elseif current_SoC <= Threshold_1
        P3_effective(t) = 0;
    end

    % Recalculate effective total demand
    E_demand_effective = (P1_effective + P2_effective + P3_effective) * delta_t;
    delta_E(t) = EPV(t) - E_demand_effective(t);

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

    % Normalize PV output (0 to 1) for price scaling
    PV_utilization = EPV / max(EPV);
    Electricity_Price = Price_max - PV_utilization * (Price_max - Price_min);

    % Compute electricity price considering storage level
    storage_factor = (1 - (SoC / SoC_max)); 
    electricity_price_marginal(t) = sLCOE + 2 * b * (EPV(t) / 1000) + c * storage_factor;
end

%% Plot Results
figure;

% Subplot 1: Energy demand for P1, P2, and P3
subplot(6, 1, 1);
plot(1:time_steps, P1_demand, '-r', 'LineWidth', 1.5);
hold on;
plot(1:time_steps, P2_demand, '-g', 'LineWidth', 1.5);
plot(1:time_steps, P3_demand, '-b', 'LineWidth', 1.5);
xlabel('Time (hours)');
ylabel('Power (W)');
title('Power Demand for P1, P2, and P3');
legend('P1 Demand', 'P2 Demand', 'P3 Demand');
grid on;

subplot(6, 1, 2);
plot(1:time_steps, P1_effective, '-r', 'LineWidth', 1.5);
hold on;
plot(1:time_steps, P2_effective, '-g', 'LineWidth', 1.5);
plot(1:time_steps, P3_effective, '-b', 'LineWidth', 1.5);
xlabel('Time (hours)'); ylabel('Power (W)');
title('Effective Power Demand');
legend('P1', 'P2', 'P3');
grid on;

subplot(6, 1, 3);
plot(1:time_steps, EPV, '-g', 'LineWidth', 1.5);
hold on;
plot(1:time_steps, E_demand_effective, '-r', 'LineWidth', 1.5);
xlabel('Time (hours)'); ylabel('Energy (Wh)');
title('Energy Demand vs PV Generation');
legend('PV Generation', 'Effective Demand');
grid on;

subplot(6, 1, 4);
plot(1:time_steps, E_storage, '-b', 'LineWidth', 1.5);
xlabel('Time (hours)'); ylabel('Energy (Wh)');
title('Storage Energy Flow');
grid on;

subplot(6, 1, 5);
plot(1:time_steps, SoC_array, '-m', 'LineWidth', 1.5);
xlabel('Time (hours)'); ylabel('SoC (%)');
title('Battery SoC');
grid on;

subplot(6, 1, 6);
plot(1:time_steps, electricity_price_marginal, '-r', 'LineWidth', 1.5);
xlabel('Time (hours)'); ylabel('Price ($/kWh)');
title('Electricity Price');
grid on;
