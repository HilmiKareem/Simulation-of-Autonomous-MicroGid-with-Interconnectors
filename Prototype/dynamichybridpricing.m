clear;

% Initialization of parameters
PV_Module_Price = 490; % $/kW
Battery_Cost = 90; % $/kWh
Rs = 1.1; % Sizing ratio
P_max_inverter = 300e3; % Maximum inverter power in Watts
PV_peak = Rs * P_max_inverter ; % Peak PV generation in watts
n_inverter = 0.9; % Efficiency of the inverter (90%)
n_charger = 0.9; % Efficiency of the charger (90%)
BESS_size = 500000; % Battery Energy Storage System size in Wh
SoC_max = 100; % Maximum State of Charge (%)
SoC_min = 0; % Minimum State of Charge (%)
SoC = 50; % Initial State of Charge (%)
P_max_charger = 500000; % Maximum charger power in Watts
time_steps = 168; % Simulation time for 7 days (hours)
delta_t = 1; % Time interval in hours

% Define SoC thresholds for load prioritization
Threshold_1 = 35; % Threshold 1 (%)
Threshold_2 = 30; % Threshold 2 (%)
DOD = 10;    % Depth of Discharge set to minimum SoC (%)

% Time settings
hours_per_day = 24;
days = 7;
total_hours = hours_per_day * days;

%%
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

% Initialize effective demand variables
P1_effective = PDomestic;
P2_effective = PCommerical;
P3_effective = PPublic;

% Total energy demand
E_demand = (P1_effective + P2_effective + P3_effective) * delta_t;

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
Annual_Energy = 8760 * Capacity_Factor * (PV_peak / 1000); % kWh/year

% Corrected sLCOE formula
sLCOE = ((Total_Capital * CRF + OnM_Cost * (PV_peak / 1000)) / Annual_Energy) + Variable_OnM_Cost;

%% Dynamic Electricity Pricing Based on PV Output
Price_min = sLCOE * 0.8; % Minimum price (low demand)
Price_max = sLCOE * 1.5; % Maximum price (high demand)

%% Define Thresholds for Hybrid Pricing
high_generation_threshold = 0.8 * max(EPV); % 80% of peak PV output
high_demand_threshold = 0.8 * max(E_demand); % 80% of peak demand

%% Compute Electricity Price with Storage Consideration
b = 0.000001; % Nonlinear coefficient for marginal cost variation
c = 0.05;  % Sensitivity factor for storage impact on pricing
electricity_price = zeros(1, time_steps);

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

    % Hybrid Pricing Logic
    if EPV(t) >= high_generation_threshold
        if E_demand_effective(t) < high_demand_threshold
            % High generation, low demand: lower prices
            electricity_price(t) = sLCOE * 0.8;
        else
            % High generation, high demand: moderate prices
            electricity_price(t) = sLCOE * 1.0;
        end
    else
        if E_demand_effective(t) >= high_demand_threshold
            % Low generation, high demand: higher prices
            electricity_price(t) = sLCOE * 1.5;
        else
            % Low generation, low demand: base price
            electricity_price(t) = sLCOE;
        end
    end

    % Incorporate storage effects
    storage_factor = (1 - (SoC / SoC_max)); 
    electricity_price(t) = electricity_price(t) + c * storage_factor;
end

%% Plot Results
figure;

% Subplot 1: Energy demand for P1, P2, and P3 (in kW)
subplot(6, 1, 1);
plot(1:time_steps, PDomestic / 1000, '-r', 'LineWidth', 1.5); % Convert to kW
hold on;
plot(1:time_steps, PCommerical / 1000, '-g', 'LineWidth', 1.5); % Convert to kW
plot(1:time_steps, PPublic / 1000, '-b', 'LineWidth', 1.5); % Convert to kW
xlabel('Time (hours)');
ylabel('Power (kW)');
title('Power Demand for P1, P2, and P3');
legend('Domestic', 'Commercial', 'Public');
grid on;

subplot(6, 1, 2);
plot(1:time_steps, P1_effective / 1000, '-r', 'LineWidth', 1.5); % Convert to kW
hold on;
plot(1:time_steps, P2_effective / 1000, '-g', 'LineWidth', 1.5); % Convert to kW
plot(1:time_steps, P3_effective / 1000, '-b', 'LineWidth', 1.5); % Convert to kW
xlabel('Time (hours)'); ylabel('Power (kW)');
title('Effective Power Demand');
legend('P1', 'P2', 'P3');
grid on;

subplot(6, 1, 3);
plot(1:time_steps, EPV / 1000, '-g', 'LineWidth', 1.5); % Convert to kW
hold on;
plot(1:time_steps, E_demand_effective / 1000, '-r', 'LineWidth', 1.5); % Convert to kW
xlabel('Time (hours)'); ylabel('Energy (kWh)');
title('Energy Demand vs PV Generation');
legend('PV Generation', 'Effective Demand');
grid on;

subplot(6, 1, 4);
plot(1:time_steps, E_storage / 1000, '-b', 'LineWidth', 1.5); % Convert to kW
xlabel('Time (hours)'); ylabel('Energy (kWh)');
title('Storage Energy Flow');
grid on;

subplot(6, 1, 5);
plot(1:time_steps, SoC_array, '-m', 'LineWidth', 1.5);
xlabel('Time (hours)'); ylabel('SoC (%)');
title('Battery SoC');
grid on;

subplot(6, 1, 6);
plot(1:time_steps, electricity_price, '-r', 'LineWidth', 1.5);
xlabel('Time (hours)'); ylabel('Price ($/kWh)');
title('Electricity Price');
grid on;