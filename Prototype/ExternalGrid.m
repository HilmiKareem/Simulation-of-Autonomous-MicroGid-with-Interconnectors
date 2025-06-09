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
Battery_Cost      = 2160.43;        % $/kWh
Rs                = 1.2;       % Sizing ratio
P_max_inverter    = 500e3;     % Maximum inverter power (W)
PV_peak           = Rs * P_max_inverter; % Peak PV generation (W)
n_inverter        = 0.9;       % Inverter efficiency
n_charger         = 0.9;       % Charger efficiency
BESS_size         = 100e3;       % Battery capacity (Wh)
SoC_max           = 100;       % Maximum State of Charge (%)
SoC_min           = 0;         % Minimum State of Charge (%)
SoC               = 50;        % Initial State of Charge (%)
P_max_charger     = 500e3;     % Maximum charger power (W)
time_steps        = 168;       % Simulation duration (hours; 7 days)
delta_t           = 1;         % Time step (hours)
min_Trading_SoC   = 30;        % (%)

%% External Grid Parameters
% External Grid 1
ext1_total_energy = 100; % kWh over 7 days (positive for selling)
ext1_energy = ext1_total_energy * ones(1, time_steps);
%ext1_energy = ext1_energy / sum(ext1_energy) * ext1_total_energy;
ext1_price = 0.5*ones(1, time_steps) * 0.5; % $0.5 to $1.0 per kWh
ext1_trans_cap = 100; % kWh per hour

% External Grid 2
ext2_total_energy = -100; % kWh over 7 days (negative for buying)
ext2_energy = ext2_total_energy * ones(1, time_steps);
%ext2_energy = ext2_energy / sum(ext2_energy) * ext2_total_energy;
ext2_price = 0.8*ones(1, time_steps) * 0.4; % $0.8 to $1.2 per kWh
ext2_trans_cap = 100; % kWh per hour

% Initialize trading and shedding variables
energy_sold_ext1 = zeros(1, time_steps);
energy_sold_ext2 = zeros(1, time_steps);
energy_bought_ext1 = zeros(1, time_steps);
energy_bought_ext2 = zeros(1, time_steps);
trading_cost = zeros(1, time_steps);
trading_profit = zeros(1, time_steps);
energy_surplus = zeros(1, time_steps);
energy_deficit = zeros(1, time_steps);
unserved_P1 = zeros(1, time_steps);
unserved_P2 = zeros(1, time_steps);
unserved_P3 = zeros(1, time_steps);
unserved_P4 = zeros(1, time_steps);
unserved_P5 = zeros(1, time_steps);
unserved_P6 = zeros(1, time_steps);
unserved_P7 = zeros(1, time_steps);

%% Load Prioritization Thresholds
Threshold_1 = 45;  % (%)
Threshold_2 = 35;  % (%)
Threshold_3 = 25;  % (%)
Threshold_4 = 20;  % (%)
Threshold_5 = 15;  % (%)
Threshold_6 = 10;  % (%; new threshold for Healthcare)
DOD         = 5;   % Depth of Discharge (%)

%% Time Settings
hours_per_day = 24;
days          = 7;
total_hours   = hours_per_day * days;

%% Load Profiles from Excel (Crest & Clover)
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
Effective_Coverage    = Coverage_per_tower / Overlap_factor;     % m² per tower
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

%% Initialize Effective Demand (Load Prioritization Order)
P1_effective = PTelco;                   % Telecommunication (Highest Priority)
P2_effective = PHealthcare;              % Healthcare
P3_effective = PCommerical_HighImpact;     % High Impact Businesses
P4_effective = PDomestic;                % Household
P5_effective = PCommerical_LowImpact;      % Low Impact Businesses
P6_effective = PEducation;               % Education
P7_effective = PPublic;                  % Public (Lowest Priority)

% Total energy demand (Wh)
E_demand = (P1_effective + P2_effective + P3_effective + ...
            P4_effective + P5_effective + P6_effective + P7_effective) * delta_t;

% Enforce inverter capacity limit
E_max_inverter = P_max_inverter * delta_t;
E_demand(E_demand > E_max_inverter) = E_max_inverter;

%% Read and Process Irradiance Data
filename = 'irradiance_data.xlsx'; 
sheet = 1; 
irradiance_data = readmatrix(filename, 'Sheet', sheet);
irradiance_data = irradiance_data(:)';  % Ensure row vector

num_hours = length(irradiance_data);
if num_hours > time_steps
    irradiance_data = irradiance_data(1:time_steps);
elseif num_hours < time_steps
    error('The irradiance data is shorter than the required simulation time.');
end
irradiance_data(irradiance_data < 0) = 0;
irradiance_data(isnan(irradiance_data)) = 0;

%% Compute PV Generation
EPV = irradiance_data * PV_peak;
EPV = EPV * n_inverter; 
EPV = min(EPV, P_max_inverter);

%% Cost Calculation
Solar_Panel_Capital = (PV_peak / 1000) * PV_Module_Price;  % PV capital cost
Battery_Capital     = Battery_Cost * (BESS_size / 1000);   % BESS capital cost
Total_Capital       = Solar_Panel_Capital + Battery_Capital;

% Financial & Operational Parameters
i                 = 0.03;  % Interest rate (3%)
OnM_Cost = 0.015 * Total_Capital;  % 1.5% of total capital cost
Variable_OnM_Cost = 0;   % Variable O&M ($/kWh)
Operational_Years = 20;    % PV lifetime
Capacity_Factor   = 0.17;  % PV capacity factor

% Capital Recovery Factor (CRF)
CRF = (i * (1 + i)^Operational_Years) / ((1 + i)^Operational_Years - 1);

% Annual Energy Output (kWh/year)
E_PV = 8760 * Capacity_Factor * (PV_peak / 1000);  % PV-generated energy
Annual_Energy = E_PV ;  % Total usable energy

% sLCOE Calculation ($/kWh)
sLCOE = ((Total_Capital * CRF + OnM_Cost) / Annual_Energy) + Variable_OnM_Cost;

% Profit margin (5%)
Price = 1.05 * sLCOE;

%% Dynamic Electricity Pricing Based on Storage
c = 0.2;       % Sensitivity factor for storage impact
b = 1e-7;       % Proportionality constant for imbalance pricing (adjust as needed)

electricity_price = zeros(1, time_steps);

%% Prepare Arrays for Storing Results
E_storage_array         = zeros(1, time_steps);
SoC_array               = zeros(1, time_steps);
SoC_array(1)            = SoC;
delta_E_array           = zeros(1, time_steps);
Cost_of_Unserved_array  = zeros(1, time_steps);
E_demand_effective_array= zeros(1, time_steps);  % For effective demand

%% Storage Energy Flow, SoC Update, and Unserved Energy Cost Calculation
for t = 1:time_steps
    % Update current SoC
    if t == 1
        current_SoC = SoC;
    else
        current_SoC = SoC_array(t-1);
    end

    % Store actual demands before load shedding
    actual_P1_Demand = P1_effective(t);
    actual_P2_Demand = P2_effective(t);
    actual_P3_Demand = P3_effective(t);
    actual_P4_Demand = P4_effective(t);
    actual_P5_Demand = P5_effective(t);
    actual_P6_Demand = P6_effective(t);
    actual_P7_Demand = P7_effective(t);

    % Compute the storage factor
    storage_factor = 1 - (current_SoC / SoC_max);

    actual_total_demand = PTelco(t) + PHealthcare(t) + PCommerical_HighImpact(t) + ...
                      PDomestic(t) + PCommerical_LowImpact(t) + PEducation(t) + PPublic(t);

    % Compute the signed imbalance: positive if demand > supply, negative if supply > demand
    imbalance = actual_total_demand - EPV(t);

    % Compute the normalized price adjustment so that the adjustment saturates at ±k
    price_adjustment = b * imbalance ;

    % Compute the instantaneous electricity price:
    electricity_price(t) = Price + c * storage_factor + price_adjustment;

    if imbalance <= 0
        % Excess PV: Store and possibly sell
        excess_energy = EPV(t) - actual_total_demand;
    
        % Limit charging rate to maximum charger power
        delta_E = min(excess_energy, P_max_charger);
    
        % Ensure delta_E does not exceed battery capacity limits
        DeltaAE = SoC_max - current_SoC;
        if delta_E > DeltaAE
            delta_E = DeltaAE;
        end
        if delta_E < 0
            AE = current_SoC - SoC_min;
            if abs(delta_E) > AE
                delta_E = -AE;
            end
        end
    
        % Update storage and SoC
        E_storage_array(t) = delta_E; % Positive value (charging)
        SoC = current_SoC + (n_charger * delta_E / BESS_size) * 100;
        SoC_array(t) = SoC;

        % Compute effective demand (original, no shedding)
        E_demand_effective = actual_total_demand * delta_t;
        E_demand_effective_array(t) = E_demand_effective;
    
        if SoC_array(t) > min_Trading_SoC
            % Check if we can sell to external grids
            available_energy = (SoC_array(t) - SoC_min) / 100 * BESS_size;
    
            % Determine which grid to sell to first (prioritize higher price)
            if ext1_price(t) > ext2_price(t)
                % Sell to External Grid 1 first
                if ext1_energy(t) < 0 && ext1_price(t) > electricity_price(t)
                    max_sell = min([ext1_trans_cap, available_energy, abs(ext1_energy(t))]);
                    if max_sell > 0
                        SoC_array(t) = SoC_array(t) - (max_sell / BESS_size) * 100;
                        profit = (ext1_price(t) - electricity_price(t)) * max_sell;
                        trading_profit(t) = trading_profit(t) + profit;
                        energy_sold_ext1(t) = max_sell;
                        % Update E_storage_array to reflect energy sold (negative flow)
                        E_storage_array(t) = E_storage_array(t) - max_sell;
                        % Reduce available energy after selling to Grid 1
                        available_energy = available_energy - max_sell;
                    end
                end
    
                % Sell to External Grid 2 (if there's still available energy)
                if ext2_energy(t) < 0 && ext2_price(t) > electricity_price(t)
                    max_sell = min([ext2_trans_cap, available_energy, abs(ext2_energy(t))]);
                    if max_sell > 0
                        SoC_array(t) = SoC_array(t) - (max_sell / BESS_size) * 100;
                        profit = (ext2_price(t) - electricity_price(t)) * max_sell;
                        trading_profit(t) = trading_profit(t) + profit;
                        energy_sold_ext2(t) = max_sell;
                        % Update E_storage_array to reflect energy sold (negative flow)
                        E_storage_array(t) = E_storage_array(t) - max_sell;
                    end
                end
            else
                % Sell to External Grid 2 first
                if ext2_energy(t) < 0 && ext2_price(t) > electricity_price(t)
                    max_sell = min([ext2_trans_cap, available_energy, abs(ext2_energy(t))]);
                    if max_sell > 0
                        SoC_array(t) = SoC_array(t) - (max_sell / BESS_size) * 100;
                        profit = (ext2_price(t) - electricity_price(t)) * max_sell;
                        trading_profit(t) = trading_profit(t) + profit;
                        energy_sold_ext2(t) = max_sell;
                        % Update E_storage_array to reflect energy sold (negative flow)
                        E_storage_array(t) = E_storage_array(t) - max_sell;
                        % Reduce available energy after selling to Grid 2
                        available_energy = available_energy - max_sell;
                    end
                end
    
                % Sell to External Grid 1 (if there's still available energy)
                if ext1_energy(t) < 0 && ext1_price(t) > electricity_price(t)
                    max_sell = min([ext1_trans_cap, available_energy, abs(ext1_energy(t))]);
                    if max_sell > 0
                        SoC_array(t) = SoC_array(t) - (max_sell / BESS_size) * 100;
                        profit = (ext1_price(t) - electricity_price(t)) * max_sell;
                        trading_profit(t) = trading_profit(t) + profit;
                        energy_sold_ext1(t) = max_sell;
                        % Update E_storage_array to reflect energy sold (negative flow)
                        E_storage_array(t) = E_storage_array(t) - max_sell;
                    end
                end
            end
        end
    else 
    
        % Apply load prioritization logic based on current SoC
        if current_SoC <= DOD
            % Shed all loads
            P1_effective(t) = 0;
            P2_effective(t) = 0;
            P3_effective(t) = 0;
            P4_effective(t) = 0;
            P5_effective(t) = 0;
            P6_effective(t) = 0;
            P7_effective(t) = 0;
        elseif current_SoC <= Threshold_6
            P2_effective(t) = 0;
            P3_effective(t) = 0;
            P4_effective(t) = 0;
            P5_effective(t) = 0;
            P6_effective(t) = 0;
            P7_effective(t) = 0;
        elseif current_SoC <= Threshold_5
            P3_effective(t) = 0;
            P4_effective(t) = 0;
            P5_effective(t) = 0;
            P6_effective(t) = 0;
            P7_effective(t) = 0;
        elseif current_SoC <= Threshold_4
            P4_effective(t) = 0;
            P5_effective(t) = 0;
            P6_effective(t) = 0;
            P7_effective(t) = 0;
        elseif current_SoC <= Threshold_3
            P5_effective(t) = 0;
            P6_effective(t) = 0;
            P7_effective(t) = 0;
        elseif current_SoC <= Threshold_2
            P6_effective(t) = 0;
            P7_effective(t) = 0;
        elseif current_SoC <= Threshold_1
            P7_effective(t) = 0;
        end
        
        % Calculate effective total demand (Wh) for current time step
        E_demand_effective = (P1_effective(t) + P2_effective(t) + P3_effective(t) + ...
                              P4_effective(t) + P5_effective(t) + P6_effective(t) + P7_effective(t)) * delta_t;
        E_demand_effective_array(t) = E_demand_effective;  % Store for plotting
        
        % Compute difference between PV generation and effective demand
        delta_E = EPV(t) - E_demand_effective;
        delta_E_array(t) = delta_E;
        
        % Calculate unserved energy (Wh) for each load
        unserved_P1 = (actual_P1_Demand - P1_effective(t)) * delta_t;
        unserved_P2 = (actual_P2_Demand - P2_effective(t)) * delta_t;
        unserved_P3 = (actual_P3_Demand - P3_effective(t)) * delta_t;
        unserved_P4 = (actual_P4_Demand - P4_effective(t)) * delta_t;
        unserved_P5 = (actual_P5_Demand - P5_effective(t)) * delta_t;
        unserved_P6 = (actual_P6_Demand - P6_effective(t)) * delta_t;
        unserved_P7 = (actual_P7_Demand - P7_effective(t)) * delta_t;
        
        % Cost of unserved energy (convert Wh to kWh)
        Cost_of_Unserved_array(t) = electricity_price(t) * ...
            (unserved_P1*VoLL_Telecommunications + unserved_P2*VoLL_Healthcare + ...
             unserved_P3*VoLL_HighImpactCommercial + unserved_P4*VoLL_Domestic + ...
             unserved_P5*VoLL_LowImpactCommercial + unserved_P6*VoLL_Education + ...
             unserved_P7*VoLL_Public) / 1000;
        
        % Limit charging/discharging rates
        if delta_E > 0
            delta_E = min(delta_E, P_max_charger);
        end
        
        DeltaAE = SoC_max - current_SoC;
        if delta_E > DeltaAE
            delta_E = DeltaAE;
        end
        if delta_E < 0
            AE = current_SoC - SoC_min;
            if abs(delta_E) > AE
                delta_E = -AE;
            end
        end
        
        % Update storage and SoC
        E_storage_array(t) = delta_E;
        SoC = current_SoC + (n_charger * delta_E / BESS_size) * 100;
        SoC_array(t) = SoC;
    end
end

total_unserved_energy_cost = sum(Cost_of_Unserved_array);

%% Plot Results
figure;

% Subplot 1: Original Power Demand (kW)
subplot(7, 1, 1);
plot(1:time_steps, PTelco / 1000, '-r', 'LineWidth', 1.5);
hold on;
plot(1:time_steps, PHealthcare / 1000, '-g', 'LineWidth', 1.5);
plot(1:time_steps, PCommerical_HighImpact / 1000, '-b', 'LineWidth', 1.5);
plot(1:time_steps, PDomestic / 1000, '-m', 'LineWidth', 1.5);
plot(1:time_steps, PCommerical_LowImpact / 1000, '-c', 'LineWidth', 1.5);
plot(1:time_steps, PEducation / 1000, '-k', 'LineWidth', 1.5);
plot(1:time_steps, PPublic / 1000, '-y', 'LineWidth', 1.5);
xlabel('Time (hours)'); ylabel('Power (kW)');
title('Original Power Demand (No Load Shedding)');
legend('Telecommunication', 'Healthcare', 'High Impact', 'Household', 'Low Impact', 'Education', 'Public');
grid on;

% Subplot 2: Effective Power Demand with Load Shedding (kW)
subplot(7, 1, 2);
plot(1:time_steps, P1_effective / 1000, '-r', 'LineWidth', 1.5);
hold on;
plot(1:time_steps, P2_effective / 1000, '-g', 'LineWidth', 1.5);
plot(1:time_steps, P3_effective / 1000, '-b', 'LineWidth', 1.5);
plot(1:time_steps, P4_effective / 1000, '-m', 'LineWidth', 1.5);
plot(1:time_steps, P5_effective / 1000, '-c', 'LineWidth', 1.5);
plot(1:time_steps, P6_effective / 1000, '-k', 'LineWidth', 1.5);
plot(1:time_steps, P7_effective / 1000, '-y', 'LineWidth', 1.5);
xlabel('Time (hours)'); ylabel('Power (kW)');
title('Effective Power Demand with Load Shedding');
legend('Telecommunication', 'Healthcare', 'High Impact', 'Household', 'Low Impact', 'Education', 'Public');
grid on;

% Subplot 3: PV Generation vs. Effective Demand (kWh)
subplot(7, 1, 3);
plot(1:time_steps, EPV / 1000, '-g', 'LineWidth', 1.5);
hold on;
plot(1:time_steps, E_demand_effective_array / 1000, '-r', 'LineWidth', 1.5);
xlabel('Time (hours)'); ylabel('Energy (kWh)');
title('PV Generation vs. Effective Demand');
legend('PV Generation', 'Effective Demand');
grid on;

% Subplot 4: Storage Energy Flow (kWh)
subplot(7, 1, 4);
plot(1:time_steps, E_storage_array / 1000, '-b', 'LineWidth', 1.5);
xlabel('Time (hours)'); ylabel('Energy (kWh)');
title('Storage Energy Flow');
grid on;

% Subplot 5: Battery SoC (%)
subplot(7, 1, 5);
plot(1:time_steps, SoC_array, '-m', 'LineWidth', 1.5);
xlabel('Time (hours)'); ylabel('SoC (%)');
title('Battery SoC');
grid on;

% Subplot 6: Electricity Price ($/kWh)
subplot(7, 1, 6);
plot(1:time_steps, electricity_price, '-r', 'LineWidth', 1.5);
xlabel('Time (hours)'); ylabel('Price ($/kWh)');
title('Electricity Price');
grid on;

% Subplot 7: Cumulative Unserved Energy Cost ($)
cumulative_Cost = cumsum(Cost_of_Unserved_array);
subplot(7, 1, 7);
plot(1:time_steps, cumulative_Cost, '-b', 'LineWidth', 1.5);
xlabel('Time (hours)'); ylabel('Cumulative Cost ($)');
title('Cumulative Unserved Energy Cost');
grid on;

%% Plot Trading Data
figure('Name', 'Trading Analysis', 'Position', [100, 100, 1200, 800]);

% External Grid Energy Profiles
subplot(3, 2, 1);
hold on
plot(ext1_energy, 'b'); 
plot(ext2_energy, 'r'); 
legend('Grid 1 Energy', 'Grid 2 Energy');
xlabel('Time (hours)'); ylabel('kWh');
title('External Grid Energy Profile'); 
grid on;

% External Grid Energy Price
subplot(3, 2, 2);
hold on
plot(ext1_price, 'b'); 
plot(ext2_price, 'r'); 
title('External Grid Energy Price');
legend('Grid 1 Price', 'Grid 2 Price');
xlabel('Time (hours)'); ylabel('$/kWh');
grid on;

% Traded Energy (Energy Sold to External Grids)
subplot(3, 2, 3);
hold on
plot(energy_sold_ext1, 'b', 'LineWidth', 1.5);
plot(energy_sold_ext2, 'r', 'LineWidth', 1.5);
legend('Energy Sold to Grid 1', 'Energy Sold to Grid 2');
xlabel('Time (hours)'); ylabel('Energy (kWh)');
title('Traded Energy (Energy Sold to External Grids)');
grid on;

% Cumulative Trade profit
cumulative_trade_profit = cumsum(trading_profit);
subplot(3, 2, 4);
plot(cumulative_trade_profit, 'k', 'LineWidth', 1.5);
xlabel('Time (hours)'); ylabel('Cumulative profit ($)');
title('Cumulative Trade profit');
grid on;


disp(['Total Unserved Energy Cost: $', num2str(total_unserved_energy_cost)]);
