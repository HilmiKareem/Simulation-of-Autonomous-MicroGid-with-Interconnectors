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
min_Trading_SoC   = 30;            % Minimum SoC for trading (%)
DOD               = 5;             % Depth of Discharge (%)
%% Time Settings
hours_per_day = 24;
days          = 7;
total_hours   = hours_per_day * days;
time_steps    = total_hours;   % Simulation duration (hours; 7 days)
delta_t       = 1;             % Time step (hours)

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
max_deviation_Wh = max(deviations)* delta_t;

% Convert to kWh for reporting:
max_deviation_kWh = max_deviation_Wh / 1000;

%% Calculate the Minimum SoC Reserve Required
% BESS_size is in Wh. The reserve required (in percentage) is the maximum deviation 
% divided by the total battery capacity, times 100.
slack_SoC_required = (max_deviation_Wh / BESS_size) * 100;
min_allowed_SoC = DOD + slack_SoC_required; % both in percentage

fprintf('Maximum deviation: %.2f kWh\n', max_deviation_kWh);
fprintf('Minimum SoC reserve required: %.2f%% of BESS capacity\n', min_allowed_SoC);

% Initialize arrays to store served loads and non-served loads (in Wh)
served_load_P1_array = zeros(1, time_steps);
served_load_P2_array = zeros(1, time_steps);
served_load_P3_array = zeros(1, time_steps);
served_load_P4_array = zeros(1, time_steps);
served_load_P5_array = zeros(1, time_steps);
served_load_P6_array = zeros(1, time_steps);
served_load_P7_array = zeros(1, time_steps);

NSL_P1_array = zeros(1, time_steps);
NSL_P2_array = zeros(1, time_steps);
NSL_P3_array = zeros(1, time_steps);
NSL_P4_array = zeros(1, time_steps);
NSL_P5_array = zeros(1, time_steps);
NSL_P6_array = zeros(1, time_steps);
NSL_P7_array = zeros(1, time_steps);

imbalance_array = zeros(1, time_steps);

%% Read and Process Irradiance Data (Replace with actual data)
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
b = 1e-7;       % Proportionality constant for imbalance pricing

electricity_price = zeros(1, time_steps);

%% Main Simulation Loop
for t = 1:time_steps
    % Update current SoC
    if t == 1
        current_SoC = SoC;
    else
        current_SoC = SoC_array(t-1);
    end


    % Compute total demand and imbalance
    actual_total_demand = PTelco(t) + PHealthcare(t) + PCommerical_HighImpact(t) + ...
                          PDomestic(t) + PCommerical_LowImpact(t) + PEducation(t) + PPublic(t);
    imbalance = actual_total_demand - EPV(t);

    % Compute electricity price
    storage_factor = 1 - (current_SoC / SoC_max);
    price_adjustment = b * imbalance;
    electricity_price(t) = Price + c * storage_factor + price_adjustment;
    
    imbalance_array(t) = imbalance;

    if imbalance <= 0
        %% Excess Energy: Charge Battery or Sell
        excess_energy = EPV(t) - actual_total_demand;
        
        % Charge battery (DC energy)
        delta_E = min([excess_energy * n_charger, ...        % Apply charger efficiency
                      (SoC_max - current_SoC)/100 * BESS_size, ...
                      P_max_charger]);
        
        SoC_array(t) = current_SoC + (delta_E / BESS_size) * 100;
        E_storage_array(t) = delta_E;  % Track DC charging
        
        % Update effective demand (no shedding)
        E_demand_effective_array(t) = actual_total_demand * delta_t;
        
        % Sell excess energy if SoC > min_Trading_SoC
        if SoC_array(t) > min_Trading_SoC
            available_energy = (SoC_array(t) - DOD)/100 * BESS_size;  % Respect DOD
            
            % Prioritize higher-priced grid
            if ext1_price(t) > ext2_price(t)
                % Sell to Grid 1
                if ext1_energy(t) < 0  % Grid 1 is a buyer
                    max_sell = min([ext1_trans_cap, ...
                                    available_energy * n_inverter, ...  % AC energy
                                   abs(ext1_energy(t))]);
                    energy_taken = max_sell / n_inverter;  % DC energy
                    if max_sell > 0
                        SoC_array(t) = SoC_array(t) - (energy_taken / BESS_size) * 100;
                        trading_profit(t) = trading_profit(t) + (ext1_price(t) - electricity_price(t)) * max_sell;
                        energy_sold_ext1(t) = max_sell;
                    end
                end
            else
                % Sell to Grid 2
                if ext2_energy(t) < 0  % Grid 2 is a buyer
                    max_sell = min([ext2_trans_cap, ...
                                   available_energy * n_inverter, ...  % AC energy
                                   abs(ext2_energy(t))]);
                    energy_taken = max_sell / n_inverter;  % DC energy
                    if max_sell > 0
                        SoC_array(t) = SoC_array(t) - (energy_taken / BESS_size) * 100;
                        trading_profit(t) = trading_profit(t) + (ext2_price(t) - electricity_price(t)) * max_sell;
                        energy_sold_ext2(t) = max_sell;
                    end
                end
            end
        end
        E_demand_effective_array(t) = actual_total_demand * delta_t;
    else
        %% Optimization-Based Deficit Handling
        % Define VoLL for each load (sorted by priority)
        VoLL = [VoLL_Telecommunications;    % P1 (highest priority)
                VoLL_Healthcare;
                VoLL_HighImpactCommercial;
                VoLL_Domestic;
                VoLL_LowImpactCommercial;
                VoLL_Education;
                VoLL_Public];               % P7 (lowest priority)

        % Load demands at timestep t (convert to kWh)
        load_demands = [PTelco(t); PHealthcare(t); PCommerical_HighImpact(t);
                        PDomestic(t); PCommerical_LowImpact(t); PEducation(t); 
                        PPublic(t)] / 1000; 

        % --- Check Seller Grids ---
        seller_grids = [];
        if ext1_total_energy > 0  % Grid 1 is a seller
            seller_grids = [seller_grids; struct('price', ext1_price(t), ...
                                               'trans_cap', ext1_trans_cap, ...
                                               'energy', ext1_energy(t), ...
                                               'grid_num', 1)];
        end
        if ext2_total_energy > 0  % Grid 2 is a seller
            seller_grids = [seller_grids; struct('price', ext2_price(t), ...
                                               'trans_cap', ext2_trans_cap, ...
                                               'energy', ext2_energy(t), ...
                                               'grid_num', 2)];
        end

        % Number of seller grids
        num_seller_grids = length(seller_grids);

        % Decision variables: [battery_usage_DC, grid_import_1, ..., grid_import_N, non_served_loads_1, ..., non_served_loads_7]
        num_loads = length(VoLL);
        num_vars = 1 + num_seller_grids + num_loads; 

        % --- Objective Function ---
        % Modify this section in your code
        f = [1.5*electricity_price(t) * n_inverter, ...    % Battery usage (DC->AC efficiency)
             [seller_grids.price], ...                 % Grid import costs (row vector)
             (VoLL * electricity_price(t))'];          % VoLL for each load (scalar multiplication)

        % --- Inequality Constraints (Upper Bounds) ---
        % 1. Battery usage <= min(available energy, max discharge rate)
        available_battery_energy = max(0, ((current_SoC - min_allowed_SoC)/100) * BESS_size / 1000); % in kWh (DC)
        max_discharge_power = (P_max_inverter / 1000) * delta_t;                 % kWh (AC)
        max_battery_usage_DC = max_discharge_power / n_inverter;                 % kWh (DC)
        battery_upper = min(available_battery_energy, max_battery_usage_DC);

        % 2. Grid imports <= transmission capacity and available energy
        grid_upper = [seller_grids.trans_cap];

        % 3. Non-served loads <= actual load demand (kWh)
        non_served_upper = load_demands;

        % Combine upper bounds
        lb = zeros(1, num_vars);       % All variables >= 0
        ub = [battery_upper, grid_upper, non_served_upper'];

        % --- Equality Constraint (Energy Balance) ---
        % n_inverter * battery_usage + sum(grid_imports) + sum(non_served) = imbalance
        Aeq = [n_inverter, ones(1, num_seller_grids), ones(1, num_loads)];
        beq = imbalance / 1000;     % Convert imbalance to kWh

        % First, calculate available battery energy before optimization
        available_battery_energy = max(0, ((current_SoC - min_allowed_SoC)/100) * BESS_size / 1000);
        
        if available_battery_energy <= 0 && isempty(seller_grids)
            % If no battery energy available and no seller grids, must shed load
            x = [];
            exitflag = -1;  % Force load shedding path
        else
            % Run optimization with updated constraints
            % Update battery upper bound to reflect actual available energy
            battery_upper = min(available_battery_energy, max_battery_usage_DC);
            ub = [battery_upper, grid_upper, non_served_upper'];
            
            % Run optimization
            options = optimset('Display', 'off');
            [x, fval, exitflag] = linprog(f, [], [], Aeq, beq, lb, ub, options);
        end
        
        % Handle optimization results
        if ~isempty(x) && exitflag == 1
            % Extract optimization results
            battery_usage_DC = x(1);                     % kWh (DC)
            grid_imports = x(2:1+num_seller_grids);        % kWh (AC)
            non_served_loads = x(2+num_seller_grids:end);  % kWh (AC), one per load profile
        
            battery_energy_used_Wh = battery_usage_DC * 1000;  
            grid_imports_Wh = sum(grid_imports) * 1000;       
        
            % Retrieve the actual individual demands (in Wh)
            actual_P1_Demand = P1_effective(t);  % Already in W (and for delta_t = 1, W = Wh)
            actual_P2_Demand = P2_effective(t);
            actual_P3_Demand = P3_effective(t);
            actual_P4_Demand = P4_effective(t);
            actual_P5_Demand = P5_effective(t);
            actual_P6_Demand = P6_effective(t);
            actual_P7_Demand = P7_effective(t);
        
            % Compute the served load per profile by subtracting non-served load.
            % non_served_loads are in kWh so multiply by 1000 to get Wh.
            served_load_P1 = actual_P1_Demand * delta_t - non_served_loads(1) * 1000;
            served_load_P2 = actual_P2_Demand * delta_t - non_served_loads(2) * 1000;
            served_load_P3 = actual_P3_Demand * delta_t - non_served_loads(3) * 1000;
            served_load_P4 = actual_P4_Demand * delta_t - non_served_loads(4) * 1000;
            served_load_P5 = actual_P5_Demand * delta_t - non_served_loads(5) * 1000;
            served_load_P6 = actual_P6_Demand * delta_t - non_served_loads(6) * 1000;
            served_load_P7 = actual_P7_Demand * delta_t - non_served_loads(7) * 1000;
            
            % The effective served demand is the sum of the served loads
            served_demand = served_load_P1 + served_load_P2 + served_load_P3 + ...
                            served_load_P4 + served_load_P5 + served_load_P6 + served_load_P7;
                        
            % Calculate total available power (Wh) from PV, battery (after inverter losses), and grid imports
            total_available_power = EPV(t) + battery_energy_used_Wh * n_inverter + grid_imports_Wh;
            
            % Verify energy balance based on served (effective) demand
            energy_deficit = served_demand - total_available_power;
            if energy_deficit > 1  % Tolerance of 1 Wh
                fprintf('Energy balance violated at time %d: Deficit = %.2f Wh\n', t, energy_deficit);
                % Optionally, adjust load shedding further if desired.
            end
        
            % Update Battery State-of-Charge and storage flow
            SoC_array(t) = current_SoC - (battery_energy_used_Wh / BESS_size) * 100;
            E_storage_array(t) = -battery_energy_used_Wh;  % Negative sign for discharge
        
            % Record the effective served demand (Wh)
            E_demand_effective_array(t) = served_demand;
            
            % Store individual served loads and non-served loads (convert non-served to Wh)
            served_load_P1_array(t) = served_load_P1;
            served_load_P2_array(t) = served_load_P2;
            served_load_P3_array(t) = served_load_P3;
            served_load_P4_array(t) = served_load_P4;
            served_load_P5_array(t) = served_load_P5;
            served_load_P6_array(t) = served_load_P6;
            served_load_P7_array(t) = served_load_P7;
            
            NSL_P1_array(t) = non_served_loads(1) * 1000;
            NSL_P2_array(t) = non_served_loads(2) * 1000;
            NSL_P3_array(t) = non_served_loads(3) * 1000;
            NSL_P4_array(t) = non_served_loads(4) * 1000;
            NSL_P5_array(t) = non_served_loads(5) * 1000;
            NSL_P6_array(t) = non_served_loads(6) * 1000;
            NSL_P7_array(t) = non_served_loads(7) * 1000;
        
            % Update trading cost and grid purchase information (remains unchanged)
            for g = 1:num_seller_grids
                if seller_grids(g).grid_num == 1
                    energy_bought_ext1(t) = grid_imports(g) * 1000;
                    trading_cost(t) = trading_cost(t) + grid_imports(g) * seller_grids(g).price;
                elseif seller_grids(g).grid_num == 2
                    energy_bought_ext2(t) = grid_imports(g) * 1000;
                    trading_cost(t) = trading_cost(t) + grid_imports(g) * seller_grids(g).price;
                end
            end
        
            % Update unserved energy cost using the individual non-served loads
            Cost_of_Unserved_array(t) = sum(non_served_loads .* VoLL .* electricity_price(t));
            
        else
            % Fallback: No feasible solution, implement load shedding directly.
            % In this branch, only PV is available.
            SoC_array(t) = current_SoC;
            E_storage_array(t) = 0;
            total_available_power = EPV(t);
            energy_deficit = actual_total_demand * delta_t - total_available_power;
            non_served_loads = distribute_deficit(load_demands, VoLL, energy_deficit/1000);
            
            % Compute served load individually for each profile
            served_load_P1 = P1_effective(t) * delta_t - non_served_loads(1) * 1000;
            served_load_P2 = P2_effective(t) * delta_t - non_served_loads(2) * 1000;
            served_load_P3 = P3_effective(t) * delta_t - non_served_loads(3) * 1000;
            served_load_P4 = P4_effective(t) * delta_t - non_served_loads(4) * 1000;
            served_load_P5 = P5_effective(t) * delta_t - non_served_loads(5) * 1000;
            served_load_P6 = P6_effective(t) * delta_t - non_served_loads(6) * 1000;
            served_load_P7 = P7_effective(t) * delta_t - non_served_loads(7) * 1000;
            
            served_demand = served_load_P1 + served_load_P2 + served_load_P3 + ...
                            served_load_P4 + served_load_P5 + served_load_P6 + served_load_P7;
            energy_deficit = served_demand - total_available_power;
            
            E_demand_effective_array(t) = served_demand;
            
            served_load_P1_array(t) = served_load_P1;
            served_load_P2_array(t) = served_load_P2;
            served_load_P3_array(t) = served_load_P3;
            served_load_P4_array(t) = served_load_P4;
            served_load_P5_array(t) = served_load_P5;
            served_load_P6_array(t) = served_load_P6;
            served_load_P7_array(t) = served_load_P7;
            
            NSL_P1_array(t) = non_served_loads(1) * 1000;
            NSL_P2_array(t) = non_served_loads(2) * 1000;
            NSL_P3_array(t) = non_served_loads(3) * 1000;
            NSL_P4_array(t) = non_served_loads(4) * 1000;
            NSL_P5_array(t) = non_served_loads(5) * 1000;
            NSL_P6_array(t) = non_served_loads(6) * 1000;
            NSL_P7_array(t) = non_served_loads(7) * 1000;
            
            Cost_of_Unserved_array(t) = sum(non_served_loads .* VoLL .* electricity_price(t));
        end
    end     
end

%% Calculate Total Unserved Energy Cost
total_unserved_energy_cost = sum(Cost_of_Unserved_array);

%% Plot Results
% Then modify your plotting section:
figure('Name', 'System Performance', 'Position', [100, 100, 1200, 1000]);

% Subplot 1: Original Power Demand (kW)
subplot(6, 1, 1);  % Changed from 7,1,1 to 6,1,1
plot(1:time_steps, PTelco / 1000, '-r', 'LineWidth', 1.5);
hold on;
plot(1:time_steps, PHealthcare / 1000, '-g', 'LineWidth', 1.5);
plot(1:time_steps, PCommerical_HighImpact / 1000, '-b', 'LineWidth', 1.5);
plot(1:time_steps, PDomestic / 1000, '-m', 'LineWidth', 1.5);
plot(1:time_steps, PCommerical_LowImpact / 1000, '-c', 'LineWidth', 1.5);
plot(1:time_steps, PEducation / 1000, '-k', 'LineWidth', 1.5);
plot(1:time_steps, PPublic / 1000, '-y', 'LineWidth', 1.5);
xlabel('Time (hours)'); ylabel('Power (kW)');
title('Original Power Demand');
legend('Telecommunication', 'Healthcare', 'High Impact', 'Household', 'Low Impact', 'Education', 'Public');
grid on;

% Remove the second subplot (Effective Power Demand with Load Shedding)

% Subplot 2: PV Generation vs. Effective Demand (kWh)
subplot(6, 1, 2);  % Changed from 7,1,3 to 6,1,2
plot(1:time_steps, EPV / 1000, '-g', 'LineWidth', 1.5);
hold on;
plot(1:time_steps, E_demand_effective_array / 1000, '-r', 'LineWidth', 1.5);
xlabel('Time (hours)'); ylabel('Energy (kWh)');
title('PV Generation vs. Effective Demand');
legend('PV Generation', 'Effective Demand');
grid on;

% Adjust the remaining subplots accordingly
subplot(6, 1, 3);  % Storage Energy Flow
plot(1:time_steps, E_storage_array / 1000, '-b', 'LineWidth', 1.5);
xlabel('Time (hours)'); ylabel('Energy (kWh)');
title('Storage Energy Flow');
grid on;

subplot(6, 1, 4);  % Battery SoC
plot(1:time_steps, SoC_array, '-m', 'LineWidth', 1.5);
xlabel('Time (hours)'); ylabel('SoC (%)');
title('Battery SoC');
grid on;

subplot(6, 1, 5);  % Electricity Price
plot(1:time_steps, electricity_price, '-r', 'LineWidth', 1.5);
xlabel('Time (hours)'); ylabel('Price ($/kWh)');
title('Electricity Price');
grid on;

subplot(6, 1, 6);  % Cumulative Unserved Energy Cost
cumulative_Cost = cumsum(Cost_of_Unserved_array);
plot(1:time_steps, cumulative_Cost, '-b', 'LineWidth', 1.5);
xlabel('Time (hours)'); ylabel('Cumulative Cost ($)');
title('Cumulative Unserved Energy Cost');
grid on;

%% Plot Individual Load Shedding per Profile
figure('Name', 'Load Shedding per Profile', 'Position', [100, 100, 1200, 800]);

subplot(4,2,1);
plot(1:time_steps, NSL_P1_array, '-r', 'LineWidth', 1.5);
xlabel('Time (hours)');
ylabel('Load Shed (Wh)');
title('Telecommunication (P1)');
grid on;

subplot(4,2,2);
plot(1:time_steps, NSL_P2_array, '-g', 'LineWidth', 1.5);
xlabel('Time (hours)');
ylabel('Load Shed (Wh)');
title('Healthcare (P2)');
grid on;

subplot(4,2,3);
plot(1:time_steps, NSL_P3_array, '-b', 'LineWidth', 1.5);
xlabel('Time (hours)');
ylabel('Load Shed (Wh)');
title('High Impact Commercial (P3)');
grid on;

subplot(4,2,4);
plot(1:time_steps, NSL_P4_array, '-m', 'LineWidth', 1.5);
xlabel('Time (hours)');
ylabel('Load Shed (Wh)');
title('Domestic (P4)');
grid on;

subplot(4,2,5);
plot(1:time_steps, NSL_P5_array, '-c', 'LineWidth', 1.5);
xlabel('Time (hours)');
ylabel('Load Shed (Wh)');
title('Low Impact Commercial (P5)');
grid on;

subplot(4,2,6);
plot(1:time_steps, NSL_P6_array, '-k', 'LineWidth', 1.5);
xlabel('Time (hours)');
ylabel('Load Shed (Wh)');
title('Education (P6)');
grid on;

subplot(4,2,7);
plot(1:time_steps, NSL_P7_array, '-y', 'LineWidth', 1.5);
xlabel('Time (hours)');
ylabel('Load Shed (Wh)');
title('Public (P7)');
grid on;

%% Combined Load Shedding Plot
figure('Name', 'Combined Load Shedding per Profile', 'Position', [100, 100, 1200, 600]);
plot(1:time_steps, NSL_P1_array, '-r', 'LineWidth', 1.5);
hold on;
plot(1:time_steps, NSL_P2_array, '-g', 'LineWidth', 1.5);
plot(1:time_steps, NSL_P3_array, '-b', 'LineWidth', 1.5);
plot(1:time_steps, NSL_P4_array, '-m', 'LineWidth', 1.5);
plot(1:time_steps, NSL_P5_array, '-c', 'LineWidth', 1.5);
plot(1:time_steps, NSL_P6_array, '-k', 'LineWidth', 1.5);
plot(1:time_steps, NSL_P7_array, '-y', 'LineWidth', 1.5);
xlabel('Time (hours)');
ylabel('Load Shed (Wh)');
title('Load Shedding for All Profiles');
legend('Telecommunication', 'Healthcare', 'High Impact', 'Domestic', 'Low Impact', 'Education', 'Public');
grid on;


%% Plot Trading Data (Updated)
figure('Name', 'Trading Analysis', 'Position', [100, 100, 1200, 1000]);

% 1. External Grid Energy Profiles
subplot(3, 2, 1);
hold on;
plot(ext1_energy, 'b'); 
plot(ext2_energy, 'r'); 
legend('Grid 1 Energy', 'Grid 2 Energy');
xlabel('Time (hours)'); ylabel('kWh');
title('External Grid Energy Profile'); 
grid on;

% 2. External Grid Energy Price
subplot(3, 2, 2);
hold on;
plot(ext1_price, 'b'); 
plot(ext2_price, 'r'); 
title('External Grid Energy Price');
legend('Grid 1 Price', 'Grid 2 Price');
xlabel('Time (hours)'); ylabel('$/kWh');
grid on;

% 3. Energy Sold to External Grids
subplot(3, 2, 3);
hold on;
plot(energy_sold_ext1, 'b', 'LineWidth', 1.5);
plot(energy_sold_ext2, 'r', 'LineWidth', 1.5);
legend('Energy Sold to Grid 1', 'Energy Sold to Grid 2');
xlabel('Time (hours)'); ylabel('Energy (kWh)');
title('Energy Sold to External Grids');
grid on;

% 4. Cumulative Trade Profit
subplot(3, 2, 4);
cumulative_trade_profit = cumsum(trading_profit);
plot(cumulative_trade_profit, 'k', 'LineWidth', 1.5);
xlabel('Time (hours)'); ylabel('Cumulative Profit ($)');
title('Cumulative Trade Profit');
grid on;

% 5. Power Import from External Grids (New)
subplot(3, 2, 5);
hold on;
plot(energy_bought_ext1, 'b', 'LineWidth', 1.5);
plot(energy_bought_ext2, 'r', 'LineWidth', 1.5);
legend('Energy Bought from Grid 1', 'Energy Bought from Grid 2');
xlabel('Time (hours)'); ylabel('Energy (kWh)');
title('Power Import from External Grids');
grid on;

% 6. Cumulative Cost of Power Import (New)
subplot(3, 2, 6);
cumulative_import_cost = cumsum(trading_cost);
plot(cumulative_import_cost, 'k', 'LineWidth', 1.5);
xlabel('Time (hours)'); ylabel('Cumulative Cost ($)');
title('Cumulative Cost of Power Import');
grid on;

%% Plot Energy Imbalance Over Time
figure('Name', 'Energy Imbalance Over Time', 'Position', [100, 100, 1200, 600]);
plot(1:time_steps, imbalance_array, '-r', 'LineWidth', 1.5);
xlabel('Time (hours)');
ylabel('Energy Imbalance (Wh)');
title('Energy Imbalance Over Time');
grid on;


%%
disp(['Total Unserved Energy Cost: $', num2str(total_unserved_energy_cost)]);
%%
function non_served = distribute_deficit(loads, VoLL, deficit_kWh)
    % Convert loads to kWh if they're in Wh
    loads = loads;
    
    % Initialize non-served loads
    non_served = zeros(size(loads));
    remaining_deficit = deficit_kWh;
    
    % Sort loads by VoLL (lowest priority first)
    [sorted_VoLL, idx] = sort(VoLL);
    sorted_loads = loads(idx);
    
    % Distribute deficit starting from lowest priority
    for i = 1:length(sorted_loads)
        if remaining_deficit <= 0
            break;
        end
        % Amount to shed from this load
        to_shed = min(sorted_loads(i), remaining_deficit);
        non_served(idx(i)) = to_shed;
        remaining_deficit = remaining_deficit - to_shed;
    end
end