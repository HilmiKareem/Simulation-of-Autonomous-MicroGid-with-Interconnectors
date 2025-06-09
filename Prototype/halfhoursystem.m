clear all

%% Define Value of Lost Load (VoLL) Parameters (in unitless)
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
P_max_inverter    = 400;         % Maximum inverter power (kW)
PV_peak           = Rs * P_max_inverter; % Peak PV generation (kW)
n_inverter        = 0.9;           % Inverter efficiency (AC output)
n_charger         = 0.9;           % Charger efficiency (DC charging)
BESS_size         = 1600;         % Battery capacity (kWh)
SoC_max           = 100;           % Maximum State of Charge (%)
SoC_min           = 0;             % Minimum State of Charge (%)
SoC               = 50;            % Initial State of Charge (%)
P_max_charger     = 500;         % Maximum charger power (kW)
min_Trading_SoC   = 30;            % Minimum SoC for trading (%)
DOD               = 5;             % Depth of Discharge (%)
%% Time Settings
hours_per_day = 24;
days          = 7;
total_hours   = hours_per_day * days;
time_steps    = total_hours*2;   % Half-hourly steps (336 steps)
delta_t       = 0.5;               % Time step (hours)

%% External Grid Parameters
% External Grid 1 (Seller)
ext1_trans_cap = 300;              % kW
% External Grid 2 (Buyer)
ext2_trans_cap = 300;              % kW

block = time_steps/4;  % = 84

ext1_power = zeros(1,time_steps);
ext2_power = zeros(1,time_steps);

% block 1
ext1_power(1:block)            =  500;
ext2_power(1:block)            = -500;

% block 2
ext1_power(block+1:2*block)    =    0;
ext2_power(block+1:2*block)    = -500;

% block 3
ext1_power(2*block+1:3*block)  =  500;
ext2_power(2*block+1:3*block)  =    0;

% block 4
ext1_power(3*block+1:end)      =    0;
ext2_power(3*block+1:end)      =    0;

% prices stay flat
ext1_price = 9e9* ones(size(ext1_power));
ext2_price = 0.8 * ones(size(ext2_power));

%% Initialize Variables
power_sold_ext1 = zeros(1, time_steps);
power_sold_ext2 = zeros(1, time_steps);
power_bought_ext1 = zeros(1, time_steps);
power_bought_ext2 = zeros(1, time_steps);
trading_cost = zeros(1, time_steps);
trading_profit = zeros(1, time_steps);
Cost_of_Unserved_array = zeros(1, time_steps);
SoC_array = zeros(1, time_steps);
SoC_array(1) = SoC;
P_storage_array = zeros(1, time_steps);


%% Load Profiles from Excel (Replace with actual data)
filename = 'C:\Users\ASUS\OneDrive - Imperial College London\Documents\MATLAB Drive\FYP\Clover Load 7days.xlsx'; 

% Read and reshape data from Excel sheets
Domestic   = reshape(readmatrix(filename, 'Sheet', 'domestic'), [], 1); %W
PDomestic = Domestic./1000; %kW
Commerical = reshape(readmatrix(filename, 'Sheet', 'commercial'), [], 1); %W
PCommerical = Commerical./1000; %kW
Public     = reshape(readmatrix(filename, 'Sheet', 'public'), [], 1); %W
PPublic = Public./1000; %kW

% Split Commercial Load into High and Low Impact
PCommerical_HighImpact = 0.2 .* PCommerical;  % 20% %kW
PCommerical_LowImpact  = 0.8 .* PCommerical;  % 80% %kW

%% Telecommunication Load Calculation
Total_Area            = 56000;        % m²
Range_per_tower       = 40000;        % m (coverage radius)
Consumption_per_tower = 5;          % kW per tower
Overlap_factor        = 1.2;          % accounts for overlap

Coverage_per_tower    = pi * Range_per_tower^2;                % m² per tower
Effective_Coverage    = Coverage_per_tower / Overlap_factor;   % m² per tower
Number_Tower_Needed   = ceil(Total_Area / Effective_Coverage);
Total_Consumption     = (Number_Tower_Needed + 2) * Consumption_per_tower;

% Create constant telecommunication load for 168 hours
PTelco = Total_Consumption * ones(time_steps/2, 1); %kW

%% Healthcare Load Profile
daily_profile        = [250, 250, 250, 250, 250, 700, 750, 810, 910, 1200, 1200, 1600, ...
                        1700, 1600, 950, 900, 700, 700, 500, 500, 400, 400, 300, 300];  % in W
daily_energy_target  = 11.5 * 1000;  % Wh (11.5 kWh)
day_variability      = 0.10;          % 10%
hour_variability     = 0.15;          % 15%
scaling_factor       = daily_energy_target / sum(daily_profile);
scaled_daily_profile = daily_profile * scaling_factor;

Healthcare = zeros(1, days * 24);
for day = 1:days
    day_factor   = 1 + (rand * 2 - 1) * day_variability;
    daily_load   = scaled_daily_profile * day_factor;
    hourly_loads = daily_load .* (1 + (rand(1, 24) * 2 - 1) * hour_variability);
    Healthcare((day-1)*24 + 1 : day*24) = hourly_loads; % W
end

PHealthcare = Healthcare./1000; %kW

%% Education Facility Load Profile
total_daily_energy_target = 13.56 ;  % kWh
peak_load                 = 3.050;           % kW (during school hours)
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

%% Create Time Vectors
t_hourly = 1:168;                  % 168 hourly time points
% To obtain 336 half-hourly data points, extend the range to 168.5
t_halfHourly = 1:0.5:168.5;          % This yields 336 points

%% Perform Linear Interpolation for each load profile

% Define parameters for each profile
sin_amplitude_PTelco = 0;      
sin_frequency_PTelco = 0;     

sin_frequency_PHealthcare = 0;     % 3 cycles/hour
sin_amplitude_PHealthcare = 0;      % kW

sin_frequency_PCommerical_HighImpact = 0;     % 3 cycles/hour
sin_amplitude_PCommerical_HighImpact = 0;      % kW

sin_frequency_PDomestic = 0;     % 3 cycles/hour
sin_amplitude_PDomestic = 0;      % kW

sin_frequency_PCommerical_LowImpact = 0;     % 3 cycles/hour
sin_amplitude_PCommerical_LowImpact = 0;      % kW

sin_frequency_PEducation = 0;     % 3 cycles/hour
sin_amplitude_PEducation = 0;      % kW

sin_frequency_PPublic = 0;     
sin_amplitude_PPublic = 0;      

% Generate sinusoidal components
sin_PTelco = sin_amplitude_PTelco * sin(2*pi*sin_frequency_PTelco*t_halfHourly);
sin_PHealthcare = sin_amplitude_PHealthcare * sin(2*pi*sin_frequency_PHealthcare*t_halfHourly);
sin_PCommerical_HighImpact = sin_amplitude_PCommerical_HighImpact * sin(2*pi*sin_frequency_PCommerical_HighImpact*t_halfHourly);
sin_PDomestic = sin_amplitude_PDomestic * sin(2*pi*sin_frequency_PDomestic*t_halfHourly);
sin_PCommerical_LowImpact = sin_amplitude_PCommerical_LowImpact * sin(2*pi*sin_frequency_PCommerical_LowImpact*t_halfHourly);
sin_PEducation = sin_amplitude_PEducation * sin(2*pi*sin_frequency_PEducation*t_halfHourly);
sin_PPublic = sin_amplitude_PPublic * sin(2*pi*sin_frequency_PPublic*t_halfHourly);

% Apply to profiles
PTelco_30mins = interp1(t_hourly, PTelco, t_halfHourly, 'linear') + sin_PTelco;
PTelco_30mins(isnan(PTelco_30mins)) = 0;

PHealthcare_30mins = interp1(t_hourly, PHealthcare, t_halfHourly, 'linear') + sin_PHealthcare;
PHealthcare_30mins(isnan(PHealthcare_30mins)) = 0;

PCommerical_HighImpact_30mins = interp1(t_hourly, PCommerical_HighImpact, t_halfHourly, 'linear') + sin_PCommerical_HighImpact;
PCommerical_HighImpact_30mins(isnan(PCommerical_HighImpact_30mins)) = 0;

PDomestic_30mins            = interp1(t_hourly, PDomestic, t_halfHourly, 'linear') + sin_PDomestic;
PDomestic_30mins(isnan(PDomestic_30mins)) = 0;

PCommerical_LowImpact_30mins = interp1(t_hourly, PCommerical_LowImpact, t_halfHourly, 'linear') + sin_PCommerical_LowImpact;
PCommerical_LowImpact_30mins(isnan(PCommerical_LowImpact_30mins)) = 0;

PEducation_30mins           = interp1(t_hourly, PEducation, t_halfHourly, 'linear') + sin_PEducation;
PEducation_30mins(isnan(PEducation_30mins)) = 0;

PPublic_30mins              = interp1(t_hourly, PPublic, t_halfHourly, 'linear') + sin_PPublic;
PPublic_30mins(isnan(PPublic_30mins)) = 0;

P1_effective = PTelco_30mins;                   % Telecommunication (Highest Priority)
P2_effective = PHealthcare_30mins;              % Healthcare
P3_effective = PCommerical_HighImpact_30mins;     % High Impact Businesses
P4_effective = PDomestic_30mins;                % Household
P5_effective = PCommerical_LowImpact_30mins;      % Low Impact Businesses
P6_effective = PEducation_30mins;               % Education
P7_effective = PPublic_30mins;                  % Public (Lowest Priority)

total_P = P1_effective + P2_effective + P3_effective + P4_effective + P5_effective + P6_effective + P7_effective;

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
max_deviation_kWh = max(deviations)* delta_t;


%% Calculate the Minimum SoC Reserve Required
% BESS_size is in Wh. The reserve required (in percentage) is the maximum deviation 
% divided by the total battery capacity, times 100.
slack_SoC_required = (max_deviation_kWh / BESS_size) * 100;
min_allowed_SoC = DOD + slack_SoC_required; % both in percentage

fprintf('Maximum deviation: %.2f kWh\n', max_deviation_kWh);
fprintf('Minimum SoC reserve required: %.2f%% of BESS capacity\n', min_allowed_SoC);

% Initialize arrays non-served loads (in kW)

NSL_P1_array = zeros(1, time_steps);
NSL_P2_array = zeros(1, time_steps);
NSL_P3_array = zeros(1, time_steps);
NSL_P4_array = zeros(1, time_steps);
NSL_P5_array = zeros(1, time_steps);
NSL_P6_array = zeros(1, time_steps);
NSL_P7_array = zeros(1, time_steps);

P_imbalance_array = zeros(1, time_steps);

%% Read and Process Irradiance Data (Replace with actual data)
filename = 'irradiance_data.xlsx'; 
sheet = 1; 
irradiance_data = readmatrix(filename, 'Sheet', sheet);
irradiance_data = irradiance_data(:)';  % Ensure row vector

%% Read and Process Irradiance Data (Interpolate to Half-Hourly)
irradiance_data_halfhourly = interp1(1:168, irradiance_data, linspace(1, 168, 336), 'linear');
irradiance_data_halfhourly(irradiance_data_halfhourly < 0) = 0;
irradiance_data_halfhourly(isnan(irradiance_data_halfhourly)) = 0;  % Handle NaN

num_hours = length(irradiance_data_halfhourly);
if num_hours > time_steps
    irradiance_data_halfhourly = irradiance_data_halfhourly(1:time_steps);
elseif num_hours < time_steps
    error('The irradiance data is shorter than the required simulation time.');
end
irradiance_data(irradiance_data < 0) = 0;
irradiance_data(isnan(irradiance_data)) = 0;

%% Compute PV Power Generation (Half-Hourly)
Power_PV = irradiance_data_halfhourly * PV_peak ; % Adjust for half-hourly where irradiance is normalized clearness index (0–1)
Power_PV = Power_PV * n_inverter; 
Power_PV = min(Power_PV, P_max_inverter ); % (kW)

%% Cost Calculation
Solar_Panel_Capital = (PV_peak) * PV_Module_Price;  % PV capital cost $
Battery_Capital     = Battery_Cost * (BESS_size);   % BESS capital cost $
Total_Capital       = Solar_Panel_Capital + Battery_Capital; %$

% Financial & Operational Parameters
i                 = 0.03;  % Interest rate (3%)
OnM_Cost = 0.015 * Total_Capital;  % 1.5% of total capital cost
Variable_OnM_Cost = 0;   % Variable O&M ($/kWh)
Operational_Years = 20;    % PV lifetime
Capacity_Factor   = 0.17;  % PV capacity factor

% Capital Recovery Factor (CRF)
CRF = (i * (1 + i)^Operational_Years) / ((1 + i)^Operational_Years - 1);

% Annual Energy Output (kWh/year)
E_PV = 8760 * Capacity_Factor * (PV_peak);  % PV-generated energy
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
    actual_total_demand = (PTelco_30mins(t) + PHealthcare_30mins(t) + PCommerical_HighImpact_30mins(t) + ...
                          PDomestic_30mins(t) + PCommerical_LowImpact_30mins(t) + PEducation_30mins(t) + PPublic_30mins(t)); %kW
    P_imbalance = actual_total_demand - Power_PV(t); %kW
    E_imbalance = P_imbalance*delta_t; %kWh

    % Compute electricity price
    storage_factor = 1 - (current_SoC / SoC_max);
    price_adjustment = b * E_imbalance;
    electricity_price(t) = Price + c * storage_factor + price_adjustment;
    
    P_imbalance_array(t) = P_imbalance; %kW

    if P_imbalance <= 0
        %% Excess Power: Charge Battery or Sell
        excess_power = abs(P_imbalance); %kW
        
        % Charge battery (DC Power)energy_sold_ext1
        Charge_Power = min([excess_power * n_charger, ...        % Apply charger efficiency
                      P_max_charger]); %kW
        
        SoC_array(t) = min(100, current_SoC + (Charge_Power*delta_t / BESS_size) * 100);

        P_storage_array(t) = Charge_Power;  % Track DC charging %kW

        % Sell excess energy if SoC > min_Trading_SoC
        if SoC_array(t) > min_Trading_SoC
            available_energy = (SoC_array(t) - min_allowed_SoC)/100 * BESS_size;  % Respect DOD
            
            % Prioritize higher-priced grid
            if ext1_price(t) > ext2_price(t)
                % Sell to Grid 1
                if ext1_power(t) < 0  % Grid 1 is a buyer
                    max_sell = min([ext1_trans_cap, ...
                                    available_energy/delta_t * n_inverter, ...  
                                   abs(ext1_power(t))]);
                    power_taken = max_sell;  
                    if max_sell > 0
                        SoC_array(t) = SoC_array(t) - (power_taken*delta_t / BESS_size) * 100;
                        trading_profit(t) = trading_profit(t) + ext1_price(t) * max_sell*delta_t;
                        power_sold_ext1(t) = max_sell;
                        P_storage_array(t) = P_storage_array(t) - power_taken;
                    end
                end
            else
                % Sell to Grid 2
                if ext2_power(t) < 0  % Grid 2 is a buyer
                    max_sell = min([ext2_trans_cap, ...
                                    available_energy/delta_t * n_inverter, ...  
                                   abs(ext2_power(t))]);
                    power_taken = max_sell;  
                    if max_sell > 0
                        SoC_array(t) = SoC_array(t) - (power_taken*delta_t / BESS_size) * 100;
                        trading_profit(t) = trading_profit(t) + ext2_price(t) * max_sell*delta_t;
                        power_sold_ext2(t) = max_sell;
                        P_storage_array(t) = P_storage_array(t) - power_taken;
                    end
                end
            end
        end

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

        % Load demands at timestep t (kW)
        load_demands = [PTelco_30mins(t); PHealthcare_30mins(t); PCommerical_HighImpact_30mins(t);
                        PDomestic_30mins(t); PCommerical_LowImpact_30mins(t); PEducation_30mins(t); 
                        PPublic_30mins(t)] ; %kW 

        % --- Check Seller Grids ---
        seller_grids = struct('price',{},'trans_cap',{},'power',{},'grid_num',{});
        if ext1_power(t) > 0  % Grid 1 is a seller
            seller_grids = [seller_grids; struct('price', ext1_price(t), ...
                                               'trans_cap', ext1_trans_cap, ...
                                               'power', ext1_power(t), ...
                                               'grid_num', 1)];
        end
        if ext2_power(t) > 0  % Grid 2 is a seller
            seller_grids = [seller_grids; struct('price', ext2_price(t), ...
                                               'trans_cap', ext2_trans_cap, ...
                                               'power', ext2_power(t), ...
                                               'grid_num', 2)];
        end

        % Number of seller grids
        num_seller_grids = length(seller_grids);

        % Decision variables: [battery_usage_DC, grid_import_1, ..., grid_import_N, non_served_loads_1, ..., non_served_loads_7]
        num_loads = length(VoLL);
        num_vars = 1 + num_seller_grids + num_loads; 

        % --- Objective Function ---
        f = [electricity_price(t), ...             % Battery usage
             [seller_grids.price], ...                 % Grid import costs (row vector)
             (VoLL * electricity_price(t))'];          % VoLL for each load (scalar multiplication)

        % --- Inequality Constraints (Upper Bounds) ---
        % 1. Battery usage <= min(available energy, max discharge rate)
        available_battery_energy = max(0, ((current_SoC - min_allowed_SoC)/100) * BESS_size* delta_t); % in kWh 
        max_discharge_power = (P_max_inverter);                 % kW 
        max_battery_usage_DC = max_discharge_power * n_inverter;                 % kW 
        battery_upper = min(available_battery_energy/delta_t, max_battery_usage_DC); % kW

        % 2. Grid imports <= transmission capacity or power available
        if isempty(seller_grids)
            grid_upper = [];
        else
            caps  = [seller_grids.trans_cap];   % row vector of each grid’s trans_cap
            pwrs  = [seller_grids.power];       % row vector of each grid’s available power
            grid_upper = min(caps, pwrs);       % element‐wise minimum
        end

        % Modify the upper bounds for non-served loads to ensure they can't exceed demand
        non_served_upper = min(load_demands, load_demands);  % Each can't exceed its own demand

        % Combine upper bounds
        lb = zeros(1, num_vars);       % All variables >= 0
        ub = [battery_upper, grid_upper, non_served_upper'];

        % --- Equality Constraint (Power Balance) ---
        % n_inverter * battery_usage + sum(grid_imports) + sum(non_served) = imbalance
        Aeq = [n_inverter, ones(1, num_seller_grids), ones(1, num_loads)];
        beq = P_imbalance;     %  kW

        % First, calculate available battery energy before optimization
        available_battery_energy = max(0, ((current_SoC - min_allowed_SoC)/100) * BESS_size* delta_t); %kWh
        
        if available_battery_energy <= 0 && isempty(seller_grids)
            % If no battery energy available and no seller grids, must shed load
            x = [];
            exitflag = -1;  % Force load shedding path
        else
            % Run optimization
            options = optimset('Display', 'off');
            [x, fval, exitflag] = linprog(f, [], [], Aeq, beq, lb, ub, options);
        end
        
        % Handle optimization results
        if ~isempty(x) && exitflag == 1
            % Extract optimization results
            battery_usage_DC = x(1);                     % kW 
            grid_imports = x(2:1+num_seller_grids);        % kW 
            non_served_loads = x(2+num_seller_grids:end);  % kW, one per load profile
       
            grid_imports_kW = sum(grid_imports);       
        
            % Retrieve the actual individual demands (in kW)
            actual_P1_Demand = P1_effective(t);  
            actual_P2_Demand = P2_effective(t);
            actual_P3_Demand = P3_effective(t);
            actual_P4_Demand = P4_effective(t);
            actual_P5_Demand = P5_effective(t);
            actual_P6_Demand = P6_effective(t);
            actual_P7_Demand = P7_effective(t);
        
            % Compute the served load per profile by subtracting non-served
            % load. (in kW)
            served_load_P1 = actual_P1_Demand  - non_served_loads(1) ;
            served_load_P2 = actual_P2_Demand  - non_served_loads(2) ;
            served_load_P3 = actual_P3_Demand  - non_served_loads(3) ;
            served_load_P4 = actual_P4_Demand  - non_served_loads(4) ;
            served_load_P5 = actual_P5_Demand  - non_served_loads(5) ;
            served_load_P6 = actual_P6_Demand  - non_served_loads(6) ;
            served_load_P7 = actual_P7_Demand  - non_served_loads(7) ;
            
            % The effective served demand is the sum of the served loads (in kW)
            served_demand = served_load_P1 + served_load_P2 + served_load_P3 + ...
                            served_load_P4 + served_load_P5 + served_load_P6 + served_load_P7;
                        
            % Calculate total available power (kW) from PV, battery (after inverter losses), and grid imports
            total_available_power = Power_PV(t) + battery_usage_DC * n_inverter + grid_imports_kW;
            
            % Verify power balance based on served (effective) demand
            power_deficit = served_demand - total_available_power;
            if power_deficit > 1e-3  % Tolerance of 1 W
                fprintf('power balance violated at time %d: Deficit = %.2f kW\n', t, power_deficit);
                % Optionally, adjust load shedding further if desired.
            end
        
            % Update Battery State-of-Charge and storage flow
            SoC_array(t) = current_SoC - (battery_usage_DC*delta_t / BESS_size) * 100;
            P_storage_array(t) = -battery_usage_DC;  % Negative sign for discharge
        

            % Store individual served loads and non-served loads (convert non-served to kW)
            NSL_P1_array(t) = non_served_loads(1) ;
            NSL_P2_array(t) = non_served_loads(2) ;
            NSL_P3_array(t) = non_served_loads(3) ;
            NSL_P4_array(t) = non_served_loads(4) ;
            NSL_P5_array(t) = non_served_loads(5) ;
            NSL_P6_array(t) = non_served_loads(6) ;
            NSL_P7_array(t) = non_served_loads(7) ;
        
            % Update trading cost and grid purchase information (remains unchanged)
            for g = 1:num_seller_grids
                if seller_grids(g).grid_num == 1
                    power_bought_ext1(t) = grid_imports(g) ;
                    trading_cost(t) = trading_cost(t) + grid_imports(g)*delta_t * seller_grids(g).price;
                elseif seller_grids(g).grid_num == 2
                    power_bought_ext2(t) = grid_imports(g) ;
                    trading_cost(t) = trading_cost(t) + grid_imports(g)*delta_t * seller_grids(g).price;
                end
            end
        
            % Update unserved energy cost using the individual non-served loads
            Cost_of_Unserved_array(t) = sum(non_served_loads*delta_t .* VoLL .* electricity_price(t));
            
        else
            % Fallback: No feasible solution, implement load shedding directly.
            % In this branch, only PV is available.
            SoC_array(t) = current_SoC;
            P_storage_array(t) = 0;
            total_available_power = Power_PV(t);
            power_deficit = P_imbalance;
            non_served_loads = distribute_deficit(load_demands, VoLL, power_deficit);
            
            % Compute served load individually for each profile
            served_load_P1 = P1_effective(t)  - non_served_loads(1) ;
            served_load_P2 = P2_effective(t)  - non_served_loads(2) ;
            served_load_P3 = P3_effective(t)  - non_served_loads(3) ;
            served_load_P4 = P4_effective(t)  - non_served_loads(4) ;
            served_load_P5 = P5_effective(t)  - non_served_loads(5) ;
            served_load_P6 = P6_effective(t)  - non_served_loads(6) ;
            served_load_P7 = P7_effective(t)  - non_served_loads(7) ;
            
            served_demand = served_load_P1 + served_load_P2 + served_load_P3 + ...
                            served_load_P4 + served_load_P5 + served_load_P6 + served_load_P7;
            
            NSL_P1_array(t) = non_served_loads(1) ;
            NSL_P2_array(t) = non_served_loads(2) ;
            NSL_P3_array(t) = non_served_loads(3) ;
            NSL_P4_array(t) = non_served_loads(4) ;
            NSL_P5_array(t) = non_served_loads(5) ;
            NSL_P6_array(t) = non_served_loads(6) ;
            NSL_P7_array(t) = non_served_loads(7) ;
            
            Cost_of_Unserved_array(t) = sum(non_served_loads .* VoLL .* electricity_price(t));
        end
    end     
end

%% Calculate Total Unserved Energy Cost
total_unserved_energy_cost = sum(Cost_of_Unserved_array);

%% Plot Results
% Then modify your plotting section:
figure('Name', 'System Performance', 'Position', [100, 100, 1200, 1000]);

subplot(6, 1, 1);  % Using 6 rows in the subplot
plot(1:time_steps, PTelco_30mins , '-r', 'LineWidth', 1.5);
hold on;
plot(1:time_steps, PHealthcare_30mins , '-g', 'LineWidth', 1.5);
plot(1:time_steps, PCommerical_HighImpact_30mins , '-b', 'LineWidth', 1.5);
plot(1:time_steps, PDomestic_30mins , '-m', 'LineWidth', 1.5);
plot(1:time_steps, PCommerical_LowImpact_30mins , '-c', 'LineWidth', 1.5);
plot(1:time_steps, PEducation_30mins , '-k', 'LineWidth', 1.5);
plot(1:time_steps, PPublic_30mins , '-y', 'LineWidth', 1.5);
xlabel('Time (hours)');
ylabel('Power (kW)');
title('Original Power Demand (30-min data)');
legend('Telecommunication', 'Healthcare', 'High Impact', 'Household', 'Low Impact', 'Education', 'Public');
grid on;

% Subplot 2: PV Generation vs. Effective Demand (kW)
subplot(6, 1, 2);  % Changed from 7,1,3 to 6,1,2
plot(1:time_steps, Power_PV , '-g', 'LineWidth', 1.5);
hold on;
plot(1:time_steps, total_P , '-r', 'LineWidth', 1.5);
xlabel('Time (hours)'); ylabel('Power (kW)');
title('PV Generation vs. Total Demand');
legend('PV Generation', 'Total Demand');
grid on;

% Adjust the remaining subplots accordingly
subplot(6, 1, 3);  % Storage Energy Flow
plot(1:time_steps, P_storage_array , '-b', 'LineWidth', 1.5);
xlabel('Time (hours)'); ylabel('Power (kW)');
title('Storage Power Flow');
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
ylabel('Load Shed (kW)');
title('Telecommunication (P1)');
grid on;

subplot(4,2,2);
plot(1:time_steps, NSL_P2_array, '-g', 'LineWidth', 1.5);
xlabel('Time (hours)');
ylabel('Load Shed (kW)');
title('Healthcare (P2)');
grid on;

subplot(4,2,3);
plot(1:time_steps, NSL_P3_array, '-b', 'LineWidth', 1.5);
xlabel('Time (hours)');
ylabel('Load Shed (kW)');
title('High Impact Commercial (P3)');
grid on;

subplot(4,2,4);
plot(1:time_steps, NSL_P4_array, '-m', 'LineWidth', 1.5);
xlabel('Time (hours)');
ylabel('Load Shed (kW)');
title('Domestic (P4)');
grid on;

subplot(4,2,5);
plot(1:time_steps, NSL_P5_array, '-c', 'LineWidth', 1.5);
xlabel('Time (hours)');
ylabel('Load Shed (kW)');
title('Low Impact Commercial (P5)');
grid on;

subplot(4,2,6);
plot(1:time_steps, NSL_P6_array, '-k', 'LineWidth', 1.5);
xlabel('Time (hours)');
ylabel('Load Shed (kW)');
title('Education (P6)');
grid on;

subplot(4,2,7);
plot(1:time_steps, NSL_P7_array, '-y', 'LineWidth', 1.5);
xlabel('Time (hours)');
ylabel('Load Shed (kW)');
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

total_shedding = NSL_P1_array + NSL_P2_array + NSL_P3_array + ...
                 NSL_P4_array + NSL_P5_array + NSL_P6_array + NSL_P7_array;

plot(1:time_steps, total_shedding, '--', 'Color', [0.3, 0.3, 0.3], 'LineWidth', 2);  % dashed dark gray line

xlabel('Time (hours)');
ylabel('Load Shed (kW)');
title('Load Shedding for All Profiles');
legend('Telecommunication', 'Healthcare', 'High Impact', 'Domestic', 'Low Impact', 'Education', 'Public', 'Total Shedding');
grid on;


%% Plot Trading Data (Updated)
figure('Name', 'Trading Analysis', 'Position', [100, 100, 1200, 1000]);

% 1. External Grid Power Profiles
subplot(3, 2, 1);
hold on;
plot(ext1_power, 'b'); 
plot(ext2_power, 'r'); 
legend('Grid 1', 'Grid 2');
xlabel('Time (hours)'); ylabel('kW');
title('External Grid Power Profile'); 
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

% 3. Power Sold to External Grids
subplot(3, 2, 3);
hold on;
plot(power_sold_ext1, 'b', 'LineWidth', 1.5);
plot(power_sold_ext2, 'r', 'LineWidth', 1.5);
legend('Energy Sold to Grid 1', 'Energy Sold to Grid 2');
xlabel('Time (hours)'); ylabel('Power (kW)');
title('Power Sold to External Grids');
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
plot(power_bought_ext1, 'b', 'LineWidth', 1.5);
plot(power_bought_ext2, 'r', 'LineWidth', 1.5);
legend('Power Bought from Grid 1', 'Power Bought from Grid 2');
xlabel('Time (hours)'); ylabel('Power (kW)');
title('Power Import from External Grids');
grid on;

% 6. Cumulative Cost of Power Import (New)
subplot(3, 2, 6);
cumulative_import_cost = cumsum(trading_cost);
plot(cumulative_import_cost, 'k', 'LineWidth', 1.5);
xlabel('Time (hours)'); ylabel('Cumulative Cost ($)');
title('Cumulative Cost of Power Import');
grid on;

%% Plot Power Imbalance Over Time
figure('Name', 'Power Imbalance Over Time', 'Position', [100, 100, 1200, 600]);
plot(1:time_steps, P_imbalance_array, '-r', 'LineWidth', 1.5);
xlabel('Time (hours)');
ylabel('Power Deficit (kW)');
title('Power Imbalance Over Time');
grid on;


%%
disp(['Total Unserved Energy Cost: $', num2str(total_unserved_energy_cost)]);
%%
function non_served = distribute_deficit(loads, VoLL, deficit_kW)
    % Initialize non-served loads
    non_served = zeros(size(loads));
    remaining_deficit = min(deficit_kW, sum(loads));  % Limit deficit to total load
    
    % Sort loads by VoLL (lowest priority first)
    [sorted_VoLL, idx] = sort(VoLL);
    sorted_loads = loads(idx);
    
    % Distribute deficit starting from lowest priority
    for i = 1:length(sorted_loads)
        if remaining_deficit <= 0
            break;
        end
        % Amount to shed from this load - strictly between 0 and actual load
        to_shed = min(sorted_loads(i), remaining_deficit);
        to_shed = max(0, to_shed); % Ensure no negative load shedding
        non_served(idx(i)) = to_shed;
        remaining_deficit = remaining_deficit - to_shed;
    end
end