clear

% Initialization of parameters
PV_peak = 600; % Peak PV generation in watts
P_max_inverter = 700; % Maximum inverter power in Watts
n_inverter = 0.9; % Efficiency of the inverter (90%)
BESS_size = 8500; % Battery Energy Storage System size in Wh
SoC_max = 100; % Maximum State of Charge (%)
SoC_min = 20; % Minimum State of Charge (%)
SoC = 50; % Initial State of Charge (%)
P_max_charger = 1000; % Maximum charger power in Watts
time_steps = 168; % Simulation time for 7 days (hours)
delta_t = 1; % Time interval in hours


% Time settings
hours_per_day = 24; % 24 hours in a day
days = 7; % Total days
total_hours = hours_per_day * days; % Total hours for simulation
%% 

% Daily demand profiles (24 hours each)
P1_daily = [20 20 20 20 30 30 40 40 50 50 50 60 40 40 30 30 50 60 70 80 80 50 40 30]; % Highest Priority (W)
P2_daily = [0 0 0 0 0 0 50 50 100 150 200 200 100 100 50 50 50 150 200 300 200 100 50 0]; % Medium Priority (W)
P3_daily = [0 0 0 0 0 0 0 0 0 200 300 400 300 200 0 0 0 200 400 500 400 300 100 0]; % Lowest Priority (W)

% Replicate daily demand for 7 days
P1_demand = repmat(P1_daily, 1, days);
P2_demand = repmat(P2_daily, 1, days);
P3_demand = repmat(P3_daily, 1, days);

% Combine into a matrix for visualization
demand_matrix_7days = [P1_demand; P2_demand; P3_demand];

% Time vector for 7 days
time_vector = 1:total_hours;

% Total energy demand
E_demand = (P1_demand + P2_demand + P3_demand) * delta_t; % Total demand per hour

% Inverter capacity limit
E_max_inverter = P_max_inverter * delta_t; % Maximum energy through inverter per hour
%E_demand(E_demand > E_max_inverter) = E_max_inverter; % Limit demand by inverter capacity

% %% Read and Process Irradiance Data
% filename = 'irradiance_data.xlsx'; % Update this with the actual filename
% sheet = 1; % Specify the sheet number if necessary
% irradiance_data = readmatrix(filename, 'Sheet', sheet);
% 
% % Ensure correct shape (convert column to row if needed)
% irradiance_data = irradiance_data(:)';
% 
% % Check data length
% num_hours = length(irradiance_data);
% if num_hours > time_steps
%     irradiance_data = irradiance_data(1:time_steps);
% elseif num_hours < time_steps
%     error('The irradiance data is shorter than the required simulation time.');
% end
% 
% % Remove negative values or NaNs
% irradiance_data(irradiance_data < 0) = 0;
% irradiance_data(isnan(irradiance_data)) = 0;
% 
% %% Compute PV Generation
% % Convert irradiance to PV power
% EPV = irradiance_data*PV_peak;
% 
% % Apply inverter efficiency
% EPV = EPV * n_inverter; 
% 
% % Limit by inverter maximum power
% EPV = min(EPV, P_max_inverter);


% Define PV generation profile
PV_peak = 600; % Peak PV generation in watts
day_start = 6; % Start of daytime (6 AM)
day_end = 18;  % End of daytime (6 PM)

EPV = zeros(1, time_steps); % Initialize PV generation array

for t = 1:time_steps
    hour_of_day = mod(t-1, 24); % Hour of the day (0-23)
    if hour_of_day >= day_start && hour_of_day <= day_end
        % Generate PV power as a sinusoidal function peaking at noon
        EPV(t) = PV_peak * sin(pi * (hour_of_day - day_start) / (day_end - day_start));
    end

    % Limit PV generation by the inverter capacity
    EPV(t) = min(EPV(t)*n_inverter, P_max_inverter*n_inverter);
end
%% 

% Storage energy flow and SoC update
E_storage = zeros(1, time_steps); % Energy flow from/to storage
SoC_array = zeros(1, time_steps); % State of Charge log
delta_E = zeros(1, time_steps); % Net energy flow (generation - demand)

for t = 1:time_steps
    % Calculate net energy (generation - demand)
    delta_E(t) = EPV(t) - E_demand(t);

    % Apply maximum charger power limit during charging
    if delta_E(t) > 0 % Charging case
        delta_E(t) = min(delta_E(t), P_max_charger); % Limit charging to max charger power
    end

    % Ensure charging does not exceed SoC_max
    DeltaAE_t = SoC_max - SoC;
    if delta_E(t) > DeltaAE_t
        delta_E(t) = DeltaAE_t; % Limit energy flow by available capacity
    end

    % Ensure discharging does not fall below SoC_min
    if delta_E(t) < 0 % Discharging case
        AE_t = SoC - SoC_min;
        if abs(delta_E(t)) > AE_t
            delta_E(t) = -AE_t; % Limit discharging energy
        end
    end

    % Update storage energy and SoC
    E_storage(t) = delta_E(t);
    SoC = SoC + (delta_E(t) / BESS_size) * 100; % Update SoC (%)
    SoC_array(t) = SoC; % Log SoC
end

%% 

% Plot results
figure;

% Subplot 1: Energy demand for P1, P2, and P3
subplot(4, 1, 1);
plot(time_vector, P1_demand, '-r', 'LineWidth', 1.5);
hold on;
plot(time_vector, P2_demand, '-g', 'LineWidth', 1.5);
plot(time_vector, P3_demand, '-b', 'LineWidth', 1.5);
xlabel('Time (hours)');
ylabel('Power (W)');
title('Power Demand for P1, P2, and P3');
legend('P1 Demand', 'P2 Demand', 'P3 Demand');
grid on;

% Subplot 2: Total energy demand vs PV generation
subplot(4, 1, 2);
plot(time_vector, E_demand, '-r', 'LineWidth', 1.5);
hold on;
plot(time_vector, EPV, '-g', 'LineWidth', 1.5);
xlabel('Time (hours)');
ylabel('Energy (Wh)');
title('Energy Demand vs PV Generation (After Inverter) (7 Days)');
legend('Total Energy Demand', 'PV Generation (After Inverter)');
grid on;

% Subplot 3: Energy flow to/from storage
subplot(4, 1, 3);
plot(time_vector, E_storage, '-b', 'LineWidth', 1.5);
xlabel('Time (hours)');
ylabel('Energy (Wh)');
title('Energy Flow to/from Storage (7 Days)');
legend('Storage Energy Flow');
grid on;

% Subplot 4: Battery State of Charge (SoC)
subplot(4, 1, 4);
plot(time_vector, SoC_array, '-m', 'LineWidth', 1.5);
xlabel('Time (hours)');
ylabel('State of Charge (%)');
title('Battery State of Charge (SoC) (7 Days)');
legend('State of Charge');
grid on;
