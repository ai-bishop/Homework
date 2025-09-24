% Thermocouple transient heating with convection and radiation
clear; clc;

% Constants
D = 0.0007;               % Diameter (m)
r = D/2;
rho = 8500;               % Density (kg/m^3)
c = 400;                  % Specific heat (J/kg-K)
k = 400;                  % Thermal conductivity (W/m-K)
epsilon = 0.85;           % Emissivity
h = 10;                   % Convective heat transfer coefficient (W/m^2-K)
sigma = 5.67e-8;          % Stefan-Boltzmann constant (W/m^2-K^4)

T_inf = 50;               % Ambient temperature (deg C)
T_init = 20;              % Initial temperature (deg C)
T_target = 49.5;          % Target temperature for ±0.5°C accuracy

% Geometry
V = (4/3)*pi*r^3;         % Volume
A = 4*pi*r^2;             % Surface area

% ODE: dT/dt = (Q_in) / (rho * V * c)
% Q_in = h*A*(T_inf - T) + epsilon*sigma*A*((T_inf+273)^4 - (T+273)^4)

% Define ODE function
dTdt = @(t, T) (h*A*(T_inf - T) + ...
    epsilon*sigma*A*((T_inf+273)^4 - (T+273)^4)) / (rho * V * c);

% Time span (in seconds)
tspan = [0 500];

% Solve ODE
[t, T] = ode45(dTdt, tspan, T_init);

% Find time when temperature reaches 49.5°C
idx = find(T >= T_target, 1);
time_to_accuracy = t(idx);

% Plot
figure;
plot(t, T, 'b', 'LineWidth', 2);
hold on;
yline(49.5, 'r--', 'LineWidth', 1.5, 'Label', '49.5°C Threshold');
plot(time_to_accuracy, T(idx), 'ko', 'MarkerFaceColor', 'g');
text(time_to_accuracy+5, T(idx), ...
    sprintf('Time to 49.5°C: %.2f sec', time_to_accuracy), ...
    'FontSize', 12);

xlabel('Time (s)');
ylabel('Temperature (°C)');
title('Thermocouple Junction Temperature vs Time');
grid on;
legend('Temperature', 'Threshold (49.5°C)', 'Target Reached');