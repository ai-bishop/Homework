%% Initialize Parameters
clear; clc;

% Thermocouple Properites
d = 0.0007; % m, of the thermocouple
r = d/2; % m, of the thermocouple
T_0 = 20; % *C, initial temp of thermocouple
T_inf = 50; % *C, temperature of surrounding air
A = 4 * pi * r^2; % m2, area of the sphere

% Air Properties
rho = 8500; % kg/m3, density of the air
C = 400; % J/kg-K, specific heat
k = 400; % w/m-k, thermal conductivity
epsilon = 0.85; % surface emissivity
h = 10; % W/m2-k, heat transfer coefficient
sigma = 5.67 * 10^(-8); % w/m2-k4, stefan-boltzmann constant
V = (4/3) * pi * r^3; % m3, volume of thermocouple

%% Assumptions
% Neglect heat conduction throught the shaft of the thermocouple junction
% Biot number < 0.1 (proven in class), therefore can used Lumped Capacitance method.

%% Solve Problem
% How long does it take the thermocouple junction to have a temperature
% reading within 0.5*C accuracy?
T_target = 49.5; % *C, target temperature

% ODE: dT/dt = (Q_in) / (rho * V * C)
% Q_in = h*A*(T_inf - T) + epsilon*sigma*A*((T_inf+273)^4 - (T+273)^4)

% define ODE function
dTdt = @(t, T) (h*A*(T_inf - T) + epsilon*sigma*A*((T_inf+273)^4 - (T+273)^4)) / (rho * V * C);

% time span (s), how long to calc for
tspan = [0 500];

% solve ODE
[t, T] = ode45(dTdt, tspan, T_0);

% find when T = T_target
idx = find(T >= T_target, 1);
time_to_accuracy = t(idx);

fprintf('The time for the thermocouple to reach the target temperature of 49.5°C is %.2f seconds \n', time_to_accuracy)

%% Graph Results

figure
plot(t, T, 'b', 'LineWidth', 2);
hold on
yline(50)
yline(49.5, 'r--', 'LineWidth', 1.5, 'Label', '49.5°C Threshold') % T_target
plot(time_to_accuracy, T(idx), 'ko', 'MarkerFaceColor', 'g')
text(time_to_accuracy+5, T(idx), sprintf('Time to 49.5°C: %.2f sec', time_to_accuracy), 'FontSize', 12, 'HorizontalAlignment', 'left', 'VerticalAlignment', 'top');
yline(20) % T_0
ylim([15,55])
xlabel('Time (s)')
ylabel('Temperature (C)')
title('Thermocouple Junction Temperature vs Time');
grid on
legend('Temperature', 'Threshold (49.5°C)', 'Target Reached', 'Location', 'east');
