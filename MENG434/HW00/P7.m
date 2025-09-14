% This file is used for of HW00P7

syms y

yprime = @(t,y) -y - cos(9 * y) + sin(t); % ODE solved for y'
y0 = 0.5; % initial condition


[t,y] = ode45(yprime, [0 12], y0); % solves for y(t), given as a series of points

% plotting
figure 
hold on
grid on

% plot function
%t = 0:0.01:12;
plot(t, y)
xlim([0,12])