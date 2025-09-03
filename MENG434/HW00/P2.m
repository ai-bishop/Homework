% This file is used for the graphing of HW00P2

y = @(t) (2 / 401) * (0.5 * cos(10 * t) + 10 * sin(10 * t)) + (400/401) * exp(-t/2);


% plotting
figure 
hold on
grid on

% plot function
t = 0:0.01:10;
plot(t, y(t))
xlim([0,10])