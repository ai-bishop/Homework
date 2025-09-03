% This file is used for the graphing of HW00P7

y = @(t) (1/2) * sin(t) - (9/82) * sin(9*t) - (1/2) * cos(t) - (1/82) * cos(9*t) + (83/82) * exp(-t);


% plotting
figure 
hold on
grid on

% plot function
t = 0:0.01:12;
plot(t, y(t))
xlim([0,12])