%% Problem 3


% Shigley’s 14.1 A steel spur pinion has a pitch of 6 teeth/in, 22 full-depth teeth,
% and a 20° pressure angle. The pinion runs at a speed of 1200 rev/min and transmits 15
% hp to a 60-tooth gear. If the face width is 2 in, estimate the bending stress.

N = 22;
P = 6;
d = N / P;
Y = 0.331;
n = 1200;
V = pi * d * n / 12; % ft/min
K_v = (1200 + V) / 1200;
H = 15;
W_t = 33000 * H / V; % lbf
F = 2;
sigma = K_v * W_t * P / (F * Y); % answer
