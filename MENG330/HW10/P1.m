%% Problem 1

% A certain application requires a ball bearing with
% the inner ring rotating, with a design life of 25 kh at a speed of 350 rev/min. The radial load
% is 2.5 kN and an application factor of 1.2 is appropriate. The reliability goal is 0.90. Find the
% multiple of rating life required, xD, and the catalog rating C10 with which to enter a bearing
% table. Choose a 02-series deep-groove ball bearing from Table 11–2, and estimate the
% reliability in use.

% initiate vars

n_d = 350;
l_d2 = 25000;
l_d = 60 * n_d * l_d2;
l_r = 10^6;


%% Part 1

x_d = l_d / l_r;

%% Part 2

F_d = 2.5 * 1.2; % kN
C_10 = 3.0 * (525 / (0.02 + (4.459 - 0.02) * (ln(1/0.9))^(1/1.483))) ^ (1/3); % approx 24.3Kn


%% Part 3

% choose bearing w c10 >= C_10
% this results in 02-35mm w/ C10 = 25.5 kN

%% Part 4

R = exp(-((525 * (3/25.5)^3 - 0.02) / (4.459 - 0.02))^1.483);
