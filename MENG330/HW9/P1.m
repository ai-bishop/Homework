%% Problem 1
% Shigley’s 10.3 (Same in 10th Edition) A helical compression spring is
% wound using 2.5-mm-diameter music wire. The spring has an outside
% diameter of 31 mm with plain ground end, and 14 total coils.
% 1. Estimate the spring rate (spring constant)
% 2. What force is needed to compress this spring to closure if we assume that force
% to induce the maximum stress at the solid compressed state? (use a safety factor
% of 1 in the comparison rather than the recommended design factor of 1.2)
% 3. What should the free length be based on the result from part (2)?
% 4. Is there a possibility that the spring might buckle in service?






d_w = 2.5; % mm, diameter of wire
d_w_in = 0.0984; % inches
d_tot = 31; % mm, outside diameter of coil
D = d_tot - d_w;
n_tot = 14; % total number of coils
n_active = n_tot - 1;
C = D / d_w; % spring index
N = 14; % number of coils
G = 81000; % MPA aka N /mm2 - originally 81GPA
n_safety = 1;


%% Part 1
% Estimate the Spring Rate aka spring constant
% k ~= d^4 * G / 8 * D^3 * N
k_est = (d_w^4 * G) / (8 * D^3 * n_active);


%% Part 2
% 2. What force is needed to compress this spring to closure if we assume that force
% to induce the maximum stress at the solid compressed state? (use a safety factor
% of 1 in the comparison rather than the recommended design factor of 1.2)

% Table 10-4
m = 0.145;
A = 2211; % MPa * mm^m
S_ut = A / (d_w^m);
S_sy = 0.45 * S_ut;
F_s = pi * d_w^3 * S_sy / (8 * D);


%% Part 3
% 3. What should the free length be based on the result from part (2)?
L_s = n_tot * d_w;

L_0 = F_s / k_est + L_s;



%% Part 4
% 4. Is there a possibility that the spring might buckle in service?

alpha = 2.63 * D / L_0;
% alpha too large
% will fail
