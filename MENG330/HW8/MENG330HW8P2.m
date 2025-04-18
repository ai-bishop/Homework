%% Shigley's 8.33


%% Initialize values
% non permanent bolt
E_steel = 207; % GPA
E_ci = 100; % GPA
% in mm:
L_A = 20;
L_B = 20;
L_C = 100;
L_D = 150;
L_E = 200;
L_F = 300;
diam = 12;
l = 40;

A_g = pi * (L_D / 2) ^ 2; % area of head of gasket

N = 10; % Number bolts

p_g = 6; % MPA

% Bolt grade: ISO 9.8
% Bolt Spec: M12 x 1.75

% Find:
% L
% n_p
% n_l
% n_0

%% Start work
% T A-31:
H = 10.8;
% L = l + H = 50.8
L = 60 % rounded up

L_t = 2 * diam + 6;
l_d = L - L_t;
l_t = l - l_d;
A_d = pi * diam^2 / 4;
A_t = 84.3; % T8-1

% bolt
k_b = (A_t * A_d * E_steel) / (A_d * l_t + A_t * l_d);

% members
% steel
t = 20;
d = 12;
D = 18;

k_1 = 0.5774 * pi * E_steel * d / log(((1.155 * t + D - d) * (D + d) / ((1.155 * t + D + d) * (D - d))));

% cast iron - only change is material
k_2 = 0.5774 * pi * E_ci * d / log(((1.155 * t + D - d) * (D + d) / ((1.155 * t + D + d) * (D - d))));

k_m = (1 / k_1 + 1 / k_2)^-1;

C = k_b / (k_b + k_m);

S_p = 650; % MPA, T8-11

% assume nonpermanent connections
F_i = 0.75 * A_t * S_p * 10^-3; % kN

P = p_g * A_g * 10^-3 / N; % kN/bolt


n_p = S_p * A_t * 10^-3 / (C * P + F_i)

n_l = (S_p * A_t  * 10^-3 - F_i) / (C * P)

n_0 = F_i / (P * (1 - C))




