%% Shigley's P8.14

%% initialize values


E_steel = 30; % M psi
h_steel = 2; % in
E_ci = 14.5; % M psi
h_ci = 1; % in
d = 0.5; % in
l = 3; % in
H = 7/16;

% one bolt
% one nut

% bolt is steel 1/2 13 UNC


%% Part 1
% determine suitable length for the bolt, rounded up to the nearest quarter inch

% reference Table A-31
% nut height = 7/16in
% L >= l + H = 3 7/16 in, need to round up
L = 3.5; % in


%% Part 2:
% determine bolt stiffness
% spring equivalent thingy

L_t = 2 * d + 1/4;
l_d = L - L_t;
l_t = l - l_d;
A_d = pi * d^2 / 4;
A_t = 0.1419; % Table 8-2
k_b = A_d * A_t * E_steel / (A_d * l_t + A_t * l_d) % Mlbf / in


%% Part 3
% determine stiffness of the members
% cumulative equivalent

% Top Steel Frustum:
t = 1.5;
D = 0.75;
k_1 = 0.5774 * pi * E_steel * d / log((1.155 * t + D - d)*(D + d) / ((1.155 * t + D + d) * (D - d)));

% Bottom Steel Frustum
t = 0.5;
D = 0.75 + 2 * tand(30);
k_2 = 0.5774 * pi * E_steel * d / log((1.155 * t + D - d)*(D + d) / ((1.155 * t + D + d) * (D - d)));

% Cast Iron Frustum
t = 1;
D = 0.75;
k_3 = 0.5774 * pi * E_ci * d / log((1.155 * t + D - d)*(D + d) / ((1.155 * t + D + d) * (D - d)));

k_m = ( 1 / k_1 + 1 / k_2 + 1 / k_3)^(-1)
