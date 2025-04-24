%% Problem 2
% An extension spring with full-twisted loop ends is made from AISI 1065 OQ&T 
% wire. The spring has 84 coils and is close-wound with a preload of 16 lbf. The wire is 
% 0.162 in in diameter. The coil outside diameter is 1 ½ in. The end hooks radii are given as
% r1 is the radius of the mean diameter of the coil and r2=¿ ½ in.
% 1. Find the closed length of the spring
% 2. Find the torsional stress in the spring body corresponding to the preload
% 3. Estimate the spring rate (spring constant)
% 4. What load would cause permanent deformation? (check hooks and body)
% 5. What is the spring deflection corresponding to the load found in part (5)?
% Assume no deflection in the end hooks.

n_tot = 84;
f_i = 16; % lbf
d_wire = 0.162; % in
d_coil = 1.5; % in
D = d_coil - d_wire;
r1 = 0.669;
r2 = 0.25 + d_wire / 2;


%% Part 1
% 1. Find the closed length of the spring

L_s = 2 * (D - d_wire) + (n_tot + 1) * d_wire;


%% Part 2
% 2. Find the torsional stress in the spring body corresponding to the preload

C = D / d_wire;

K_ = (4 * C + 2) / (4 * C - 3);

tau = K_ * 8 * f_i * D / (pi * d_wire^3);


%% Part 3
% 3. Estimate the spring rate (spring constant)

% Table 10-5
G = 11.4 * 10^6; % psi
E = 28.5 * 10^6; % psi

n_a = n_tot + G / E;

k_est = d_wire^4 * G / (8 * D^3 * n_a); % lbf/in



%% Part 4
% 4. What load would cause permanent deformation? (check hooks and body)

% Table 10-4
m = 0.187;
A = 147; % psi * in^m
S_ut = A / d_wire^m;
S_y = 0.75 * S_ut;
S_sy = 0.5 * S_ut;

% body
F_body = 1000 * pi * d_wire^3 * S_sy / (8 * K_ * D); % need to adjust OAM

% tors @ B
C2 = 2 * r2 / d_wire;
K_b = (4 * C2 - 1) / (4 * C2 - 4);
F_b = 1000 * pi * d_wire^3 * S_sy / (8 * K_b * D); % adjust OAM

% norm @ A
C1 = 2 * r1 / d_wire;
K_a = (4 * C1^2 - C1 - 1) / ( 4 * C1 * (C1 - 1));
F_a = 1000 * S_y / ((16 * K_a * D / (pi * d_wire^3)) + 4 / (pi * d_wire^2));

F_min = min([F_body, F_a, F_b])



%% Part 5
% 5. What is the spring deflection corresponding to the load found in part (5)?

defl = (F_min - f_i) / k_est
