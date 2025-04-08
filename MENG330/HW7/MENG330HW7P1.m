%% Initialize problem
% create variables
Sy = 420; % N/mm^2
Sut = 560; % N/mm^2
SafetyFactor = 2.5; % safety factor
F_angle = 20; % degrees
Torque_a = 340 * 10^3; % N*mm
BC = 250; % mm
CD = 100; % mm
D_diam = 150; % mm

% D.E.: n = Sy/sigma
% MGP sigma_a_0 / Se + sigma_m_0 / Sut = 1/n


%% Initial Properties

% find force in shaft
% torque = r * F * sin theta
% solve for F

Force = Torque_a / (D_diam / 2 * cosd(F_angle));

% calculate support reactions
syms R_b R_c % N
% forces
sup_eqns = [ R_b + R_c == Force ...
    R_c * BC == Force * (BC + CD)];

sup_eqns_soln = solve(sup_eqns);
R_b = double(sup_eqns_soln.R_b);
R_c = double(sup_eqns_soln.R_c);


% Calculate shear and bending moments
syms V(x) M(x) % N, N*mm 
% x dist into BC from B
syms x2 % to stand in for x
V(x) = R_b * heaviside(x) + R_c * heaviside(x - BC) - Force * heaviside(x - (BC - CD));
M(x) = int(subs(V(x), x, x2), x2, 0, x);



%% Part A: Distortion Energy

% can ignore stress concentrations because steel is sufficiently ductile and this is static failure
% assume diameter is within [25 75] mm
% worst case scenario is at B because that is where the greatest shear and moment magnitudes are

syms diam % functions of diameter

% define values at B
x_max = 250; % mm
V_max = double(subs(V, x, x_max));
M_max = double(subs(M, x, x_max));

polar_inertia = pi * diam^4 / 32;
inertia = polar_inertia / 2;
area = pi * (diam / 2 )^2;

ttorsion = Torque_a * (diam / 2) / polar_inertia; % uniform for both shear and moment

% check at max shear and max moment
% max shear
tshear = 4 * V_max / (3 * area);
tau_tot = ttorsion + tshear;

vm_shear = sqrt(3 * tau_tot ^ 2);

% Max moment
sigma_de = M_max * (diam / 2) / inertia;

vm_moment = sqrt(sigma_de^2 + 3 * tau_tot^2);

% solve
diam_s = double(solve(Sy / vm_shear == SafetyFactor));
diam_m = double(solve(Sy / vm_moment == SafetyFactor));

diam_static = max(diam_s, diam_m)









%% Part B: Mod-Goodman

