%% Initialize problem
% create variables
Sy = 420; % N/mm^2 aka MPA
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
syms t % to stand in for x
V(x) = R_b * heaviside(x) + R_c * heaviside(x - BC) - Force * heaviside(x - (BC + CD));
M(x) = int(subs(V(x), x, t), t, 0, x);



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
sigma_mom = M_max * (diam / 2) / inertia;

vm_moment = sqrt(sigma_mom^2 + 3 * tau_tot^2);

% solve
diam_s = double(vpasolve(Sy / vm_shear == SafetyFactor, diam, [0 inf]));
diam_m = double(vpasolve(Sy / vm_moment == SafetyFactor, diam, [0 inf]));

diam_static = max(diam_s, diam_m)


%% Part B: Mod-Goodman

syms diam % everything is a function of diameter here

% assumption: infinitie life desired
% moment gives a larger force than shear as per Part A

% calculate without stress concs
sigma_m = sqrt(3*ttorsion^2);
sigma_a_0 = sqrt(sigma_mom^2);
% factor in stress concs later once have

% stress conc @ shoulders + gear, T7-1
K_t = 2.7;

% conversion from metric to english as the q equations are in english
mm_in = 1/25.4; % mm to inches conversion factor
mpa_ksi = 0.145038; % ksi per mpa

% assume r/d = 0.02 - justifiable as sharp radii value
% find radius in inches
rad_in = diam * 0.02 / mm_in; % in

root_a = 0.246 - 3.08e-3 * (Sut * mpa_ksi) + 1.51e-5 * (Sut * mpa_ksi)^2 - 2.67e-8 * (Sut * mpa_ksi)^3; % in^0.5

% calculate q, given equation
q = 1 / (1 + root_a / sqrt(rad_in));

% find K_f
K_f = 1 + q * (K_t - 1);

% find alternating stress at concentration
sigma_a = K_f * sigma_a_0;


% shaft is steel, so Se' = Sut / 2 if Sut < 1400 mpa
Se_0 = Sut / 2;

% find k factors
% k factors are taken from textbook tables, also on slides
% assume shaft is machined as not specified in problem
k_a = 3.04 * Sut ^ (-0.217);

% diameter ~ 32mm so falls in first range for k_b
k_b = 1.51 * diam^-0.157;
% k_b = (diam/7.62) ^ -0.107;

% using von mises stress so k_c = 1
k_c = 1;

% no temperature stated - assume operating at room temperature
k_d = 1;

% no reliability mentioned - assume 50%
k_e = 1;

% no misc vals
k_f = 1;

k_tot = k_a * k_b * k_c * k_d * k_e * k_f;
Se = k_tot * Se_0;

diam_fatigue = vpasolve(sigma_a / Se + sigma_m / Sut == 1 / SafetyFactor, diam, [0 inf])


























