%% Problem 2
clear
clc

%% initialize values
mm_in = 1/25.4; % mm to inches conversion factor


% distances of the center fo the parts of interest from the leftmost point, in
loc_keyfan = 1;
loc_keyrgear = 9.35;
loc_notchrgear = 9.995;
loc_filletrgear = 10.67;

% locations
l_1 = 0;
l_2 = 2;
l_3 = 2.75;
l_4 = 2.9;
l_5 = 8.5;
l_6 = 9.985;
l_7 = 10.085;
l_8 = 10.67;
l_9 = 12.12;
l_10 = 12.87;


% diameters and their locations, in
d_1 = 1;
d_2 = 1.181;
d_3 = 1.7;
d_4 = 1.750;
d_5 = 2;
d_6 = 1.4;

d_loc_1 = l_1;
d_loc_2 = l_2;
d_loc_3 = l_4;
d_loc_4 = l_5;
d_loc_5 = l_7;
d_loc_6 = l_8;
d_loc_7 = l_9;

d = [d_1 d_2 d_3 d_4 d_5 d_6 d_1];
d_loc = [d_loc_1 d_loc_2 d_loc_3 d_loc_4 d_loc_5 d_loc_6 d_loc_7];

% keyways
key_1 = [1/4 1/8];
key_2 = [3/8 3/16];

% radii of the chamfers
rad_1 = 1/16;
rad_2 = 1/32;
rad_3 = 1/8;
rad_4 = 0.1;
% rad_5 = rad_3
% rad_6 = rad_2

pitch_diam = 8; % in

% forces from gear lbf
F_rad = 230;
F_tang = 633;

% combine into once force
F_comb = sqrt(F_rad^2 + F_tang^2);

% find torque
Torque = F_tang * pitch_diam / 2;

% material properties, ksi
Sut = 68;
Sy = 57;
Se_prime = 0.5 * Sut; % Sut is within the range to have Se' be 1/2 Sut
% generate k_a factor - material property
k_a = 2.7 * Sut ^ -0.265; % 0.833

%% Support Reactions
% two supports (ball bearings), one force (gear)
% force is at gear, torque
R_1_loc = mean([l_2 l_3]);
R_2_loc = mean([l_9 l_10]);

% moment equation centered @ R_1_loc
R_2 = (F_comb * (loc_keyrgear - R_1_loc)) / (R_2_loc - R_1_loc);

% find R_1
R_1 = F_comb - R_2;

% find moments - shear and bending moment diagrams on paper
% units in lbf in
Mom_keyfan = 0;
Mom_keyrgear = 1459;
Mom_notchrgear = 1115;
Mom_filletrgear = 845;


% due to shaft rotation
% torsion is sigma_m
% and bending is sigma_a



%% 1. Keyway of Fan Blade, n
% tau = Tr/J
J_keyfan = pi * d_1^4 / 32;
tau_torsion = Torque * (d_1 / 2) / J_keyfan;

vm_sigma_fan = tau_torsion * sqrt(3) / 1000; % convert to kips

n_de_fan = Sy / vm_sigma_fan;


% find sigma_a and sigma_m
sigma_m_fan = vm_sigma_fan;

I_keyfan = J_keyfan/2;

sigma_a_fan_0 = Mom_keyfan * (d_1 / 2) / I_keyfan; % no moment here

% assume r/d = 0.02 - typical keyway size
rad_in = d_1 * mm_in * 0.02;

K_t = 2.2;
root_a = 0.246 - 3.08e-3 * (Sut) + 1.51e-5 * (Sut)^2 - 2.67e-8 * (Sut)^3; % in^0.5
q = 1 / (1 + root_a / sqrt(rad_in));

K_f = 1 + q * (K_t - 1);
sigma_a_fan = sigma_a_fan_0 * K_f;

% define k_b to k_f
k_b = (1/0.3) ^ -0.107; % d = 1in
k_c = 1; % bending
k_d = 1; % no temp delta spec'd
k_e = 1; % no rel spec'd
k_f = 1; % no extra factors spec'd
Se_fan = k_a * k_b * Se_prime;


n_g_fan = 0.5 * (Sut / sigma_m_fan)^2 * (sigma_a_fan / Se_fan) * (-1 + sqrt(1+(2 * sigma_m_fan * Se_fan / (Sut * sigma_a_fan))));

n_1 = min([n_de_fan n_g_fan]);


%% 2. Notch right of right gear, n

% Fig A15-14 applicable here
% r = 0.1, D = 1.75, d = 1.55
d_notch = d_4 - 2 * rad_4;
r_d = rad_4/d_notch;
D_d = d_4/d_notch;

K_t = 2.1; % ref Fig A15-14
q = 0.76; % Fig 6-20
K_f = 1 + q * (K_t - 1); 

% k_b factor - shape dependent
k_b = (d_notch / 0.3) ^ -0.107;

Se_notch = k_a * k_b * Se_prime;

% now compute e
% 
% 
% 
% quation for safety factor

Torque;

I_notch = pi * d_notch^4 / 64;


sigma_a_notch_0 = Mom_notchrgear * (d_notch / 2) / I_notch;

sigma_a_notch = K_f * sigma_a_notch_0;



sigma_m_notch = tau_torsion * sqrt(3) / 1000; % convert from psi to ksi to keep units consistent

n_de_notch = Sy / sigma_m_notch;

n_g_notch = 0.5 * (Sut / sigma_m_notch)^2 * (sigma_a_notch / Se_notch) * (-1 + sqrt(1+(2 * sigma_m_notch * Se_notch / (Sut * sigma_a_notch))));

n_3 = min([n_de_notch n_g_notch])








































%% 3. Find Deflections


















