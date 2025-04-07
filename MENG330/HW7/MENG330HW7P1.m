%% Initialize problem
% create variables
Sy = 420; % N/mm^2
Sut = 560; % N/mm^2
n = 2.5; % safety factor
F_angle = 20; % degrees
Torque_a = 340 * 10^3; % N*mm
BC = 250; % mm
CD = 100; % mm
D_diam = 150; % mm

% D.E.: n = Sy/sigma
% MGP sigma_a_0 / Se + sigma_m_0 / Sut = 1/n


%% Part A: Distortion Energy

sigma = Sy/n; % sigma required

% find force in shaft
% torque = r * F * sin theta
% solve for F

Force_BC = Torque_a / (D_diam * sind(F_angle))






