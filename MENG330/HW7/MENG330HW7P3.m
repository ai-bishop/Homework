%% Requirements & problem statement
% 18lbf gear left
% 32lbf gear right
% estimate: [rpm]
% first critical speed due to the loads
% critical speed without the loads
% critical speed of combination
% operating speed


% locations, in, from the leftmost point aka origin
l_0 = 0;
l_1 = 1;
l_2 = 2; % location of 18lbf force
l_3 = 9;
l_4 = 14; % location of 32lbf
l_5 = 15;
l_6 = 16;

% forces, lbf
F_1 = 18;
F_2 = 32;

% radii, in
r_1 = 2;
r_2 = 2.472;
r_3 = 2.763;

% grav
grav = 32.174 * 12; % in/s2

% supports at l_0, l_6
% NEGLECT MASS OF SHAFT
% use moment equation centered at l_0
R_2 = (F_1 * l_2 + F_2 * l_4) / l_6; % support at l_6
R_1 = F_1 + F_2 - R_2; % support at l_0 aka origin

%% Deflection




















