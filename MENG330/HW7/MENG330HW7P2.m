%% Problem 2


%% initialize values

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

%% start actual work!


































