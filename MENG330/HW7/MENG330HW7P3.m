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

% dkameters, in
d_1 = 2;
d_2 = 2.472;
d_3 = 2.763;

% grav
grav = 32.174 * 12; % in/s2

% supports at l_0, l_6
% NEGLECT MASS OF SHAFT
% use moment equation centered at l_0
R_2 = (F_1 * l_2 + F_2 * l_4) / l_6; % support at l_6
R_1 = F_1 + F_2 - R_2; % support at l_0 aka origin

% conversion of rad/s --> rpm
rad_rpm = 60 / (2 * pi);


%% Deflection

% parms for funct
F = [F_1 F_2];
F_loc = [l_2 l_4];
d = [d_1 d_2 d_3 d_1];
d_loc = [l_0 l_1 l_3 l_5];
R_loc = [l_0 l_6];
L = l_6;

% exectute function given by prof on canvas
[x,y,dydx, M, MdEI, R, diam, EI] = ShaftDeflectionEnglish(F,F_loc,d,d_loc,R_loc,L);




%% Critical Speeds of attached elements

[zzz, y_idx_1] = min(abs(l_2 - x)); % wrt F_1
[zzz, y_idx_2] = min(abs(l_4 - x)); % wrt F_2

y_1 = y(y_idx_1); % wrt F_1
y_2 = y(y_idx_2); % wrt F_2

omega_attached = rad_rpm * sqrt(grav * (F_1 * y_1 + F_2 * y_2) / (F_1 * y_1^2 + F_2 * y_2^2));



%% Critical Speed of Shaft

% steel density
steel_density = 0.284; % lbf/in3

% use more accurate method
% divide shaft into four sections
d; d_loc; % initialized earlier. radiuses of shaft and where they become that radius

% find weight forces
f_steel_loc = [d_loc(1) / 2, (d_loc(1) + d_loc(2)) / 2, (d_loc(2) + d_loc(3)) / 2, (d_loc(3) + d_loc(4)) / 2];
f_steel = steel_density * pi / 4 * [d(1)^2 * d_loc(1), d(2)^2 * (d_loc(2) - d_loc(1)), d(3)^2 * (d_loc(3) - d_loc(2)), d(4)^2 * (d_loc(4) - d_loc(3))];

R_loc; L; % defined earlier

% run deflection equation again
[x,y,dydx, M, MdEI, R, diam, EI] = ShaftDeflectionEnglish(f_steel,f_steel_loc,d,d_loc,R_loc,L);


% initialize numerator and denominator values
numerator = 0; denominator = 0;


for i = 1:length(f_steel_loc)
    [zzz, y_idx] = min(abs(f_steel_loc(i) - x));
    y_val = y(y_idx);
    numerator = numerator + f_steel(i) * y_val;
    denominator = denominator + f_steel(i) * y_val^2;
end

omega_shaft = rad_rpm * sqrt(grav * numerator / denominator);


%% Combine shaft and attached critical speeds

% find total critical speed
omega_total = sqrt((1/omega_attached^2 + 1/omega_shaft^2)^-1);

% find operating speed
operating_speed = 0.5 * omega_total



