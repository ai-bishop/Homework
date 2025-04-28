%% Problem 2

% Shigley’s 13-34 A compound reverted gear train is to be designed as a speed
% increaser to provide a total increase of speed of exactly 40 to 1. Use spur gears with a
% 20° pressure angle and a diametral pitch of 6 teeth/in. Have gear set 2-3 have the larger
% gear ratio of the two gear sets. Have the gear ratios of the sets be integer values.
% 1. Specify appropriate numbers of teeth to minimize the gearbox size while
% avoiding the interference problem in the teeth.
% 2. Determine the gearbox outside dimension Y, assuming the wall thickness of the
% box is 0.75 in, and allowing 0.5 in clearances between the tips of the gear teeth
% and the walls.


speed_change = 40; % ratio
pressure_angle = 20; % degrees
teeth_min = 17; % minimum teeth without interference for pressure angle
p_d = 6; % teeth/in
% gears 2-3 have larger ratio, integers


%% Part 1
% spec app. # of teeth to min size

% N4 / N5 > N2 / N3
% (N4/N5) * (N2/N3) = 40
% use ratios of 8 & 5
% sum to 40 but are not too far apart
% N2 = 5N3
% N4 = 8N5
% N2 + N3 = N4 + N5
% therefore N5 = 7/4 N3
% choose N3 = 4
% but this causes interference! 4 & 7 dont work!
% scale up by 3
N3 = 12; % choose
N5 = 21;
N2 = 105;
N4 = 96;

n_teeth = N2 + N3 + N4 + N5; % answer

%% Part 2
% find gearbox outside dimension Y. t = 0.75in, 0.5in clearancee

d_2 = N2 / p_d;
d_3 = N3 / p_d;
d_4 = N4 / p_d;
d_5 = N5 / p_d;

cl_ = 1; % clearance sum top + bottom
th_ = 1.5; % side thicknesses
G1_rad = d_2 / 2;
G2_rad = d_5 / 2;

Y = cl_ + th_ + G1_rad + G2_rad;








