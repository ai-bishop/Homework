%% Shigley's 8.76


%% Initialize values
% vertical channel 152x76
% cantilever beam
% channel: AISI 1006 HR
% bar: AISI 1015 HR
% shoulder bolts: M10x1.5 ISO 5.8
% assume bolt threads dont extend into the joint
% n = 2
% find F_min, lowest safe force that can be applied to the cantilever

% distances, mm
AO = 50;
OB = 50;
BC = 26;
CF = 125;
width = 50;
thickness = 12;
diam = 10;
n = 2;
t = 6.4;

Sy_bolt = 420; % MPA, T8-11
Sy_channel  = 170; % TA-20
t_channel = 6.4; % TA-7
Sy_cantilever = 190;





%% Find:
% F_min
% shear on bolt
% bearing on bolt
% bearing on channel
% bearing cantilver
% bending of cantilever


%% Solution:

% force divided over three bolts = F/3
% Fa_prime = Fb_prime = Fc_prime = F/3
% M = (50 + 26 + 125)F = 201F
% Fa_prime2 = Fc_prime2 = 201F / 2*50 = 2.01F
% max force: 
% Fc = Fc_prime + Fc_prime2 = 2.343F

% shear on bolts:
A_s = pi * diam^2 / 4;
Ssy_bolt = 0.577 * Sy_bolt;
% tau = Fc/A_s = Ssy / n
% Fc = 2.343F
F_shear_bolt = (Ssy_bolt / n) * (A_s / 2.343)

% bearing on bolts
A_b1 = t * diam;
F_bearing_bolt = (Sy_bolt / n) * (A_b1 / 2.343)

% bearing on channel
A_b2 = t * diam;
F_bearing_channel = (Sy_channel / n) * (A_b2 / 2.343)

% bearing on cantilever
A_b3 = thickness * diam;
F_bearing_cantilever = (Sy_cantilever / n) * (A_b3 / 2.343)

% Bending of Cantilever @ C
I = (1/12) * thickness * (width^3 - diam^3);
F_bending_cantilever = (Sy_cantilever / n) * (I / (151 * (width / 2)))



F_min = min([F_shear_bolt F_bearing_bolt F_bearing_channel F_bearing_cantilever F_bending_cantilever])


