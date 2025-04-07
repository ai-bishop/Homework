%% Shaft Deflection Driver

clc
clear

%% checking xy plane deflections/slopes
%Forces
F=[197, 885]; %magnitudes of forces are applied
Floc=[2.75, 8.5]; %location of where forces are applied
d=[1,1.4,1.625,2,1.625,1.4,1]; %diameters of the shaft
dloc=[0,1.25,2.0,3.5,7.5,9.5,10.25]; %locations of start of diameter change
Rloc=[0.75,10.75]; %locations of bearings
L=11.5; %length of the shaft

[x1,y1,dydx1, M1, MdEI1, Ry, diam1, EI1]=ShaftDeflectionEnglish(F,Floc,d,dloc,Rloc,L);

%% checking xz plane deflections/slopes
%Forces
F=[540, -2431];
Floc=[2.75, 8.5];

[x2,y2,dydx2, M2, MdEI2, Rz, diam2, EI2]=ShaftDeflectionEnglish(F,Floc,d,dloc,Rloc,L);

%% vector summing to find total slopes and deflections
x=x1; %should be same x value for both dimensions
y=sqrt(y1.^2+y2.^2); %resultant transverse deflection
dydx=sqrt(dydx1.^2+dydx2.^2); %resultant angular deflection
M=sqrt(M1.^2+M2.^2); %resultant moment

figure
subplot(3,1,1)
plot(x,diam1/2,'r')
title('Half-Shaft geometry')
xlim([0,L])
ylim([0,1.2*max(d)/2])
ylabel('Radius')

subplot(3,1,2)
plot(x,y,'g')
xlim([0,L])
title('Magnitude of total deflection')
ylabel('Deflection')

subplot(3,1,3)
plot(x,dydx,'b')
xlim([0,L])
title('Magnitude of total Slope')
ylabel('Slope (rad)')