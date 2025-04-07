%% HW 04
% 4.4 & 4.5
clear 
close all
clc


%% 4.4 1D Interpolation on the Natural Line
% xi ∈ [-1,1]

%% 4.4.1 From x to xi
% Derive the linear interpolation function between a value of x ∈ [x1, x2]
% and ξ ∈ [−1, 1]. Show we can write it as

% derive
% xi = slope * x + intercept

%

% parameters
xbound = [0, 5]; % [x1,x2], variable
xibound = [-1,1]; % defined by problem
x = xbound(1):1:xbound(2);

% point one
p1 = [xbound(1), xibound(1)];

% point two
p2 = [xbound(2), xibound(2)];

% derive slope
% slope is the average rise per unit
% as two points, 2 in the numerator
% and the difference between the bounds in the denominator
slope = 2 / (xbound(2)-xbound(1));


% derive intercept
% intercept is average value between x1 and x2
% so sum of x1 and x2, divided by the difference between x1 and x2
intercept = (xbound(1) + xbound(2))/(xbound(1)-xbound(2));

% combine the two terms
xi1 = slope * x + intercept;

xi2 = (2/(xbound(2) - xbound(1))) * x + (xbound(1) + xbound(2))/(xbound(1) - xbound(2));




%% 4.4.2 Linear Shape Functions

% known:
% xi = [-1, 1];

% derive u(xi), w parms a & b
% u(xi) = a*xi + b

% because xi is two values, must create two equations
% size(a) = 2 x 1; size(b) = 1 x 1
% sysm xi
% u(xi1) = a * xi(1) + b;
% u(xi2) = a * xi(2) + b;


% define u values
% u1 = -a + b;
% u2 = a + b;

% add u1, u2 and solve for b
% u1 + u2 = a + b - a + b;
% b = (u1 + u2) / 2


% now solve for a
% u2 - u1 = a + b - (-a + b)
% 2a = u2 - u1
% a = (u2 - u1) / 2

% substitute u values back in
% u(xi) = ((u2-u1) / 2) * xi + (u1 + u2) / 2

% rearrainge equation
% u(xi) = u1 * (1-xi) / 2 + u2 (1+xi) / 2

% shift back into matrix form
% final result
% u(xi) = [(1-xi)/2 (1+xi)/2] [u1; u2]



%% 4.4.3 Quadratic Shape Functions

% If this time we know 3 points, we can determine a quadratic polynomial to
% interpolate over the field. By convention the first two known points are
% at the ends ξ1 = −1, ξ2 = 1 while the 3rd point is at the center
% ξ3 = 0. As before, assume that the field we are interpolating u has known
% points ui corresponding to the location ξi. Derive the quadratic
% interpolation function for a scalar function u(ξ) by determining the 
% constants a, b, and c

% u(xi) = a * xi^2 + b * xi + c

% Show that we can write the interpolation as

% u(xi) = [xi(xi - 1)/2 xi(xi + 1)/2 1-xi^2] [u1; u2; u3]

% if u = transpose([3,2,4]), plot u(xi)

% derivation

xi = [-1 1 0]; % define xi, as given in problem statement
u = transpose([3,2,4]); % define u, as given in problem



% define ui as u(xi(i))

% u1 = a * (-1)^2 + b * -1 + c = a - b + c
% u2 = a * (1)^2 + b * 1 + c = a + b + c
% u3 = a * (0)^2 + b * 0 + c = c

% take these and format into a system of equations

% u1 = a - b + c
% u2 = a + b + c
% u3 = c

% we can substitue u3 in for c, such that
% u1 = a - b + u3 ==> a - b = u1 - u3
% u2 = a + b + u3 ==> a + b = u2 - u3

% if we sum the equations
% (a + b) + (a - b) = u1 + u2 + 2 * u3
% 2a = u1 + u2 - 2 * u3
% a = (u1 + u2 - 2 * u3) / 2

% if we subtract the equations
% (a + b) - (a - b) = (u2 - u3) - (u1 - u3)
% 2b = u2 - u1
% b = (u2 - ui) / 2

% thus,
% a = (u1 + u2 - 2 * u3) / 2
% b = (u2 - u1) / 2
% c = u3

% if we were to sub the values into u(xi) = a * xi^2 + b * xi + c
% u(xi) = ((u1 + u2 - 2 * u3) / 2) * xi^2 + ((u2 - u1) / 2) * xi + u3

% you can then factor out this equation

% u(xi) = [xi(xi - 1)/2 xi(xi + 1)/2 1-xi^2] [u1; u2; u3]

u_xi = @(xi) ( ( u(1) + u(2) - 2 * u(3) ) / 2 ) * xi.^2 + ((u(2) - u(1)) / 2) * xi + u(3);

xi = -1:0.01:1;

% plot results
figure
plot(xi, u_xi(xi))
xlabel('xi')
ylabel('u(xi)')
title('Quadratic Interpolation Function')
grid on



%% 4.5 2D Interpolation on the Natural Square


%% 4.5.1 Bilinear interpolation
% Consider the square in Figure 4.2 where data is known at each of the
% four corners. We want to be able to linearly interpolate a scalar u(ζ, eta)
% across the square and we need to determine the coefficients in

% EQUATION 4.8
% u(xi,eta) = a * xi + b * eta + c + d * xi * eta

% Derive the interpolation function and show that we can write the results as

% EQUATION 4.9
% u(xi, eta) = [(1-xi)(1-eta)/4 (1-xi)(1+eta)/4 (1+xi)(1-eta)/4
% (1+xi)(1+eta)/4] * [u1; u2; u3; u4]

% xi = [-1, 1]
% eta = [-1, 1]

% Define u-vals (xi, eta)
% u1 = (-1,-1)
% u2 = (+1,-1)
% u3 = (+1,+1)
% u4 = (-1,+1)

% define u-values in terms of coefficients by plugging into 4.8
% u1 = - a - b + c + d; #1
% u2 = - a + b + c - d; #2
% u3 = + a - b + c - d; #3
% u4 = + a + b + c + d; #4

% these four u equations are a system of equations

% solve for a & c:

% sum #1 & #2
% u1 + u2 = (− a − b + c + d) + (− a + b + c − d)
% u1 + u2 = -2 * a + 2 * c ==>
% a = (u1 + u2 - 2 * c) / -2 OR
% -a + c = (u1 + u2) / 2; #5

% sum #3 & #4
% u3 ​+ u4​ = a − b + c − d) + (a + b + c + d)
% u3 ​+ u4​ = 2 * a + 2 * c ==> 
% a = (u3 ​+ u4 ​− 2 * c) / 2​ OR
% a + c = (u3 + u4) / 2; #6

% sum #5 & #6
% c = (u1 + u2 + u3 + u4) / 4; #7

% sub #7 into #6
% a + (u1 + u2 ​+ u3​ + u4​) / 4​ = (u3 ​+ u4) / 2
% a = (-u1 - u2 + u3 + u4) / 4; #8


% solve for b & d:

% subtract #2 from #1
% u2 - u1 = (-a - b + c + d) - (-a + b - c - d)
% u2 - u1 = -2 * b + 2 * d ==>
% -b + d = (u1 - u2) / 2​; #9

% subtract #3 from #4
% u4 ​− u3​ = (a + b + c + d) − (a − b + c − d)
% u4 - u3 = 2 * b + 2 * d ==>
% b + d = (u4 - u3) / 2; #10

% sum #9 & #10
% d = (u1 ​− u2 ​- u3 ​+ u4) / 4​​; #11

% sub #11 into #10
% b + (u1 ​− u2​ - u3 + u4​) / 4 ​​=(-u3​ + u4) / 2
% b = (-u1 + u2​​ - u3 + u4) / 4; #12

% sub coeffs a, b, c, d (#7, #8, #11, #12) into EQ4.8
% not writing this out for sake of brevity

% can simplify into EQ4.9
% u(xi, eta) = [(1-xi)(1-eta)/4 (1-xi)(1+eta)/4 (1+xi)(1-eta)/4 
% (1+xi)(1+eta)/4] * [u1; u2; u3; u4]




%% 4.5.2 Quadratic interpolation

% Let’s extend the process of §4.5.1 to quadratic interpolation, as
% numbered in Fig 4.3 where data is known at each of the four corners, the
% center of each edge, and the center. We want to be able to quadratically
% interpolate a scalar u(ζ, η) across the square and we need to determine
% the coefficients in
% u(xi, eta) = a1 + a2 * xi + a3 * eta + a4 * xi * eta + a5 * eta^2 +
% a6 * xi^2 + a7 * xi^2 * eta + a8 * xi * eta^2 + a9 * xi^2 * eta^2 (4.10)

% Figure 4.3: A square in natural coordinates (ξ, η) with 9 nodes.
% Derive the interpolation function and show that we can write the results as
% u(zeta, eta) = N(xi, eta) * u (4.11)
% where
% u = [ u1 u2 u3 u4 u5 u6 u7 u8 u9]T (4.12)
% and
% N(xi, eta) = [ N1 N2 N3 N4 N5 N6 N7 N8 N9] (4.13a)
% N1 = eta * xi (eta − 1) (xi − 1) / 4 (4.13b)
% N2 = eta * xi (eta − 1) (xi + 1) / 4 (4.13c)
% N3 = eta * xi (eta + 1) (xi + 1) / 4 (4.13d)
% N4 = eta * xi (eta + 1) (xi − 1) / 4 (4.13e)
% N5 = − eta (eta − 1) (xi − 1) (xi + 1) / 2 (4.13f)
% N6 = − xi (eta − 1) (eta + 1) (xi + 1) / 2 (4.13g)
% N7 = − eta (eta + 1) (xi − 1) (xi + 1) / 2 (4.13h)
% N8 = − xi (eta − 1) (eta + 1) (xi − 1) / 2 (4.13i)
% N9 = (eta − 1) (eta + 1) (xi − 1) (xi + 1) (4.13j)


% WRITE OUT DERIVATION
% 3 SOE'S w 3 unknowns each?

% Summary
% The N values are derived such that they are 0 elsewhere to their node.
% these N values are then put into their array
% the u coefficients are solved the same way the u coefficients were solved
% in 4.5.1. The u values were assigned points on the grid and the xi,eta
% values were then inserted to create a series of equations with 9
% equations and nine unknowns [a:i]
% solving these equations gave the coefficients to insert intp the initial
% equation, 4.10. These can be factored out into a matrix of xi and eta
% factors, also known as N / shape factors, and the u values. This creats
% the function 4.11
% the shape factor array in 4.5.2 is the same as the first array in 4.5.1,
% aka EQ4.9. It contains more values so it is shown as a matrix of
% equations due to length, rather than showing the complete equations.
% u_func = @(zeta, eta, xi) N(xi, eta) * u;

% xi = [-1, 1]
% eta = [-1, 1]

% Equation 4.10
% u(ξ, η) = a1 + a2ξ + a3η + a4ξη + a5ξ2 + a6η2 + a7ξ2η + a8ξη2 + a9ξ2η2

% Define u-vals (xi, eta)
% u1 = (-1,-1)
% u2 = (+1,-1)
% u3 = (+1,+1)
% u4 = (-1,+1)
% u5 = ( 0,-1)
% u6 = (+1, 0)
% u7 = ( 0,+1)
% u8 = (-1, 0)
% u9 = ( 0, 0)


% define u-values in terms of coefficients by plugging into 4.10
% u1 = a1 - a2 - a3 + a4 + a5 + a6 - a7 - a8 + a9; #1
% u2 = a1 + a2 - a3 - a4 + a5 + a6 - a7 + a8 + a9; #2
% u3 = a1 + a2 + a3 + a4 + a5 + a6 + a7 + a8 + a9; #3
% u4 = a1 - a2 + a3 - a4 + a5 + a6 + a7 - a8 + a9; #4
% u5 = a1 + a3 + a6; #5
% u6 = a1 + a2 + a5; #6
% u7 = a1 - a2 + a5; #7
% u8 = a1 - a3 + a6; #8
% u9 = a1; #9

% these values effectively give us 4 systems of equations with two
% equations each

% a1 already completely defined
% CAN SUB IN u9 FOR a1 AT ANY TIME - #9
% a1 = u9





% sum #6 & #7
% u6 + u7 = 2a1 + 2a5
% a5 = (u6 + u7) / 2 + u9 #10

% sum #5 + #8
% u5 + u8 = 2a1+ 2a6
% a6 = (u5 + u8) / 2 - u9 #11

% subtract #7 from #6
% u6 - u7 = 2a2
% a2 = (u6 - u7) / 2 #12

% subtract #8 from #5
% u5 - u8 = 2a3
% a3 = (u5 - u8) / 2 #13

% sum #2 & #3
% u2 + u3 = 2a1 + 2a2 + 2a5 + 2a6 + 2a8 + 2a9
% a8 + a9 = -u9 - a2 - a5 - a6 + (u2 + u3) / 2
% a8 + a9 = -(u6 - u7) / 2 - (u6 + u7) / 2 - (u5 - u8) / 2 + u9 - u9 + (u2
% + u3) / 2 - u9
% a8 + a9 = -u6 - (u5 - u8) / 2 + (u2 + u3) / 2 #14

% sum #1 & #4
% u1 + u4 = 2a1 - 2a2 + 2a5 + 2a6 - 2a8 + 2a9
% -a8 + a9 = -a1 + a2 - a5 - a6 + (u1 + u4) / 2
% -a8 + a9 = -u9 + (u6 - u7) / 2 - (u6 + u7) / 2 - (u5 + u8) / 2 + (u1 +
% u4) / 2 
% -a8 + a9 = -u7 - u9 - (u5 + u8) / 2 + (u1 + u4) / 2 #15

% sum #14 & #15
% 2a9 = -u6 - (u5 - u8) / 2 + (u2 + u3) / 2 - u9 - u7 - u9 - u6 - (u5 - u8) / 2 + (u2 + u3) / 2
% 2a9 = -2u6 - u7 - 2u9 + (u2 + u3)
% a9 = -u6 - u7 - u9 + (u2 + u3) / 2 #16

% sub #16 into #14
% a8 - u7 - u9 = -(u5 - u8) / 2
% a8 = (-u5 + 2u7 + u8 + 2u9) / 2 # 17


% subtract #4 from #1
% u1 - u4 = -2a3 + 2a4 - 2a7
% -a4 + a7 = -a3 - (u1 - u4) / 2 #18

% subtract #3 from #2
% u2 - u3 = -2a3 -2a4 - 2a7
% a4 + a7 = -a3 - (u2 - u3) / 2 #19

% sum #18 & #19
% 2a7 = -a3 - (u1 - u4) / 2 - a3 - (u2 - u3) / 2
% 2a7 = -2a3 - (u1 + u2 - u3 - u4) / 2
% a7 = -a3 - (u1 + u2 - u3 - u4) / 4 #20
% a7 = (-u1 - u2 + u3 + u4 + 2u5 - 2u8) / 4 #20

% sub #20 into #19
% a4 - a3 - (u1 + u2 - u3 - u4) / 4 = -a3 - (u2 - u3) / 2
% a4 = -(u2 - u3) / 2 + (u1 + u2 - u3 - u4) / 4
% a4 = (u1 -u2 - 3u3 - u4) / 4 # 21

% we can then plug these nine a-values back into the initial equation

% a1 = u9
% a2 = (u6 - u7) / 2 #12
% a3 = (u5 - u8) / 2 #13
% a4 = (u1 -u2 - 3u3 - u4) / 4 # 21
% a5 = (u6 + u7) / 2 + u9 #10
% a6 = (u5 + u8) / 2 - u9 #11
% a7 = (-u1 - u2 + u3 + u4 + 2u5 - 2u8) / 4 #20
% a8 = (-u5 + 2u7 + u8 + 2u9) / 2 # 17
% a9 = -u6 - u7 - u9 + (u2 + u3) / 2 #16

% u(ξ, η) = u9 + ((u6 - u7) / 2)ξ + ((u5 - u8) / 2)η + ((u1 -u2 - 3u3 - u4) / 4)ξη
% + ((u6 + u7) / 2 + u9)ξ2 + ((u5 + u8) / 2 - u9)η2 
% + ((-u1 - u2 + u3 + u4 + 2u5 - 2u8) / 4)ξ2η
% + ((-u5 + 2u7 + u8 + 2u9) / 2)ξη2
% + (-u6 - u7 - u9 + (u2 + u3) / 2)ξ2η2

% We can factor out the u-values via a matrix, along with the xi and eta
% values. This operation will result in 2 matrices - our shape factor
% matrix, N, and our scalar matrix, u.




%% 4.5.3 Example

% Consider the information of a 4-node element in Table 4.2. Using the
% relations derived in the previous problems and sections, determine the following.
% 1. The value of u at xi = 0.1, eta = −0.2
% 2. The value of u at x = 0, y = 0
% 3. Visualize u(x, y) in a contour plot.
 

% Table 4.2:
% Node  x       y       u
% 1     −0.5    −0.3    9
% 2     1.0     −0.1    15
% 3     0.7     1.5     6
% 4     −0.8    0.9     −1

% Use equation from 4.5.1
u1 = 9;
u2 = 15;
u3 = 6;
u4 = -1;

% goes x = [-0.5, 1] & y = [-0.3,1.5]
u_func = @(xi, eta) [(1-xi) * (1-eta)/4 (1-xi) * (1+eta)/4 (1+xi) * (1-eta)/4  (1+xi) * (1+eta)/4] * [u1; u2; u3; u4];

x = [-0.5, 1];
y = [-0.3,1.5];
xi = x * (4/3) - (1/3);
eta = y * (10/9) -(2/3);


% 1: 
u_func(0.1,-0.2)


% 2: want (x,y) = (0,0)
% therefore (xi,eta) = (-1/3, -2/3)
u_func(-1/3,-2/3)


% 3:
figure
hold on
h = fcontour(u_func);
xlabel('xi')
ylabel('eta')
xlim([-1,1])
ylim([-1,1])
% clabel(h.LevelList) % attempt at getting contour labels. did not succeed
% + internet & matlab function definitions did not help