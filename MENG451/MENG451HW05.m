%% HW 05
% 5.4, 5.6, 5.7, and 5.8

%% 5.4
% Consider the mechanism shown in Figure 5.1, where a spring k holds the
% follower A against a cam. The equation describing the profile of
% the simple cam is
% R^2 + 6Rcosθ − 27 = 0 (5.13)
% where R is in inches. The weight of the follower is w, and the
% coefficient of dry sliding friction between the cam and follower is μk.
% The known parameter values are given in Table 5.2

%% 5.4.1 Kinematics

% known parameters
% m = 0.5; % lbs
% uk = 0.2; % unitless

% Let’s consider the steady-state motion of the cam follower system when
% the angular velocity ˙θ(t) = ω0 is constant.


% 1. Determine the position, velocity, and acceleration of the follower as
% a function of the cam angle θ.
% angular velocity ˙θ(t) = ω0 is constant.



% position = (Rcos(theta), Rsin(theta)) (x,y)
% velocity (m/s) = angular vel (rad/s) * distance from R (m)
% accel. (m/s2) = d/dt (velocity) = angular vel (rad/s) * 
% rate of change of dist from center (m/s)

% predefine t
t = 0;

% angular velocity; deg/s
omega = 1; % arbitrary value, does not matter.

% angle with respect to center of camshaft; 
theta = @(t) omega * t; % angle

% diameter of circle; in
Diam = 6;

% distance from center of camshaft; in
R = @(theta) 3 + 9 * sind(theta / 2);
R0 = 3;

% calc deriv of R
Rprime = @(theta) (9/2) * cosd(theta/2);

% define values for theta
theta_values = linspace(0,360,1000);

% function for position
Position = @(theta) [Diam * cosd(theta), Diam * sind(theta)];

% function fo0r velocity
Velocity = @(theta) R(theta) * omega;

% function for acceleration
Acceleration = @(theta) Rprime(theta) * omega;







% precalc the position, velocity, and acceleration
positions = arrayfun(@(theta) Position(theta), theta_values, 'UniformOutput', false);
positions = cell2mat(positions')';  % Convert cell array to matrix

velocities = arrayfun(@(theta) Velocity(theta), theta_values);
accelerations = arrayfun(@(theta) Acceleration(theta), theta_values);


% 2. Make plots of the above results, suitably normalized, for one cycle
% of the motion. Hint, the acceleration curve should look like Fig 5.2
% theta [0,360]

% create figure
figure
hold on

% plot position
plot(theta_values, positions(1, :)/R0, 'b', 'DisplayName', 'Position (X)', 'LineWidth', 1.5)
plot(theta_values, positions(2, :)/R0, 'c', 'DisplayName', 'Position (Y)', 'LineWidth', 1.5)

% plot velocity
plot(theta_values, velocities/R0, 'r', 'DisplayName', 'Velocity', 'LineWidth', 1.5)

% plot acceleration
plot(theta_values, accelerations/R0, 'g', 'DisplayName', 'Acceleration', 'LineWidth', 1.5)

% graph labels
grid on
xlim([0,360])
xlabel('Angle (Degrees)')
ylabel('Magnitude over R0')
legend
title('Position, Velocity, and Acceleration of the Cam')




%% 5.4.2 Dynamics

% known parameters
m = 0.5; % lbs
uk = 0.2; % unitless

%Let’s now consider the forces and losses experienced by this system
% for a particular running condition.


% 1. If the constant angular speed is ω0 = 2500 rpm, determine the lower
% bound of the spring constant that ensures the follower remains attached
% to the cam. Assume that the spring is compressed 0.25 inches at the
% smallest value of R.

omega0 = 2500/60; % rps

% force of spring out will be 0.25k (lbf) @ R = 3 (aka R0)
% Fspring = k (lb/ft) * x (ft) 
% Fspring,max = 0.25k

Fcam = @(theta) uk * m * (omega0)^2 * cosd(theta) * pi; % (distance per rotation)
Fcam_max = Fcam(0)

% set these two equal
k_min = 4 * Fcam_max / 12 % lbf/in
% fitz got ~70


% 2. How does ω0 compare to the natural frequency √k/m of the follower
% modeled as a mass-spring system?

% natural frequency = sqrt(k/m)
NaturalFrequency_min = sqrt(12 * 12 * k_min/m) % adjust for in ==> ft; (lbf/in / lbm)  [=] 1/s

% omega0 shares the same units as the natural frequency of the follower
% modeled as a mass spring system.
% if the natural frequency is less than omega0, then there is a resonance
% and changes to the system need to be made. as NaturalFrequency_min >
% omega0, there is no resonance and no need to adjust the system in its
% current state.



% 3. Determine the energy lost to friction per cycle for a reasonable
% spring value above the minimum found in the previous part.
% Hint: Work is ∫ P dt, where P is power.

% v = velocity of the follower
v = @(theta) Rprime(theta) * omega0;

% power of friction
P = @(theta) Fcam(theta) .* v(theta);

% energy lost to friction in a single cycle
E_Friction = integral(P, 0, 360)


% 4. What would change in the friction analysis if we used a lubricated
% model instead of dry-friction

% If we were to use a lubricated model, the coefficient of friction would
% greatly decrease. This would result in a spring with a decreased spring
% constant, a lesser Natural Frequency, and less energy lost to friction
% per cycle.




%% 5.6 The period of a simple pendulum
% Consider an ideal pendulum with its mass m lumped a length l from the
% frictionless hinge. The angle from the vertical is θ, and the starting
% state of the pendulum is θ(0) = θ0, ˙θ(0) = 0. The equation of motion
% is well-known to be
% θ'' + g/l sin θ = 0 (5.21)
% The linear approximation for small angles, gives the natural frequency
% of the oscillation as
% ω0 = √ g/l = 2π/T0 (5.22)
% However, as the angle increases we know that approximation will break
% down. Through some clever manipulation of the equation of motion, it
% can be shown that the period of oscillation is (Nayfeh and Mook, 1995)
% T = 4√l/g ∫(0, π/2) dθ/√(1 − k^2 sin2θ) (5.23)
% where k = sin(θ0/2).
% This integral is known as an elliptic integral of the first kind.
% Make a plot showing how much T /T0 changes as a function of θ0.
% Discuss the range of θ0 for which the linear approximation is useful.


% initalize constants
l = 1;
g = 9.81; % gravity

% approximation for small angles
% initial angular frequency
omega0 = sqrt(g/l);
% inital period
T0 = 2 * pi / omega0;

iter = 0;
% define theta0 - initial angle - main variable
for theta0 = 0:0.01:pi/2

    iter = iter+1;
    % define k
    k = sin(theta0/2);

    % create equation for integral
    dtheta = @(theta) 1 ./ sqrt(1 - k^2 * (sin(theta)) .^ 2);

    T = 4 * (1/omega0) * integral(dtheta,0,(pi/2));
    % can sub in omega0 for sqrt(l/g) bc eqn 5.22



    Tnorm(iter) = T / T0;



end

rads = 0:0.01:pi/2;

figure
plot(rads,Tnorm)
xlabel('Rads')
ylabel('Normalization of T (T/T0)')
title('Normalized T as a function of radians')


% At Theta0 = 0, the natural frequency is the same as ideal. At a Theta0 of
% pi/2, the natural frequency is 18% above the ideal.
% Normally, small angle approximation is valid up to 0.26 rads or 15
% degrees. An initial angle of 0.26 rads of this gives a normalized period
% of 1.0039, or 0.37% greater than an initial angle of zero, based on the
% period of the pendulum. Based on this, we get an acceptable margin of
% error of 0.4% on the normalized period of the pendulum. This gives us a
% range of theta0 of [0,0.26] rads for small angle approximation.




%% 5.7 Center of Mass
% close all

% Recall that the definition of the center of mass of a body is
% rg = 1/m ∫(Ω) r dm (5.24)
% where the total mass can be computed from the density ρ
% m =∫(Ω) dm = ∫(Ω) ρ dV . (5.25)

% Let’s consider a single element from a two-dimensional finite
% element mesh. We have the shape-functions (interpolation functions)
% N(ξ, η) so we can build the position of a point inside the element with
% r = xˆi + yˆj → r = [ N(ξ, η)x; N(ξ, η)y] (5.26)

% In order to compute either integral above, we need to do a
% change of variables
% dV = dx dy = J dξ dη (5.27)
% where J is the determinant of the Jacobian matrix
% J = [ x,ξ x,η; y,ξ y,η] = [ N,ξ x N,η x; N,ξ y N,η y] (5.28)

% Now we can recast the integrals as
% m = ∫(-1,1) ∫(-1,1) ρJ dξ dη (5.29)
% rg = 1/m ∫(-1,1) ∫(-1,1) ρ [ N(ξ, η)x; N(ξ, η)y] J dξ dη (5.30)

% which is computable with Gauss-Legendre quadrature.


% 5.7.1 Code
% Build a function that takes the following inputs for a 4-node quadrilateral
% 1. x : the list of the x node locations
% 2. y : the list of the y node locations
% 3. ρ: the function of the density ρ(x, y)
% 4. N : the number of integration points to use in each direction for
% the Gauss-Legendre rule.
% and returns the mass and the location of the center of mass.


% shape functions?
% have final integration be from -N to N and divide jacobian J by N?



function [Mass, CenterOfMass] = COM(x, y, rho, N)

    % x, y: Coordinates of the 4 nodes [x1, x2, x3, x4]
    % rho: Function handle for the density, rho(x, y)
    % N: Number of Gauss-Legendre integration points in each direction
    
    % find Gauss-Legendre points + weights
    [xz, wi] = lgwt(N,-1,1);

    % initialize mass + moment values
    mass = 0;
    mx = 0;
    my = 0;

    % loop over GL points
    
    for i = 1:N
        for j = 1:N
        
            xi = xz(i);
            eta = xz(j);


            % shape functions
            N1 = 0.25 * (1 - xi) * (1 - eta);
            N2 = 0.25 * (1 + xi) * (1 - eta);
            N3 = 0.25 * (1 + xi) * (1 + eta);
            N4 = 0.25 * (1 - xi) * (1 + eta);

            % shape function derivs
            dN1_dxi = -0.25 * (1 - eta);
            dN2_dxi =  0.25 * (1 - eta);
            dN3_dxi =  0.25 * (1 + eta);
            dN4_dxi = -0.25 * (1 + eta);
            
            dN1_deta = -0.25 * (1 - xi);
            dN2_deta = -0.25 * (1 + xi);
            dN3_deta =  0.25 * (1 + xi);
            dN4_deta =  0.25 * (1 - xi);
            
            % jac matrix [dxdxi dxdeta; dydxi dydeta]
            dxdxi = dN1_dxi * x(1) + dN2_dxi * x(2) + dN3_dxi * x(3) + dN4_dxi * x(4);
            dxdeta = dN1_deta * x(1) + dN2_deta * x(2) + dN3_deta * x(3) + dN4_deta * x(4);
            dydxi = dN1_dxi * y(1) + dN2_dxi * y(2) + dN3_dxi * y(3) + dN4_dxi * y(4);
            dydeta = dN1_deta * y(1) + dN2_deta * y(2) + dN3_deta * y(3) + dN4_deta * y(4);
            
            % jac
            Jac = [dxdxi dxdeta; dydxi dydeta];

            % det jac
            detJ = det(Jac);
            
            % eval density
            x_current = N1 * x(1) + N2 * x(2) + N3 * x(3) + N4 * x(4);
            y_current = N1 * y(1) + N2 * y(2) + N3 * y(3) + N4 * y(4);
            rho_current = rho(x_current, y_current);

            % find weight and mass
            weight = wi(i) * wi(j) * detJ;
            mass = mass + rho_current * weight;
            mx = mx + x_current * rho_current * weight;
            my = my + y_current * rho_current * weight;

        end

    end

    % label outputs
    Mass = mass;

    % Center of Mass in Gauss-Legendre
    CenterOfMassGL = [mx / mass, my / mass];

    % Convert to Center of Mass in x,y by multiplying by the jacobian
    CenterOfMass = CenterOfMassGL * Jac;

end





% 5.7.2 Example case
% To practice on, consider the node positions given in Table 5.3 for a
% 4-node quadrilateral. The (admittedly goofy) density function is
% ρ(x, y) = (3x + 7)(sin(3y) sin(3x) + 1) . (5.31)
% Compute the mass of the element and determine the location of the
% center of mass. Visualize the results

% Table 5.3: Node locations for the 4-node quadrilateral.
% Node x y
% 1 −0.5 −0.3
% 2 1.0 −0.1
% 3 0.7 1.5
% 4 −0.8 0.9

% define values
xvals = [-0.5 1.0 0.7 -0.8];
yvals = [-0.3 -0.1 1.5 0.9];

% define density function
rho = @(x,y) (3 .* x + 7) .* (sin(3 * y) .* sin(3 * x) + 1);

% Number points
N = 10;

% initialize trigger
trig = true;
Mass = 0; % initalize mass such that it forces next iteration

% run iteration to find necessary number of integration points
while trig == true

    N = N + 1;
    
    Mass1 = Mass; % capture old Mass value

    % run function
    [Mass, CenterOfMass] = COM(xvals, yvals, rho, N);

    % Check how much mass value has changed

    MassDiff = abs((Mass - Mass1) / Mass1);

    if MassDiff < 0.1 % 1% threshold

        trig = false;

    end


    if N > 100

        trig = false;
        fprintf('Maximum iterations reached: %f', N)

    end


end

Mass
CenterOfMass


% madelyn got ~12 for mass
% CoM should be pos pos

% creat figure
figure
hold on
grid on
h = fcontour(rho,[min(xvals) max(xvals) min(yvals) max(yvals)]);
scatter([xvals CenterOfMass(1)], [yvals CenterOfMass(2)])



%% 5.8 Inertia Tensor

% Let’s extend the discussion of center of mass to compute the mass
% moment of inertia tensor. Recall the definition of the inertia tensor
% in 3d about a point O can be stated as
% IO = ∫(Ω) ρ ((r · r)E3 − r ⊗ r) dV (5.32)
% where E3 is the 3x3 identity. The symbol is changed here to not
% conflict with I for the inertia. The ‘o times’
% is called the outer-product and we can compute it as
% r ⊗ r → r ⊗ r = rr^T (5.33)
% where the result is a 3x3 matrix. The vector r in the integrand is
% the position of the point in the body relative to the point O.

% 3x3 is the sigma x, sigma y, and tau xy for all 3 planes (xy, xz, yz)


% 1. Using the discussion of §5.7 recast the mass moment of inertia in
% terms of finite element interpolations and a change of variables so we
% can use Gauss-Legendre quadrature.

% you have to replace dV = dxdy with dxideta
% convert all points to gl space
% create pointer vector based on new coordinates




% 2. Build a function for a 4-node quadrilateral and returns the mass
% moment of inertia as a 3x3 matrix. The function inputs are the same as
% in §5.7 plus the location of reference point O.


function [MMOI] = MMOI3d(x, y, rho, int, O)
    % inputs:
    % x,y: positions of points in xy plane.
    % rho: function to calculate height / z-value of point
    % int: number of Gauss-Legendre points in each direction
    % O: point we take the moment of inertia with respect to

    MMOI = zeros(3,3); % initialize variable
    E3 = eye(3); % identity matrix

    % find Gauss-Legendre points + weights
    [xz, wi] = lgwt(int,-1,1);
    % number of points
    % weight of points

    xi = xz;
    eta = xz;

    for m = 1:int
        for n = 1:int

            % create shape functions
            [N, dN_dxi, dN_deta] = shapefunctions(xi(m), eta(n));

            % create jacobian
            jac = [dN_dxi(1), dN_dxi(2); dN_deta(1), dN_deta(2)] * [x(1), x(2); y(1), y(2)];
        
            % det jacobian
            detJ = det(jac);

            % creates pointer to current point in the Gauss-Legendre space

            xcurr = N(1) * x(1) + N(2) * x(2) + N(3) * x(3) + N(4) * x(4);
            ycurr = N(1) * y(1) + N(2) * y(2) + N(3) * y(3) + N(4) * y(4);

            r = [xcurr,ycurr] - O; 

            % do integration
            for i = 1:3
                for j = 1:3
        
                % MMOI = integral over domain of:
                % rho * (dot(r,r) * E3 - r * transpose(r))
                % + prev MMOI val
                
                MMOIf = @(i,j) rho(i,j) * (dot(r, r) * E3(i,j) - (r * r')) * detJ * wi(i) * wi(j);
                MMOI(i,j) = MMOI(i,j) + MMOIf(i,j);

                end

            end

        end

    end
   
end



% 3. Compute the inertia with respect to the center of mass IG using the
% values from §5.7.2. Show how your results converge as the number of
% integration points is increased. What is the needed number of
% integrations points for this problem?


% Table 5.3: Node locations for the 4-node quadrilateral.
% Node x y
% 1 −0.5 −0.3
% 2 1.0 −0.1
% 3 0.7 1.5
% 4 −0.8 0.9

[MassOld, CenterOfMassOld] = COM(xvals, yvals, rho, N);

% define values
xvals = [-0.5 1.0 0.7 -0.8];
yvals = [-0.3 -0.1 1.5 0.9];

% define density function
rho = @(x,y) (3 .* x + 7) .* (sin(3 * y) .* sin(3 * x) + 1);

% Number points
N = 5;

[MMOI] = MMOI3d(xvals, yvals, rho, N, CenterOfMassOld)


% 4. Compute the principle values of the inertia and their directions.
% Visualize their directions and their relative values on the element.

% eigenvalues?


[K, lambda] = eig(MMOI)

figure
hold on
grid on

% plot the contour background
h = fcontour(rho,[min(xvals) max(xvals) min(yvals) max(yvals)]);
scatter([xvals CenterOfMassOld(1)], [yvals CenterOfMassOld(2)])

% extract eigenvectors and make 2D
K1 = [K(1,1) K(2,1)];
K2 = [K(1,2) K(2,2)];

% plot eigenvectors from center of mass
quiver(CenterOfMassOld(1), CenterOfMassOld(2), K1(1), K1(2))
quiver(CenterOfMassOld(1), CenterOfMassOld(2), K2(1), K2(2))


%% lgwt.m

function [x,w]=lgwt(N,a,b)

% lgwt.m
%
% This script is for computing definite integrals using Legendre-Gauss 
% Quadrature. Computes the Legendre-Gauss nodes and weights  on an interval
% [a,b] with truncation order N
%
% Suppose you have a continuous function f(x) which is defined on [a,b]
% which you can evaluate at any x in [a,b]. Simply evaluate it at all of
% the values contained in the x vector to obtain a vector f. Then compute
% the definite integral using sum(f.*w);
%
% Written by Greg von Winckel - 02/25/2004
N=N-1;
N1=N+1; N2=N+2;

xu=linspace(-1,1,N1)';

% Initial guess
y=cos((2*(0:N)'+1)*pi/(2*N+2))+(0.27/N1)*sin(pi*xu*N/N2);

% Legendre-Gauss Vandermonde Matrix
L=zeros(N1,N2);

% Derivative of LGVM
Lp=zeros(N1,N2);

% Compute the zeros of the N+1 Legendre Polynomial
% using the recursion relation and the Newton-Raphson method

y0=2;

% Iterate until new points are uniformly within epsilon of old points
while max(abs(y-y0))>eps
    
    L(:,1)=1;
    Lp(:,1)=0;
    
    L(:,2)=y;
    Lp(:,2)=1;
    
    for k=2:N1
        L(:,k+1)=( (2*k-1)*y.*L(:,k)-(k-1)*L(:,k-1) )/k;
    end
 
    Lp=(N2)*( L(:,N1)-y.*L(:,N2) )./(1-y.^2);   
    
    y0=y;
    y=y0-L(:,N2)./Lp;
    
end

% Linear map from[-1,1] to [a,b]
x=(a*(1-y)+b*(1+y))/2;      

% Compute the weights
w=(b-a)./((1-y.^2).*Lp.^2)*(N2/N1)^2;

end




%% shapefunctions.m

function [N, dN_dxi, dN_deta] = shapefunctions(xi, eta)

    % Shape functions for 4-node quadrilateral element (Lagrange interpolation)
    N = 1/4 * [(1 - xi) * (1 - eta); (1 + xi) * (1 - eta); (1 + xi) * (1 + eta); (1 - xi) * (1 + eta)];
    dN_dxi = 1/4 * [-(1 - eta); (1 - eta); (1 + eta); -(1 + eta)];
    dN_deta = 1/4 * [-(1 - xi); -(1 + xi); (1 + xi); (1 - xi)];

end