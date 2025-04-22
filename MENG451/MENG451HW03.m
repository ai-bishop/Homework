%% MENG451HW03
% PROBLEM 3.1
clear all
close all
clc

%% Problem 3.1.1

%% Part 1
% Solve (3.6) using a built-in solver. Such as Matlab’s or Julia’s \
% operator. Plot the results on a contour lot

% create temperature matrix
Temperature = zeros(5,5);

% initialize edge values
% not initializing corners because they are not used in calculation
Temperature(2:4, 1) = 75; % left of matrix
Temperature(2:4, 5) = 50; % right of matrix
Temperature(1, 2:4) = 100; % top of matrix
Temperature(5, 2:4) = 0; % bottom of matrix

% create coefficients matrix for large solve
Coeffs = [4 -1 0, -1 0 0, 0 0 0;
        -1 4 -1, 0 -1 0, 0 0 0;
        0 -1 4, 0 0 -1, 0 0 0;
        -1 0 0, 4 -1 0, -1 0 0;
        0 -1 0, -1 4 -1, 0 -1 0;
        0 0 -1, 0 -1 4, 0 0 -1;
        0 0 0, -1 0 0, 4 -1 0;
        0 0 0, 0 -1 0, -1 4 -1;
        0 0 0, 0 0 -1, 0 -1 4];

% create solutions for large solve
Results = transpose([75 0 50 75 0 50 175 100 150]);

% solve - matrix given is 9x1 
DT = transpose(Coeffs\Results);

% NOTE: temps matrix  DT(9,1) gives result as
%  = [T11 T21 T31 T12 T22 T32 T13 T23 T33] on diagram
% This equates to [T31 T32 T33 T21 T22 T23 T11 T12 T13] in TempCore

TempCore = [DT(7) DT(8) DT(9);
            DT(4) DT(5) DT(6);
            DT(1) DT(2) DT(3)];

Temperature(2:4, 2:4) = TempCore


%% Part 2
% Determine the condition number (Moler, 2017; Wikipedia contributors,
% 2018) based on the ratio of the singular values σ of the system matrix
% you solved in the previous part of the problem.
% κ = σmax(A)/ σmin(A) (3.7)
% Compare κ to the results of a built-in condition number function.

% calc k = singmax/singmin

% Calculate the singular values. They are outputted along the diagonal in S
% S(1) is the maximum, S(5) is the minimum
[U,S,V] = svd(Temperature);

SingMax = S(1);
SingMin = S(5);
Khand = SingMax/SingMin;

% Built in ocndition number function
Kcommand = cond(Temperature);

Khand
Kcommand




%% Part 3
% Implement the unpreconditioned biconjugate gradient stabilized
% method, usually written as BiCGSTAB, of van der Vorst (1992). Solve
% (3.6) using your BiCGSTAB function. During the iterations, keep track
% of the residuals, and make a plot showing how your solver converged.

% A = Coeffs
% B = Results
% x0 = guess

xi = zeros(9,1); % create initial guess, x(i)

% NOTE: the i indices are advanced one as MATLAB starts in a matrix at 1,
% instead of 0

% DO NOT NEED TO CARRY OLD VALUES
r = Results - Coeffs*xi; % initialize r(1,:)
rtilde = 10*rand(9,1); % initialize rtilde
rho = dot(r, rtilde); % initialize rho(1)
p = r; % initialize p(:,1)
i = 0; % initialize i
trigger = 0; % initalize trigger for function
margin1 = 0.01; % toilerance for s
margin2 = 0.01; % tolerance for r
safetyfactor = 1000;
[m,n] = size(Temperature);


% do BiCGSTAB - proceedure from wikipedia
while trigger == 0
        i = i + 1; % safety factor to prevent recursion
        v = Coeffs * p;
        alpa = rho/dot(rtilde, v);
        h = xi + alpa * p;
        s = r - alpa * v;

        if norm(s) < margin1
            trigger = 1;
        else
            t = Coeffs*s;
            omega = dot(t, s)/dot(t, t);
            xi = h + omega*s;
            r = s - omega*t;
        end

        if norm(r) < margin2
            trigger = 1;
        else
            rho1 = dot(rtilde, r);
            beta = (rho1/rho)*(alpa/omega);
            p = r + beta*(p - omega*v);
            rho = rho1; % update rhos
        end

       

        
        % create temps matrix from BiCGSTAB
        DT1 = transpose((Results-r)\Coeffs);

        BiCGSTABTempCore = [DT1(7) DT1(8) DT1(9);
                            DT1(4) DT1(5) DT1(6);
                            DT1(1) DT1(2) DT1(3)];
        
        % initialize edge values
        % not initializing corners because they are not used in calculation
        BiCGSTABTemp(2:4, 1) = 75; % left of matrix
        BiCGSTABTemp(2:4, 5) = 50; % right of matrix
        BiCGSTABTemp(1, 2:4) = 100; % top of matrix
        BiCGSTABTemp(5, 2:4) = 0; % bottom of matrix



        BiCGSTABTemp(2:4, 2:4) = BiCGSTABTempCore;


        % create residual for iteration
        Residual = zeros(m,n); % create residual matrix; reset every iteration; based off BiCGSTABTemp


        for j = 2:(m-1)
            for k = 2:(n-1)
                Residual(j,k) = -4*BiCGSTABTemp(j,k) + BiCGSTABTemp(j+1,k) + BiCGSTABTemp(j-1,k) + BiCGSTABTemp(j,k+1) + Temperature (j,k-1);
            end
        end

        R(i) = norm(Residual); % put the resiudal value into R matrix

        if i > safetyfactor
            break
        end

end



% create plot of residuals
figure
plot(R)
ylabel('Residual')
xlabel('Iteration')




%% Problem 3.1.2

%% Part 1
% Implement a Successive Over Relaxation (SOR) solver for this system
% of n equally-spaced grid points in each direction

tol = 1; % tolerance values
trigger = 0;
iter = 0;
R = zeros(1,1);
relerror = 1000;

sized = 7; % determines size of matrix aka n

% generate P - input
% size of (size-2, size-2)
P = randi(100,[sized-2 sized-2]);

Told = zeros(sized, sized); % initalize size parameter

% initialize edge values
% not initializing corners because they are not used in calculation
Told(2:(sized-1), 1) = 75; % left of matrix
Told(2:(sized-1), sized) = 50; % right of matrix
Told(1, 2:(sized-1)) = 100; % top of matrix
Told(sized, 2:(sized-1)) = 0; % bottom of matrix

Told(2:(sized-1),2:(sized-1)) = P;
Tnew = Told;
Tnew(2:sized-1,2:sized-1) = 0;

omega = 1.00111; % parameter for SOR, can change

while relerror > tol

    iter = iter + 1;

    % perform SOR
    for i=2:(sized-1)
        for j=2:(sized-1)
    
            phi = (1/4) * (Told(i+1,j)+Told(i-1,j)+Told(i,j+1)+Told(i,j-1));
    
            Tnew(i,j) = (1-omega)*Told(i,j) + omega * phi;
    
        end
    end

    % calc error
    relerror = norm((Tnew - Told));

    % calculate residual
    Residual = zeros(sized,sized); % initialize residuals matrix
    
    for i = 2:(sized-1)
        for j = 2:(sized-1)
            Residual(i,j) = 4*Tnew(i,j) - (Tnew(i+1,j) + Tnew(i-1,j) + Tnew(i,j+1) + Tnew(i,j-1));
        end
    end
    
    R(iter) = relerror;

    Told = Tnew;

    if iter > 1000
        break
    end

end

% create plot of residuals
figure
plot(R)
ylabel('Residual')
xlabel('Iteration')





%% Part 2
% Plot the convergence as a function of grid points. Also make a
% contour-type plot of the finest grid.

% ok yes im probably doing a lot of extra work but idc :p


for sized = 7:15
    tic
    iter = 0;

    N = sized - 2;
    
    % generate P - interior values
    % size of NxN
    P = randi(100,[N,N]);
    
    x = zeros(N,N);
    
    % generate corners
    x(1,1) = -4 * P(1,1) + P(2,1) + P(1,2); % NW Corner
    x(1,N) = -4 * P(1,N) + P(2,N) + P(1,N-1); % NE Corner
    x(N,1) = -4 * P(N,1) + P(N,2) + P(N-1,1); % SW Corner
    x(N,N) = -4 * P(N,N) + P(N-1,N) + P(N,N-1); % NW Corner
    
    % generate north edge - i = 1
    for j = 2:N-1
        x(1,j) = -4 * P(1,j) + P(2,j) + P(1,j+1) + P(1,j-1);
    end

    % generate south edge - i = N
    for j = 2:N-1
        x(N,j) = -4 * P(N,j) + P(N-1,j) + P(N,j+1) + P(N,j-1);
    end

    % generate west edge - j = 1
    for i = 2:N-1
        x(i,1) = -4 * P(i,1) + P(i,2) + P(i+1,1) + P(i-1,1);
    end

    % generate east edge - j = N
    for i = 2:N-1
        x(i,N) = -4 * P(i,N) + P(i,N-1) + P(i+1,N) + P(i-1,N);
    end

    % generate interior
    for i = 2:N-1
        for j = 2:N-1

            x(i,j) = -4 * P(i,j) + P(i+1,j) + P(i-1,j) + P(i,j+1) + P(i,j-1);

        end
    end

    xold = zeros(sized,sized);
    xold(2:sized-1, 2:sized-1) = x;


    % initialize edge values
    % not initializing corners because they are not used in calculation
    xold(2:sized-1, 1) = 75; % left of matrix
    xold(2:sized-1, sized) = 50; % right of matrix
    xold(1, 2:sized-1) = 100; % top of matrix
    xold(sized, 2:sized-1) = 0; % bottom of matrix


    tol = 1;
    xnew = xold;
    err = 10;
    omega = 1.01;

    while err > tol
        
        iter = iter + 1;
        
        % perform SOR
        for i=2:sized-1
            for j=2:sized-1
    
                phi = (1/4) * (xold(i+1,j)+xold(i-1,j)+xold(i,j+1)+xold(i,j-1));
    
                xnew(i,j) = (1-omega)*xold(i,j) + omega * phi;
    
            end
        end

        err = norm(xnew-xold);
    
        xold = xnew;

        if iter > 1000
            break
        end

    end
resds(sized-6) = err;
iterations312(sized-6) = iter;
time312(sized-6) = toc;

end

% create plot of iterations/convergence
figure
plot(iterations312)
ylabel('iterations')
xlabel('number of points to a side')


figure
contour(xnew)
xlabel('Temp')
ylabel('Temp')



%% Part 3
% Implement the BiCGSTAB method without building the entire matrix A
% as a single dense array. Hint: think carefully about what operations
% are needed in the method and how to do them in terms of (3.5).
%Compare the following things between SOR and BiCGSTAB:
%(a) Number of iterations vs grid size
%(b) Wall time vs grid size
%(c) For a large grid case, compare how the residuals converge.

% dont create A b/c extremely large and barely filled (~1%)

% "implement" means cant use bicgstab() ^~^

% BiCGSTAB as from WIKIPEDIA
    % v = Api−1
    % α = ρi−1/(r̂0, v)
    % h = xi−1 + αpi−1
    % s = ri−1 − αv
    % If h is accurate enough, i.e., if s is small enough, then set xi = h and quit
    % t = As
    % ω = (t, s)/(t, t)
    % xi = h + ωs
    % ri = s − ωt
    % If xi is accurate enough, i.e., if ri is small enough, then quit
    % ρi = (r̂0, ri)
    % β = (ρi/ρi−1)(α/ω)
    % pi = ri + β(pi−1 − ωv)




% Define grid size and parameters - assuming square

for n = 7:15
    tic
    [Resids(:,n-6), iterations313(n-6)] = uncon_noA_BiCGSTAB(n);
    time313(n-6) = toc;

end

Resids;

n = 7:15;

% a) # iterations v grid size
figure
hold on
plot(iterations313,n,'DisplayName','BiCGSTAB')
plot(iterations312,n,'DisplayName', 'SOR'); 
xlabel('Iterations')
ylabel('Grid size (n * n)')
legend

% b) wall time v grid size
n = 7:15;
figure
hold on
plot(time313,n,'DisplayName','BiCGSTAB');
plot(time312,n,'DisplayName','SOR');
ylabel('Grid size (n * n)')
xlabel('Wall Time')
legend()

% c) for grid = large, how residuals converge?

figure
hold on
plot(Resids,n,'DisplayName','BiCGSTAB');
plot(resds,n,'DisplayName','SOR');
ylabel('Grid size (n * n)')
xlabel('Residuals')
legend()




function [Res, iter] = uncon_noA_BiCGSTAB(n)

    
    tol = 1e-6;  % Convergence tolerance
    
    % Initialize temperature grid and boundary conditions
    b = zeros(n * n, 1);  % Boundary condition vector
    
    % Assign boundary values
    b(1:n) = 100;  % Top boundary
    b(end-(n+1):end) = 0;  % Bottom boundary
    b(1:n:end) = 75;  % Left boundary
    b(n:n:end) = 50;  % Right boundary
    
    
    xtest = randi(100,n);
    x = xtest(:);  % Initial guess for the solution
    r = b - InteriorCalc(x, n);  % Initial residual
    rhat = r;  % Set r = r0
    rho_old = dot(r, rhat);  
    p = r;
    iter = 0; % starts iteration tracker
    trigger = 0;
    maxiter = 1000;
    
    while trigger == 0
        iter = iter + 1;
    
        v2 = InteriorCalc(x, n)\x * p;
        alpha = rho_old\dot(rhat, v2);
        h = x + alpha * rho_old;
        s = r - alpha * v2;
    
        % If h is accurate enough, i.e., if s is small enough, then set xi = h and quit
    
        t2 = InteriorCalc(x, n)\x * s;
        omega = dot(t2, s)/dot(t2, t2);
        xi = h + omega * s;
        ri = s - omega * t2;
    
        % If xi is accurate enough, i.e., if ri is small enough, then quit
    
        rho_new = dot(rhat, ri);
        beta = (rho_new / rho_old) * (alpha / omega);
        
        pi = ri + beta * (p - omega * v2);
    
        % update values for next iteration
        p = pi;
        r = ri;
        x = xi;
        rhat = r;
        rho_old = rho_new;
    
    
        re = norm(rhat);
        Res(iter) = re;
    
    
        if re < tol
            trigger = 1;
        end
    
        if iter >= maxiter
            trigger = 1;
        end
    
    end
end





function T_new = InteriorCalc(T, n)
    % Solve initial inside of matrix
    T_new = zeros(n * n, 1);  % Initialize new temperature array
    for i = 2:n-1
        for j = 2:n-1
            index = sub2ind([n, n], i, j);  % Convert 2D indices to 1D index
            % Apply the finite difference stencil (central difference)
            T_new(index) = (T(index+1) + T(index-1) + T(index+n) + T(index-n) - 4*T(index));
        end
    end
end
