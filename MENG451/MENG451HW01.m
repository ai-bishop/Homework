%% Problems:
% 1.2 Moody Diagram
% 1.3.2 Sine
% 1.3.3 Pow


%% Problem 1.2: Moody Diagram
% Using fixed point iteration to solve for fD, Construct the Moody Diagram
% for the turbulent region (4*10^3 to 10^8)
% define terms
    % f_D = Darcy Friction Factor
    % eD = Relative Roughness
    % Re = Reynolds Number
% moody diagram has fD yaxis Re xaxis
close all
clear
clc
format long

% create Reynolds Number Range
Re = logspace(log10(4e3), log10(1e8), 100);  % Reynolds number range

% create relative roughness eD values
eD = [0.0001 0.001 0.01 0.05]; % to give multiple lines on moody diagram

% create Friction Factor Values
f_D = logspace(log10(0.01), log10(0.1), 100);

% calc friction factor
friction_factor = @(Re, eD) (1 ./ (-2 * log10((eD ./ 3.7) + (2.51 ./ (Re .* sqrt(f_D)))))).^2;

% plotting
figure 
hold on
grid on

% loop through each eD and plot
for i = 1:length(eD)
    f_values = friction_factor(Re, eD(i));  % compute f_D for each Re and eD
    plot(Re, f_values, 'DisplayName', ['eD = ' num2str(eD(i))]);
end

% configure plot
title('Moody Diagram')
xlabel('Reynolds Number (Re)')
ylabel('Friction Factor (f_D)')
set(gca, 'XScale', 'log', 'YScale', 'log');


%% Problem 1.3.2 Taylor Series Sine:
% compares sinx to arbitrary relative error
% compare time used for range of inputs
% take advantage of cyclical nature of sine to shift x and make calc easier
format long
% Taylor Series of Sin(x):
    % = x^1/1! - x^3/3! + x^5/5!+ ...
    % (-1)^(n-1) * x^n/n!
close all
clear
clc
format long

% run tests - do 3

[mysin,iter] = mySine(3*pi/2)

% test #1: x = 1
tic;
[mysin1, iter1] = mySine(1) % x = 1
test1 = toc

% funct 1: x = 1
tic
sin1 = sin(1)
sin1time = toc

% test #2: x = 5
tic
[mysin2, iter2] = mySine(5) % x = 5
test2 = toc

% funct 2: x= 5
tic
sin2 = sin(5)
sin2time = toc

% test #3: x = -3
tic
[mysin3, iter3] = mySine(-3) % x = -3
test3 = toc

% funct 3: x = -3
tic
sin3 = sin(-3)
sin3time = toc


function [mysin, iter] = mySine(x1)
    % evaluate x and see if outside of 2pi
    x = abs(x1); 
    x = mod(x, 2*pi); % reduce below 2pi

    % initialize values
    tol = 1e-8; % tolerance
    iter = 0; % number of iterations
    maxiter = 1000; % max number of iterations
    check = 0; % create check
    mysin = x; % initalize mysin w first sine value
    newterm = x; % create newterm function
    
    while check == 0
        iter = iter + 1;

        % update mysine
        newterm = -1 * newterm * x * x / ((2*iter) * (2*iter + 1)); % increase power of x by 2 and factorial initializer by 2
        mysin = mysin + newterm;

        % calculate error
        err = abs(newterm);

        if iter >= maxiter
            check = 1;
            error("Failed to converge in %d steps", iter)
        elseif err < tol
            check = 2;
        end
    end

    if x1 < 0
        
        mysin = -mysin;

    end

end


%% Problem 1.3.3 Taylor Series Power
% calc x^y using simple arithmatic
% use rel error
% compare time for range of inputs
% compare y=0.5 to square root function
% x^y = e^(yln(x)) and can decompose lnx into taylor series and then use
% myexp
% taylor series for lnx is (-1)^(n-1) * ((x-1)^n)/n
close all
clear
clc
format long

% tests

tic
[pow1 iterlnx1 iterex1] = MyPow(3, 2)
toc

tic
[pow2 iterlnx2 iterex2] = MyPow(4, 0.5)
toc

tic
[pow3 iterlnx3 iterex3] = MyPow(9, 0.5)
toc

tic
[pow4 iterlnx4 iterex4] = MyPow(0, 2)
toc

tic
[pow5 iterlnx5 iterex5] = MyPow(2, 0)
toc

function [mypower, iterlnx, iterex] = MyPow(x1,y)
    
    checklnx = 0;
    checkex = 0;
    
    % error checks - special cases

    if y == 0

        checklnx = 1;
        checkex = 1;
        iterlnx = 0;
        iterex = 0;
        mypower = 1;
        return
    
    end

    if x1 == 0
    
        checklnx = 1;
        checkex = 1;
        iterlnx = 0;
        iterex = 0;
        mypower = 0;
        return

    end

    if x1 == 1

        checklnx = 1;
        checkex = 1;
        iterlnx = 0;
        iterex = 0;
        mypower = 1;
        return
    
    end


    % check if x1 neg
    if x1 < 0

        x = abs(x1); % define x
        negval = 0; % reference later

    else 

        x = x1;
        negval = 1;

    end

    % initialize values
    tollnx = 1e-12;
    tolex = 1e-8;
    iterlnx = 0;
    iterex = 0;
    maxiter = 1000;
    mylnx = 0;
    newtermlnx = -1; % to account for (-1)^(n-1)
    newtermex = 1;
    myex = 1;

    % begin decomposition calculation w lnx
    while checklnx == 0
        
        if x > 1

            xex = 1/x;
            flaginverse = 0;

        else

            xex = x;
            flaginverse = 1;

        end

        iterlnx = iterlnx + 1;

        % taylor series for lnx
        newtermlnx = newtermlnx * (-1) * (xex-1);
        mylnx = mylnx + newtermlnx/iterlnx;

        % calc error
        err = abs(newtermlnx / iterlnx);

        % testcheck = [iterlnx mylnx] % was inserted as tester when was
        % throwing errors. can ignore

        if iterlnx >= maxiter

            checklnx = 1;
            error("Ln(x) failed to converge in %d steps", maxiter)

        elseif err < tollnx

            checklnx = 2;

        end
        
    end

    if flaginverse == 0

        mylnx =  -mylnx;

    end
    

    % do y * lnx to get the power that e goes to
    epower = y * mylnx;

    while checkex == 0

    iterex = iterex + 1;

    % update ex
    newtermex = newtermex * epower / iterex;
    myex = myex + newtermex;

    % calc err
    err = abs(newtermex);

        if iterex >= maxiter

            checkex = 1;
            error("e^x failed to converge in %d steps", maxiter)

        elseif err < tolex

            checkex = 2;

        end

    end

% put mypower as answer
mypower = myex;

    % account for negative x value
    if negval == 0

        mypower = -mypower;
    
    end 

end

