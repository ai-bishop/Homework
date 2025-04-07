%% Problem 2.3.2: Matrix Multiplication
% Implement a naïve (simple loop-based) function that multiplies two rectangular matrices.
% Check the performance against the built-in matrix multiply on a range of random matrices of increasing size
close all
clear
clc

% determine number of runs
n = 5;

for a = 1:n

    % create random bounds
    bound1 = randi([2, 6]);
    bound2 = randi([3, 7]);
    bound3 = randi([1, 9]);

    % generate A & B with multiplicable bounds
    A = rand(bound1, bound2);
    B = rand(bound2, bound3);

    % time created function
    tic
    MatrixMultiplication(A, B);
    test(a) = toc;

    % time innate function
    tic
    A * B;
    funct(a) = toc;

end

test
funct

function C = MatrixMultiplication(A, B)
    
    % extract bounds from input matrices
    [n m] = size(A); % A is n x m
    [m2 p] = size(B); % B is m2 x p

    % check for invalid dimensions
    if m ~= m2

        error('Matrix Dimensions Invalid For Matrix Multiplication')

        return
    
    end

    % Initialize matrix C w zeros
    C = zeros(n, p);

    % sum each step
    % rows
    for i = 1:n

        % columns
        for j = i:p
        
            sum = 0;

            % multiplying each location
            for k = 1:m
            
                sum = sum + A(i, k) * B(k, j);

            end
    
            % inputting into matrix C
            C(i,j) = sum;

        end

    end

end