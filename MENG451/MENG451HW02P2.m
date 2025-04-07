%% HW03 P2.7
close all
clear
clc

%% 2.7.1


%% 2.7.2


%% 2.7.3
sx = [(-2 -2i), -2 + 2i, -4 - 4i, -4 + 4i];
dvals = poly(sx);


%% 2.7.4
% a
B = transpose([0 0 1 0]);

A = [0 0 1 0; 0 0 0 1; -30 10 -0.6 0.2; 10 -10 0.2 -0.2];

AB = A*B;

A2B = A*AB;

A3B = A*A2B;
 
Cm = [B AB A2B A3B];

det(Cm);

% b
[Q,R] = qr(Cm);

% c
[U,S,V] = svd(Cm);
Sdiag = diag(S);
nnz(Sdiag);

% d

lambda = eig(Cm);

% e

rank(Cm);


%% 2.7.5

C = [0 1 0 0];
CA = C*A;
CA2 = CA*A;
CA3 = CA2 * A;
Om = [C; CA; CA2; CA3];
rank(Om);
det(Om);

%% 2.7.6

lambda2 = eig(A);
coeffs2 = poly(lambda2);
Atilde = [0 1 0 0; 0 0 1 0; 0 0 0 1; -200 -8 -40.08 -0.8];
Btilde = transpose([0 0 0 1]);
AtildeBtilde =  Atilde * Btilde;
Atilde2Btilde = Atilde * AtildeBtilde;
Atilde3Btilde = Atilde * Atilde2Btilde;
Cz = [Btilde AtildeBtilde Atilde2Btilde Atilde3Btilde];
P = Cm * inv(Cz);


avals = poly(eig(A));
tempkz = fliplr(dvals - avals);
Kz = tempkz(1:4);
Kx = -Kz * inv(P);

%% 2.7.8

lambda3 = eig(transpose(A));
avals2 = poly(lambda3);
ct = transpose(C);
atc = Atilde * ct;
at2c = Atilde * atc;
at3c = Atilde * at2c;

Cl = [ct atc at2c at3c];

newroots = 10*sx;
dvals2 = poly(newroots);

tempkz2 = fliplr(dvals2 - avals);
Kz2 = tempkz2(1:4);
Kx2 = -Kz2 * inv(P);
L = transpose(Kx2);


%% 2.7.9

P2 = A*A*A*A + dvals(2) * A*A*A + dvals(3)* A*A + dvals(4) * A + dvals(5);
K2 = [0 0 0 1] * inv(Cm) * P2;



%% 2.7.10



u = @(t,x) -Kx*x;

dx = @(t,x) A*x + B*u(t,x);

x0 = [0; 0; 0; 1];

sol = ode45(dx, [0,1],x0);

figure
plot(sol.x, sol.y(1:2,:), "LineWidth",3)










% u = @(t,x) -Kx*x;
% 
% dx = @(t,x) A*x + B*u(t,x); 
% 
% x0 = [0; 0; 0; 0];
% 
% 
% 
% y = [0, 0, 0, 1, 0, 0, 0, 0];
% 
% 
% x = y(1:4);
% xhat = y(5:end);
% u1 = -Kx * transpose(xhat);
% dx = A*transpose(x) + B*u1;
% dxhat = A*transpose(xhat) + B*u1 + L*C*transpose(x-xhat);
% dw = [dx; dxhat];
% 
% dw1 = transpose(dw)
% 
% sol = ode45(dw1, [0,1], y);
% 
% figure
% plot(sol.x, sol.y(1:2,:), "LineWidth",3)
% 
% 
% Ox = [ C; 
%        C*A;
%        C*A^2;
%        C*A^3];
% rank(Ox)
% 
% 


%%
% function dw = mysys(w)
% 
%     x = w(1:4);
%     xhat = w(5:end);
%     u1 = -Kx * xhat;
%     dx = A*x + B*u1;
%     dxhat = A*xhat + B*u1 + L*C*(x-xhat);
%     dw = [dx; dxhat];
% 
% end
% 
% 
% 
