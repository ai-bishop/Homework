function [N, dN_dxi, dN_deta] = shape_functions2(xi, eta)
    % Shape functions for 4-node quadrilateral element (Lagrange interpolation)
    
    % Define shape functions
    N1 = 1/4 * (1 - xi) * (1 - eta);
    N2 = 1/4 * (1 + xi) * (1 - eta);
    N3 = 1/4 * (1 + xi) * (1 + eta);
    N4 = 1/4 * (1 - xi) * (1 + eta);
    
    % Combine shape functions in a column vector
    N = [N1; N2; N3; N4];
    
    % Derivatives of the shape functions with respect to xi and eta
    dN1_dxi = -1/4 * (1 - eta);
    dN2_dxi =  1/4 * (1 - eta);
    dN3_dxi =  1/4 * (1 + eta);
    dN4_dxi = -1/4 * (1 + eta);
    
    dN1_deta = -1/4 * (1 - xi);
    dN2_deta = -1/4 * (1 + xi);
    dN3_deta =  1/4 * (1 + xi);
    dN4_deta =  1/4 * (1 - xi);
    
    % Combine the derivatives into column vectors
    dN_dxi = [dN1_dxi; dN2_dxi; dN3_dxi; dN4_dxi];
    dN_deta = [dN1_deta; dN2_deta; dN3_deta; dN4_deta];
end