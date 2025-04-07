function [N, dN_dxi, dN_deta] = shapefunctions(xi, eta)

    % Shape functions for 4-node quadrilateral element (Lagrange interpolation)
    N = 1/4 * [(1 - xi) * (1 - eta); (1 + xi) * (1 - eta); (1 + xi) * (1 + eta); (1 - xi) * (1 + eta)];
    dN_dxi = 1/4 * [-(1 - eta); (1 - eta); (1 + eta); -(1 + eta)];
    dN_deta = 1/4 * [-(1 - xi); -(1 + xi); (1 + xi); (1 - xi)];

end