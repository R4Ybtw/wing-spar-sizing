function [w, theta, x] = beam_fem(L, EI, q, ne, P_tip)
%BEAM_FEM  Euler-Bernoulli cantilever by the finite-element method.
%
%   [w, theta, x] = beam_fem(L, EI, q, ne)
%   [w, theta, x] = beam_fem(L, EI, q, ne, P_tip)
%
%   L      length (m), clamped at x = 0, free at x = L
%   EI     bending stiffness: a scalar, a vector (one value per element), or a
%          function handle EI(x) evaluated at each element's midpoint
%   q      distributed load as a function handle q(x) in N/m (upward positive)
%   ne     number of elements
%   P_tip  optional point load at the tip (N)
%   Returns deflection w and slope theta at the ne+1 nodes x.
%
%   Method
%   * Each node has two unknowns: deflection w and slope theta. Within an
%     element of length h, w(s) is a cubic built from four Hermite shape
%     functions (s = 0..1 along the element):
%         N1 = 1 - 3s^2 + 2s^3      N2 = h (s - 2s^2 + s^3)
%         N3 = 3s^2 - 2s^3          N4 = h (-s^2 + s^3)
%   * Element stiffness: k_ij = integral of EI N_i'' N_j'' dx over the element.
%   * Element load vector: f_i = integral of q(x) N_i(x) dx. Computed by 6-point
%     Gauss-Legendre quadrature. A smooth load's "consistent" nodal moments
%     are not zero.
%   * Assemble: element e couples global unknowns 2e-1, 2e, 2e+1, 2e+2.
%   * Boundary condition: w = theta = 0 at the root, so delete unknowns 1 and 2
%     before solving K u = F.
%
%   Verified by tests/test_wing.m: uniform load (exact with one element), tip
%   point load, elliptic load against a virtual-work integral, and second-order
%   convergence on a tapered beam.

    if nargin < 5, P_tip = 0; end
    nn = ne + 1;
    x  = linspace(0, L, nn)';
    if isa(EI, 'function_handle'), EIe = EI(0.5*(x(1:end-1) + x(2:end)));
    elseif isscalar(EI),           EIe = EI*ones(ne, 1);
    else,                          EIe = EI(:);
    end
    K = zeros(2*nn);
    F = zeros(2*nn, 1);
    [gp, gw] = gauss_legendre(6);        % points on [-1, 1] and their weights

    for e = 1:ne
        a = x(e); h = x(e+1) - a;

        % ---- TODO 1: element stiffness matrix k (4x4) for stiffness EIe(e) ----
        k = [];
        error('beam_fem:todo', 'beam_fem.m TODO 1: write the element stiffness matrix.');

        % ---- TODO 2: consistent load vector f (4x1) by Gauss quadrature ----
        %   map gp(g) from [-1,1] to s in [0,1]: s = 0.5*(gp(g) + 1), dx = 0.5*h*ds-weight
        f = zeros(4, 1);

        % ---- TODO 3: add k and f into K and F at the right global indices ----

    end

    F(2*nn - 1) = F(2*nn - 1) + P_tip;   % tip point load acts on the last deflection unknown
    free = 3:2*nn;                       % root clamped: unknowns 1 (w) and 2 (theta) removed
    u = zeros(2*nn, 1);
    u(free) = K(free, free)\F(free);
    w = u(1:2:end);
    theta = u(2:2:end);
end

function [x, w] = gauss_legendre(n)
    % Gauss-Legendre points and weights on [-1, 1] (Golub-Welsch method)
    b = (1:n-1)./sqrt(4*(1:n-1).^2 - 1);
    [V, Dg] = eig(diag(b, 1) + diag(b, -1));
    x = diag(Dg); [x, i] = sort(x);
    w = 2*V(1, i)'.^2;
end
