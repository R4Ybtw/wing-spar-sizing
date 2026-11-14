function [V, M] = cantilever_statics(x, q)
%CANTILEVER_STATICS  Shear force and bending moment along a cantilever from equilibrium.
%
%   [V, M] = cantilever_statics(x, q)
%   x  stations from root (0) to tip (L), increasing
%   q  distributed load at those stations (N/m)
%
%   A cantilever with a free tip is statically determinate: cut it at x, and
%   the part outboard of the cut must be in equilibrium on its own, so
%       V(x) = integral from x to L of q(s) ds
%       M(x) = integral from x to L of q(s) (s - x) ds
%   No stiffness needed, so these are exact up to the integration rule. Stress
%   therefore comes from statics; the FEM is needed for DEFLECTION (and for
%   structures that aren't determinate, such as a strut-braced wing).

    x = x(:); q = q(:);
    n = numel(x);
    V = zeros(n, 1); M = zeros(n, 1);
    for k = 1:n-1
        s = x(k:end); w = q(k:end);
        V(k) = trapz(s, w);
        M(k) = trapz(s, w.*(s - x(k)));
    end
end
