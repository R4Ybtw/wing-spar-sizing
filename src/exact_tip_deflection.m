function w_tip = exact_tip_deflection(L, EI_fun, q_fun)
%EXACT_TIP_DEFLECTION  Tip deflection of a cantilever by the unit-load method (test oracle).
%
%   w_tip = exact_tip_deflection(L, EI_fun, q_fun)
%
%   Virtual work: put a unit load at the tip. Its bending moment at x is (L - x).
%   Then w_tip = integral from 0 to L of M(x) (L - x) / EI(x) dx, with M(x)
%   from statics. Completely independent of the finite-element method, which
%   is what makes it a good check on it.
%   Computed with adaptive quadrature, so it is accurate to ~1e-10 for smooth
%   loads; the elliptic load's infinite slope at the tip is handled fine.

    Mfun = @(x) arrayfun(@(xx) integral(@(s) q_fun(s).*(s - xx), xx, L, 'RelTol', 1e-12), x);
    w_tip = integral(@(x) Mfun(x).*(L - x)./EI_fun(x), 0, L, 'RelTol', 1e-10);
end
