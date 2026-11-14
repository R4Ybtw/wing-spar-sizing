function L = wing_loads(ac, model)
%WING_LOADS  Spanwise load on one half-wing at limit and ultimate load factor.
%
%   L = wing_loads(ac)                     lifting line (default)
%   L = wing_loads(ac, 'elliptic')         Schrenk's classic textbook shortcut
%   L = wing_loads(ac, 'uniform')          crude upper bound on root moment
%
%   L.x                stations from root (0) to tip (semispan), m
%   L.q_lim, L.q_ult   lift per unit span at limit / ultimate load (N/m)
%   L.V_*, L.M_*       shear force and bending moment from statics
%   L.M_root_lim/ult   root bending moment
%   L.alpha_deg        root angle of attack that gives CL_req (lifting line only)
%
%   The load factor fixes the TOTAL lift: n W / 2 per half-wing. Aerodynamics
%   decides only its SHAPE along the span, and the shape sets the root bending
%   moment. A rectangular wing carries relatively more lift near the tip than
%   an elliptic one, so its root moment is higher.
%   Not included: the wing's own weight, which relieves the bending
%   (conservative to leave out), and the lift carried across the fuselage.

    if nargin < 2, model = 'lifting_line'; end
    s = ac.semispan;
    x = linspace(0, s, 401)';
    switch model
        case 'lifting_line'
            cfun = @(y) ac.chord*ones(size(y));
            aL0 = ac.alpha_L0_deg*pi/180;
            LL1 = lifting_line(ac.span, cfun, ac.a0_per_rad, aL0, aL0 + 0.1, 25);
            CLa = LL1.CL/0.1;                               % wing lift slope (per rad)
            alpha = aL0 + ac.CL_req/CLa;
            LL = lifting_line(ac.span, cfun, ac.a0_per_rad, aL0, alpha, 25);
            shape = interp1(LL.y, LL.lift_shape, x);
            L.alpha_deg = alpha*180/pi; L.CL_alpha = CLa; L.e = LL.e; L.CL = LL.CL;
        case 'elliptic'
            shape = sqrt(max(0, 1 - (x/s).^2))*4/(pi*s);
        case 'uniform'
            shape = ones(size(x))/s;
        otherwise
            error('model must be lifting_line, elliptic or uniform');
    end
    shape = shape/trapz(x, shape);                          % integrate to exactly 1
    L.x = x; L.model = model; L.shape = shape;
    L.q_lim = ac.L_half_lim*shape;
    L.q_ult = ac.L_half_ult*shape;
    [L.V_lim, L.M_lim] = cantilever_statics(x, L.q_lim);
    [L.V_ult, L.M_ult] = cantilever_statics(x, L.q_ult);
    L.M_root_lim = L.M_lim(1); L.M_root_ult = L.M_ult(1);
    L.centroid = trapz(x, x.*shape);                        % spanwise centre of lift
end
