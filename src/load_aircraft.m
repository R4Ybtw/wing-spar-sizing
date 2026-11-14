function ac = load_aircraft(file)
%LOAD_AIRCRAFT  Read data/aircraft.json and derive the load case.
%
%   ac = load_aircraft()
%
%   Derived fields:
%     W, S, AR, semispan              weight (N), wing area, aspect ratio, half span
%     V_stall, V_A                    1 g stall speed and maneuvering speed
%     n_ult                           ultimate load factor = n_limit x factor of safety
%     L_half_lim, L_half_ult          lift on ONE half-wing at limit and ultimate load
%     CL_req                          lift coefficient needed for n_limit at V_design
%     spar_depth, D_max               airfoil thickness at the spar, and the largest
%                                     tube that fits inside with skin clearance
%
%   Why these numbers: the structure must carry the LIMIT load (the most the
%   aircraft should ever see) with no permanent deformation, and the ULTIMATE
%   load (1.5 x limit) without breaking. The load factor fixes the TOTAL lift;
%   aerodynamics only decides how it is spread along the span.

    if nargin < 1
        file = fullfile(fileparts(mfilename('fullpath')), '..', 'data', 'aircraft.json');
    end
    ac = jsondecode(fileread(file));
    g = 9.80665;
    ac.W = ac.mass*g;
    ac.S = ac.span*ac.chord;
    ac.AR = ac.span^2/ac.S;
    ac.semispan = ac.span/2;
    ac.V_stall = sqrt(2*ac.W/(ac.rho*ac.S*ac.CL_max));
    ac.V_A = ac.V_stall*sqrt(ac.n_limit);
    ac.n_ult = ac.n_limit*ac.factor_of_safety;
    ac.L_half_lim = ac.n_limit*ac.W/2;
    ac.L_half_ult = ac.n_ult*ac.W/2;
    q = 0.5*ac.rho*ac.V_design^2;
    ac.q = q;
    ac.CL_req = ac.n_limit*ac.W/(q*ac.S);

    % NACA 4-digit thickness at the spar station (thickness is independent of camber)
    t = str2double(ac.airfoil(3:4))/100; x = ac.spar_x_over_c;
    yt = 5*t*(0.2969*sqrt(x) - 0.1260*x - 0.3516*x^2 + 0.2843*x^3 - 0.1015*x^4);
    ac.spar_depth = 2*yt*ac.chord;
    ac.D_max = ac.spar_depth - 2*ac.skin_clearance;
end
