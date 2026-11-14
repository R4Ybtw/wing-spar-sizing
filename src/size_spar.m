%% size_spar.m — lightest aluminium tube spar that meets every requirement
%
% Requirements (each one is a separate column in the printout, so you can see
% which one decides the answer):
%   1. Strength at ULTIMATE load: root bending stress <= Ftu       (no failure)
%   2. Strength at LIMIT load:    root bending stress <= Fty       (no permanent bend)
%   3. Stiffness at LIMIT load:   tip deflection <= 5% of the semispan (beam FEM)
%   4. Local buckling of the thin wall: stress <= knockdown x 0.605 E t/R
%      (classical cylinder buckling stress, cut to 30% because real tubes
%      buckle far below the classical value — crude but conservative)
%   5. Fit: tube diameter <= airfoil depth at the spar minus skin clearance
%
% Needs beam_fem.m to be written (see its TODO blocks).

clear; clc;
here = fileparts(mfilename('fullpath')); addpath(here);
ac = load_aircraft();
Lw = wing_loads(ac);
Le = wing_loads(ac, 'elliptic');
mat = ac.material;

fprintf('%s: W = %.2f N, S = %.3f m^2, AR = %.1f\n', ac.name, ac.W, ac.S, ac.AR);
fprintf('Stall %.2f m/s, maneuvering speed %.2f m/s; design point %.0f m/s needs CL = %.3f at n = %g\n', ...
        ac.V_stall, ac.V_A, ac.V_design, ac.CL_req, ac.n_limit);
fprintf('Lifting line: alpha = %.2f deg, wing lift slope %.3f /rad, e = %.3f\n', Lw.alpha_deg, Lw.CL_alpha, Lw.e);
fprintf('Root bending moment (ultimate): lifting line %.2f N m, elliptic %.2f N m (%+.1f%%)\n', ...
        Lw.M_root_ult, Le.M_root_ult, 100*(Lw.M_root_ult - Le.M_root_ult)/Le.M_root_ult);
fprintf('Spar depth available %.1f mm -> largest tube %.1f mm\n\n', 1000*ac.spar_depth, 1000*ac.D_max);

Ds = (12:1:30)/1000;
ts = [0.5 0.7 0.8 1.0 1.2 1.5 2.0]/1000;
qfun = @(xx) interp1(Lw.x, Lw.q_lim, xx, 'pchip');
res = [];
for D = Ds
    for t = ts
        if t >= D/2, continue; end
        T = tube_props(D, t);
        s_ult = Lw.M_root_ult*T.c/T.I;
        s_lim = Lw.M_root_lim*T.c/T.I;
        w = beam_fem(ac.semispan, mat.E*T.I, qfun, 40);
        tip = w(end);
        s_cr = ac.buckling_knockdown*0.605*mat.E*t/(D/2);
        mass = 2*T.A*ac.semispan*mat.density;              % both half-spars
        ok = [s_ult <= mat.Ftu, s_lim <= mat.Fty, tip <= ac.tip_deflection_limit*ac.semispan, ...
              s_ult <= s_cr, D <= ac.D_max];
        res(end+1, :) = [D t mass s_ult s_lim tip s_cr ok]; %#ok<SAGROW>
    end
end
feas = all(res(:, 8:12), 2);
[~, order] = sort(res(:,3));
fprintf('Lightest designs overall, then the lightest that pass (stress MPa, deflection mm, mass g; 1 = pass):\n');
fprintf('%6s %5s %7s %8s %8s %7s %8s | %s\n', 'D', 't', 'mass', 'sig_ult', 'sig_lim', 'tip', 'sig_cr', 'ult lim stiff buck fit');
fo = order(feas(order));
show = [order(1:min(6, end)); fo(1:min(6, end))];
for k = show'
    r = res(k, :);
    tag = ''; if all(r(8:12)), tag = '<- feasible'; end
    fprintf('%6.1f %5.2f %7.0f %8.0f %8.0f %7.1f %8.0f | %3d %3d %4d %4d %3d %s\n', 1000*r(1), 1000*r(2), ...
            1000*r(3), r(4)/1e6, r(5)/1e6, 1000*r(6), r(7)/1e6, r(8:12), tag);
end
best = order(find(feas(order), 1));
if isempty(best)
    fprintf('\nNo tube in the grid meets every requirement — rethink the airfoil depth or the load case.\n');
else
    b = res(best, :);
    fprintf('\nLightest feasible: %.0f x %.2f mm, %.0f g for both spars, tip deflection %.1f mm at limit load\n', ...
            1000*b(1), 1000*b(2), 1000*b(3), 1000*b(6));
    names = {'ultimate strength', 'yield', 'stiffness', 'buckling', 'fit'};
    margins = [mat.Ftu/b(4), mat.Fty/b(5), ac.tip_deflection_limit*ac.semispan/b(6), b(7)/b(4), ac.D_max/b(1)] - 1;
    [~, i] = min(margins);
    fprintf('Closest constraint: %s (margin %.0f%%)\n', names{i}, 100*margins(i));
end

figure; hold on; grid on;
plot(1000*res(~feas,1), 1000*res(~feas,3), 'x', 'Color', [0.7 0.7 0.7]);
plot(1000*res(feas,1), 1000*res(feas,3), 'o');
if ~isempty(best), plot(1000*res(best,1), 1000*res(best,3), 'p', 'MarkerSize', 14, 'LineWidth', 2); end
xlabel('tube diameter [mm]'); ylabel('mass of both spars [g]');
title('Spar design space (o feasible, x fails a requirement)');

