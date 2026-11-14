%% fea_loads_table.m — segment forces for SolidWorks Simulation, and the answer to expect
%
% The half-span is split into n_segments strips (the CAD tube has its surface
% split the same way). Each strip gets the lift that falls on it. Apply each
% force to its strip face in SolidWorks Simulation, in +Y.
%
% Prints the tip deflection the beam FEM predicts for exactly this stepped loading:
% that, not the smooth-load value, is what SolidWorks should reproduce.
%
% Tube under analysis:
D = 0.025; t = 0.0007;     % m  (lightest feasible from size_spar.m with the default inputs)

here = fileparts(mfilename('fullpath')); addpath(here);
ac = load_aircraft(); Lw = wing_loads(ac);
n = ac.n_segments; s = ac.semispan;
edges = linspace(0, s, n + 1)';
F_lim = zeros(n, 1); F_ult = zeros(n, 1);
for k = 1:n
    in = Lw.x >= edges(k) & Lw.x <= edges(k+1);
    F_lim(k) = trapz(Lw.x(in), Lw.q_lim(in));
    F_ult(k) = trapz(Lw.x(in), Lw.q_ult(in));
end
F_lim = F_lim*ac.L_half_lim/sum(F_lim);                     % remove trapezoid round-off
F_ult = F_ult*ac.L_half_ult/sum(F_ult);

out = fullfile(here, '..', 'data', 'results', 'fea_segment_loads.csv');
fid = fopen(out, 'w');
fprintf(fid, 'segment,z_start_mm,z_end_mm,F_limit_N,F_ultimate_N\n');
fprintf(fid, '%d,%.1f,%.1f,%.3f,%.3f\n', [(1:n)' 1000*edges(1:end-1) 1000*edges(2:end) F_lim F_ult]');
fclose(fid);
fprintf('Segment loads written to %s\n', out);
fprintf('%4s %10s %10s %10s %10s\n', 'seg', 'z0 [mm]', 'z1 [mm]', 'F_lim [N]', 'F_ult [N]');
fprintf('%4d %10.1f %10.1f %10.3f %10.3f\n', [(1:n)' 1000*edges(1:end-1) 1000*edges(2:end) F_lim F_ult]');

T = tube_props(D, t); EI = ac.material.E*T.I;
qstep = @(xx) step_load(xx, edges, F_lim);
w_step = beam_fem(s, EI, qstep, 20*n);
qsm = @(xx) interp1(Lw.x, Lw.q_lim, xx, 'pchip');
w_sm = beam_fem(s, EI, qsm, 200);
[~, Mst] = cantilever_statics(Lw.x, Lw.q_ult);
fprintf('\nTube %.0f x %.2f mm, EI = %.1f N m^2\n', 1000*D, 1000*t, EI);
fprintf('Expected tip deflection at LIMIT load: %.2f mm (stepped strips), %.2f mm (smooth load)\n', ...
        1000*w_step(end), 1000*w_sm(end));
fprintf('Expected bending stress at the root, ULTIMATE: %.0f MPa (beam theory).\n', Mst(1)*T.c/T.I/1e6);
M100 = interp1(Lw.x, Mst, 0.100);
fprintf('Expected bending stress 100 mm from the root, ULTIMATE: %.0f MPa  <- compare the FEA here\n', M100*T.c/T.I/1e6);
fprintf('Read FEA stress at least one diameter out from the fixed end: the fixture creates a local peak.\n');
