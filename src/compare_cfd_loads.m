%% compare_cfd_loads.m — spanwise lift: lifting line against Flow Simulation strips
%
% Input data/results/flowsim_strips.csv, one row per strip of the HALF wing
% (root strip first), forces straight from the strip goals:
%   strip,z_start_mm,z_end_mm,Fx,Fy
%   1,0,120,0.31,9.85
% ALPHA_DEG is the angle of the Flow Simulation run (lift is perpendicular to the inflow).
%
% What to compare: the SHAPE (lift per unit span divided by the half-wing
% total), because the structure is sized with the total fixed at n W / 2.
% Compare the total too: it tests the lift slope.

clear; clc;
ALPHA_DEG = 4.1;                 % deg, the n = 4 design point
here = fileparts(mfilename('fullpath')); addpath(here);
ac = load_aircraft();
S = dlmread(fullfile(here, '..', 'data', 'results', 'flowsim_strips.csv'), ',', 1, 0);
z0 = S(:,2)/1000; z1 = S(:,3)/1000; w = z1 - z0; zc = 0.5*(z0 + z1);
a = ALPHA_DEG*pi/180;
Lstrip = S(:,5)*cos(a) - S(:,4)*sin(a);           % lift per strip, N
lp = Lstrip./w;                                  % lift per unit span, N/m
shape_cfd = lp/sum(Lstrip);

cfun = @(y) ac.chord*ones(size(y));
LL = lifting_line(ac.span, cfun, ac.a0_per_rad, ac.alpha_L0_deg*pi/180, a, 25);
CL_cfd = 2*sum(Lstrip)/(ac.q*ac.S);              % half model -> whole wing
fprintf('CL at %.1f deg: lifting line %.3f, CFD %.3f (%+.1f%%)\n', ALPHA_DEG, LL.CL, CL_cfd, ...
        100*(CL_cfd - LL.CL)/LL.CL);
fprintf('(q = %.1f Pa uses V_design and rho from aircraft.json — run the CFD at those values)\n', ac.q);

% Root bending moment for the same total lift, from each shape
M_cfd = sum(Lstrip.*zc)/sum(Lstrip);             % centre of lift, m from root
M_ll  = trapz(LL.y, LL.y.*LL.lift_shape);
fprintf('Spanwise centre of lift: lifting line %.3f m, CFD %.3f m (%+.1f%%) -> root moment changes by the same %%\n', ...
        M_ll, M_cfd, 100*(M_cfd - M_ll)/M_ll);

figure; hold on; grid on;
plot(LL.y, LL.lift_shape, 'k-', 'LineWidth', 1.5, 'DisplayName', 'Lifting line');
stairs([z0; z1(end)], [shape_cfd; shape_cfd(end)], 'r-', 'LineWidth', 1.5, 'DisplayName', 'Flow Sim strips');
el = sqrt(max(0, 1 - (LL.y/ac.semispan).^2))*4/(pi*ac.semispan);
plot(LL.y, el, 'b--', 'DisplayName', 'Elliptic');
xlabel('distance from root [m]'); ylabel('lift per span / half-wing lift [1/m]');
legend('Location', 'southwest');
