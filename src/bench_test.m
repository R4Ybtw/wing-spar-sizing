%% bench_test.m — measured bending stiffness of a real tube, and the stiffness of the clamp
%
% Clamp the tube to a table, hang known masses at a measured distance L from
% the clamp edge, and measure the tip deflection (dial gauge, or a ruler and
% a camera). Record in data/results/bench_test.csv:
%   length_m,mass_kg,deflection_mm
%   0.800,0.0,0.0
%   0.800,0.5,3.1
% Ideally measure at TWO lengths (e.g. 0.6 and 0.9 m). Then the fit separates
% the tube's bending stiffness EI from the clamp's rotational flexibility k:
%   deflection / load  =  L^3/(3 EI)  +  L^2/k
% (a clamp that lets the root rotate slightly adds the second term). With one
% length the clamp's flexibility gets blamed on the tube, and EI comes out low.

clear; clc;
D = 0.025; t = 0.0007;          % m, measured outer diameter and wall of the test tube
E = 68.9e9;                     % Pa, any aluminium alloy (alloy changes strength, not stiffness)
g = 9.80665;

here = fileparts(mfilename('fullpath')); addpath(here);
B = dlmread(fullfile(here, '..', 'data', 'results', 'bench_test.csv'), ',', 1, 0);
Ls = unique(B(:,1));
C = zeros(numel(Ls), 1);
for k = 1:numel(Ls)
    s = B(:,1) == Ls(k);
    P = B(s,2)*g; d = B(s,3)/1000;
    p = polyfit(P, d, 1);                       % slope = compliance (m/N); intercept = zero offset
    C(k) = p(1);
    res = d - polyval(p, P);
    fprintf('L = %.3f m: compliance %.4e m/N, fit residual max %.2f mm\n', Ls(k), C(k), 1000*max(abs(res)));
end
T = tube_props(D, t); EI_nom = E*T.I;
if numel(Ls) == 1
    EI = Ls^3/(3*C);
    fprintf('\nEI measured %.1f N m^2 (one length: includes clamp flexibility)\n', EI);
else
    X = [Ls.^3/3 Ls.^2];                        % unknowns [1/EI ; 1/k]
    sol = X\C;
    EI = 1/sol(1); kc = 1/sol(2);
    fprintf('\nEI measured %.1f N m^2, clamp rotational stiffness %.0f N m/rad\n', EI, kc);
end
fprintf('EI nominal  %.1f N m^2 (E = %.1f GPa, D = %.1f mm, t = %.2f mm)\n', EI_nom, E/1e9, 1000*D, 1000*t);
fprintf('Difference %+.1f%%. Wall thickness enters I almost linearly: measure it carefully.\n', ...
        100*(EI - EI_nom)/EI_nom);

figure; hold on; grid on;
for k = 1:numel(Ls)
    s = B(:,1) == Ls(k);
    plot(B(s,2)*g, B(s,3), 'o-', 'DisplayName', sprintf('L = %.2f m', Ls(k)));
    PP = linspace(0, max(B(s,2))*g, 10);
    plot(PP, 1000*PP*Ls(k)^3/(3*EI_nom), '--', 'DisplayName', 'beam theory, nominal EI');
end
xlabel('tip load [N]'); ylabel('tip deflection [mm]'); legend('Location', 'northwest');
