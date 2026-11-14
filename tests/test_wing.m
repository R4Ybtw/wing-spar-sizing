%% test_wing.m — verification against exact answers
%
% Run from the project folder:   addpath tests; test_wing
% Section 4 tests src/beam_fem.m (reported as NOT IMPLEMENTED while it has open TODOs).

function test_wing()
here = fileparts(mfilename('fullpath'));
addpath(fullfile(here, '..', 'src'));
ok = true;

%% 1. Lifting line: an elliptic wing has the exact answer
%  CL = a0 (alpha - aL0) / (1 + a0/(pi AR)), e = 1, and a constant downwash
b = 2.4; c0 = 0.3; a0 = 2*pi; aL0 = 0; al = 5*pi/180;
cfun = @(y) c0*sqrt(max(0, 1 - (2*y/b).^2));
LL = lifting_line(b, cfun, a0, aL0, al, 25);
AR = b^2/(pi*c0*b/4);
ok = chk('Elliptic wing CL (exact)', LL.CL, a0*al/(1 + a0/(pi*AR)), 1e-3) && ok;
ok = chk('Elliptic wing span efficiency e = 1', LL.e, 1, 1e-3) && ok;

%  Rectangular wing: converged in the number of terms, e < 1
cr = @(y) 0.25*ones(size(y));
L1 = lifting_line(2.4, cr, a0, aL0, al, 15);
L2 = lifting_line(2.4, cr, a0, aL0, al, 45);
ok = chk('Rectangular wing CL converged in N', L1.CL, L2.CL, 2e-3) && ok;
ok = chk('Rectangular wing e between 0.9 and 1', double(L2.e > 0.9 && L2.e < 1), 1, 0) && ok;
%  Helmbold's approximation is independent and accurate to a few % at this AR
Ah = 2.4^2/0.6; CLa_h = a0/(sqrt(1 + (a0/(pi*Ah))^2) + a0/(pi*Ah));
ok = chk('Rectangular CL slope ~ Helmbold (within 3%)', L2.CL/al, CLa_h, 3e-2) && ok;

%% 2. Statics: root moment of standard load shapes
x = linspace(0, 1.2, 2001)'; s = 1.2; P = 100;
[V, M] = cantilever_statics(x, P/s*ones(size(x)));
ok = chk('Uniform load: root moment = P s/2', M(1), P*s/2, 1e-6) && ok;
ok = chk('Uniform load: root shear = P', V(1), P, 1e-9) && ok;
[~, M] = cantilever_statics(x, P*4/(pi*s)*sqrt(max(0, 1 - (x/s).^2)));
ok = chk('Elliptic load: root moment = P 4s/(3 pi)', M(1), P*4*s/(3*pi), 2e-3) && ok;

%% 3. Section properties and the load case
T = tube_props(0.025, 0.001);
ok = chk('Tube 25x1 I [mm^4]', T.I*1e12, pi*(25^4 - 23^4)/64, 1e-12) && ok;
ac = load_aircraft();
ok = chk('Load case: CL needed at n = 4, 25 m/s', ac.CL_req, 4*ac.W/(0.5*1.225*25^2*0.6), 1e-12) && ok;
ok = chk('NACA 4412 depth at 25% chord [mm]', 1000*ac.spar_depth, 29.70, 1e-3) && ok;
Lw = wing_loads(ac);
ok = chk('Half-wing lift integrates to n W / 2', trapz(Lw.x, Lw.q_lim), ac.L_half_lim, 1e-9) && ok;

%% 4. Finite-element beam
try
    EI = 1; L = 1;
    w = beam_fem(L, EI, @(xx) 1 + 0*xx, 1);
    ok = chk('FEM uniform load, ONE element: qL^4/8EI', w(end), 1/8, 1e-10) && ok;
    w = beam_fem(L, EI, @(xx) 0*xx, 4, 1);
    ok = chk('FEM tip point load: PL^3/3EI', w(end), 1/3, 1e-10) && ok;
    q = @(xx) sqrt(max(0, 1 - xx.^2));
    ex = exact_tip_deflection(1, @(xx) 1 + 0*xx, q);
    w = beam_fem(1, 1, q, 32);
    ok = chk('FEM elliptic load vs unit-load oracle', w(end), ex, 2e-4) && ok;
    EIt = @(xx) 1 - 0.7*xx;                    % tapered: stiffness falls to 30% at tip
    ex = exact_tip_deflection(1, EIt, q);
    w8 = beam_fem(1, EIt, q, 8); w16 = beam_fem(1, EIt, q, 16); w32 = beam_fem(1, EIt, q, 32);
    e8 = abs(w8(end) - ex); e16 = abs(w16(end) - ex); e32 = abs(w32(end) - ex);
    fprintf('   tapered beam error: %.2e, %.2e, %.2e (8, 16, 32 elements)\n', e8, e16, e32);
    ok = chk('FEM convergence order on tapered beam ~ 2', log2(e16/e32), 2, 0.1) && ok;
catch err
    if strcmp(err.identifier, 'beam_fem:todo')
        fprintf('%-48s NOT IMPLEMENTED (src/beam_fem.m has open TODOs)\n', 'FEM tests');
        ok = false;
    else
        rethrow(err);
    end
end

%% 5. Bench-test fit separates tube stiffness from clamp flexibility
EI0 = 1500; kc = 400; Ls = [0.6 0.9];                    % synthetic test rig
C = Ls.^3/(3*EI0) + Ls.^2/kc;
X = [Ls(:).^3/3 Ls(:).^2]; sol = X\C(:);
ok = chk('Two-length fit recovers EI', 1/sol(1), EI0, 1e-9) && ok;

if ok, fprintf('\nALL TESTS PASS\n'); else, fprintf('\nNOT ALL TESTS PASS (see above)\n'); end
end

function ok = chk(name, got, want, reltol)
    if want == 0, e = abs(got); else, e = abs(got - want)/abs(want); end
    ok = e <= reltol;
    s = 'FAIL'; if ok, s = 'PASS'; end
    fprintf('%-48s got %12.6g  want %12.6g  -> %s\n', name, got, want, s);
end
