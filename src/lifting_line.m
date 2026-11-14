function LL = lifting_line(span, chord_fun, a0, alpha_L0, alpha, N, twist_fun)
%LIFTING_LINE  Prandtl lifting-line theory, solved by Glauert's Fourier-series method.
%
%   LL = lifting_line(span, chord_fun, a0, alpha_L0, alpha, N)
%   LL = lifting_line(..., twist_fun)
%
%   span       full span b (m)
%   chord_fun  chord c(y) as a function handle, y = 0 at the root, b/2 at the tip
%   a0         2D lift slope of the section (per radian)
%   alpha_L0   zero-lift angle of the section (radians)
%   alpha      root angle of attack (radians)
%   N          number of Fourier terms (odd terms only, by symmetry), e.g. 25
%   twist_fun  optional geometric twist (radians) vs y, added to alpha
%
%   Output LL:
%     .CL, .CDi, .e       wing lift, induced drag, span efficiency
%     .y                  stations on the half-span, root -> tip
%     .Gamma              circulation per unit (V = 1): multiply by V for the real value
%     .cl                 local section lift coefficient
%     .lift_shape         spanwise lift per unit span divided by the half-wing
%                         total, so it integrates to 1 over 0..b/2
%     .A                  Fourier coefficients
%
%   Idea: a finite wing sheds trailing vortices. They induce a downwash that
%   lowers each section's effective angle of attack, more strongly toward the
%   tips. Writing the circulation as Gamma = 2 b V sum(A_n sin(n theta)),
%   with y = (b/2) cos(theta), turns the integral equation into one linear
%   equation per station:
%       sum_n A_n sin(n theta) (mu n + sin theta) = mu (alpha_eff - alpha_L0) sin theta,
%       mu = a0 c / (4 b)
%   Then CL = pi AR A_1 and CDi = pi AR sum(n A_n^2).
%   Valid for straight wings of moderate-to-high aspect ratio (above ~4)
%   in attached flow. It doesn't model sweep or stall.

    if nargin < 7 || isempty(twist_fun), twist_fun = @(y) 0*y; end
    n  = 1:2:(2*N - 1);                          % odd terms (symmetric loading)
    th = (1:N)'*pi/(2*N);                        % stations on (0, pi/2], tip excluded
    y  = (span/2)*cos(th);
    c  = chord_fun(y);
    S  = 2*integral_half(chord_fun, span/2);
    AR = span^2/S;
    mu = a0*c/(4*span);
    Aeff = alpha + twist_fun(y);
    M = sin(th*n).*(mu*n + repmat(sin(th), 1, N));
    rhs = mu.*(Aeff - alpha_L0).*sin(th);
    A = M\rhs;

    LL.A = A; LL.n = n;
    LL.CL  = pi*AR*A(1);
    LL.CDi = pi*AR*sum(n(:).*A.^2);
    LL.e   = LL.CL^2/(pi*AR*LL.CDi);
    LL.AR = AR; LL.S = S;

    % Dense output on the half span, root -> tip
    yy = linspace(0, span/2, 401)';
    tt = acos(min(yy/(span/2), 1));
    G  = 2*span*sin(tt*n)*A;                     % Gamma / V
    LL.y = yy;
    LL.Gamma = G;
    LL.cl = 2*G./chord_fun(yy);                  % l' = rho V Gamma = 0.5 rho V^2 c cl
    LL.cl(end) = NaN;                            % 0/0 at the tip
    half = trapz(yy, G);
    LL.lift_shape = G/half;                      % integrates to 1 over the half span
end

function A = integral_half(f, s)
    yy = linspace(0, s, 2001);
    A = trapz(yy, f(yy));
end
