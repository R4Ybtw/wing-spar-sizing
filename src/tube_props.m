function T = tube_props(D, t)
%TUBE_PROPS  Section properties of a round tube (outer diameter D, wall t, metres).
%
%   T.I   second moment of area (m^4)     pi (D^4 - d^4)/64
%   T.A   cross-section area (m^2)
%   T.c   distance to the outer fibre (D/2), so bending stress = M c / I
%   For thin walls I ~ pi R^3 t: diameter cubed, thickness only to the first
%   power. That's why a bigger, thinner tube is usually the lighter answer.

    d = D - 2*t;
    T.I = pi*(D.^4 - d.^4)/64;
    T.A = pi*(D.^2 - d.^2)/4;
    T.c = D/2;
    T.D = D; T.t = t;
end
