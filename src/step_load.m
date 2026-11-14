function q = step_load(x, edges, F)
%STEP_LOAD  Piecewise-constant distributed load: strip k carries F(k) spread
%   evenly over [edges(k), edges(k+1)). Matches applying each strip force to
%   its split face in SolidWorks Simulation.
    q = zeros(size(x));
    w = diff(edges(:));
    for k = 1:numel(F)
        in = x >= edges(k) & x < edges(k+1);
        if k == numel(F), in = in | x == edges(end); end
        q(in) = F(k)/w(k);
    end
end
