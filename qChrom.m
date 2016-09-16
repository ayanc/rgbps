% function rho = qChrom(Q)
%
% Return a set of QxQ chromaticities uniformly sampled on the
% positive eigth sphere.
%
% Copyright (C) 2016, Ayan Chakrabarti <ayanc@ttic.edu>
function rho = qChrom(Q)

g = linspace(0,1,Q+2); g = g(2:end-1);
th = linspace(0,pi/2,Q+2); th = th(2:end-1);

[g,th] = meshgrid(g,th); g = g(:); th = th(:);

rho = repmat(g,[1 3]);

g = sqrt(1-g.^2);

rho(:,1) = g.*cos(th);
rho(:,3) = g.*sin(th);