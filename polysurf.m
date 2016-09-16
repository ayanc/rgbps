% [z,nx,ny] = polysurf(sz,deg)
%
%    Construct z, nx, ny matrices so that each matrix * coeff gives
%    surface and normals respectively, where coeff is the
%    coefficient vector. For sz x sz patch, polynomial of degree dg.
%
% Copyright (C) 2016, Ayan Chakrabarti <ayanc@ttic.edu>
function [z,nx,ny] = polysurf(sz,deg)

y = [1:sz]-(sz+1)/2;
[x,y] = meshgrid(y,y);
x = x(:); y = y(:);

z = []; nx = []; ny = [];

for xdeg = 0:deg
  for ydeg = 0:(deg-xdeg)
    
    % No base offset
    if xdeg + ydeg == 0
      continue;
    end;

    z = [z (x.^xdeg).*(y.^ydeg)];
    nx = [nx xdeg.*(x.^(max(0,xdeg-1))).*(y.^ydeg)];
    ny = [ny ydeg.*(y.^(max(0,ydeg-1))).*(x.^xdeg)];
  end;
end;