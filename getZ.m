% Z = getZ(nrm)
%
% Integrate normals nrm to get depth map Z.
%
% Copyright (C) 2016, Ayan Chakrabarti <ayanc@ttic.edu>
function Z = getZ(nrm)

NUMIT=400; % Number of iterations
fac = 1;   % Up-sampling factor (set to 1 by default)

% Derivative filters (data term)
fx = [0 0 0; 1 0 -1; 0 0 0]/2/fac; fy = fx';

% Regularization filters
rfilt = [-1 2 -1]*1e-2;
rf1 = [0 0 0; rfilt; 0 0 0]; rf2 = diag(rfilt);
regf = cat(3,rf1,rf1',rf2,fliplr(rf2));

% Set up mask and normals
msk0 = single(nrm(:,:,3) > 0) .* max(eps,nrm(:,:,3).^2); 
h = size(msk0,1); w = size(msk0,2);
msk = zeros(([h w]-1)*fac+1); gx=msk; gy=msk;
msk(1:fac:end,1:fac:end) = msk0;
gx0 = nrm(:,:,1) ./ max(eps,nrm(:,:,3));
gy0 = nrm(:,:,2) ./ max(eps,nrm(:,:,3));
gx(1:fac:end,1:fac:end) = gx0;
gy(1:fac:end,1:fac:end) = gy0;


Z = zeros(size(msk)+2);

% Move to GPU
Z = gpuArray(Z); gx = gpuArray(gx); gy = gpuArray(gy);
msk = gpuArray(msk); regf = gpuArray(regf); 
fx = gpuArray(fx); fy = gpuArray(fy);

r = conv2(gx.*msk,fliplr(fx),'full') + ...
    conv2(gy.*msk,flipud(fy),'full');
p = r;

rsum = sum(r(:).^2);
for it = 1:NUMIT
  fprintf('\r Iteration %04d of %04d',it,NUMIT);
  Ap = conv2(conv2(p,fx,'valid').*msk,fliplr(fx),'full') + ...
       conv2(conv2(p,fy,'valid').*msk,flipud(fy),'full');
  for j = 1:size(regf,3)
    Ap = Ap + conv2(conv2(p,regf(:,:,j),'valid'), ...
		    fliplr(flipud(regf(:,:,j))),'full');
  end;
  
  alpha = rsum / sum(p(:).*Ap(:));
  Z = Z + alpha * p;
  r = r - alpha * Ap;

  beta = rsum; rsum = sum(r(:).^2);
  beta = rsum / beta;
  p = r + beta*p;
end;
fprintf('\n');

Z = gather(Z);
Z = Z(2:fac:end-1,2:fac:end-1);

Z(nrm(:,:,3) == 0) = nan;