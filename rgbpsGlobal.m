% nrm = rgbpsGlobal(cf,score,psz,ropts) 
%
%   cf, score: Cell of matrices Returned by rgbpsRestr
%   psz: Patch Sizes
%   ropts: Structure with option values.
%
% Returns: 
%   nrm: Normal map
%
% Copyright (C) 2016, Ayan Chakrabarti <ayanc@ttic.edu>
function nrm = rgbpsGlobal(cfs,scores,psizes,ropts) 

lambdas = ropts.gLS;
THR = ropts.gTHR;
MODEL=ropts.MODEL; 
NMAX=ropts.NMAX;

tic;

% Set up estimates
imh = size(cfs{1},1) + psizes{1}-1;
imw = size(cfs{1},2) + psizes{1}-1;
nx = zeros([imh imw],'single','gpuArray');
ny = zeros([imh imw],'single','gpuArray');

cnv_mask = zeros([imh imw],'single','gpuArray');

oszs = {}; nxms = {}; nyms = {};
for t = 1:length(psizes)
  psz = psizes{t}; cf = cfs{t}; score = scores{t};
  
  
  % Set up model stuff
  [zm,nxm,nym] = polysurf(psz,MODEL);
  MSZ = size(zm,2);
  A = [nxm; nym]; AA = A'*A; A = A*inv(AA); 
  nxms{t} = reshape(nxm,[psz psz MSZ]);
  nyms{t} = reshape(nym,[psz psz MSZ]);

  % Set up im2col stuff
  im2p = reshape([1:imh*imw],[imh imw]);
  im2p = im2col(im2p,[psz psz],'sliding')';  % Indices

  % convmask
  cnv_mask = cnv_mask + ...
      conv2(gpuArray(single(isfinite(score(:,:,1,1)))),ones(psz,psz),'full');

  % Move other stuff to GPU
  osz = size(cf); osz = osz(1:3); oszs{t} = osz;
  NPROP=size(scores{t},4);

  
  cf=reshape(cf,[size(im2p,1) size(cf,3) NPROP]);
  score = reshape(score,[size(im2p,1) NPROP]);
  valid = find(isfinite(score(:,1)));

  cf = cf(valid,:,:); score = score(valid,:,:);
  im2p = gpuArray(single(im2p(valid,:)));
  cf = single(cf); %cf = gpuArray(cf);
  score = single(score); %score = gpuArray(score);
  valid = gpuArray(single(valid));

  A = gpuArray(single(A)); AA = gpuArray(single(AA));
  A1 = A(1:psz^2,:); A2 = A(psz^2+1:end,:);
  
  % Put back in cell array
  im2ps{t} = im2p; cfs{t} = cf; scores{t} = score;
  AAs{t} = AA; A1s{t} = A1; A2s{t} = A2;
  valids{t} = valid;
end;
mask = gather(cnv_mask > eps);
cnv_mask = max(eps,cnv_mask);

NUMIT=length(lambdas);
for i = 1:NUMIT
  lambda = lambdas(i);
  fprintf('\r Iteration %03d of %03d    ',i,NUMIT);
  
  nx0 = nx; ny0 = ny; nx(:) = 0; ny(:) = 0;
  
  
  for t = 1:length(psizes)
  
    cf = gpuArray(cfs{t}); score = gpuArray(scores{t});
    %cf = cfs{t}; score = scores{t};
    
    %%%%%% Upward
  
    % Compute projections based on current map
    nxi = nx0(im2ps{t}); nyi = ny0(im2ps{t});
    cf_z0 = nxi*A1s{t} + nyi*A2s{t};
  
  
    % Compute error with respect to local estimates and find the best
    % one (or declare outlier)
    cf_z = cf_z0;
    best_score = THR*ones(length(valids{t}),1,'single','gpuArray');
    for j = 1:size(score,2)
      cf_d = lambda/(1+lambda) * (cf_z-cf(:,:,j));
      sc_z = score(:,j) + sum(cf_d .* (cf_d*AAs{t}),2);
      idx = find(sc_z < best_score);
      best_score(idx) = sc_z(idx);
      cf_z(idx,:) = (cf(idx,:,j) + lambda*cf_z0(idx,:)) / (1+lambda);
    end;
  
    %%%%%%% Downward
    osz = oszs{t};
    cf_z0 = zeros(osz(1)*osz(2),osz(3),'single','gpuArray');
    cf_z0(valids{t},:) = cf_z;
    cf_z = reshape(cf_z0,osz);
    for j = 1:MSZ
      nx = nx + conv2(cf_z(:,:,j),nxms{t}(:,:,j),'full');
      ny = ny + conv2(cf_z(:,:,j),nyms{t}(:,:,j),'full');
    end;
  end;
  
  nx = nx ./ cnv_mask;  ny = ny ./ cnv_mask;
  
  fac = sqrt(nx.^2+ny.^2); fac = min(NMAX,fac) ./ max(1e-8,fac);
  nx = nx.* fac; ny = ny.* fac;
  
end;

nz = 1./sqrt(1+nx.^2+ny.^2);
nx = nx.*nz; ny = ny.*nz;

nx = gather(nx);ny = gather(ny);nz = gather(nz);
nrm = cat(3,nx,ny,nz.*mask);

fprintf('\n');