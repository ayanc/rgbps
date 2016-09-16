% hist = rgbpsHist(img,mask,l,patch_size,ropts) 
%
% Get albedo histogram of scores to help identify small number of
% albedos in scene.
%
%   img = Observed image
%   mask = Mask
%   l = 3x3 Light matrix, [lr lg lb]
%   patch_size = patch size
%   ropts = struct with marapeters
%
% Returns: hist = size #ChromaBins x #LuminanceBins of scores.
%
% Copyright (C) 2016, Ayan Chakrabarti <ayanc@ttic.edu>
function hist = rgbpsHist(img,mask,l,psz,ropts)

tic;

% Set up luminance-chromaticity discretization
Q = ropts.Q; LQ = ropts.LQ; LMAX = ropts.LMAX;
rho = qChrom(Q);

% Set up polynomial model matrices
MODEL=ropts.MODEL; 
[zm,nxm,nym] = polysurf(psz,MODEL); msz = size(zm,2);
A = [nxm; nym]; A = A*inv(A'*A); 
nxm = nxm'; nym = nym';
% Split coeff matrix into *nx and *ny components
A1 = A(1:(psz^2),:); A2 = A((psz^2)+1:end,:);

% Set up im2col stuff
imh = size(img,1); imw = size(img,2);
im2p = reshape([1:imh*imw],[imh imw]);
im2p = im2col(im2p,[psz psz],'sliding');  % Indices

% Valid patches (all pixels valid)
msk = mask(im2p); valid = uint32(find(all(msk,1)));
nValid = length(valid); im2p = im2p(:,valid)';

% Move stuff to GPU in single precision
rho = single(rho); img = single(img);
A1 = single(A1); A2 = single(A2);
linv = single(inv(l)); l = single(l);
nxm = single(nxm); nym = single(nym);

rho = gpuArray(rho); img = gpuArray(img);
A1 = gpuArray(A1); A2 = gpuArray(A2); 
linv = gpuArray(linv); l = gpuArray(l);
nxm = gpuArray(nxm); nym = gpuArray(nym);
im2p = gpuArray(im2p);
%%%%

% Pixel-level image
img = reshape(img,[imh*imw 3]);
% RGB Image patches
imR = img(:,1); imR = imR(im2p);
imG = img(:,2); imG = imG(im2p);
imB = img(:,3); imB = imB(im2p);

% Allocate memory for stuff
cp = zeros(size(im2p),'single','gpuArray');
nx = zeros(size(im2p),'single','gpuArray');
ny = zeros(size(im2p),'single','gpuArray');
sc_rh = zeros(nValid,1,'single','gpuArray');
LQidx = gpuArray(single([1:LQ]));

% Final results on CPU
hist = zeros(Q^2,LQ,'single');


% Normalized rendering error
snrm = sum(imR.^2,2) + sum(imG.^2,2) + sum(imB.^2,2) + 1e-8;
NMAX = ropts.NMAX;
thresh = ropts.h_thresh;

stime = toc;
fprintf('Setup Time: %.4fs ',stime);
fprintf('Starting optimization for %d valid patches.\n',nValid);
% Go through all chromaticities
for rh = 1:size(rho,1)

  fprintf('\r Testing chromaticity %d of %d         [%.4f s]    ',rh,size(rho,1),(toc-stime)/(rh-1));

  % Solve per-pixel normals
  linv2 = diag(1./rho(rh,:))*linv;
  
  nxp = img*linv2(:,1);  nyp = img*linv2(:,2);  nzp = img*linv2(:,3);
  
  cpp = sqrt(nxp.^2+nyp.^2+nzp.^2);
  nxp = nxp ./ nzp; nyp = nyp ./ nzp;
  
  % Clip to +- NMAX
  fac = sqrt(nxp.^2+nyp.^2); fac = min(NMAX,fac) ./ max(1e-8,fac);
  nxp = nxp.* fac; nyp = nyp.* fac;

  % Pixel to patches
  cp = cpp(im2p); nx = nxp(im2p); ny = nyp(im2p);
  
  c = mean(cp,2);              % Fit to equal albedo
  coeff = nx*A1 + ny*A2;       % Fit to poly surface

  % Compute re-rendering score
  nxF = coeff*nxm; nyF = coeff*nym;
  l2 = l*diag(rho(rh,:));
  getSSD;

  % Normalize, threshold, and add to histogram
  sc_rh = sc_rh ./ snrm;
  sc_rh = max(0,thresh-sc_rh);
  
  c = c / LMAX * (LQ-1) + 1;
  rhHist = bsxfun(@eq,floor(c),LQidx) + bsxfun(@eq,ceil(c),LQidx);
  rhHist = sum( bsxfun(@times,rhHist,sc_rh), 1);
  hist(rh,:) = gather(rhHist);
end;

fprintf('\nTotal Time: %.4f s\n',toc);
