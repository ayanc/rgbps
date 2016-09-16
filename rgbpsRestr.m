% [coeff,score] = rgbpsRestr(img,mask,l,patch_size,h_q,l_q,ropts) 
%
%   img = Observed image
%   mask = Mask
%   l = 3x3 Light matrix, [lr lg lb]
%   patch_size = patch size
%   h_q, h_lq = output of hMax (selected albedos)
%   ropts = struct with marapeters
%
% Returns: all matrices are of size 
%          (img_height-psz+1)x(img_width-psz+1)xNxlength(q)
%   coeff: N = 20, coefficients of a 5 degree polynomial model
%   score: N = 1. Normalized Rendering error.
%
% Copyright (C) 2016, Ayan Chakrabarti <ayanc@ttic.edu>
function [ocoeff,score] = rgbpsRestr(img,mask,l,psz,h_q,h_lq,ropts)

tic;

% Set up luminance-chromaticity discretization
Q = ropts.Q; LQ = ropts.LQ; LMAX = ropts.LMAX;
rho = qChrom(Q); rho = rho(h_q,:);
c_lb = (h_lq-2)/(LQ-1)*LMAX;
c_ub = (h_lq)/(LQ-1)*LMAX;


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

% Allocate memory for outputs on cpu
oh = (imh-psz+1); ow = (imw-psz+1);
osz = ow*oh;

ocoeff = zeros([osz msz length(h_q)],'single');
score = Inf*ones([osz 1 length(h_q)],'single');


% Normalized rendering error
snrm = sum(imR.^2,2) + sum(imG.^2,2) + sum(imB.^2,2) + 1e-8;
NMAX = ropts.NMAX;

stime = toc;
fprintf('Setup Time: %.4fs ',stime);
fprintf('Starting optimization for %d valid patches.\n',nValid);
% Go through all chromaticities
for rh = 1:size(rho,1)

  fprintf('\r Testing albedo %d of %d         [%.4f s]    ',rh,size(rho,1),(toc-stime)/(rh-1));

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
  c = max(c_lb(rh),min(c_ub(rh),c)); % Restrict to bin
  coeff = nx*A1 + ny*A2;       % Fit to poly surface

  % Compute re-rendering score
  nxF = coeff*nxm; nyF = coeff*nym;
  l2 = l*diag(rho(rh,:));
  getSSD;

  % Gather scores
  sc_rh = sc_rh ./ snrm;

  score(valid,1,rh) = gather(sc_rh);
  ocoeff(valid,:,rh) = gather(coeff);
end;

score = reshape(score,[oh ow 1 length(h_q)]);
ocoeff = reshape(ocoeff,[oh ow msz length(h_q)]);
fprintf('\nTotal Time: %.4f s\n',toc);
