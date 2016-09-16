% nrm = doRGBPS(img,mask,lights,ropts)
%
% Compute surface normals from an image img (with mask m) captured
% in an RGB-PS setup. Here, lights is the 3x3 matrix that describes
% the lighting environment, i.e., lights = [lr lg lb].
%
% The structure ropts is an optional argument that can be used to
% customize the parameters used by the estimation algorithm. If
% ommitted, the defaults (as described in the paper) are used.
%
% Use 'help defOpts' to get a list of options and default options.
%
% The returned matrix nrm contains the estimated normals---each
% nrm(i,j,:) will be a unit norm vector, except at places that are
% masked out where all three values will be zero. (Note that the
% output mask will be slightly smaller than the input mask).
%
% Copyright (C) 2016, Ayan Chakrabarti <ayanc@ttic.edu>
function nrm = fullRGBPS(img,mask,lights,ropts)

if ~exist('ropts')
  ropts=struct;
end;
if length(ropts) == 0
  ropts=struct;
end;
ropts = defOpts(ropts);


hpsz = ropts.hpsz;
psizes = ropts.psizes;


hist = rgbpsHist(img,mask,lights,hpsz,ropts);
[q,lq] = hMax(hist,ropts);

cfs = {}; scs = {};
for i = 1:length(psizes)
  [cf,sc] = rgbpsRestr(img,mask,lights,psizes{i},q,lq,ropts);
  cfs{i} = cf; scs{i} = sc;
end;
nrm = rgbpsGlobal(cfs,scs,psizes,ropts);

% Slightly erode mask
mask = single(imerode(mask > 0,strel('disk',5)));
nrm = bsxfun(@times,nrm,mask);