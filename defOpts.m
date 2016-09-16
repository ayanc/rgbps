% The following is a list of the available options that can be
% passed as the last argument to the doRGBPS function, along with
% their default values. Note that you only need to define a
% structure with members corresponding to the options that you want
% to change.
%
% For example, to only change the parameter nms_num, define:
% 
% >> ropts = struct; ropts.nms_num = 200;
%
% and pass ropts as the last argument to doRGBPS.
%
% Options
% -------
%
% ropts.MODEL:     Polynomial degree for shape model, D in paper 
%                  (default: 5)
%
% ropts.Q:          Number of bins to quantize each of the two
%                   chromaticity dimensions, i.e. number of bins
%                   will be ropts.Q^2  
%                   (default: 64)
%
% ropts.LQ:         Number of bins to quantize luminance 
%                   (default: 100)
%
% ropts.LMAX:       Maximum value of luminance to consider,
%                   \tau_{max} in paper 
%                   (default: 3)
%
% ropts.h_thresh:   Clipping threshold for histogram scores,
%                   h_{max} in paper.
%                   (default: 1e-2)
% 
% ropts.nms_num:    Size of global albedo set, K in paper.
%                   (default: 100)
%
% ropts.gTHR:       Outlier threshold in globalization, \gamma in
%                   paper.
%                   (default: 4)
%
% ropts.gLS:        This is a vector whose length is the number of
%                   iterations that you want globalization to run
%                   for, with ropts.gLS(i) specifying the value of
%                   lambda in iteration i.
%                   (default: 2.^[-64:0.5:8])
%
% ropts.hpsz,
% ropts.psizes:     These two specify the patch sizes to be used
%                   for estimation. hpsz (a scalar) that specifies
%                   the patch-size to use to build the histogram.
%                   psizes is a *list* of patch-sizes over which to
%                   compute local shape distributions---even though
%                   in practice we only compute distributions for
%                   one patch size, the code supports incorporating
%                   patches from multiple scales.
%                   (default: 8 & {8} to use 8x8 patches for everything)
%
% Copyright (C) 2016, Ayan Chakrabarti <ayanc@ttic.edu>
function ropts = defOpts(ropts)

%%% Set default values
r = struct;

% Common
r.MODEL=5; r.LMAX=3; r.LQ=100; r.Q=64;

r.NMAX=1000; % Clip normals for stability

% rgbpsHist
r.h_thresh = 1e-2;
% hMax
r.nms_num = 100;
% globalization
r.gLS=2.^[-64:0.5:8]; r.gTHR=4;

% Main
r.hpsz = 8; r.psizes = {8};

%%% Replace undefined values in ropts with defaults
flds = fieldnames(r);
for i = 1:length(flds)
  if ~isfield(ropts,flds{i})
    ropts = setfield(ropts,flds{i},getfield(r,flds{i}));
  end;
end;