# RGB-PS Image Data for Real Objects
Copyright (C) 2016, Adobe Systems, Inc.

## LICENSE

This directory contains datasets associated with the paper
"Single-image RGB Photometric Stereo With Spatially-varying 
Albedo". This data was captured by Adobe Systems, Inc. and
is being released under the following license:

**This work is licensed under the Creative Commons Attribution 
4.0 International License. To view a copy of this license, visit 
http://creativecommons.org/licenses/by/4.0/ or send a letter to 
Creative Commons, PO Box 1866, Mountain View, CA 94042, USA.**

## Format

The data for each object is provided as a separate `.mat` file. Each
file has the following variables:

| Variable Name | Description                                                                            |
|:--------------|:---------------------------------------------------------------------------------------|
| `img_rgbps`   | Image of object captured under RGB-PS setup.                                           |
| `mask`        | Mask of valid pixels (i.e., contained in object).                                      |
| `lights`      | 3x3 lighting matrix `[lr lg lb]`                                                       |
| `nrm_gt`      | Ground truth normals (captured with multi-shot PS)                                     |
| `albedos_gt`  | Ground truth albedos (captured with multi-shot PS)                                     |
| `img_3s`      | Simulated 3 full color images captured from the same lighting directions (in `lights`) |
| `nrm_est`     | Normals estimated from `img_rgbps` using our approach                                  |
| `Z_est`       | Depth map from integrating `nrm_est`                                                   |

