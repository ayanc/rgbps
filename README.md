# Single-shot RGB-PS With Spatially-varying Albedo

This directory contains a reference implementation of the algorithm described in the paper:

Ayan Chakrabarti and Kalyan Sunkavalli, "**[Single-image RGB Photometric Stereo With Spatially-varying Albedo](https://arxiv.org/abs/1609.04079)**", Proc. of the IEEE Intl. Conf. on 3D Vision (3DV) 2016. 

It is being released for non-commercial research use only. If you find this code useful in your research, we request that you cite the above paper. Along with the implementation source code, we are also releasing real datasets captured under the RGB-PS setup in the `data/` sub-directory of this repository/distribution. These datasets are owned by Adobe Systems, and are being released under a CC-by-attribution license: please see the `data/LICSENSE.txt` and `data/README.md` files for details.

## Usage

Our implementation requires a modern version of MATLAB with support for the `gpuArray` class, a CUDA-capable GPU (our experiments were run on a Titan X), as well as a distribution of CUDA (with the `nvcc` compiler).

You will first need to compile the mex function in `getSSD.cu`. Please see the documentation for your version of MATLAB to figure out how to compile mex files with CUDA code. If you have the latest version of MATLAB and all paths set up, this should be as simple as:

```MATLAB
>> mexcuda getSSD.cu
```

The main function to do estimation is `doRGBPS`. Given an input RGB image, a mask, and the lighting matrix of the RGB-PS setup, this function returns estimated normals for the surface. Please use `help doRGBPS` in MATLAB to get detailed instructions on usage. We also provide a simple utility function `getZ` that can be used to integrate these normals to get a depth map.

Please see the [project page](http://www.ttic.edu/chakrabarti/rgbps/) and contact <ayanc@ttic.edu> with any questions.
