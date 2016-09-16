/*
  Fast CUDA version that implements the following method. Use caller
  space variable pointers in an un-official way to prevent memory
  duplication on the GPU.

  Expects the following local variables to be present in the workspace
  of the calling function: psz,nValid,nxF,nyF,l2,c,imR,imG,imB,sc_rh
  
    nlF = bsxfun(@rdivide,c,sqrt(nxF.^2 + nyF.^2 + 1));
    nxF = nxF .* nlF; nyF = nyF .* nlF;

    sc_rh = sum( (l2(1,1)*nxF + l2(2,1)*nyF + l2(3,1)*nlF - imR).^2, 2);
    sc_rh = sc_rh + ...
	    sum( (l2(1,2)*nxF + l2(2,2)*nyF + l2(3,2)*nlF - imG).^2, 2);
    sc_rh = sc_rh + ...
	    sum( (l2(1,3)*nxF + l2(2,3)*nyF + l2(3,3)*nlF - imB).^2, 2);

  If any of the above variables don't exist, or aren't on the gpu when
  they're expected to be, this function might crash very badly!
	    
  Copyright (C) 2016, Ayan Chakrabarti <ayanc@ttic.edu>
*/

#include "mex.h"
#include "gpu/mxGPUArray.h"
#include <stdint.h>

#define F float

#define NUMT 1024

void __global__ getSSD(F * l2, F * nxf, F * nyf, F * cf,
		       F * imr, F * img, F * imb,
		       F * sc_rh,
		       int psz, int nV) {


	__shared__ F ell[9];
	F nxv, nyv, nlv, cv, val, sum;
	int i, j, nloc;
	

	i = threadIdx.x;

	if(i < 9) ell[i] = l2[i];
	__syncthreads();

	i += blockDim.x * blockIdx.x;
	if(i < nV) {
		cv = cf[i]; sum = 0; nloc = i;
		for(j = 0; j < psz; j++) {
			nxv = nxf[nloc]; nyv = nyf[nloc]; 
			nlv = cv / sqrtf(1.0+nxv*nxv+nyv*nyv); 
			nxv = nxv*nlv; nyv = nyv*nlv;

			val = ell[0]*nxv + ell[1]*nyv + ell[2]*nlv - imr[nloc];
			sum += val*val;

			val = ell[3]*nxv + ell[4]*nyv + ell[5]*nlv - img[nloc];
			sum += val*val;

			val = ell[6]*nxv + ell[7]*nyv + ell[8]*nlv - imb[nloc];
			sum += val*val;

			nloc += nV;
		}
		sc_rh[i] = sum;
	}

}


F * getGPUmem(const char * name) {

	const mxGPUArray * tmp;
	F * dptr;

	if(!mxIsGPUArray(mexGetVariablePtr("caller",name)))
		mexPrintf("%s is not on gpu!\n",name);

	tmp = mxGPUCreateFromMxArray(mexGetVariablePtr("caller",name));
	dptr = (F*) mxGPUGetDataReadOnly(tmp);
	mxGPUDestroyGPUArray(tmp);

	return (F*) dptr;
}

/* function getSSD */
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {

	F * l2, * nxf, * nyf, * cf, * imr, * img, * imb, * sc_rh;
	int psz, nV, nB;

	psz = mxGetScalar(mexGetVariablePtr("caller","psz")); psz *= psz;
	nV = mxGetScalar(mexGetVariablePtr("caller","nValid"));

	nB = (nV+NUMT-1)/NUMT;


	l2 = getGPUmem("l2"); sc_rh = getGPUmem("sc_rh");
	nxf = getGPUmem("nxF");nyf = getGPUmem("nyF");cf = getGPUmem("c");
	imr = getGPUmem("imR");img = getGPUmem("imG");imb = getGPUmem("imB");

	getSSD<<<nB,NUMT>>>(l2,nxf,nyf,cf,imr,img,imb,sc_rh,psz,nV);
}
