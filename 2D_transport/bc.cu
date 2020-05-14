#include "parameters.h"
#include "2D_bij.h"
#include "bc.h"


__global__ void bc(float *f)
{
	int tidx = threadIdx.x + blockIdx.x*blockDim.x;
	int tidy = threadIdx.y + blockIdx.y*blockDim.y;

	if(tidy == 0) //use only nf threads for 2D BC
	{
		f[ij(0   , tidx)] = f[ij(1   , tidx)];//bottom
		f[ij(nf-1, tidx)] = f[ij(nf-2, tidx)];//top
	}

	if(tidx == 0) //use only nf threads for 2D BC
	{
		f[ij(tidy, 0   )] = f[ij(tidy, 1   )];//left
		f[ij(tidy, lf-1)] = f[ij(tidy, lf-2)];//right
	}

}
