#include "parameters.h"
#include "2D_bij.h"
#include "init.h"

__global__ void initialize(float f[], float x[], float y[], const float dx, const float dy)
{
	int tidx = threadIdx.x + blockIdx.x*blockDim.x;
	int tidy = threadIdx.y + blockIdx.y*blockDim.y;

	float xt = -0.5f + (float(tidx)+0.5f)*dx; //centers of cells
	float yt = -0.5f + (float(tidy)+0.5f)*dy;

	if(tidx < lf)
	{
		x[tidx] = xt;
	}

	if(tidy < nf)
	{
		y[tidy] = yt;
	}

	if(tidx < lf){
		if(tidy < nf){
			f[ij(tidy,tidx)] = exp(-0.5f*(xt*xt + yt*yt));
		}
	}

}
