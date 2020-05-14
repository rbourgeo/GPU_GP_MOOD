#include "parameters.h"
#include "2D_bij.h"
#include "init.h"

// FCTS kernel, result is put in fout
__global__ void ftcs(const float f[], float fout[], const float dx, const float dy, const float k, const float dt)
{
	int tidx = threadIdx.x + blockIdx.x*blockDim.x;
	int tidy = threadIdx.y + blockIdx.y*blockDim.y;

	if(tidx > 0 && tidx < lf-1) // Skip boundaries!
	{
		if(tidy > 0 && tidy < nf-1)
		{
			{
				float temp =  f[ij(tidy,tidx)] + k*dt/(dx*dx)*(f[ij(tidy  ,tidx-1)] - 2*f[ij(tidy,tidx)] + f[ij(tidy  ,tidx+1)]);
				temp       =              temp + k*dt/(dy*dy)*(f[ij(tidy-1,tidx  )] - 2*f[ij(tidy,tidx)] + f[ij(tidy+1,tidx  )]);
				fout[ij(tidy,tidx)] = temp;
			}
		}
	}

}
