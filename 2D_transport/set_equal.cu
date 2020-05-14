#include "parameters.h"
#include "2D_bij.h"
#include "set_equal.h"
// put fout's terms in f, usefull to avoid race cdt ( the out
//  notation is confusing but it fits the variables names)
__global__ void set_equal( float f[], const float fout[])
{
	int tidx = threadIdx.x + blockIdx.x*blockDim.x;
	int tidy = threadIdx.y + blockIdx.y*blockDim.y;


	if(tidx > 0 && tidx < lf-1) // Skip boundaries!
	{
		if(tidy > 0 && tidy < nf-1)
		{
			{
				f[ij(tidy,tidx)] = fout[ij(tidy,tidx)];
			}
		}
	}

}
