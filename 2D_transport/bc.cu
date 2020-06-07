#include "bc.h"
#include "parameters.h"

__global__ void bc(double f[])
{
    int tidx = c2f(threadIdx.x + blockIdx.x * blockDim.x);
    int tidy = c2f(threadIdx.y + blockIdx.y * blockDim.y);
    //

    // if(tidy == 0) //use only nf threads for 2D BC
    // {
    // 	f[ij(0   , tidx)] = f[ij(1   , tidx)];//bottom
    // 	f[ij(nf-1, tidx)] = f[ij(nf-2, tidx)];//top
    // }
    //
    // if(tidx == 0) //use only nf threads for 2D BC
    // {
    // 	f[ij(tidy, 0   )] = f[ij(tidy, 1   )];//left
    // 	f[ij(tidy, lf-1)] = f[ij(tidy, lf-2)];//right
    // }

    for (int k = 1; k <= ngc; k++) {
        if (tidy == k) //use only nf threads for 2D BC
        {
            f[ij(1 - k, tidx)] = f[ij(nf - k + 1, tidx)]; //bottom
            f[ij(nf + k, tidx)] = f[ij(1 + k - 1, tidx)]; //top
        }

        if (tidx == k) //use only nf threads for 2D BC
        {
            f[ij(tidy, 1 - k)] = f[ij(tidy, lf - k + 1)]; //left
            f[ij(tidy, lf + k)] = f[ij(tidy, 1 + k - 1)]; //right
        }
    }
}
