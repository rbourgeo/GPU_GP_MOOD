#include "parameters.h"
#include "set_equal.h"

__global__ void set_equal(double out[], const double in[])
{
    int tidx = c2f(threadIdx.x + blockIdx.x * blockDim.x);
    int tidy = c2f(threadIdx.y + blockIdx.y * blockDim.y);

    if (tidx >= 1 - ngc && tidx <= lf + ngc) {
        if (tidy >= 1 - ngc && tidy <= nf + ngc) {
            {
                out[ij(tidy, tidx)] = in[ij(tidy, tidx)];
            }
        }
    }
}

__global__ void set_comb_lin(double out[], const double in1[], const double in2[], const double coef1, const double coef2)
{
    int tidx = c2f(threadIdx.x + blockIdx.x * blockDim.x);
    int tidy = c2f(threadIdx.y + blockIdx.y * blockDim.y);

    if (tidx >= 1 - ngc && tidx <= lf + ngc) {
        if (tidy >= 1 - ngc && tidy <= nf + ngc) {
            {
                out[ij(tidy, tidx)] = coef1 * in1[ij(tidy, tidx)] + coef2 * in2[ij(tidy, tidx)];
            }
        }
    }
}
