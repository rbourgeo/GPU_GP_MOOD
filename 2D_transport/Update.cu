#include "Update.h"
#include "parameters.h"

__global__ void Update_A(double f[], const double fluxes_x[], const double fluxes_y[])
{
    int tidx = c2f(threadIdx.x + blockIdx.x * blockDim.x);
    int tidy = c2f(threadIdx.y + blockIdx.y * blockDim.y);

    if (tidx >= 0 && tidx <= lf) // Skip boundaries!
    {
        if (tidy >= 0 && tidy <= nf) {
            {
                f[ij(tidy, tidx)] += -fluxes_x[ij(tidy, tidx)] - fluxes_y[ij(tidy, tidx)];
            }
        }
    }
}

__global__ void Update_B(double f[], const double fluxes_x[], const double fluxes_y[])
{
    int tidx = c2f(threadIdx.x + blockIdx.x * blockDim.x);
    int tidy = c2f(threadIdx.y + blockIdx.y * blockDim.y);

    if (tidx >= 0 && tidx <= lf) // Skip boundaries!
    {
        if (tidy >= 0 && tidy <= nf) {
            {
                f[ij(tidy + 1, tidx)] += +fluxes_y[ij(tidy, tidx)];
            }
        }
    }
}

__global__ void Update_C(double f[], const double fluxes_x[], const double fluxes_y[])
{
    int tidx = c2f(threadIdx.x + blockIdx.x * blockDim.x);
    int tidy = c2f(threadIdx.y + blockIdx.y * blockDim.y);

    if (tidx >= 0 && tidx <= lf) // Skip boundaries!
    {
        if (tidy >= 0 && tidy <= nf) {
            {
                f[ij(tidy, tidx + 1)] += +fluxes_x[ij(tidy, tidx)];
            }
        }
    }
}
