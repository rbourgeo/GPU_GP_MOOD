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
              for (int i_cons = i_rho; i_cons <= i_ener; i_cons++) {

                f[ij_sol(tidy, tidx, i_cons)] += -fluxes_x[ij_sol(tidy, tidx, i_cons)] - fluxes_y[ij_sol(tidy, tidx, i_cons)];

              }
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
              for (int i_cons = i_rho; i_cons <= i_ener; i_cons++) {

                f[ij_sol(tidy + 1, tidx, i_cons)] += +fluxes_y[ij_sol(tidy, tidx, i_cons)];

              }
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
              for (int i_cons = i_rho; i_cons <= i_ener; i_cons++) {

                f[ij_sol(tidy, tidx + 1, i_cons)] += +fluxes_x[ij_sol(tidy, tidx, i_cons)];

              }
            }
        }
    }
}
