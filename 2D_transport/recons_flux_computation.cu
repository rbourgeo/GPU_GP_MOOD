#include "Godunov_Flux.h"
#include "parameters.h"
#include "reconstruction.h"

__global__ void ftcs(double f[]
                        , double fluxes_x[]
                        , double fluxes_y[]
                        , const double dx
                        , const double dy
                        , const double dt
                        , const double zT[]
                        , const int index[]
                        , const double gauss_weight[])
{
    int tidx = c2f(threadIdx.x + blockIdx.x * blockDim.x);
    int tidy = c2f(threadIdx.y + blockIdx.y * blockDim.y);

    double fx, fy, ul, ur, ub, ut;

    if (tidx >= 0 && tidx <= lf) // from 1's left intf to lf's right interface
    {
        if (tidy >= 0 && tidy <= nf) {
            {
                fluxes_x[ij(tidy, tidx)] = 0.0;
                fluxes_y[ij(tidy, tidx)] = 0.0;

                for (int r = 0; r <= ngp - 1; r++) {

                    /* high order reconstruction*/
                    reconstruction(ul, ur, tidy, tidx, f, dir_x, zT, index, r);
                    reconstruction(ub, ut, tidy, tidx, f, dir_y, zT, index, r);

                    /*Fluxes  computation*/
                    fx = Godunov_Flux(ul, ur, dir_x);
                    fy = Godunov_Flux(ub, ut, dir_y);

                    /*Fluxes storage*/
                    fluxes_x[ij(tidy, tidx)] += fx * gauss_weight[a2l_ngp(ngp - 1, r)];
                    fluxes_y[ij(tidy, tidx)] += fy * gauss_weight[a2l_ngp(ngp - 1, r)];
                }
                fluxes_x[ij(tidy, tidx)] *= dt / dx;
                fluxes_y[ij(tidy, tidx)] *= dt / dy;
            }

        }
    }
}
