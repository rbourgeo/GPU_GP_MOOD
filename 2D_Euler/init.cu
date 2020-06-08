#include "exact_sol.h"
#include "init.h"
#include "parameters.h"

__global__ void initialize(double f[], double x[], double y[], const double dx, const double dy)
{
    int tidx = c2f(threadIdx.x + blockIdx.x * blockDim.x);
    int tidy = c2f(threadIdx.y + blockIdx.y * blockDim.y);

    double xt = (double(tidx) - 0.5f) * dx; //centers of cells
    double yt = (double(tidy) - 0.5f) * dy;

    if (tidx <= lf + ngc) {
        x[f2c(tidx)] = xt;
    }

    if (tidy <= nf + ngc) {
        y[f2c(tidy)] = yt;
    }

    if (tidx <= lf + ngc) {
        if (tidy <= nf + ngc) {

            f[ij_sol(tidy, tidx, i_rho)] = exact_soln(tidx, tidy, dx, dy, x, y);

            for (int i_cons = i_momx; i_cons <= i_ener; i_cons++) {
              f[ij_sol(tidy, tidx, i_cons)] = 0.0;
            }
        }
    }
}
