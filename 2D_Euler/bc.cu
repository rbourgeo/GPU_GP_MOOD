#include "bc.h"
#include "parameters.h"

__global__ void bc(double f[])
{
    int tidx = c2f(threadIdx.x + blockIdx.x * blockDim.x);
    int tidy = c2f(threadIdx.y + blockIdx.y * blockDim.y);

    for (int k = 1; k <= ngc; k++) {

      for (int i_cons = i_rho; i_cons <= i_ener; i_cons++) {

        if (tidy == k) //use only nf threads for 2D BC
        {
            f[ij_sol( 1 - k, tidx, i_cons)] = f[ij_sol(nf - k + 1, tidx, i_cons)]; //bottom
            f[ij_sol(nf + k, tidx, i_cons)] = f[ij_sol( 1 + k - 1, tidx, i_cons)]; //top
        }

        if (tidx == k) //use only nf threads for 2D BC
        {
            f[ij_sol(tidy,  1 - k, i_cons)] = f[ij_sol(tidy, lf - k + 1, i_cons)]; //left
            f[ij_sol(tidy, lf + k, i_cons)] = f[ij_sol(tidy,  1 + k - 1, i_cons)]; //right
        }
      }
      
    }

}
