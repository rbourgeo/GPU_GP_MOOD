#include "parameters.h"
#include "set_equal.h"
#include "constants.h"


__global__ void set_equal(double out[], const double in[])
{
  int tidx = c2f(threadIdx.x + blockIdx.x * blockDim.x);
  int tidy = c2f(threadIdx.y + blockIdx.y * blockDim.y);

  if (tidx >= 1 - ngc && tidx <= lf + ngc) {
    if (tidy >= 1 - ngc && tidy <= nf + ngc) {
      {
        for (int i_cons = i_rho; i_cons <= i_ener; i_cons++) {
          out[ij_sol(tidy, tidx, i_cons)] = in[ij_sol(tidy, tidx, i_cons)];
        }
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
      for (int i_cons = i_rho; i_cons <= i_ener; i_cons++) {
        {
          out[ij_sol(tidy, tidx, i_cons)] = coef1 * in1[ij_sol(tidy, tidx, i_cons)] + coef2 * in2[ij_sol(tidy, tidx, i_cons)];
        }
      }
    }
  }
}
