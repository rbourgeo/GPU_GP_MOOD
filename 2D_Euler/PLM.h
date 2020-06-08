#include "constants.h"
#include "minmod.h"
#include "parameters.h"

/* PLM / MUSCL 2nd order methods
outputs: ul and ur, high order pointwise face centered values
inputs : idx idy, coordinates of the interface
f: solution
dir: direction of reconstructon */
__host__ __device__ inline void PLM_(double ul[], double ur[], const int idy, const int idx, const double f[], const int dir)
{

  double slope;
  if (dir == dir_x) {

    for (int i_cons = i_rho; i_cons <= i_ener; i_cons++) {


      slope = minmod(f[ij_sol(idy, idx, i_cons)] - f[ij_sol(idy, idx - 1, i_cons)], f[ij_sol(idy, idx + 1, i_cons)] - f[ij_sol(idy, idx, i_cons)]);

      ul[i_cons] = f[ij_sol(idy, idx, i_cons)] + 0.5 * slope;

      slope = minmod(f[ij_sol(idy, idx + 1, i_cons)] - f[ij_sol(idy, idx, i_cons)], f[ij_sol(idy, idx + 2, i_cons)] - f[ij_sol(idy, idx + 1, i_cons)]);

      ur[i_cons] = f[ij_sol(idy, idx + 1, i_cons)] - 0.5 * slope;
    }

  } else if (dir == dir_y) {
    for (int i_cons = i_rho; i_cons <= i_ener; i_cons++) {

      slope = minmod(f[ij_sol(idy, idx, i_cons)] - f[ij_sol(idy - 1, idx, i_cons)], f[ij_sol(idy + 1, idx, i_cons)] - f[ij_sol(idy, idx, i_cons)]);

      ul[i_cons] = f[ij_sol(idy, idx, i_cons)] + 0.5 * slope;

      slope = minmod(f[ij_sol(idy + 1, idx, i_cons)] - f[ij_sol(idy, idx, i_cons)], f[ij_sol(idy + 2, idx, i_cons)] - f[ij_sol(idy + 1, idx, i_cons)]);

      ur[i_cons] = f[ij_sol(idy + 1, idx, i_cons)] - 0.5 * slope;
    }
  }
}
