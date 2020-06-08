#include "constants.h"
#include "parameters.h"

/* First order method, we just copy the cell averaged values
* outputs: ul and ur, high order pointwise face centered values
* inputs : idx idy, coordinates of the interface
* f: solution
* dir: direction of reconstructon */
__host__ __device__ inline void FOG_(double ul[], double ur[], const int idy, const int idx, const double f[], const int dir)
{

    if (dir == dir_x) {
      for (int i_cons = i_rho; i_cons <= i_ener; i_cons++) {

        ul[i_cons] = f[ij_sol(idy, idx    , i_cons)];
        ur[i_cons] = f[ij_sol(idy, idx + 1, i_cons)];
      }
    } else if (dir == dir_y) {
      for (int i_cons = i_rho; i_cons <= i_ener; i_cons++) {

        ul[i_cons] = f[ij_sol(idy,     idx, i_cons)];
        ur[i_cons] = f[ij_sol(idy + 1, idx, i_cons)];
      }
    }
}
