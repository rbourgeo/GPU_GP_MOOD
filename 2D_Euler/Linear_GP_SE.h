#include "constants.h"
#include "gp_stencil.h"
#include "parameters.h"

/* Linear GP method
outputs: ul and ur, high order pointwise face centered values
inputs : idx idy, coordinates of the interface
f: solution
dir: direction of reconstructon
GP : object of the class GP_stencil*/

/*routine for extrapolation */
__device__ __host__ double Extrapolate(const int idy,
  const int idx,
  const int iFace,
  const double zT[],
  const int index[],
  const int gaussian_pt,
  const double f[],
  const int i_cons)
  {

    double result = 0.0;

    for (int k = 0; k <= nop - 1; k++) {
      result += zT[iFace*nop*ngp + gaussian_pt*nop + k] * f[ij_sol(idy + index[dir_y*nop + k], idx + index[dir_x*nop + k], i_cons)];
    }

    return result;

  }

  // }

  __host__ __device__ inline void Linear_GP_SE_(double ul[],
    double ur[],
    const int idy,
    const int idx,
    const double f[],
    const int dir,
    const double zT[],
    const int index[],
    const int gaussian_point)
    {


      if (dir == dir_x) {
        for (int i_cons = i_rho; i_cons <= i_ener; i_cons++) {

          ul[i_cons] = Extrapolate(idy, idx, iR, zT, index, gaussian_point, f, i_cons);
        }

        for (int i_cons = i_rho; i_cons <= i_ener; i_cons++) {

          ur[i_cons] = Extrapolate(idy, idx + 1, iL, zT, index, gaussian_point, f, i_cons);
        }

      } else if (dir == dir_y) {
        for (int i_cons = i_rho; i_cons <= i_ener; i_cons++) {

          ul[i_cons] = Extrapolate(idy, idx, iT, zT, index, gaussian_point, f, i_cons);
        }

        for (int i_cons = i_rho; i_cons <= i_ener; i_cons++) {

          ur[i_cons] = Extrapolate(idy + 1, idx, iB, zT, index, gaussian_point, f, i_cons);
        }
      }

    }
