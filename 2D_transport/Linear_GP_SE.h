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
                                        const double f[])
  {

    double result;

    result = 0.0;

    for (int k = 0; k <= nop - 1; k++) {
      result += zT[iFace*nop*ngp + gaussian_pt*nop + k] * f[ij(idy + index[dir_y*nop + k], idx + index[dir_x*nop + k])];
    }
    return result;

  }

  // }

  __host__ __device__ inline void Linear_GP_SE_(double& ul,
    double& ur,
    const int idy,
    const int idx,
    const double f[],
    const int dir,
    const double zT[],
    const int index[],
    const int gaussian_point)
    {

      if (dir == dir_x) {

        ul = Extrapolate(idy, idx, iR, zT, index, gaussian_point, f);

        ur = Extrapolate(idy, idx + 1, iL, zT, index, gaussian_point, f);

      } else if (dir == dir_y) {
        ul = Extrapolate(idy, idx, iT, zT, index, gaussian_point, f);

        ur = Extrapolate(idy + 1, idx, iB, zT, index, gaussian_point, f);
      }
    }
