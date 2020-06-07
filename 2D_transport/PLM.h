#include "constants.h"
#include "minmod.h"
#include "parameters.h"

/* PLM / MUSCL 2nd order methods
 outputs: ul and ur, high order pointwise face centered values
 inputs : idx idy, coordinates of the interface
 f: solution
 dir: direction of reconstructon */
__host__ __device__ inline void PLM_(double& ul, double& ur, const int idy, const int idx, const double f[], const int dir)
{

    double slope;
    if (dir == dir_x) {

        slope = minmod(f[ij(idy, idx)] - f[ij(idy, idx - 1)], f[ij(idy, idx + 1)] - f[ij(idy, idx)]);

        ul = f[ij(idy, idx)] + 0.5 * slope;

        slope = minmod(f[ij(idy, idx + 1)] - f[ij(idy, idx)], f[ij(idy, idx + 2)] - f[ij(idy, idx + 1)]);

        ur = f[ij(idy, idx + 1)] - 0.5 * slope;

    } else if (dir == dir_y) {
        slope = minmod(f[ij(idy, idx)] - f[ij(idy - 1, idx)], f[ij(idy + 1, idx)] - f[ij(idy, idx)]);

        ul = f[ij(idy, idx)] + 0.5 * slope;

        slope = minmod(f[ij(idy + 1, idx)] - f[ij(idy, idx)], f[ij(idy + 2, idx)] - f[ij(idy + 1, idx)]);

        ur = f[ij(idy + 1, idx)] - 0.5 * slope;
    }
}
