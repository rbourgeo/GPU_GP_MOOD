#include "constants.h"
#include "parameters.h"

/* First order method, we just copy the cell averaged values
* outputs: ul and ur, high order pointwise face centered values
* inputs : idx idy, coordinates of the interface
* f: solution
* dir: direction of reconstructon */
__host__ __device__ inline void FOG_(double& ul, double& ur, const int idy, const int idx, const double f[], const int dir)
{

    if (dir == dir_x) {
        ul = f[ij(idy, idx)];
        ur = f[ij(idy, idx + 1)];
    } else if (dir == dir_y) {
        ul = f[ij(idy, idx)];
        ur = f[ij(idy + 1, idx)];
    }
}
