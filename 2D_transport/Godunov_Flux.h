#include "constants.h"
#include "parameters.h"

/*
* Godunov_Flux : Computes the Numerical flux from left and righ states
* ul : left state
* ur : right state
*/
__host__ __device__ inline double Godunov_Flux(const double ul, const double ur, const int dir)
{
    if (dir == dir_x) {

        if (ax >= 0.) {
            return ax * ul;
        } else {
            return ax * ur;
        }

    } else if (dir == dir_y) {

        if (ay >= 0.) {
            return ay * ul;
        } else {
            return ay * ur;
        }
    }
}
