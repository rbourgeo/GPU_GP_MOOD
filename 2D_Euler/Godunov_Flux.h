#include "constants.h"
#include "parameters.h"

/*
* Godunov_Flux : Computes the Numerical flux from left and righ states
* ul : left state
* ur : right state
*/
__host__ __device__ inline void Godunov_Flux(const double ul[], const double ur[], const int dir, double flux[])
{

  for (int i_cons = i_rho; i_cons <= i_ener; i_cons++) {

    if (dir == dir_x) {

        if (ax >= 0.) {
            flux[i_cons] = ax * ul[i_cons];
        } else {
            flux[i_cons] = ax * ur[i_cons];
        }

    } else if (dir == dir_y) {

        if (ay >= 0.) {
            flux[i_cons] = ay * ul[i_cons];
        } else {
            flux[i_cons] = ay * ur[i_cons];
        }
    }
    
  }

}
