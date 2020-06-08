#include "FOG.h"
#include "Linear_GP_SE.h"
#include "PLM.h"
#include "constants.h"
#include "gp_stencil.h"
#include "parameters.h"

/* Fill ul and ur with pointwise high order reconstruction of u
 ur and ul are outputs
 idy, idx are the coordinates
 f[] is the solution
 dir is the direction of reconstruction, x or y */

__host__ __device__ inline void reconstruction(double ul[],
                                               double ur[],
                                              const int idy,
                                              const int idx,
                                              const double f[],
                                              const int dir,
                                              const double zT[],
                                              const int index[],
                                              const int gaussian_pt)
{

    if (space_method == FOG) {

        FOG_(ul, ur, idy, idx, f, dir);

    } else if (space_method == PLM) {

        PLM_(ul, ur, idy, idx, f, dir);

    } else if (space_method == Linear_GP_SE) {

        Linear_GP_SE_(ul, ur, idy, idx, f, dir, zT, index, gaussian_pt);

    } else {
        //std::cout<< "space method not programmed"<<std::endl;
    }
}
