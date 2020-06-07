/*ftcs : computes the fluxes at the interfaces
*f : solutions
*f: fluxes_x/y : fluxes on the x/y interfaces
*dx, dy, dt : space/time steps */
#include "gp_stencil.h"
__global__ void ftcs(double d_f[]
                        , double d_fluxes_x[]
                        , double d_fluxes_y[]
                        , const double dx
                        , const double dy
                        , const double dt
                        , const double zT[]
                        , const int index[]
                        , const double gauss_weight[]);
