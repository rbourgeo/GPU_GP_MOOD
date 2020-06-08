#include "Godunov_Flux.h"
#include "parameters.h"
#include "constants.h"
#include "reconstruction.h"

__global__ void ftcs(const double d_f[]
  , double d_fluxes_x[]
  , double d_fluxes_y[]
  , const double dx
  , const double dy
  , const double dt
  , const double d_zT[]
  , const int d_index[]
  , const double d_gauss_weight[])
  {
    int tidx = c2f(threadIdx.x + blockIdx.x * blockDim.x);
    int tidy = c2f(threadIdx.y + blockIdx.y * blockDim.y);

      double ul[4], ur[4], ub[4], ut[4];
      double fx[4], fy[4];


    //
    //
    if (tidx >= 0 && tidx <= lf) // from 1's left intf to lf's right interface
    {
      if (tidy >= 0 && tidy <= nf) {
        {

          for (int i_cons = i_rho; i_cons <= i_ener; i_cons++) {
            d_fluxes_x[ij_sol(tidy, tidx, i_cons)] = 0.0;
            d_fluxes_y[ij_sol(tidy, tidx, i_cons)] = 0.0;
          }

          for (int r = 0; r <= ngp - 1; r++) {

            /* high order reconstruction*/
             reconstruction(ul, ur, tidy, tidx, d_f, dir_x, d_zT, d_index, r);
             reconstruction(ub, ut, tidy, tidx, d_f, dir_y, d_zT, d_index, r);

            /*Fluxes  computation*/

              Godunov_Flux(ul, ur, dir_x, fx);
              Godunov_Flux(ub, ut, dir_y, fy);


            /*Fluxes storage*/
            for (int i_cons = i_rho; i_cons <= i_ener; i_cons++) {

               d_fluxes_x[ij_sol(tidy, tidx, i_cons)] += fx[i_cons] * d_gauss_weight[a2l_ngp(ngp - 1, r)];
               d_fluxes_y[ij_sol(tidy, tidx, i_cons)] += fy[i_cons] * d_gauss_weight[a2l_ngp(ngp - 1, r)];

            }
          }

          for (int i_cons = i_rho; i_cons <= i_ener; i_cons++) {

             d_fluxes_x[ij_sol(tidy, tidx, i_cons)] *= dt / dx;
             d_fluxes_y[ij_sol(tidy, tidx, i_cons)] *= dt / dy;

            //d_fluxes_x[ij_sol(tidy, tidx, i_cons)] = 0.0;
            //d_fluxes_y[ij_sol(tidy, tidx, i_cons)] = 0.0;
          }
        }

      }
    }

  }
