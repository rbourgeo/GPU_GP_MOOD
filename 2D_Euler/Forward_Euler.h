#ifndef FE_H

#include "Update.h"
#include "gp_stencil.h"
#include "parameters.h"
#include "set_equal.h"
#include "bc.h"
#include "recons_flux_computation.h"



#define FE_H

inline void Forward_Euler(
    double d_f[]
  , double d_fout[]
  , double d_fluxes_x[]
  , double d_fluxes_y[]
  , const double dx
  , const double dy
  , const double dt
  , const double zT[]
  , const int index[]
  , const double gauss_weight[])
  {

    //Call BC
    bc<<<dimGrid, dimBlock>>>(d_f);
    cudaDeviceSynchronize();
    
    //initialize fout at fin
    set_equal<<<dimGrid, dimBlock>>>(d_fout, d_f); //fout = f
    cudaDeviceSynchronize();

    /*Call the stencil routine -> fluxes*/
    ftcs<<<dimGrid, dimBlock>>>(d_f, d_fluxes_x, d_fluxes_y, dx, dy, dt, zT, index, gauss_weight);
    cudaDeviceSynchronize();

    /*Update U with the fluxes (avoid race condition)*/
     Update_A<<<dimGrid, dimBlock>>>(d_fout, d_fluxes_x, d_fluxes_y);
     cudaDeviceSynchronize();
     Update_B<<<dimGrid, dimBlock>>>(d_fout, d_fluxes_x, d_fluxes_y);
     cudaDeviceSynchronize();
     Update_C<<<dimGrid, dimBlock>>>(d_fout, d_fluxes_x, d_fluxes_y);
     cudaDeviceSynchronize();
  }

  #endif
