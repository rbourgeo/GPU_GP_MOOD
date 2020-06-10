#ifndef FE_H

#include "Update.h"
#include "gp_stencil.h"
#include "parameters.h"
#include "set_equal.h"
#include "bc.h"
#include "recons_flux_computation.h"



#define FE_H

__global__ void fill_int(int out[], const int i)
{
  int tidx = c2f(threadIdx.x + blockIdx.x * blockDim.x);
  int tidy = c2f(threadIdx.y + blockIdx.y * blockDim.y);

  if (tidx >= 1 - ngc && tidx <= lf + ngc) {
    if (tidy >= 1 - ngc && tidy <= nf + ngc) {
      out[ij(tidy,tidx)] = i;
    }
  }

}

__global__ void fill_bool(bool out[], const bool i)
{
  int tidx = c2f(threadIdx.x + blockIdx.x * blockDim.x);
  int tidy = c2f(threadIdx.y + blockIdx.y * blockDim.y);

  if (tidx >= 1 - ngc && tidx <= lf + ngc) {
    if (tidy >= 1 - ngc && tidy <= nf + ngc) {
      out[ij(tidy,tidx)] = i;
    }
  }

}

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
    // 
    // int  CellGPO[le*ne];
    // bool DetCell[le*ne], DetFace_x[le*ne], DetCell_y[le*ne];
    //
    // fill_bool<<<dimGrid, dimBlock>>>(DetCell,1);
    // cudaDeviceSynchronize();
    //
    // fill_bool<<<dimGrid, dimBlock>>>(DetFace_x,1);
    // cudaDeviceSynchronize();
    //
    // fill_bool<<<dimGrid, dimBlock>>>(DetFace_y,1);
    // cudaDeviceSynchronize();







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
