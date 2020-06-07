#include "Forward_Euler.h"
#include "Update.h"
#include "gp_stencil.h"
#include "parameters.h"
#include "set_equal.h"

inline void SSP_RK2_(double d_f[]
                        , double d_fout[]
                        , double d_f1[]
                        , double d_fluxes_x[]
                        , double d_fluxes_y[]
                        , const double dx
                        , const double dy
                        , const double dt
                        , const double zT[]
                        , const int index[]
                        , const double gauss_weight[])
{
    // We save initial state
    set_equal<<<dimGrid, dimBlock>>>(d_f1, d_f); //d_f1 = d_f
    cudaDeviceSynchronize();

    // Evolve 1 time
    Forward_Euler(d_f, d_fout, d_fluxes_x, d_fluxes_y, dx, dy, dt, zT, index, gauss_weight);

    // Save half step
    set_equal<<<dimGrid, dimBlock>>>(d_f, d_fout); //d_f = d_fout
    cudaDeviceSynchronize();

    // Evolve 2nd time
    Forward_Euler(d_f, d_fout, d_fluxes_x, d_fluxes_y, dx, dy, dt, zT, index, gauss_weight);

    // Take the averages
    set_comb_lin<<<dimGrid, dimBlock>>>(d_fout, d_f1, d_fout, 0.5, 0.5); //d_f = d_fout
    cudaDeviceSynchronize();
}
