#include "exact_sol.h"
#include "init.h"
#include "parameters.h"
#include "physics.h"

__global__ void initialize(double f[], double x[], double y[], const double dx, const double dy)
{
  int tidx = c2f(threadIdx.x + blockIdx.x * blockDim.x);
  int tidy = c2f(threadIdx.y + blockIdx.y * blockDim.y);

  double xt = (double(tidx) - 0.5f) * dx; //centers of cells
  double yt = (double(tidy) - 0.5f) * dy;

  double u[4];

  if (tidx <= lf + ngc) {
    x[f2c(tidx)] = xt;
  }

  if (tidy <= nf + ngc) {
    y[f2c(tidy)] = yt;
  }

  if (tidx <= lf + ngc) {
    if (tidy <= nf + ngc) {

      if (IC_type == Lin_gauss) {

        u[i_rho ] = exact_soln(tidx, tidy, dx, dy, x, y);
        u[i_momx] = ax;
        u[i_momy] = ay;
        u[i_ener] = 1.0/gr_gamma;

      }

      if (IC_type == Sod) {

        if (yt < 0.5) {
          u[i_rho ] = 1.;
          u[i_momx] = 0.;
          u[i_momy] = 0.;
          u[i_ener] = 1.;
        }
        else{
          u[i_rho ] = 0.125;
          u[i_momx] = 0.;
          u[i_momy] = 0.;
          u[i_ener] = 0.1;
        }
      }

      if (IC_type == RP_3) {


        if ((xt<=4./5)&&(yt<=4./5)) { u[i_rho] = 0.138 ; u[i_momx] = 1.206 ; u[i_momy] = 1.206  ; u[i_ener] = 0.029 ;}
        if ((xt>=4./5)&&(yt<=4./5)) { u[i_rho] = 0.5323; u[i_momx] = 0.0   ; u[i_momy] = 1.206  ; u[i_ener] = 0.3; }
        if ((xt<=4./5)&&(yt>=4./5)) { u[i_rho] = 0.5323; u[i_momx] = 1.206 ; u[i_momy] = 0.0    ; u[i_ener] = 0.3; }
        if ((xt>=4./5)&&(yt>=4./5)) { u[i_rho] = 1.5   ; u[i_momx] = 0.0   ; u[i_momy] = 0.0    ; u[i_ener] = 1.5 ;}
      }

      primitive_to_conservative(u);

      for (int i_cons = i_rho; i_cons <= i_ener; i_cons++){
        f[ij_sol(tidy, tidx, i_cons)] = u[i_cons];
      }


    }
  }
}
