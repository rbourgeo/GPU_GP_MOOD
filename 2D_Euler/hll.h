#include "constants.h"
#include "parameters.h"
#include "physics.h"

__host__ __device__ inline void HLL_flux(
  const double ul[],
  const double ur[],
  const int    dir,
  double Flux[])


  {
    bool adm = 0;
    double pressure, c, v, vml, vpl, vpr, vmr, sR, sL;
    double FL[4], FR[4];

    if (dir == dir_x) {
      v        = xvel     (          ul);
    } else {
      v        = yvel     (          ul);
    }

    pressure = pres       (          ul);
    c        = sound_speed(pressure, ul);
    vml      = v - c;
    vpl      = v + c;
    Flux_fct(ul, FL, dir, pressure, v);

    if (dir == dir_x) {
      v        = xvel     (          ur);
    } else {
      v        = yvel     (          ur);
    }

    pressure = pres       (          ur);
    c        = sound_speed(pressure, ur);
    vmr      = v - c;
    vpr      = v + c;
    Flux_fct(ur, FR, dir, pressure, v);

    sL = min(vml, vmr);
    sR = max(vpl, vpr);

    if (sL > 0.)
    {
      for (int i_cons = i_rho; i_cons <= i_ener; i_cons++)
      {
        Flux[i_cons] = FL[i_cons];
      }
    }
    else if ((sL <= 0.) && (sR >0. ))
    {
      for (int i_cons = i_rho; i_cons <= i_ener; i_cons++)
      {
        Flux[i_cons] = (1./(sR-sL))*( sR*FL[i_cons] - sL*FR[i_cons] + sL*sR*(ur[i_cons] - ul[i_cons]) );
      }
    }
    else
    {
      for (int i_cons = i_rho; i_cons <= i_ener; i_cons++)
      {
        Flux[i_cons] = FR[i_cons];
      }
    }

  }
