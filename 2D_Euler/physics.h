#ifndef PHYSICS_H
#include "constants.h"
#include "parameters.h"

#define PHYSICS_H

/* pres : computes and returns the internal pressure */
__host__ __device__ inline double pres(const double u[])
{
    return (gr_gamma - 1.0f) * (u[i_ener] - 0.5 * (u[i_momx] * u[i_momx] + u[i_momy] * u[i_momy]) / u[i_rho]);
};

/* xvel : computes and returns the x velocity */
__host__ __device__ inline double xvel(const double u[])
{
    return u[i_momx] / u[i_rho];
};

/* yvel : computes and returns the y velocity */
__host__ __device__ inline double yvel(const double u[])
{
    return u[i_momy] / u[i_rho];
};

/* Computes and returns the sound speed
 When we need SS, we usually have computed the pressure before */
__host__ __device__ inline double sound_speed(const double pressure, const double u[])
{
    return sqrt(gr_gamma * pressure / u[i_rho]);
};


/* Converts the State from p to c*/
__host__ __device__ inline void primitive_to_conservative(double u[])
{

    double vx, vy, pressure;

    vx = u[i_momx];
    vy = u[i_momy];
    pressure = u[i_ener];

    u[i_momx] = u[i_rho] * vx;
    u[i_momy] = u[i_rho] * vy;
    u[i_ener] = ((pressure / (gr_gamma - 1.0f)) + 0.5 * u[i_rho] * (vx * vx + vy * vy) );

    /* u(ener)= ( (v(4)/(y-1.)) + 0.5*v(1)*(v(2)**2+v(3)**2) ) */
}
/* Converts the State from c to p*/

__host__ __device__ inline void conservative_to_primitive(double u[])
{

    double vx, vy, pressure;
    vx = xvel(u);
    vy = yvel(u);
    pressure = pres(u);

    u[i_momx] = vx;
    u[i_momy] = vy;
    u[i_ener] = pressure;
};

__host__ __device__ inline void Flux_fct(const double u[],
                                               double f[],
                                               const int dir,
                                               const double pressure,
                                               const double v)
{
if (dir == dir_x) {
  f[i_rho ] = u[i_momx];
  f[i_momx] = u[i_rho ]*v*v + pressure;
  f[i_momy] = u[i_momy]*v;
  f[i_ener] = v*(u[i_ener]  + pressure);
} else {
  f[i_rho ] = u[i_momy];
  f[i_momx] = u[i_momx]*v;
  f[i_momy] = u[i_rho ]*v*v + pressure;
  f[i_ener] = v*(u[i_ener]  + pressure);
}

}


#endif
