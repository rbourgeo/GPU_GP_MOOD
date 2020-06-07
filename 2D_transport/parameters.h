
#ifndef PARAMETERS_H

#include "constants.h"

#define PARAMETERS_H

const int ngc = 6; // Number of ghost cells

const int lf = 256 - 2*ngc; // number of x mesh
const int nf = lf; // number of y mesh

const int ngp = 3; // number of gaussian point on edges

const double Lx = 1.0;
const double Ly = 1.0;

const double CFL = 0.4;
const double tmax = 1.0;

const int le = lf + 2 * ngc; // total number of cells including ghosts
const int ne = nf + 2 * ngc;

const int space_method = Linear_GP_SE;
const int time_method = SSP_RK3;

const int Mord = 5; /*Order of gp extrapolation*/
const int ell_over_dx = 0;
const double ell = 0.1;

const double ax = 1.0; // speed
const double ay = 1.0; // speed

const int nop = 2*Mord-1; //number of points in the gp stencil

//Kernel parameters
const dim3 dimBlock(16, 16, 1);
const dim3 dimGrid(le / 16, ne / 16, 1);

/*
* ij : Computes the list coordinate from matrix-like entries
* i : The y in f90 format: ghost cells from [1-ngc:0] solution from [1:lf], ghost cells from [lf+1:lf+ngc]
* y : The x coordinate, same
*/

inline __host__ __device__ int f2c(int i_for)
{
    return i_for + ngc - 1; // c index from fortran indexs
}

inline __host__ __device__ int c2f(int i_c)
{
    return i_c - ngc + 1; // fortran index from c indexs
}

inline __host__ __device__ int ij(int i_for, int j_for)
{
    int i_c = f2c(i_for); // c index from fortran indexs
    int j_c = f2c(j_for);

    return i_c * le + j_c;
}

#endif
