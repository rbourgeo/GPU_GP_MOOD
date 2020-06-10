// TO do :: dt fully //
#ifndef PARAMETERS_H

#include "constants.h"
#include <quadmath.h>


#define PARAMETERS_H

const int ngc = 4; // Number of ghost cells

const int lf = 512 - 2*ngc; // number of x mesh
const int nf = lf; // number of y mesh

const int ngp = 1; // number of gaussian point on edges

const __float128 Lx_16 = 1;
const __float128 Ly_16 = 1;

const double CFL = 0.4;
const double tmax = 0.8;

const int le = lf + 2 * ngc; // total number of cells including ghosts
const int ne = nf + 2 * ngc;

const int space_method = PLM;
const int time_method = SSP_RK2;

const int IC_type = RP_3;
const int BC_type = Neumann;

const int Mord = 3; /*Order of gp extrapolation*/
const int ell_over_dx = 0;
const __float128 ell = longdouble(1)/10;

const double ax = 0.0; // speed
const double ay = 0.0; // speed

const int nop = 2*Mord-1; //number of points in the gp stencil

//Kernel parameters
const dim3 dimBlock(16, 16, 1);
const dim3 dimGrid(le / 16, ne / 16, 1);

/*
* ij : Computes the list coordinate from matrix-like entries
* i : The y in f90 format: ghost cells from [1-ngc:0] solution from [1:lf], ghost cells from [lf+1:lf+ngc]
* y : The x coordinate, same
*/

inline __host__ __device__ int f2c(const int i_for)
{
    return i_for + ngc - 1; // c index from fortran indexs
}

inline __host__ __device__ int c2f(const int i_c)
{
    return i_c - ngc + 1; // fortran index from c indexs
}

inline __host__ __device__ int ij(const int i_for, const int j_for)
{
    int i_c = f2c(i_for); // c index from fortran indexs
    int j_c = f2c(j_for);

    return i_c * le + j_c;
}

inline __host__ __device__ int ij_sol(const int i_for, const int j_for, const int iCons)
{
    int i_c = f2c(i_for); // c index from fortran indexs
    int j_c = f2c(j_for);

    return i_c*(le*4) + j_c*(4) + iCons;
}

#endif
