#ifndef EXACT_SOL_H
#include "constants.h"
#include "parameters.h"

#define EXACT_SOL_H

const int MAQUEUE = 10;
/*

Give the exact cell averaged solution in cell l,n for the gaussian advection problem

*/
__host__ __device__ inline double exact_soln(const int l, const int n, const double dx, const double dy, const double x[], const double y[])
{

    double xx, yy, intg;

    xx = x[f2c(l)];
    yy = y[f2c(n)];

    intg = (-0.0886227 * erf(-5 * dx - 10 * xx + 5.) + 0.0886227 * erf(5 * dx - 10 * xx + 5.));
    intg = intg * (-0.0886227 * erf(-5 * dy - 10 * yy + 5.) + 0.0886227 * erf(5 * dy - 10 * yy + 5.));

    intg = (intg + dx * dy) / (dx * dy);

    return intg;
}

// Compute the L1 error for the gaussian advection problem
__host__ inline double error_(double f[], double x[], double y[], const double dx, const double dy)
{

    double error = 0.0;
    double solution;

    for (int l = 1; l <= lf; l++) {
        for (int n = 1; n <= nf; n++) {
            solution = exact_soln(l, n, dx, dy, x, y);
            error += abs(f[ij(n, l)] - solution) * dx * dy;
        }
    }

    return error;
}

#endif
