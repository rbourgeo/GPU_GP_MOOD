#ifndef LINALG_H
#define LINALG_H
#include "constants.h"
#include <quadmath.h>
inline int a2l_lin(const int i, const int j, const int n)
{
  return i * n + j;
}

inline void Cholesky_Decomposition(const __float128 matrix[], int n, __float128 lower[])
{
  // Decomposing a matrix into Lower Triangular
  __float128 sum = longdouble(0);

  for (int i = 0; i <= n - 1; i++) {

    sum = longdouble(0);

    for (int k = 0; k <= i - 1; k++) {

      sum += powq(lower[a2l_lin(i, k, n)], 2);
    }

    lower[a2l_lin(i, i, n)] = sqrtq(matrix[a2l_lin(i, i, n)] - sum);

    for (int j = 1 + i; j <= n - 1; j++) {

      sum = longdouble(0);

      for (int k = 0; k <= i - 1; k++) {
        sum += lower[a2l_lin(i, k, n)] * lower[a2l_lin(j, k, n)];
      }
      lower[a2l_lin(j, i, n)] = (matrix[a2l_lin(j, i, n)] - sum) / lower[a2l_lin(i, i, n)];
    }
  }
}

inline void Reschol(const __float128 matrix[], __float128 x[], const __float128 b[], const int n)
{

  __float128* y;
  __float128* lower;

  lower = new __float128[n * n];
  y = new __float128[n];

  for (int i = 0; i <= n - 1; i++) {
    for (int j = 0; j <= n - 1; j++) {
      lower[a2l_lin(i, j, n)] = 0.0;
    }
  }

  Cholesky_Decomposition(matrix, n, lower); //Cholesky_Decomposition of matrix;

  __float128 sum;

  for (int i = 0; i <= n - 1; i++) {
    sum = longdouble(0);
    for (int j = 0; j <= i - 1; j++) {
      sum = sum + lower[a2l_lin(i, j, n)] * y[j];
    }
    y[i] = (b[i] - sum) / lower[a2l_lin(i, i, n)];
  }

  for (int i = n - 1; i >= 0; i--) {
    sum = longdouble(0);
    for (int j = i + 1; j <= n - 1; j++) {
      sum = sum + lower[a2l_lin(j, i, n)] * x[j]; //transpose
    }
    x[i] = (y[i] - sum) / lower[a2l_lin(i, i, n)];
  }

  delete y;
  delete lower;
}

inline void matrix_inversion(const __float128 matrix[], const int n, __float128 matrix_inv[])
{

  __float128 *e, *x;

  e = new __float128[n];
  x = new __float128[n];

  for (int k = 0; k <= n - 1; k++) {

    for (int i = 0; i <= n - 1; i++) {
      e[i] = longdouble(0);
    }
    e[k] = longdouble(1);

    Reschol(matrix, x, e, n);

    for (int i = 0; i <= n - 1; i++) {
      matrix_inv[a2l_lin(k, i, n)] = x[i];
    }
  }

  /*test inversion */
//  int width = 46;
//  char buf[128];

  for (int i = 0; i <= n-1; i++) {
    for (int j = 0; j <= n-1; j++) {
      __float128 sum = longdouble(0);

      for (int k = 0; k <= n-1; k++) {
        sum+=matrix_inv[a2l_lin(i, k, n)]*matrix[a2l_lin(k, j, n)];
      }
      //int n = quadmath_snprintf (buf, sizeof buf, "%+-#*.20Qe", width, sum);
      //printf ("%s\n", buf);
    }
  }

  delete e;
  delete x;
}

#endif
