#ifndef LINALG_H
#define LINALG_H

inline int a2l_lin(const int i, const int j, const int n)
{
    return i * n + j;
}

inline void Cholesky_Decomposition(const long double matrix[], int n, long double lower[])
{
    // Decomposing a matrix into Lower Triangular
    long double sum = 0.0;

    for (int i = 0; i <= n - 1; i++) {

        sum = 0.0;

        for (int k = 0; k <= i - 1; k++) {

            sum += powl(lower[a2l_lin(i, k, n)], 2);
        }

        lower[a2l_lin(i, i, n)] = sqrtl(matrix[a2l_lin(i, i, n)] - sum);

        for (int j = 1 + i; j <= n - 1; j++) {

            sum = 0.0;

            for (int k = 0; k <= i - 1; k++) {
                sum += lower[a2l_lin(i, k, n)] * lower[a2l_lin(j, k, n)];
            }
            lower[a2l_lin(j, i, n)] = (matrix[a2l_lin(j, i, n)] - sum) / lower[a2l_lin(i, i, n)];
        }
    }
}

inline void Reschol(const long double matrix[], long double x[], const long double b[], const int n)
{

    long double* y;
    long double* lower;

    lower = new long double[n * n];
    y = new long double[n];

    for (int i = 0; i <= n - 1; i++) {
        for (int j = 0; j <= n - 1; j++) {
            lower[a2l_lin(i, j, n)] = 0.0;
        }
    }

    Cholesky_Decomposition(matrix, n, lower); //Cholesky_Decomposition of matrix;

    long double sum;

    for (int i = 0; i <= n - 1; i++) {
        sum = 0.;
        for (int j = 0; j <= i - 1; j++) {
            sum = sum + lower[a2l_lin(i, j, n)] * y[j];
        }
        y[i] = (b[i] - sum) / lower[a2l_lin(i, i, n)];
    }

    for (int i = n - 1; i >= 0; i--) {
        sum = 0.;
        for (int j = i + 1; j <= n - 1; j++) {
            sum = sum + lower[a2l_lin(j, i, n)] * x[j]; //transpose
        }
        x[i] = (y[i] - sum) / lower[a2l_lin(i, i, n)];
    }

    delete y;
    delete lower;
}

inline void matrix_inversion(const long double matrix[], const int n, long double matrix_inv[])
{

    long double *e, *x;

    e = new long double[n];
    x = new long double[n];

    for (int k = 0; k <= n - 1; k++) {

        for (int i = 0; i <= n - 1; i++) {
            e[i] = 0.;
        }
        e[k] = 1.;

        Reschol(matrix, x, e, n);

        for (int i = 0; i <= n - 1; i++) {
            matrix_inv[a2l_lin(k, i, n)] = x[i];
        }
    }

    delete e;
    delete x;
}

#endif
