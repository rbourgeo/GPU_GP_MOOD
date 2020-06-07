#ifndef GP_STENCIL_H
#include "constants.h"
#include "linalg.h"
#include "parameters.h"
#define GP_STENCIL_H

/* Class GP_stencil


*/

__host__ __device__ int inline a2l_ngp(int j, int i)
{
    return i * ngp + j;
}

/*array to list index*/
int inline a2l(int j, int i)
{
    return i * nop + j;
}


class GP_stencil {

private:
public:
    /* Variables sent to device*/
    double* zT_L;
    double* zT_T;
    double* zT_R;
    double* zT_B;
    int *index_x, *index_y;

    int order, number_of_points;

    long double l;
    /* Variables not used in the simulations*/

    long double* gauss_point;

    long double *dist_kh_x, *dist_kh_y;
    long double *dist_L_ks_x, *dist_L_ks_y;
    long double *dist_T_ks_x, *dist_T_ks_y;
    long double *dist_R_ks_x, *dist_R_ks_y;
    long double *dist_B_ks_x, *dist_B_ks_y;

    long double* T_pred_L;
    long double* T_pred_T;
    long double* T_pred_R;
    long double* T_pred_B;

    long double* D2x_Tpred;
    long double* D2y_Tpred;

    long double* COV;
    long double* inv_COV;


    GP_stencil(){};
    /*constructor*/
    GP_stencil(const int ord_in, const double l_in)
    {

        order = ord_in;
        number_of_points = 2 * order - 1;
        l = l_in;

        index_x = new int[number_of_points];
        index_y = new int[number_of_points];
    };

    // /*routine for extrapolation */
    // __device__ __host__ double Extrapolate(const int idy, const int idx, const int iFace, const int gaussian_pt, const double f[])
    // {
    //
    //     int r = gaussian_pt;
    //     double result;
    //
    //     result = 0.0;
    //
    //     for (int k = 0; k <= number_of_points - 1; k++) {
    //         if (iFace == iL) {
    //             result += zT_L[a2l(k, r)] * f[ij(idy + index_y[k], idx + index_x[k])];
    //         };
    //         if (iFace == iT) {
    //             result += zT_T[a2l(k, r)] * f[ij(idy + index_y[k], idx + index_x[k])];
    //         };
    //         if (iFace == iR) {
    //             result += zT_R[a2l(k, r)] * f[ij(idy + index_y[k], idx + index_x[k])];
    //         };
    //         if (iFace == iB) {
    //             result += zT_B[a2l(k, r)] * f[ij(idy + index_y[k], idx + index_x[k])];
    //         };
    //     }
    //
    //     return result;
    // }

    __host__ size_t size_gp_stencil()
    {

        return sizeof(long double) + (2 + number_of_points) * sizeof(int) + 4 * number_of_points * sizeof(double);
        /* l + order + number of points + index_x and index_y and all 4 predictions  vector*/
    };

    /* Deletes everyting that is no usefull in the simulation */
    void init_deallocate()
    {

        delete COV;
        delete inv_COV;
        delete dist_kh_x;
        delete dist_kh_y;
        delete dist_L_ks_x;
        delete dist_L_ks_y;
        delete dist_T_ks_x;
        delete dist_T_ks_y;
        delete dist_R_ks_x;
        delete dist_R_ks_y;
        delete dist_B_ks_x;
        delete dist_B_ks_y;
        delete T_pred_L;
        delete T_pred_T;
        delete T_pred_R;
        delete T_pred_B;
        delete D2x_Tpred;
        delete D2y_Tpred;
    };

    void _deallocate()
    {

        delete index_x;
        delete index_y;
        delete zT_L;
        delete zT_T;
        delete zT_R;
        delete zT_B;
    }
    void fill_index()
    {

        int irad = 1;
        int count = 0;

        index_x[0] = 0;
        index_y[0] = 0;

        for (int k = 1; k <= number_of_points - 1; k++) {

            if (count == 0) {
                index_x[k] = -irad;
                index_y[k] = 0;
            }

            if (count == 1) {
                index_x[k] = 0;
                index_y[k] = irad;
            }

            if (count == 2) {
                index_x[k] = irad;
                index_y[k] = 0;
            }
            if (count == 3) {

                index_x[k] = 0;
                index_y[k] = -irad;
            }

            count = count + 1;

            if (count == 4) {

                count = 0;
                irad = irad + 1;
            }
        };
    };

    /* Formula for 1D intg kernel SE*/
    long double intg_kernel_SE_1D(const long double dist, const long double h)
    {

        long double a11, a12, a21, a22, a31, a32, lod, result;

        lod = l / h;

        a11 = (dist + 1.) / (sqrtl(2.) * lod);
        a12 = (dist - 1.) / (sqrtl(2.) * lod);

        a21 = -powl(dist + 1., 2) / (2. * powl(lod, 2));
        a22 = -powl(dist - 1., 2) / (2. * powl(lod, 2));

        a31 = dist / (sqrtl(2.) * lod);
        a32 = -powl(dist, 2) / (2. * powl(lod, 2));

        result = a11 * erfl(a11) + a12 * erfl(a12);
        result += ospi * (expl(a21) + expl(a22));
        result += -2 * (a31 * erfl(a31) + ospi * expl(a32));
        result *= sqrtl(pi) * powl(lod, 2);

        return result;
    }

    /* Formula for 1D pred vector*/
    long double Tpred_1D_fun(const long double delta, const long double h)
    {

        long double expm, expp, result;

        expp = (delta + 0.5) / (sqrtl(2.) * (l / h));
        expm = (delta - 0.5) / (sqrtl(2.) * (l / h));

        result = sqrtl(pi / 2) * (l / h) * (erfl(expp) - erfl(expm));

        return result;
    }
    /* Formula for 1D derivative pred vector*/
    long double Tpred_d2x_1D_fun(const long double delta, const long double h)
    {

        long double expm, expp, result, fm, fp;

        fm = delta - 0.5;
        fp = delta + 0.5;

        expm = -(powl(fm, 2)) / (2 * pow(l / h, 2));
        expp = -(powl(fp, 2)) / (2 * pow(l / h, 2));

        result = (1. / pow(l, 2)) * (fm * expl(expm) - fp * expl(expp));
        return result;
    }
    /*fill the distances, cov matrix and cov^-1*/
    void fill_dist_and_COV(long double dx, long double dy)
    {
        dist_kh_x = new long double[number_of_points * number_of_points];
        dist_kh_y = new long double[number_of_points * number_of_points];
        COV = new long double[number_of_points * number_of_points];
        inv_COV = new long double[number_of_points * number_of_points];

        for (int k = 0; k <= number_of_points - 1; k++) {
            for (int h = 0; h <= number_of_points - 1; h++) {

                dist_kh_x[a2l(k, h)] = index_x[k] - index_x[h];
                dist_kh_y[a2l(k, h)] = index_y[k] - index_y[h];
                COV[a2l(k, h)] = intg_kernel_SE_1D(dist_kh_x[a2l(k, h)], dx);
                COV[a2l(k, h)] *= intg_kernel_SE_1D(dist_kh_y[a2l(k, h)], dy);
            }
        }

        matrix_inversion(COV, number_of_points, inv_COV);

        //  Cholesky_Decomposition(COV,number_of_points);

        dist_L_ks_x = new long double[number_of_points * ngp];
        dist_L_ks_y = new long double[number_of_points * ngp];
        dist_T_ks_x = new long double[number_of_points * ngp];
        dist_T_ks_y = new long double[number_of_points * ngp];
        dist_R_ks_x = new long double[number_of_points * ngp];
        dist_R_ks_y = new long double[number_of_points * ngp];
        dist_B_ks_x = new long double[number_of_points * ngp];
        dist_B_ks_y = new long double[number_of_points * ngp];

        gauss_point = new long double[ngp * ngp];

        gauss_point[a2l_ngp(0, 0)] = 0.;

        if (ngp >= 2) {
            gauss_point[a2l_ngp(1, 0)] = -0.5 / sq3;
            gauss_point[a2l_ngp(1, 1)] = 0.5 / sq3;
        }

        if (ngp >= 3) {
            gauss_point[a2l_ngp(2, 0)] = -0.5 * sq35;
            gauss_point[a2l_ngp(2, 1)] = 0.0;
            gauss_point[a2l_ngp(2, 2)] = 0.5 * sq35;
        }

        for (int k = 0; k <= number_of_points - 1; k++) {
            for (int r = 0; r <= ngp - 1; r++) {

                dist_L_ks_x[a2l(k, r)] = index_x[k] - (-1. / 2);
                dist_L_ks_y[a2l(k, r)] = index_y[k] - gauss_point[a2l_ngp(ngp - 1, r)];

                dist_T_ks_x[a2l(k, r)] = index_x[k] - gauss_point[a2l_ngp(ngp - 1, r)];
                dist_T_ks_y[a2l(k, r)] = index_y[k] - (+1. / 2);

                dist_R_ks_x[a2l(k, r)] = index_x[k] - (+1. / 2);
                dist_R_ks_y[a2l(k, r)] = index_y[k] - gauss_point[a2l_ngp(ngp - 1, r)];

                dist_B_ks_x[a2l(k, r)] = index_x[k] - gauss_point[a2l_ngp(ngp - 1, r)];
                dist_B_ks_y[a2l(k, r)] = index_y[k] - (-1. / 2);
            }
        }

        delete gauss_point;
    }

    /*Fill the prediction vector*/
    void fill_Pred_Vectors(const long double dx, const long double dy)
    {
        T_pred_L = new long double[number_of_points * ngp];
        T_pred_T = new long double[number_of_points * ngp];
        T_pred_R = new long double[number_of_points * ngp];
        T_pred_B = new long double[number_of_points * ngp];
        D2x_Tpred = new long double[number_of_points];
        D2y_Tpred = new long double[number_of_points];

        for (int k = 0; k <= number_of_points - 1; k++) {
            for (int r = 0; r <= ngp - 1; r++) {
                T_pred_L[a2l(k, r)] = Tpred_1D_fun(dist_L_ks_x[a2l(k, r)], dx);
                T_pred_L[a2l(k, r)] *= Tpred_1D_fun(dist_L_ks_y[a2l(k, r)], dy);

                T_pred_T[a2l(k, r)] = Tpred_1D_fun(dist_T_ks_x[a2l(k, r)], dx);
                T_pred_T[a2l(k, r)] *= Tpred_1D_fun(dist_T_ks_y[a2l(k, r)], dy);

                T_pred_R[a2l(k, r)] = Tpred_1D_fun(dist_R_ks_x[a2l(k, r)], dx);
                T_pred_R[a2l(k, r)] *= Tpred_1D_fun(dist_R_ks_y[a2l(k, r)], dy);

                T_pred_B[a2l(k, r)] = Tpred_1D_fun(dist_B_ks_x[a2l(k, r)], dx);
                T_pred_B[a2l(k, r)] *= Tpred_1D_fun(dist_B_ks_y[a2l(k, r)], dy);
            }

            D2x_Tpred[k] = Tpred_d2x_1D_fun(dist_kh_x[a2l(k, 0)], dx);
            D2x_Tpred[k] = D2x_Tpred[k] * Tpred_1D_fun(dist_kh_y[a2l(k, 0)], dy);

            D2y_Tpred[k] = Tpred_d2x_1D_fun(dist_kh_y[a2l(k, 0)], dy);
            D2y_Tpred[k] = D2y_Tpred[k] * Tpred_1D_fun(dist_kh_x[a2l(k, 0)], dx);
        }
    }
    /* FIll the exrapolating vectors*/
    void fill_zTs()
    {

        zT_L = new double[number_of_points * ngp];
        zT_T = new double[number_of_points * ngp];
        zT_R = new double[number_of_points * ngp];
        zT_B = new double[number_of_points * ngp];

        double sum;

        for (int k = 0; k <= number_of_points - 1; k++) {
            for (int r = 0; r <= ngp - 1; r++) {
                //computation

                sum = 0.0;
                for (int l = 0; l <= number_of_points - 1; l++) {
                    sum += inv_COV[a2l(k, l)] * T_pred_L[a2l(l, r)];
                }

                zT_L[a2l(k, r)] = sum;

                sum = 0.0;
                for (int l = 0; l <= number_of_points - 1; l++) {
                    sum += inv_COV[a2l(k, l)] * T_pred_T[a2l(l, r)];
                }

                zT_T[a2l(k, r)] = sum;

                sum = 0.0;
                for (int l = 0; l <= number_of_points - 1; l++) {
                    sum += inv_COV[a2l(k, l)] * T_pred_R[a2l(l, r)];
                }

                zT_R[a2l(k, r)] = sum;

                sum = 0.0;
                for (int l = 0; l <= number_of_points - 1; l++) {
                    sum += inv_COV[a2l(k, l)] * T_pred_B[a2l(l, r)];
                }

                zT_B[a2l(k, r)] = sum;
            }
        }

        //normalisation
        for (int r = 0; r <= ngp - 1; r++) {

            sum = 0.0;
            for (int l = 0; l <= number_of_points - 1; l++) {
                sum += zT_L[a2l(l, r)];
            }
            for (int l = 0; l <= number_of_points - 1; l++) {
                zT_L[a2l(l, r)] = zT_L[a2l(l, r)] / sum;
            }

            sum = 0.0;
            for (int l = 0; l <= number_of_points - 1; l++) {
                sum += zT_T[a2l(l, r)];
            }
            for (int l = 0; l <= number_of_points - 1; l++) {
                zT_T[a2l(l, r)] = zT_T[a2l(l, r)] / sum;
            }

            sum = 0.0;
            for (int l = 0; l <= number_of_points - 1; l++) {
                sum += zT_R[a2l(l, r)];
            }
            for (int l = 0; l <= number_of_points - 1; l++) {
                zT_R[a2l(l, r)] = zT_R[a2l(l, r)] / sum;
            }

            sum = 0.0;
            for (int l = 0; l <= number_of_points - 1; l++) {
                sum += zT_B[a2l(l, r)];
            }
            for (int l = 0; l <= number_of_points - 1; l++) {
                zT_B[a2l(l, r)] = zT_B[a2l(l, r)] / sum;
            }
        }
    }

    /* Initalize the GP extapolation vectors*/
    void init(long double dx, long double dy)
    {
        fill_index();
        fill_dist_and_COV(dx, dy);
        fill_Pred_Vectors(dx, dy);
        fill_zTs();

        init_deallocate();
    }
};

#endif
