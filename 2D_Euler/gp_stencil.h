#ifndef GP_STENCIL_H
#include "constants.h"
#include "linalg.h"
#include "parameters.h"
#include <quadmath.h>
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
  __float128* zT_L;
  __float128* zT_T;
  __float128* zT_R;
  __float128* zT_B;
  int *index_x, *index_y;

  int order, number_of_points;

  __float128 l;
  /* Variables not used in the simulations*/

  __float128* gauss_point;

  __float128 *dist_kh_x, *dist_kh_y;
  __float128 *dist_L_ks_x, *dist_L_ks_y;
  __float128 *dist_T_ks_x, *dist_T_ks_y;
  __float128 *dist_R_ks_x, *dist_R_ks_y;
  __float128 *dist_B_ks_x, *dist_B_ks_y;

  __float128* T_pred_L;
  __float128* T_pred_T;
  __float128* T_pred_R;
  __float128* T_pred_B;

  __float128* D2x_Tpred;
  __float128* D2y_Tpred;

  __float128* COV;
  __float128* inv_COV;


  GP_stencil(){};
  /*constructor*/
  GP_stencil(const int ord_in, const __float128 l_in)
  {

    order = ord_in;
    number_of_points = 2 * order - 1;
    l = l_in;

    index_x = new int[number_of_points];
    index_y = new int[number_of_points];
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
  __float128 intg_kernel_SE_1D(const __float128 dist, const __float128 h)
  {

    __float128 a11, a12, a21, a22, a31, a32, lod, result;

    lod = l / h;
    a11 = (dist + longdouble(1)) / (sqrtq(longdouble(2)) * lod);
    a12 = (dist - longdouble(1)) / (sqrtq(longdouble(2)) * lod);

    a21 = -powq(dist + longdouble(1), 2) / (longdouble(2) * powq(lod, 2));
    a22 = -powq(dist - longdouble(1), 2) / (longdouble(2) * powq(lod, 2));

    a31 = dist / (sqrtq(longdouble(2)) * lod);
    a32 = -powq(dist, 2) / (longdouble(2) * powq(lod, 2));



    result = a11 * erfq(a11) + a12 * erfq(a12);
    result = result + ospi * (expq(a21) + expq(a22));
    result = result -longdouble(2) * (a31 * erfq(a31) + ospi * expq(a32));
    result = result*sqrtq(pi) * powq(lod, 2);

    return result;
  }

  /* Formula for 1D pred vector*/
  __float128 Tpred_1D_fun(const __float128 delta, const __float128 h)
  {

    __float128 expm, expp, result;

    expp = (delta + 0.5) / (sqrtq(longdouble(2)) * (l / h));
    expm = (delta - 0.5) / (sqrtq(longdouble(2)) * (l / h));

    result = sqrtq(pi / 2) * (l / h) * (erfq(expp) - erfq(expm));

    return result;
  }
  /* Formula for 1D derivative pred vector*/
  __float128 Tpred_d2x_1D_fun(const __float128 delta, const __float128 h)
  {

    __float128 expm, expp, result, fm, fp;

    fm = delta - 0.5;
    fp = delta + 0.5;

    expm = -(powq(fm, 2)) / (2 * powq(l / h, 2));
    expp = -(powq(fp, 2)) / (2 * powq(l / h, 2));

    result = (longdouble(1) / powq(l, 2)) * (fm * expq(expm) - fp * expq(expp));
    return result;
  }

  inline __float128 longdouble(const int i)
  {
    __float128 result;

    result = i;
    return result;
  }
  /*fill the distances, cov matrix and cov^-1*/
  void fill_dist_and_COV(__float128 dx, __float128 dy)
  {
    dist_kh_x = new __float128[number_of_points * number_of_points];
    dist_kh_y = new __float128[number_of_points * number_of_points];
    COV     = new __float128[number_of_points * number_of_points];
    inv_COV = new __float128[number_of_points * number_of_points];

    for (int k = 0; k <= number_of_points - 1; k++) {
      for (int h = 0; h <= number_of_points - 1; h++) {
  
        dist_kh_x[a2l(k, h)] = longdouble(index_x[k]) - longdouble(index_x[h]);
        dist_kh_y[a2l(k, h)] = longdouble(index_y[k]) - longdouble(index_y[h]);
        COV[a2l(k, h)]  = intg_kernel_SE_1D(dist_kh_x[a2l(k, h)], dx);
        COV[a2l(k, h)] = COV[a2l(k, h)]*intg_kernel_SE_1D(dist_kh_y[a2l(k, h)], dy);


      }
    }

    matrix_inversion(COV, number_of_points, inv_COV);

    //  Cholesky_Decomposition(COV,number_of_points);

    dist_L_ks_x = new __float128[number_of_points * ngp];
    dist_L_ks_y = new __float128[number_of_points * ngp];
    dist_T_ks_x = new __float128[number_of_points * ngp];
    dist_T_ks_y = new __float128[number_of_points * ngp];
    dist_R_ks_x = new __float128[number_of_points * ngp];
    dist_R_ks_y = new __float128[number_of_points * ngp];
    dist_B_ks_x = new __float128[number_of_points * ngp];
    dist_B_ks_y = new __float128[number_of_points * ngp];

    gauss_point = new __float128[ngp * ngp];

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

        dist_L_ks_x[a2l(k, r)] = longdouble(index_x[k]) - (-longdouble(1) / 2);
        dist_L_ks_y[a2l(k, r)] = longdouble(index_y[k]) - gauss_point[a2l_ngp(ngp - 1, r)];

        dist_T_ks_x[a2l(k, r)] = longdouble(index_x[k]) - gauss_point[a2l_ngp(ngp - 1, r)];
        dist_T_ks_y[a2l(k, r)] = longdouble(index_y[k]) - (+longdouble(1) / 2);

        dist_R_ks_x[a2l(k, r)] = longdouble(index_x[k]) - (+longdouble(1) / 2);
        dist_R_ks_y[a2l(k, r)] = longdouble(index_y[k]) - gauss_point[a2l_ngp(ngp - 1, r)];

        dist_B_ks_x[a2l(k, r)] = longdouble(index_x[k]) - gauss_point[a2l_ngp(ngp - 1, r)];
        dist_B_ks_y[a2l(k, r)] = longdouble(index_y[k]) - (-longdouble(1) / 2);
      }
    }

    delete gauss_point;
  }

  /*Fill the prediction vector*/
  void fill_Pred_Vectors(const __float128 dx, const __float128 dy)
  {
    T_pred_L = new __float128[number_of_points * ngp];
    T_pred_T = new __float128[number_of_points * ngp];
    T_pred_R = new __float128[number_of_points * ngp];
    T_pred_B = new __float128[number_of_points * ngp];
    D2x_Tpred = new __float128[number_of_points];
    D2y_Tpred = new __float128[number_of_points];

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

    zT_L = new __float128[number_of_points * ngp];
    zT_T = new __float128[number_of_points * ngp];
    zT_R = new __float128[number_of_points * ngp];
    zT_B = new __float128[number_of_points * ngp];

    __float128 sum;

    for (int k = 0; k <= number_of_points - 1; k++) {
      for (int r = 0; r <= ngp - 1; r++) {
        //computation

        sum = longdouble(0);
        for (int l = 0; l <= number_of_points - 1; l++) {
          sum += inv_COV[a2l(k,l )] * T_pred_L[a2l(l, r)];
        }

        zT_L[a2l(k, r)] = sum;

        sum = longdouble(0);
        for (int l = 0; l <= number_of_points - 1; l++) {
          sum += inv_COV[a2l(k, l)] * T_pred_T[a2l(l, r)];
        }

        zT_T[a2l(k, r)] = sum;

        sum = longdouble(0);
        for (int l = 0; l <= number_of_points - 1; l++) {
          sum += inv_COV[a2l(k, l)] * T_pred_R[a2l(l, r)];
        }

        zT_R[a2l(k, r)] = sum;

        sum = longdouble(0);
        for (int l = 0; l <= number_of_points - 1; l++) {
          sum += inv_COV[a2l(k, l)] * T_pred_B[a2l(l, r)];
        }

        zT_B[a2l(k, r)] = sum;
      }
    }
    std::cout.precision(17);


    /*normalisation*/
    for (int r = 0; r <= ngp - 1; r++) {

      sum = longdouble(0);
      for (int l = 0; l <= number_of_points - 1; l++) {
        sum += zT_L[a2l(l, r)];
      }
      for (int l = 0; l <= number_of_points - 1; l++) {
        zT_L[a2l(l, r)] = zT_L[a2l(l, r)] / sum;
      }

      sum = longdouble(0);
      for (int l = 0; l <= number_of_points - 1; l++) {
        sum += zT_T[a2l(l, r)];
      }
      for (int l = 0; l <= number_of_points - 1; l++) {
        zT_T[a2l(l, r)] = zT_T[a2l(l, r)] / sum;
      }

      sum = longdouble(0);
      for (int l = 0; l <= number_of_points - 1; l++) {
        sum += zT_R[a2l(l, r)];
      }
      for (int l = 0; l <= number_of_points - 1; l++) {
        zT_R[a2l(l, r)] = zT_R[a2l(l, r)] / sum;
      }

      sum = longdouble(0);
      for (int l = 0; l <= number_of_points - 1; l++) {
        sum += zT_B[a2l(l, r)];
      }
      for (int l = 0; l <= number_of_points - 1; l++) {
        zT_B[a2l(l, r)] = zT_B[a2l(l, r)] / sum;
      }
    }
  }

  /* Initalize the GP extapolation vectors*/
  void init(__float128 dx, __float128 dy)
  {
    fill_index();
    fill_dist_and_COV(dx, dy);
    fill_Pred_Vectors(dx, dy);
    fill_zTs();

    // int width = 46;
    // char buf[128];
    // std::cout<<"L,r=1" << std::endl;
    //
    // for (int k = 0; k <= 2*3-2; k++) {
    //   int n = quadmath_snprintf (buf, sizeof buf, "%+-#*.30Qe", width, zT_L[a2l(k, 0)]);
    //   printf ("%s\n", buf);
    // }
    //
    //
    //
    // std::cout<<"T,r=1" << std::endl;
    //
    // for (int k = 0; k <= 2*3-2; k++) {
    //   int n = quadmath_snprintf (buf, sizeof buf, "%+-#*.30Qe", width, zT_T[a2l(k, 0)]);
    //   printf ("%s\n", buf);
    // }
    //
    //
    // std::cout<<"R,r=1" << std::endl;
    //
    // for (int k = 0; k <= 2*3-2; k++) {
    //   int n = quadmath_snprintf (buf, sizeof buf, "%+-#*.30Qe", width, zT_R[a2l(k, 0)]);
    //   printf ("%s\n", buf);
    //
    // }
    //
    //
    // std::cout<<"B=1" << std::endl;
    //
    // for (int k = 0; k <= 2*3-2; k++) {
    //   int n = quadmath_snprintf (buf, sizeof buf, "%+-#*.30Qe", width,  zT_B[a2l(k, 0)]);
    //   printf ("%s\n", buf);
    // }
    //
    // std::cout<<"T_predL,r=1" << std::endl;
    //
    // for (int k = 0; k <= 2*3-2; k++) {
    //   int n = quadmath_snprintf (buf, sizeof buf, "%+-#*.30Qe", width, T_pred_L[a2l(k, 0)]);
    //   printf ("%s\n", buf);
    // }
    //
    //
    //
    // std::cout<<"T_predT,r=1" << std::endl;
    //
    // for (int k = 0; k <= 2*3-2; k++) {
    //   int n = quadmath_snprintf (buf, sizeof buf, "%+-#*.30Qe", width, T_pred_T[a2l(k, 0)]);
    //   printf ("%s\n", buf);
    // }
    //
    //
    // std::cout<<"T_predR,r=1" << std::endl;
    //
    // for (int k = 0; k <= 2*3-2; k++) {
    //   int n = quadmath_snprintf (buf, sizeof buf, "%+-#*.30Qe", width, T_pred_R[a2l(k, 0)]);
    //   printf ("%s\n", buf);
    //
    // }
    //
    //
    // std::cout<<"T_predB=1" << std::endl;
    //
    // for (int k = 0; k <= 2*3-2; k++) {
    //   int n = quadmath_snprintf (buf, sizeof buf, "%+-#*.30Qe", width, T_pred_B[a2l(k, 0)]);
    //   printf ("%s\n", buf);
    // }


    init_deallocate();
  }
};






#endif
