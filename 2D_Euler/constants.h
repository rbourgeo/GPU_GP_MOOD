#ifndef CONSTANTS_H

#include <cmath>
#include <fstream>
#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <quadmath.h>
inline __float128 longdouble(const int i)
{
  __float128 result;

  result = i;
  return result;
}


#define CONSTANTS_H

const double gr_gamma = 1.4; //isentropic index

const int FOG = 1, PLM = 2, Linear_GP_SE = 3;

const int dir_x = 0, dir_y = 1;

const int SSP_RK1 = 5, SSP_RK2 = 6, SSP_RK3 = 7;

const __float128 pi = 4*atanq(longdouble(1)), ospi = longdouble(1)/sqrtq(pi);

const __float128 sq3 = sqrtl(3.), sq35 = sqrtl(3. / 5);

const int iL = 0, iT = 1, iR = 2, iB = 3;

const int ord1 = 0, ord3 = 1, ord5 = 2;

const int i_rho     = 0, i_momx = 1, i_momy = 2, i_ener = 3;
const int Lin_gauss = 3, Sod    = 4, RP_3 = 5;
const int Neumann = 1, Periodic = 2 ;

#endif
