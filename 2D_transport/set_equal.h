/*
* set_equal    : set out at in
* fluxes[] : the array in question
*/
__global__ void set_equal(double out[], const double in[]);
__global__ void set_comb_lin(double out[], const double in1[], const double in2[], const double coef1, const double coef2);
