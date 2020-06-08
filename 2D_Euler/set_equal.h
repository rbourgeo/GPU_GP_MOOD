/*
* set_equal    : set out at in
* set comb_lin : out = coef1*in1 + coef2*in2
*/
__global__ void set_equal(double out[], const double in[]);
__global__ void set_comb_lin(double out[], const double in1[], const double in2[], const double coef1, const double coef2);
