/*
* Update_A/B/C  : sequence that update the solution given the fluxes at interfaces
                it's splitted in 3 steps to avoid race condiitions
* fluxes_x/y[]  : the fluxes in question
* f[]			      : the solution
*/

__global__ void Update_A(double f[], const double fluxes_x[], const double fluxes_y[]);
__global__ void Update_B(double f[], const double fluxes_x[], const double fluxes_y[]);
__global__ void Update_C(double f[], const double fluxes_x[], const double fluxes_y[]);
