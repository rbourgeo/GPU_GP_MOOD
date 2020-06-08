#include "parameters.h"

/*
* initialize : initialize the mesh and the solution
* f[]  : pointer to the solution
* x[]  : pointer to x mesh
* y[]  : pointer to y mesh
* dx   : x step
* dy   : y step
*/
__global__ void initialize(double f[], double x[], double y[], const double dx, const double dy);
