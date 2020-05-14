#include "parameters.h"
#include "2D_bij.h"

__global__ void ftcs(const float f[], float fout[], const float dx, const float dy, const float k, const float dt);
