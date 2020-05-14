#include "parameters.h"
#include "2D_bij.h"

// bijection bewteen list to array
 __host__ __device__ int ij(int i, int j)
{
	return i*lf + j;
}
