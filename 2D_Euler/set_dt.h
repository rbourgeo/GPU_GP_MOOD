#include "constants.h"
#include "parameters.h"
#include "physics.h"

__global__ void find_min( double d_out[] , double d_in[] )
{
  int myId = threadIdx . x + blockDim . x * blockIdx . x ;
  int tid  = threadIdx . x ;
  // do reduction over global memory
  for ( unsigned int s = blockDim . x/2; s > 0; s >>= 1)
  {
    if ( tid < s )
    {
      d_in [ myId ] = min(d_in[ myId + s ],d_in [ myId ]);
    }
    __syncthreads () ;
    // Maske sure all adds at one stage
  }
  // only thread 0 writes result for this block back to global mem
  if ( tid == 0)
  {
    d_out [ blockIdx . x ] = d_in [ myId ];
  }
}




__global__ void fill_dt_array(const double f[], double dt_array[],  const double CFL,  const double dx,  const double dy)
{
  int tidx = c2f(threadIdx.x + blockIdx.x * blockDim.x);
  int tidy = c2f(threadIdx.y + blockIdx.y * blockDim.y);

  double u[4];
  double vx, vy, c, pressure, carac_time;


  if (tidx >= 1-ngc  && tidx <= lf+ngc ) {
    if (tidy >= 1-ngc  && tidy <= nf+ngc ) {
      {
        for (int i_cons = i_rho; i_cons <= i_ener; i_cons++)
        {
          u[i_cons] = f[ij_sol(tidy, tidx, i_cons)];
        }

        vx       = xvel(u);
        vy       = yvel(u);
        pressure = pres(u);
        c        = sound_speed(pressure, u);
        carac_time = CFL * min(dx / max(abs(vx - c), abs(vx + c)), dy / max(abs(vy - c), abs(vy + c)));
        dt_array[ij(tidy, tidx)] = carac_time;
      }
    }
  }
}


__host__ double find_min_host(const double array[]){

  double min = 10.0;

  for (size_t i = 0; i <= lf; i++) {
    for (size_t j = 0; j <= nf; j++) {
      if (array[ij(j,i)]<min) {
        min = array[ij(j,i)] ;
      }

    }
  }
  return min;
}



inline double set_dt(
  const double f[],
  double dt_array[],
  double dt_array_out[],
  const double CFL,
  const double dx,
  const double dy)
  {

    double result;
    double *host_array;

    host_array = new double[ne*le];

    fill_dt_array<<<dimGrid, dimBlock>>>(f, dt_array, CFL, dx, dy);
    cudaDeviceSynchronize();

    // //Kernel parameters
    // dim3 dimBlock_1D(16,1,1);
    // dim3 dimGrid_1D(le*ne/dimBlock.x, 1,1);
    //
    // find_min<<<dimGrid_1D, dimBlock_1D>>>( dt_array_out, dt_array);
    // cudaDeviceSynchronize();



    /* sequential min search */

    cudaMemcpy(host_array, dt_array, le*ne*sizeof(double), cudaMemcpyDeviceToHost);
    result = find_min_host(host_array);

    std :: cout << result << std::endl;

    return result;
  }
