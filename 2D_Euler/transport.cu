#include "Forward_Euler.h"
#include "SSP_RK2.h"
#include "SSP_RK3.h"
#include "bc.h"
#include "exact_sol.h"
#include "gp_stencil.h"
#include "init.h"
#include "output.h"
#include "parameters.h"
#include "recons_flux_computation.h"
#include "set_equal.h"
#include "state.h"

int main()
{

  /*time measurement variables*/
  float time;
  cudaEvent_t start, stop;
  cudaEventCreate(&start);
  cudaEventCreate(&stop);

  //Integration Parameters

  /* Numerical Mesh Configuration */
  double dx = Lx / double(lf);
  double dy = Ly / double(nf);

  double dtfinal;

  double ell_;

  /*init and computation of GP_stencil in the class*/
  if (ell_over_dx == 0)
  {
    ell_ = ell;
  }
  else
  {
    ell_ = dx*ell_over_dx;
  }

  GP_stencil GP(Mord,ell_);
  GP.init(dx, dy);

  int *index, *d_index;
  double *zT, *d_zT;

  zT = new double[4*nop*ngp]; //4 faces, 2*mord-1 points, ngp points per face
  index = new int[2*nop    ]; // x and y

  for (int k = 0; k <= nop-1; k++) {
    for (int r = 0; r <= ngp-1; r++) {

      zT[iL*(nop*ngp) + r*nop + k] = GP.zT_L[a2l(k, r)];
      zT[iT*(nop*ngp) + r*nop + k] = GP.zT_T[a2l(k, r)];
      zT[iR*(nop*ngp) + r*nop + k] = GP.zT_R[a2l(k, r)];
      zT[iB*(nop*ngp) + r*nop + k] = GP.zT_B[a2l(k, r)];

      index[dir_x*nop + k] = GP.index_x[k];
      index[dir_y*nop + k] = GP.index_y[k];

    }
  }


  GP._deallocate();
  /* Send the GP stencil data to the device */
  cudaMalloc(&d_index, 2*nop*    sizeof(int   ));
  cudaMalloc(&d_zT   , 4*nop*ngp*sizeof(double));

  cudaMemcpy(d_zT,    zT   , 4*nop*ngp*sizeof(double), cudaMemcpyHostToDevice);
  cudaMemcpy(d_index, index, 2*nop*    sizeof(int   ), cudaMemcpyHostToDevice);


  /* Time stepping parameters */

  //double dt = CFL * pow(1.0f / (abs(ax) / dx + abs(ay) / dy),5./3);
  double dt = CFL *1.0f / (abs(ax) / dx + abs(ay) / dy);


  double t = 0.0f, tio = 0.5f;

  //Allocate Memory

  size_t total_size     =     ne * le * sizeof(double);
  size_t total_size_sol = 4 * ne * le * sizeof(double); /* 4 conserved quantities*/

  size_t size_x  = le  *       sizeof(double);
  size_t size_y  = ne  *       sizeof(double);
  size_t size_gw = ngp * ngp * sizeof(double);

  double *f,   *x,   *y,   *gauss_weight;
  double *d_f, *d_x, *d_y, *d_fout,      *d_fluxes_x, *d_fluxes_y, *d_f1, *d_f2, *d_gauss_weight;

  f = new double[4 * le * ne];
  x = new double[le         ];
  y = new double[ne         ];

  gauss_weight = new double[ngp * ngp];

  gauss_weight[a2l_ngp(0, 0)] = 1.;

  if (ngp >= 2) {
    gauss_weight[a2l_ngp(1, 0)] = 0.5;
    gauss_weight[a2l_ngp(1, 1)] = 0.5;
  }

  if (ngp >= 3) {
    gauss_weight[a2l_ngp(2, 0)] = 0.5 * 5. / 9;
    gauss_weight[a2l_ngp(2, 1)] = 0.5 * 8. / 9;
    gauss_weight[a2l_ngp(2, 2)] = 0.5 * 5. / 9;
  }

  //Send the gaussian point to the device
  cudaMalloc(&d_gauss_weight, size_gw);
  cudaMemcpy( d_gauss_weight, gauss_weight, size_gw, cudaMemcpyHostToDevice);


  cudaMalloc(&d_f       , total_size_sol);
  cudaMalloc(&d_f1      , total_size_sol);
  cudaMalloc(&d_f2      , total_size_sol);
  cudaMalloc(&d_fout    , total_size_sol);
  cudaMalloc(&d_fluxes_x, total_size_sol);
  cudaMalloc(&d_fluxes_y, total_size_sol);
  cudaMalloc(&d_x, size_x);
  cudaMalloc(&d_y, size_y);

  //Apply Initial Condition
  initialize<<<dimGrid, dimBlock>>>(d_f, d_x, d_y, dx, dy);

  //Copy for IO operation
  cudaMemcpy(x, d_x, size_x, cudaMemcpyDeviceToHost);
  cudaMemcpy(y, d_y, size_y, cudaMemcpyDeviceToHost);

  /*====================== Perform Integration =======================*/
  std::string f2;
  int kk = 0;

  cudaEventRecord(start, 0); // We only measure the computation time

  while (t < tmax) {

    dtfinal =  tmax-t;
    if (dt>dtfinal) {
      dt= dtfinal;
    }

    if (time_method == SSP_RK1) {

      Forward_Euler(d_f, d_fout,              d_fluxes_x, d_fluxes_y, dx, dy, dt, d_zT, d_index, d_gauss_weight);

    } else if (time_method == SSP_RK2) {
      SSP_RK2_(d_f, d_fout, d_f1,        d_fluxes_x, d_fluxes_y, dx, dy, dt, d_zT, d_index, d_gauss_weight);
    } else if (time_method == SSP_RK3) {
      SSP_RK3_(d_f, d_fout, d_f1, d_f2,  d_fluxes_x, d_fluxes_y, dx, dy, dt, d_zT, d_index, d_gauss_weight);
    }

    set_equal<<<dimGrid, dimBlock>>>(d_f, d_fout); //d_f = d_fout
    cudaDeviceSynchronize();

    //std::cout<<t<<std::endl;
    t += dt;
  }

  cudaEventRecord(stop, 0);
  cudaEventSynchronize(stop);
  cudaEventElapsedTime(&time, start, stop);

  printf("Time to compute:  %3.1f ms \n", time); //faire mieux

  //Final output
  f2 = "output/final_sol.dat";
  cudaMemcpy(f, d_f, total_size_sol, cudaMemcpyDeviceToHost);

  std::cout << "Error = " << error_(f, x, y, dx, dy) << std::endl;
  io_fun(f2, x, y, f);

  //deallocate memory
  delete x;
  delete y;
  delete f;
  delete zT;
  delete index;
  cudaFree(d_f);
  cudaFree(d_f1);
  cudaFree(d_f2);
  cudaFree(d_gauss_weight);
  cudaFree(d_fluxes_x);
  cudaFree(d_fluxes_y);
  cudaFree(d_fout);
  cudaFree(d_x);
  cudaFree(d_y);
  cudaFree(d_zT);
  cudaFree(d_index);



  return 0;
}

//	if(fmod(t, tio) == 0.0f)
//	{
//	//IO function
//	f2 = "sol" + std::to_string(kk) + ".dat";
//	cudaMemcpy(f,d_f, total_size, cudaMemcpyDeviceToHost);
//	io_fun(f2, x, y, f);
//	kk++;
//	std::cout<< "output at t = "<<t<< std::endl;
//}
