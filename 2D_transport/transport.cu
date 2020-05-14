#include "parameters.h"
#include "2D_bij.h"
#include "bc.h"
#include "set_equal.h"
#include "init.h"
#include "FE.h"


__host__ void io_fun(std::string file, float *x, float *y, float *f)
{
	std::ofstream myfile_tsN;
	myfile_tsN.open(file);
	for(int i = 0; i < nf; i++)
	{
		for(int j = 0; j < lf; j++)
		{
			//	std::cout << "i" << i << "yi" << y[i]<< "j" << j << "xj" << x[j] << std::endl;
			myfile_tsN << x[j] << '\t'<< y[i] << '\t';
			myfile_tsN << f[ij(i,j)] << std::endl;
		}
	}

	myfile_tsN.close();
}





int main()
{

	//time measurement variables
	float time;
	cudaEvent_t start, stop  ;
	cudaEventCreate(&start)  ;
	cudaEventCreate(&stop)   ;

	//Integration Parameters
	float k = 1.0f; //Thermal Conductivity
	/* Numerical Mesh Configuration */
	float dx = 1.0f/float(lf);
	float dy = 1.0f/float(nf);

	float dt = 0.25f*(dx*dx)/k;

	float tmax = 2.0f;

	float t = 0.0f, tio = 0.5f;



	//Allocate Memory
	size_t sz  = lf*nf*sizeof(float);
	size_t szx = lf*sizeof(float);
	size_t szy = nf*sizeof(float);


	float *f, *x, *y;

	f    = new float[lf*nf];

	x = new float[lf];
	y = new float[nf];


	float *d_f, *d_x, *d_y, *d_fout;

	cudaMalloc(&d_f   , sz );
	cudaMalloc(&d_fout, sz );

	cudaMalloc(&d_x, szx);
	cudaMalloc(&d_y, szy);

	//Kernel parameters
	dim3 dimBlock(16,16,1);
	dim3 dimGrid(lf/dimBlock.x, nf/dimBlock.y,1);

	//Apply Initial Condition You could also create a kernel for this

	initialize<<<dimGrid, dimBlock>>>(d_f, d_x, d_y, dx, dy);

	bc<<<dimGrid, dimBlock>>>(d_f); // First BC
	cudaDeviceSynchronize();


	//Copy for IO operation
	cudaMemcpy(x, d_x, szx, cudaMemcpyDeviceToHost);
	cudaMemcpy(y, d_y, szy, cudaMemcpyDeviceToHost);


	//device x is no longer needed
	cudaFree(d_x);
	cudaFree(d_y);


	/*====================== Perform Integration =======================*/
	std::string f2;
	int kk = 0;

	cudaEventRecord(start, 0); // We only measure the computation time

	while(t<tmax)
	{

		//Call the stencil routine -> fout
		ftcs<<<dimGrid, dimBlock>>>(d_f, d_fout, dx, dy, k, dt);
		cudaDeviceSynchronize();

		// put fout in in. (avoid race condition)
		set_equal<<<dimGrid, dimBlock>>>(d_f, d_fout);
		cudaDeviceSynchronize();

		//Call BC
		bc<<<dimGrid, dimBlock>>>(d_f);
		cudaDeviceSynchronize();


	//	if(fmod(t, tio) == 0.0f)
	//	{
		//	//IO function
		//	f2 = "sol" + std::to_string(kk) + ".dat";
		//	cudaMemcpy(f,d_f, sz, cudaMemcpyDeviceToHost);
		//	io_fun(f2, x, y, f);
		//	kk++;
		//	std::cout<< "output at t = "<<t<< std::endl;
		//}

		t+=dt;
	}

	cudaEventRecord(stop, 0);
	cudaEventSynchronize(stop);
	cudaEventElapsedTime(&time, start, stop);

	printf("Time to compute:  %3.1f ms \n", time); //faire mieux


	//Final output
	f2 = "final_sol.dat";
	cudaMemcpy(f,d_f, sz, cudaMemcpyDeviceToHost);
	io_fun(f2, x, y, f);


	//deallocate memory
	delete x, y, f;
	cudaFree(d_f);
	cudaFree(d_fout);
}
