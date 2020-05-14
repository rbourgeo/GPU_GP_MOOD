#include "parameters.h"
#include "2D_bij.h"


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
