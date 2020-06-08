#include "parameters.h"

/*
* io_fun : write the outpu file
* file   : name of the file
* x     : pointer to x mesh
* y     : pointer to y mesh
* f     : solution
*/
__host__ void io_fun(std::string file, double* x, double* y, double* f)
{
    std::ofstream myfile_tsN;
    myfile_tsN.open(file);
    for (int i = 1; i <= nf; i++) {
        for (int j = 1; j <= lf; j++) {
            //	std::cout << "i" << i << "yi" << y[i]<< "j" << j << "xj" << x[j] << std::endl;
            myfile_tsN << x[f2c(j)] << '\t' << y[f2c(i)] << '\t';
            myfile_tsN << f[ij_sol(i, j, i_rho)] <<" "<<f[ij_sol(i, j, i_momy)]<< std::endl;
        }
    }

    myfile_tsN.close();
}
