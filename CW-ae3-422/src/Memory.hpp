#ifndef MEMORY_HPP_INCLUDED
#define MEMORY_HPP_INCLUDED

double **allocateDbl(int Nx, int Ny)
{	double **data  = new double*[Ny]();
	for (int i = 0; i < Ny; ++i)
	{	data[i]  = new double[Nx]();
	}
	return data;
}

double *allocateDbl(int Ny)
{ double *data;
  data = new double[Ny]();
  return data;
}

#endif // MEMORY_HPP_INCLUDED