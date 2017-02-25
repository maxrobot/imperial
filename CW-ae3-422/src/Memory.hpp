#ifndef MEMORY_HPP_INCLUDED
#define MEMORY_HPP_INCLUDED

double *allocateDbl(int Ny)
{ double *data  = 0;
  data = new double[Ny];
  return data;
}

#endif // MEMORY_HPP_INCLUDED