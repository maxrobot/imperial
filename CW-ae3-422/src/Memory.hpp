#ifndef MEMORY_HPP_INCLUDED
#define MEMORY_HPP_INCLUDED


int ***allocateInt(int Nx, int Ny, int Q)
{ int ***data  = new int**[Ny];
  for (int i = 0; i < Ny; ++i)
  { data[i]  = new int*[Nx];
    for (int j = 0; j < Nx; ++j)
    { data[i][j]  = new int[Q];
    }
  }
  return data;
}

int **allocateInt(int Nx, int Ny)
{ int **data  = new int*[Ny];
  for (int i = 0; i < Ny; ++i)
  { data[i]  = new int[Nx];
  }
  return data;
}

int *allocateInt(int Nx)
{ int *data  = new int[Nx];
  return data;
}

double ***allocateDbl(int Nx, int Ny, int Q)
{ double ***data  = new double**[Ny];
  for (int i = 0; i < Ny; ++i)
  { data[i]  = new double*[Nx];
    for (int j = 0; j < Nx; ++j)
    { data[i][j]  = new double[Q];
    }
  }
  return data;
}

double **allocateDbl(int Nx, int Ny)
{ double **data  = new double*[Ny];
  for (int i = 0; i < Ny; ++i)
  { data[i]  = new double[Nx];
  }
  return data;
}

double *allocateDbl(int Ny)
{ double *data  = new double[Ny];
  return data;
}

void deallocateInt(int Ny, int **t)
{ for (int i=0; i < Ny; i++)
  { delete[] t[i];
  }
  delete[] t;
  t = 0;
}

void deallocateDbl(int Ny, double **t)
{ for (int i=0; i < Ny; i++)
  { delete[] t[i];
  }
  delete[] t;
  t = 0;
}

void deallocateDbl(int Ny, int Nx, double ***t)
{ for(int i(0); i < Ny; i++)
  { for(int j(0); j < Nx; j++)
    { delete[] t[i][j];
    }
    delete[] t[i];
  }
  delete[] t;
  t = 0;
}

#endif // MEMORY_HPP_INCLUDED