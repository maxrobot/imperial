// #ifndef __COMMONMPI_H
// #define __COMMONMPI_H

// #include <mpi.h>
// #include <string>

// //#define  _LONG_LONG_INT_
// //#ifdef _LONG_LONG_INT_
// // typedef  long long int mint8 ;
// //#else
// // typedef  long int mint8 ;
// //#endif
// // typedef long long int mint8 ;
// //# define MPI_MINT8 MPI_LONG_LONG
// // ibermejoAdded
// typedef long long int mint8;
// #define MPI_MINT8 MPI_LONG_LONG
// // typedef int mint8 ;
// //# define MPI_MINT8 MPI_INT
// // ibermejoAddedEnd

// namespace MPI {

// // these namespace members are declared "extern" so we can include
// // the namespace definition in a header file included by
// // multiple routines

// extern int mpi_rank;
// extern int mpi_size;
// extern int ndims;
// extern int reorder;
// extern int n_rhs;
// extern int n_abv;
// extern int *mpi_coords;
// extern int *n_coords;
// extern int *n_rank;
// extern int *periodic;
// extern int *mpi_dims;

// extern MPI_Comm mpi_comm;
// extern MPI_Comm row_comm;
// extern MPI_Comm cartcomm;

// extern MPI_Request request;
// extern MPI_Request request2;

// extern MPI_Status status;

// // call this method if you really are running mpi...
// extern void initMpiStuff();
// extern void initMpiDomain();

// void splitComm(int i);

// void exchangeData(double **dat, int d1, int d2);

// void getNeighbours();

// void getRankCoords();

// void getGlobalBounds(int **x_bnd, int **y_bnd, int *Nx, int *Ny);

// void getLocalBounds(int **l_bnd, int **x_bnd, int **y_bnd);

// }
// #endif