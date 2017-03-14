#ifndef __COMMONMPI_H
#define __COMMONMPI_H

#include <mpi.h>
#include <string>

//#define  _LONG_LONG_INT_
//#ifdef _LONG_LONG_INT_
// typedef  long long int mint8 ;
//#else
// typedef  long int mint8 ;
//#endif
// typedef long long int mint8 ;
//# define MPI_MINT8 MPI_LONG_LONG
// ibermejoAdded
typedef long long int mint8;
#define MPI_MINT8 MPI_LONG_LONG
// typedef int mint8 ;
//# define MPI_MINT8 MPI_INT
// ibermejoAddedEnd

namespace MPI {

// these namespace members are declared "extern" so we can include
// the namespace definition in a header file included by
// multiple routines

extern int mpi_rank;
extern int mpi_size;
extern int ndims;
extern int reorder;
extern int n_rhs;
extern int n_abv;
extern int *mpi_coords;
extern int *n_coords;
extern int *n_rank;
extern int *periods;
extern int *mpi_dims;

// Clbas Variables
// extern int nrow;
// extern int ncol;
// extern char order;
// extern int ctx;
// extern int mype;
// extern int npe;
// extern int myrow;
// extern int mycol;


extern MPI_Comm mpi_comm;
extern MPI_Comm row_comm;
extern MPI_Comm cartcomm;

extern MPI_Request request;
extern MPI_Request request2;

extern MPI_Status status;

// call this method if you really are running mpi...
extern void initMpiStuff();
extern void initCblacsStuff();
extern void initMpiDomain();
extern void getNeighbours();
extern void copyVecConts(double *U, int Nghost_);
extern void exchangeVecConts(double *U, int Nghost_);
extern void getLeftData(double *d1, double *d2);
extern void getRightData(double *d1, double *d2);
extern void gatherData(double *output, double *U, int Nvar_g, int Nvar_);

}
#endif