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
extern int nrow;
extern int ncol;
extern char order;
extern int ctx;
extern int mype;
extern int npe;
extern int myrow;
extern int mycol;


extern MPI_Comm mpi_comm;
extern MPI_Comm row_comm;
extern MPI_Comm cartcomm;

extern MPI_Request request;
extern MPI_Request request2;

extern MPI_Status status;

// call this method if you really are running mpi...
extern void initMpiStuff();
extern void initCblacsStuff();

}
#endif