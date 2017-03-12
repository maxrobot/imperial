#include <iostream>
#include <stdio.h> // for gcc >= 4.4 compatibility
#include "CommonMPI.hpp"
#include "Common.hpp"

using namespace std;

namespace MPI
{

// give single processor initial values to these
// in practice, the user should call initMpiStuff
// near the start of there code, after MPI_Init(*,*)...
int mpi_rank;
int mpi_size = 1;
int nrow = 1;
int ncol = 2;
char order = 'R';
int ctx;
int mype;
int npe;
int myrow;
int mycol;


MPI_Comm mpi_comm = MPI_COMM_WORLD;
MPI_Comm row_comm;
MPI_Comm cartcomm;

MPI_Request request;
MPI_Request request2;

MPI_Status status;

	void initMpiStuff() 
	{	MPI_Comm_rank(mpi_comm, &mpi_rank);
	  	MPI_Comm_size(mpi_comm, &mpi_size);
	  	MPI_Comm_dup(MPI_COMM_WORLD, &mpi_comm);
	}

	void initCblacsStuff()
	{   Cblacs_pinfo(&mype, &npe);
	    Cblacs_get(0, 0, &ctx);
	    Cblacs_gridinit(&ctx, &order, 1, npe);
	    Cblacs_gridinfo(ctx, &nrow, &ncol, &myrow, &mycol);
	}
	
	void checkMPI(int retval)
	{	if(retval!=MPI_SUCCESS) 
		{	printMessage("An error occurred initialising MPI");
		}
	}
}