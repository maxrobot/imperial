 #include <iostream>
#include <stdio.h> // for gcc >= 4.4 compatibility
#include "CommonMPI.hpp"
#include "Common.hpp"

namespace MPI 
{
// give single processor initial values to these
// in practice, the user should call initMpiStuff
// near the start of there code, after MPI_Init(*,*)...
int mpi_rank;
int mpi_size = 1;

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

	void checkMPI(int retval)
	{	if(retval!=MPI_SUCCESS) 
		{	printMessage("An error occurred initialising MPI");
		}
	}
}