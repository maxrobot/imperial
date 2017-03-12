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
int row_rank;
int row_size;
int ndims = 1; // Number of dimensions in domain decomposition
int reorder = 1;
int n_rhs = 0;
int n_abv = 0;
int colour;
int *periods;
int *mpi_coords;
int *n_coords;
int *n_rank;
int *mpi_dims;

// CLbacs things
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

	void initMpiDomain()
	{	
		mpi_dims = new int[1];
		periods = new int[1];
		mpi_coords = new int[2];

		MPI_Dims_create(mpi_size, ndims, mpi_dims);
		MPI_Cart_create(mpi_comm, 1, mpi_dims, periods, reorder, &cartcomm);
		MPI_Cart_coords(cartcomm, mpi_rank, ndims, mpi_coords);
  	}

	void getNeighbours()
	{	// rhs is right hand neighbour, abv is neighbour above current cell...
		int *rhs = new int[2];
		rhs[0] = mpi_coords[0];
		rhs[1] = mpi_coords[1];

		// Find neighbours to the right...
		if (mpi_coords[0] < mpi_dims[0])
		{	if (mpi_coords[0] == (mpi_dims[0]-1))
			{	rhs[0] = 0;
				rhs[1] = mpi_coords[1];
				MPI_Cart_rank(cartcomm, rhs, &n_rhs);
			}
			else
			{	rhs[0] = mpi_coords[0]+1;
				rhs[1] = mpi_coords[1];
				MPI_Cart_rank(cartcomm, rhs, &n_rhs);
			}
		}
	}

	void exchangeLeftData(double *d1, double *d2)
	{	MPI::getNeighbours();
		// MPI_Comm_rank(row_comm, &row_rank);
		MPI_Request request;
		MPI_Status status;

		// Find the neighbours...
	    int right = n_rhs;
	    int left  = MPI::mpi_rank-1;
	    if (left<0)
	    {   left = MPI::mpi_size - 1;
	    }
	    // MPI_Sendrecv
	    MPI_Barrier(MPI_COMM_WORLD);
		MPI_Sendrecv(d1, 1, MPI_DOUBLE, left, 123, d2, 1, MPI_DOUBLE, right, 123, MPI_COMM_WORLD, &status);
	}
	
	void exchangeRightData(double *d1, double *d2)
	{	MPI::getNeighbours();
		// MPI_Comm_rank(row_comm, &row_rank);
		MPI_Request request;
		MPI_Status status;

		// Find the neighbours...
	    int right = n_rhs;
	    int left  = MPI::mpi_rank-1;
	    if (left<0)
	    {   left = MPI::mpi_size - 1;
	    }
	    // MPI_Sendrecv 
		MPI_Sendrecv(d1, 1, MPI_DOUBLE, right, 123, d2, 1, MPI_DOUBLE, left, 123, MPI_COMM_WORLD, &status);
	}
}