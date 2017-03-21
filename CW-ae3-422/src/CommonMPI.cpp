#include <iostream>
#include <stdio.h>
#include "CommonMPI.hpp"
#include "Common.hpp"
#include "Output.hpp"

using namespace std;

namespace MPI
{

// give single processor initial values to these
// in practice, the user should call initMpiStuff
// near the start of there code, after MPI_Init(*,*)...
int mpi_rank;
int mpi_size;
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

// Cblacs things
int nrow = 1;
int ncol;
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
	{   ncol = mpi_size;
		Cblacs_pinfo(&mype, &npe);
	    Cblacs_get(0, 0, &ctx);
	    Cblacs_gridinit(&ctx, &order, nrow, ncol);
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
		getNeighbours();
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

	void copyVecConts(double *U, int Nghost_)
	{	// Set exchange length
		int bnd = 3;
		int Sghost_ = 3;

		// Temp array for old input
		double *temp =  new double[Nghost_]();
		for (int i = 0; i < Nghost_; ++i)
		{	temp[i] = U[i];
		}
		for (int i = 0; i < bnd; ++i)
		{	int pnt = Sghost_ + i;
			int pnt2 = Nghost_ - bnd  + i;
			double d1 = U[pnt], d2;
			getLeftData(&d1, &d2);
			if (mpi_rank<(MPI::mpi_size-1))
			{	U[pnt2] = d2;
			}
			pnt = Nghost_ - bnd - Sghost_ + i;
			pnt2 = i;
			d1 = temp[pnt];
			getRightData(&d1, &d2);
			if (mpi_rank>0)
			{	U[pnt2] = d2;
			}
		}
	}

	void exchangeVecConts(double *U, int Nghost_, int Sghost_)
	{	// Set exchange length
		int bnd = 6;

		// Temp array for old input
		double *temp =  new double[Nghost_]();
		for (int i = 0; i < Nghost_; ++i)
		{	temp[i] = U[i];
		}
		for (int i = 0; i < bnd; ++i)
		{	int pnt = 0 + i;
			int pnt2 = Nghost_-bnd + i;
			double d1 = U[pnt], d2;
			getLeftData(&d1, &d2);
			if (mpi_rank<(MPI::mpi_size-1))
			{	U[pnt2] += d2;
			}
			d1 = temp[pnt2];
			getRightData(&d1, &d2);
			if (mpi_rank>0)
			{	U[pnt]+= d2;
			}
		}
	}

	void getLeftData(double *d1, double *d2)
	{	MPI_Request request;
		MPI_Status status;

		// Find the neighbours...
	    int right = n_rhs;
	    int left  = MPI::mpi_rank-1;
	    if (left<0)
	    {   left = MPI::mpi_size - 1;
	    }
	    // MPI_Sendrecv
		MPI_Sendrecv(d1, 1, MPI_DOUBLE, left, 123, d2, 1, MPI_DOUBLE, right, 123, MPI_COMM_WORLD, &status);
	}
	
	void getRightData(double *d1, double *d2)
	{	MPI_Request request;
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

	void gatherData(double *output, double *U, int Nvar_g, int Nvar_)
	{	MPI_Request request;
		MPI_Status status;
		int right = n_rhs, left  = MPI::mpi_rank-1;
		int pnt = 0;
		double sendArray[Nvar_];
		double d1, d2;

		// Fill in root proc
		if (MPI::mpi_rank==0)
		{	for (int i = 0; i < Nvar_; ++i)
			{	output[i] = U[i];
				pnt += 1;
			}
		}

		// Fill in following procs
		for (int i = 1; i < MPI::mpi_size; ++i)
		{	int bnd;
			if (MPI::mpi_rank==i)
			{	d1 = Nvar_;
				bnd = Nvar_;
				MPI_Send(&d1, 1, MPI_DOUBLE, 0, 123, MPI_COMM_WORLD);
			}
			else if (MPI::mpi_rank==0)
			{	MPI_Recv(&d2, 1, MPI_DOUBLE, i, 123, MPI_COMM_WORLD, &status);
				bnd = d2;
			}
			for (int j = 0; j < bnd; ++j)
			{	if (MPI::mpi_rank==i)
				{	d1 = U[j];
					MPI_Send(&d1, 1, MPI_DOUBLE, 0, 123, MPI_COMM_WORLD);
				}
				else if (MPI::mpi_rank==0)
				{	MPI_Recv(&d2, 1, MPI_DOUBLE, i, 123, MPI_COMM_WORLD, &status);
					output[pnt] = d2;
					pnt += 1;
				}
			}
			MPI_Barrier(MPI_COMM_WORLD);
		}
	}
}