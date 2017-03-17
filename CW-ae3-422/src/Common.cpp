#include "Common.hpp"
#include "CommonMPI.hpp"

using namespace std;

// Shows vector
void showVec(double *M, int N)
{	for (int i = 0; i < N; ++i)
	    cout << setprecision(5) << setw(12) << M[i] << endl;
	cout << endl;
}

// Shows parallel vector without overlap, calls showParVec(M,N,Sghost)
void showDenVec(double *M, int N, int Sghost_)
{	if (MPI::mpi_rank==0)
	{	cout << endl;
		for (int i = 0; i < N-(Sghost_+3); ++i)
		    cout << setprecision(5) << setw(12) << M[i] << endl;
		cout << endl;
    }
    for (int i = 1; i < MPI::mpi_size-1; ++i)
    {	if (MPI::mpi_rank==i)
	    {	for (int i = 0; i < N-(Sghost_+3); ++i)
			    cout << setprecision(5) << setw(12) << M[i] << endl;
			cout << endl;
	    }
    }
	if (MPI::mpi_rank==MPI::mpi_size-1)
	{   for (int i = 0; i < N-(Sghost_+3); ++i)
		    cout << setprecision(5) << setw(12) << M[i] << endl;
		cout << endl;
    }
}

// Shows parallel vector
void showParVec(double *M, int N)
{	MPI_Barrier(MPI_COMM_WORLD);
    for (int i = 0; i < MPI::mpi_size; ++i)
    {	if (MPI::mpi_rank==i)
	    {	MPI_Barrier(MPI_COMM_WORLD);
	    	showVec(M, N);
	    	MPI_Barrier(MPI_COMM_WORLD);
	    }
    	MPI_Barrier(MPI_COMM_WORLD);
    }
	MPI_Barrier(MPI_COMM_WORLD);
}

// Shows parallel vector, but doesn't output overlap
void showParVec(double *M, int N, int Sghost_)
{	MPI_Barrier(MPI_COMM_WORLD);
	if (MPI::mpi_rank==0)
	{   MPI_Barrier(MPI_COMM_WORLD);
		showDenVec(M, N, Sghost_);
	   	MPI_Barrier(MPI_COMM_WORLD);
    }
	MPI_Barrier(MPI_COMM_WORLD);
    for (int i = 1; i < MPI::mpi_size-1; ++i)
    {	if (MPI::mpi_rank==i)
	    {	MPI_Barrier(MPI_COMM_WORLD);
			showDenVec(M, N, Sghost_);
	    	MPI_Barrier(MPI_COMM_WORLD);
	    }
    	MPI_Barrier(MPI_COMM_WORLD);
    }
	MPI_Barrier(MPI_COMM_WORLD);
	if (MPI::mpi_rank==MPI::mpi_size-1)
	{   MPI_Barrier(MPI_COMM_WORLD);
		showDenVec(M, N, Sghost_);
	   	MPI_Barrier(MPI_COMM_WORLD);
    }
	MPI_Barrier(MPI_COMM_WORLD);
}

// Shows NxN, matrix M = input matrix
void showMat(double *M, int N)
{	cout << endl;
	for (int i = 0; i < N; ++i)
  	{ for (int j = 0; j < N; ++j)
	    { int pnt = i*N + j;
	      cout << setprecision(3) << setw(9)  << M[pnt] << setw(9);
	    }
	    cout << endl;
	}
	cout << endl;
}

// Shows NxO matrix M = input matrix
void showMat(double *M, int N, int O)
{	cout << endl;
	for (int i = 0; i < N; ++i)
  	{ 	for (int j = 0; j < O; ++j)
	    { 	int pnt = j*N + i;
			cout << setprecision(3) << setw(9)  << M[pnt] << setw(9);
	    }
	    cout << endl;
	}
	cout << endl;
}

// Shows parallel NxO matrix M = input matrix
void showParMat(double *M, int N, int O)
{	MPI_Barrier(MPI_COMM_WORLD);
    for (int i = 0; i < MPI::mpi_size; ++i)
    {	if (MPI::mpi_rank==i)
	    {	MPI_Barrier(MPI_COMM_WORLD);
	    	showMat(M, N, O);
	    	MPI_Barrier(MPI_COMM_WORLD);
	    }
    	MPI_Barrier(MPI_COMM_WORLD);
    }
	MPI_Barrier(MPI_COMM_WORLD);
}