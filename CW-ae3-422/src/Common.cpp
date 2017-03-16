#include "Common.hpp"
#include "CommonMPI.hpp"

using namespace std;

void showVec(double *M, int N)
{	for (int i = 0; i < N; ++i)
	    cout << setprecision(5) << setw(12) << M[i] << endl;
	cout << endl;
}

void showDenVec(double *M, int N, int Sghost_)
{	if (MPI::mpi_rank==0)
	{   for (int i = 0; i < N-(Sghost_+3); ++i)
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
	    	showVec(M, N);
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

void showMat(double *M, int N)
{	for (int i = 0; i < N; ++i)
  	{ for (int j = 0; j < N; ++j)
	    { int pnt = i*N + j;
	      cout << setprecision(3) << setw(9)  << M[pnt] << setw(9);
	    }
	    cout << endl;
	}
	cout << endl;
}

void showMat(double *M, int N, int O)
{	for (int i = 0; i < N; ++i)
  	{ 	for (int j = 0; j < O; ++j)
	    { 	int pnt = j*N + i;
			cout << setprecision(3) << setw(9)  << M[pnt] << setw(9);
	    }
	    cout << endl;
	}
	cout << endl;
}

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

void writeVec(double *M, int N, std::string test)
{	ofstream myfile;
	myfile.open ("./output/data/" + test + ".txt");
	// myfile.open ("./output/data/" + test + ".txt");
	myfile << setprecision(10) << 0. << setw(20) << 0.  <<  setw(20) << 0.  << endl;
	for (int i = 0; i < N-1; ++i)
	{	int pnt = i*3;
		myfile << setprecision(10) << M[pnt] << setw(20) << M[pnt+1] <<  setw(20) << M[pnt+2] << endl;
	}
	myfile << setprecision(10) << 0. << setw(20) << 0.  <<  setw(20) << 0.  << endl;
	myfile.close();
}

void writeVec(double *M, int N, int step, std::string test)
{	ofstream myfile;
	char numstr[21]; 
	sprintf(numstr, "%d", step);
	myfile.open ("./output/data/" + test + numstr + ".txt");
	myfile << setprecision(10) << 0. << setw(20) << 0.  <<  setw(20) << 0.  << endl;
	for (int i = 0; i < N-1; ++i)
	{	int pnt = i*3;
		myfile << setprecision(10) << M[pnt] << setw(20) << M[pnt+1] <<  setw(20) << M[pnt+2] << endl;
	}
	myfile << setprecision(10) << 0. << setw(20) << 0.  <<  setw(20) << 0.  << endl;
	myfile.close();
}

void writeParVec(double *M, int N, int Nvar_g, int Nvar_, int step, std::string test)
{	ofstream myfile;
	char numstr[21]; 
	sprintf(numstr, "%d", step);
	double *output = new double[Nvar_g]();
	MPI::gatherData(output, M, Nvar_g, Nvar_);

	if (MPI::mpi_rank==0)
	{	//showVec(output, Nvar_g);
		myfile.open ("./output/data/" + test + numstr + ".txt");
		myfile << setprecision(10) << 0. << setw(20) << 0.  <<  setw(20) << 0.  << endl;
		for (int i = 0; i < N-1; ++i)
		{	int pnt = i*3;
			myfile << setprecision(10) << output[pnt] << setw(20) << output[pnt+1] <<  setw(20) << output[pnt+2] << endl;
		}
		myfile << setprecision(10) << 0. << setw(20) << 0.  <<  setw(20) << 0.  << endl;
		myfile.close();
	}
}

void printMessage(std::string message)
{	if (MPI::mpi_rank==0)
	{	std::cout << message << std::endl;
	}
}