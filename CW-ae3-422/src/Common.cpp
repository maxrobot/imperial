#include "Common.hpp"
#include "CommonMPI.hpp"

using namespace std;

void showVec(double *M, int N)
{	for (int i = 0; i < N; ++i)
	    cout << setprecision(5) << setw(12) << M[i] << endl;
	cout << endl;
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

void printMessage(std::string message)
{	if (MPI::mpi_rank==0)
	{	std::cout << message << std::endl;
	}
}