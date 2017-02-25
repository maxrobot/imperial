#include <iostream>
#include <cmath>
#include <fstream>
#include <sstream>
#include <iomanip>
#include "Common.hpp"
#include "GlobalVars.hpp"

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
	      cout << setprecision(5) << setw(12)  << M[pnt] << setw(12);
	    }
	    cout << endl;
	}
	cout << endl;
}