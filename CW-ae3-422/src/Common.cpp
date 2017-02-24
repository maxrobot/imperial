#include <iostream>
#include <cmath>
#include <fstream>
#include <sstream>
#include <iomanip>
#include "Common.hpp"
#include "GlobalVars.hpp"

using namespace std;

void ShowMatrix(double *M, int N)
{ for (int i = 0; i < N; ++i)
  {
    for (int j = 0; j < N; ++j)
    { 
      int pnt = i*N + j;
      cout << M[pnt] << " ";
    }
    cout << endl;
  }
}