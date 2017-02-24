#include <iostream>
#include <cmath>
#include "Common.hpp"
#include "BuildFunction.hpp"
#include "GlobalVars.hpp"

using namespace std;

void buildKele(double *K, double lx_e)
{	double EI(E_*I_), a(A_*E_/lx_e), b(12*EI/pow(lx_e, 3)), 
			c(6*EI/pow(lx_e, 2)), d(12*EI/lx_e), e(2*EI/lx_e);
	int N = 6;
	{ for (int i = 0; i < N; ++i)
	  {	if (i == 0)
		{	K[i*N+0] = a;
			K[i*N+3] = -a;
			K[3*N] = -a;
		}
		if (i == 1)
		{	K[i*N+1] = b;
			K[i*N+2] = c;
			K[2*N+i] = c;
		}
		if (i == 2)
		{	K[i*N+2] = d;
			K[i*N+4] = -c;
			K[i*N+5] = e;
			K[4*N+i] = -c;
			K[5*N+i] = e;
		}
		if (i == 3)
		{	K[i*N+3] = a;
		}
		if (i == 4)
		{	K[i*N+4] = b;
			K[i*N+5] = -c;
			K[5*N+i] = -c;
		}
		if (i == 5)
		{	K[i*N+5] = c;
		}

	  }
	}
}

void buildFele(double *K, double lx_e)
{		
}