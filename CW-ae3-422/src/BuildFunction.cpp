#include <iostream>
#include <cmath>
#include "Common.hpp"
#include "BuildFunction.hpp"
#include "GlobalVars.hpp"

using namespace std;

void buildKglb(double *Kg, double *ke, int Nvar_, int Nx_g)
{	for (int i = 0; i < Nx_g; ++i)
	{	for (int j = 0; j < 6; ++j)
		{	for (int k = 0; k < 6; ++k)
			{	int pnt = j*6 + k;
				int pnt2 = j*Nvar_ + k + i*(3*Nvar_ + 3);
				// cout << Kg[pnt2] << " " << &ke[pnt] << endl;
				Kg[pnt2] = ke[pnt];
			}
		}
	}
}
void buildFglb(double *Kg, double *ke, int Nx_g)
{	for (int i = 0; i < Nx_g; ++i)
	{	for (int j = 0; j < 6; ++j)
		{	int pnt = i*3 + j;
			Kg[pnt] += ke[j];
		}
	}
}

void buildKele(double *K, double lx_e, double A_, double E_, double I_)
// {	double a(1), b(1), c(1), d(1), e(1);
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

void buildFele(double *K, double lx_e, double qx_, double qy_)
{	K[0] = K[3] = qx_/2;
	K[1] = K[4] = qy_/2;
	K[2] = K[5] = qy_*lx_e/12;
}

void zeroMat(double *K, int N)
{	fill(K, K+N, 0);
}