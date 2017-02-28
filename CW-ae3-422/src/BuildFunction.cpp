#include "Common.hpp"
#include "BuildFunction.hpp"
#include "GlobalVars.hpp"

using namespace std;

void buildKglb(double *Kg, double *ke, int Nvar_, int Nx_g)
{	// LHS Boundary
	for (int j = 3; j < 6; ++j)
	{	for (int k = 3; k < 6; ++k)
		{	int pnt = j*6 + k;
			int pnt2 = (j-3)*Nvar_ + (k-3);
			Kg[pnt2] += ke[pnt];
		}
	}
	// Central Elements
	for (int i = 1; i < Nx_g-1; ++i)
	{	for (int j = 0; j < 6; ++j)
		{	for (int k = 0; k < 6; ++k)
			{	int pnt = j*6 + k;
				int pnt2 = j*Nvar_ + (i-1)*(3*Nvar_) + k + (i-1)*3;
				Kg[pnt2] += ke[pnt];
			}
		}
	}
	// RHS Boundary
	for (int j = Nvar_-3; j < Nvar_; ++j)
	{	for (int k = Nvar_-3; k < Nvar_; ++k)
		{	int pnt = (j-(Nvar_-3))*6 + (k-(Nvar_-3)) ;
			int pnt2 = j*Nvar_ + k;
			Kg[pnt2] += ke[pnt];
		}
	}
}

void buildKglbSparse(double *Kg, double *ke, int Nvar_, int Nx_g, int buf)
{	const int jmp =  9 + buf;
	// LHS Boundary
	for (int i = 0; i < 3; ++i)
	{	// Build Diagonals
		int pnt = (i+3)*6 + (i+3);
		int pnt2 = (i*jmp) +4+buf;
	 	Kg[pnt2] += ke[pnt];
	 	if (i==1)
	 	{	pnt = (i+3)*6 + (i+3+1);
			pnt2 = (i*jmp+4)+1+buf;
		 	Kg[pnt2] += ke[pnt];	
	 	}
	 	if (i==2)
	 	{	pnt = (i+3-1)*6 + (i+3);
			pnt2 = (i*jmp+4)-1+buf;
		 	Kg[pnt2] += ke[pnt];	
	 	}
	}
	
	// Central Elements
	for (int i = 0; i < Nx_g-2; ++i)
	{	
		for (int j = 0; j < 6; ++j)
		{	// Build Diagonal
			int pnt = j*6 + j;
			int pnt2 = (j*jmp + 4+buf) + i*(jmp*3);
		 	Kg[pnt2] += ke[pnt];
		 	int max = 6-j;
		 	int bnd = 5;
		 	// Build Upper
		 	if (max < bnd)
		 	{	bnd = max;
		 	}
			for (int k = 1; k < bnd; ++k)
			{	int pos = (j*jmp + 4+buf) + i*(jmp*3) + k;
			 	Kg[pos] += ke[pnt+k];	
			}
		 	// Build Lower
			for (int k = 1; k <= j; ++k)
			{	int pos = (j*jmp + 4+buf) + i*(jmp*3) - k;
				Kg[pos] += ke[pnt-k];
			}
		}
	}
	// RHS Boundary
	for (int i = Nvar_-3; i < Nvar_; ++i)
	{	// Build Diagonals
		int pnt = (i - (Nvar_-3))*6 + (i - (Nvar_-3));
		int pnt2 = (i*jmp) +(4+buf);
	 	Kg[pnt2] += ke[pnt];	
	 	if (i==Nvar_-2)
	 	{	pnt += 1;
			pnt2 += 1;
		 	Kg[pnt2] += ke[pnt];	
	 	}
	 	if (i==Nvar_-1)
	 	{	pnt += -1;
			pnt2 += -1;
		 	Kg[pnt2] += ke[pnt];	
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
{	double EI(E_*I_), a(A_*E_/lx_e), b(12*EI/pow(lx_e, 3)), 
			c(6*EI/pow(lx_e, 2)), d(4*EI/lx_e), e(2*EI/lx_e);
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
			K[i*N+4] = -b;
			K[i*N+5] = c;
			K[2*N+i] = c;
			K[4*N+i] = -b;
			K[5*N+i] = c;
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
		{	K[i*N+5] = d;
		}
	  }
	}
}

void buildFele(double *K, double lx_e, double qx_, double qy_)
{	K[0] = K[3] = lx_e*qx_/2;
	K[1] = K[4] = lx_e*qy_/2;
	K[2] = (qy_*lx_e*lx_e)/12;
	K[5] = -(qy_*lx_e*lx_e)/12;
}

void zeroMat(double *K, int N)
{	fill(K, K+N, 0);
}