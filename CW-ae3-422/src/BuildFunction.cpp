#include "BuildFunction.hpp"
#include "Common.hpp"
#include "CommonMPI.hpp"

using namespace std;

// Builds an identity matrix
void buildEye(double *K, int Nvar_)
{	for (int i = 0; i < Nvar_; ++i)
	{	K[(i*Nvar_)+i] = 1;
	}
}

// Build global stiffness matrix K in dense format
void buildKg(double *Kg, double *ke, int Nvar_, int Nx_g)
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

// Builds global stiffness K in banded format
void buildKgband(double *Kg, double *ke, int Nvar_, int Nx_g, int buf_)
{	const int jmp = 9 + buf_;
	const int cnt = 4 + buf_;
	
	// Central Elements
	for (int i = 0; i < Nx_g-2; ++i)
	{	
		for (int j = 0; j < 6; ++j)
		{	// Build Diagonal
			int pnt = j*6 + j;
			int pnt2 = cnt + (j+(i*3))*jmp;
		 	Kg[pnt2] += ke[pnt];
		 	int max = 5-j;
		 	int bnd = 4;
		 	// Build Upper
		 	if (max < bnd)
		 	{	bnd = max;
		 	}
			for (int k = 1; k <= bnd; ++k)
			{	int pos = pnt2 + k;
			 	Kg[pos] += ke[pnt+k];	
			}
			bnd = 5;
		 	// Build Lower
		 	if (j < bnd)
		 	{	bnd = j+1;
		 	}
			for (int k = 1; k < bnd; ++k)
			{	int pos = pnt2 - k;
				Kg[pnt2-k] += ke[pnt-k];
			}
		}
		// cout << endl;
	}
	// LHS Boundary
	for (int i = 0; i < 3; ++i)
	{	// Build Diagonals
		int pnt = (i+3)*6 + (i+3);
		int pnt2 = cnt + i*jmp;
	 	Kg[pnt2] += ke[pnt];
	 	if (i==1)
	 	{	pnt = (i+4)*(6) + (i+3);
			pnt2 += 1;
		 	Kg[pnt2] += ke[pnt];	
	 	}
	 	if (i==2)
	 	{	pnt = (i+2)*6 + (i+3);
			pnt2 -= 1;
		 	Kg[pnt2] += ke[pnt];	
	 	}
	}
	// RHS Boundary
	for (int i = Nvar_-3; i < Nvar_; ++i)
	{	// Build Diagonals
		int pnt = (i - (Nvar_-3))*6 + (i - (Nvar_-3));
		int pnt2 = cnt + i*jmp;
	 	Kg[pnt2] += ke[pnt];	
	 	if (i==Nvar_-2)
	 	{	pnt += 6;
			pnt2 += 1;
		 	Kg[pnt2] += ke[pnt];
	 	}
	 	if (i==Nvar_-1)
	 	{	pnt -= 6;
			pnt2 -= 1;
		 	Kg[pnt2] += ke[pnt];	
	 	}
	}
}

// Builds global stiffness K in banded format across processes
void buildKgBandPar(double *Kg, double *ke, int Nvar_, int Nx_g, int buf_)
{	const int jmp = 9 + buf_;
	const int cnt = 4 + buf_;
	
	if (MPI::mpi_rank==0)
	{
		int ind = 0;
		for (int i = 0; i < Nvar_; ++i)
		{	// Build Diagonal
			int pnt = ind*6 + ind;
			int pnt2 = i*jmp + cnt;
			Kg[pnt2] += ke[pnt];
		 	int max = 6-ind;
		 	int bnd = 5;
		 	if (max < bnd)
		 	{	bnd = max;
		 	}
		 	// Build Lower
			for (int k = 1; k < bnd; ++k)
			{	int pos = pnt2+k;
			 	Kg[pos] += ke[pnt+k];	
			}
		 	// Build Upper
			for (int k = 1; k <= ind; ++k)
			{	int pos = pnt2-k;
				Kg[pos] += ke[pnt-k];
			}
			ind += 1;
			if (ind==6)
			{	ind = 0;
				i-=3;
			}
		}
		// Build first boundary
		for (int i = 3; i < 6; ++i)
		{	int pnt = i*6 + i;
			int pnt2 = (i-3)*jmp + cnt;
			Kg[pnt2] += ke[pnt];
			if (i==4)
			{	Kg[pnt2+1] += ke[pnt+1];
			}
			if (i==5)
			{	Kg[pnt2-1] += ke[pnt-1];
			}
		}
	}
	else if (MPI::mpi_rank==(MPI::mpi_size-1))
	{
		int ind = 3;
		for (int i = 0; i < Nvar_; ++i)
		{	// Build Diagonal
			int pnt = ind*6 + ind;
			int pnt2 = i*jmp + cnt;
			Kg[pnt2] += ke[pnt];
		 	int max = 6-ind;
		 	int bnd = 5;
		 	if (max < bnd)
		 	{	bnd = max;
		 	}
		 	// Build Lower
		 	if (i < Nvar_-3)
		 	{	for (int k = 1; k < bnd; ++k)
				{	int pos = pnt2+k;
				 	Kg[pos] += ke[pnt+k];	
				}
		 	}
		 	// Build Upper
			for (int k = 1; k <= ind; ++k)
			{	int pos = pnt2-k;
				Kg[pos] += ke[pnt-k];
			}
			ind += 1;
			if (ind==6)
			{	ind = 0;
				i-=3;
			}
		}	
	}
	else
	{	// Central Elements
		int ind = 3;
		for (int i = 0; i < Nvar_; ++i)
		{	// Build Diagonal
			int pnt = ind*6 + ind;
			int pnt2 = i*jmp + cnt;
			Kg[pnt2] += ke[pnt];
		 	int max = 6-ind;
		 	int bnd = 5;
		 	if (max < bnd)
		 	{	bnd = max;
		 	}
		 	// Build Lower
			for (int k = 1; k < bnd; ++k)
			{	int pos = pnt2+k;
			 	Kg[pos] += ke[pnt+k];	
			}
		 	// Build Upper
			for (int k = 1; k <= ind; ++k)
			{	int pos = pnt2-k;
				Kg[pos] += ke[pnt-k];
			}
			ind += 1;
			if (ind==6)
			{	ind = 0;
				i-=3;
			}
		}	
	}
}

// Builds global stiffness matrix K in banded
void buildKgBand(double *Kg, double *ke, int Nvar_, int Nx_g, int buf)
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

// Build global mass in banded format
void buildMgBand(double *Kg, double *ke, int Nvar_, int Nx_g, int buf)
{	const int jmp =  9 + buf;
	// LHS Boundary
	for (int i = 0; i < 3; ++i)
	{	// Build Diagonals
		int pnt = i*6 + i;
		int pnt2 = i;
	 	Kg[pnt2] += ke[pnt];
	}	
	// Central Elements
	for (int i = 0; i < Nx_g-2; ++i)
	{	
		for (int j = 0; j < 6; ++j)
		{	// Build Diagonal
			int pnt = j*6 + j;
			int pnt2 = j + 3*i;
		 	Kg[pnt2] += ke[pnt];
		}
	}
	// RHS Boundary
	for (int i = Nvar_-3; i < Nvar_; ++i)
	{	// Build Diagonals
		int pnt = (i - (Nvar_-3))*6 + (i - (Nvar_-3));
		int pnt2 = i;
	 	Kg[pnt2] += ke[pnt];	
	}
}

// Build global mass in banded format across all formats
void buildMgBandPar(double *Kg, double *ke, int Nvar_, int Nx_g, int buf)
{	const int jmp =  9 + buf;
	int iter = Nvar_/3;
	for (int i = 0; i < iter; ++i)
	{	for (int j = 0; j < 3; ++j)
		{	int pnt = j + 3*i;
			int pnt2 = j + j*6;
			Kg[pnt] += 2*ke[pnt2];
		}	
	}
}

// Build global force
void buildFgBand(double *Kg, double *ke, double coeff, int Nx_g, int Nvar_)
{	// Build Center	
	if (MPI::mpi_rank==0)
	{	Kg[1] = ke[1]; 
		Kg[2] = ke[2];
		for (int i = 1; i < Nvar_/3; ++i)
		{	int pnt = i*3 + 1;
			Kg[pnt] =  2*ke[1];
		}
	}
	if (MPI::mpi_rank>0)
	{	for (int i = 0; i < Nvar_/3; ++i)
		{	int pnt = i*3 + 1;
			Kg[pnt] =  2*ke[1];
		}
	}
	// Add the concentrated force to the middle node
	if (MPI::mpi_rank==((MPI::mpi_size)/2 -1) && MPI::mpi_size>1)
	{	Kg[(Nvar_-2)] += coeff*1000;
	}
	else if (MPI::mpi_size==1)
	{	Kg[(Nvar_-1)/2] += coeff*1000;
	}
}

void buildFg(double *Kg, double *ke, int Nx_g, int Nvar_)
{	// Build Center	
	for (int i = 0; i < Nx_g-1; ++i)
	{	for (int j = 0; j < 6; ++j)
		{	if (i==Nx_g-2 && j>=4)
			{	break;
			}
			else
			{	int pnt = i*3 + j;
				Kg[pnt] += ke[j];
			}
		}
	}
	Kg[(Nvar_-1)/2] +=1000;
}

void buildMe(double *K, double A_, double rho_, double lx_e, double Al_, double dt_)
// Brought the time integration into the mass matrix
{	const double p =  (rho_*A_*lx_e)/(dt_*dt_);
	int N = 6;
	for (int i = 0; i < N; ++i)
	{	if (i == 0)
		{	K[i*N+0] = p*.5;
		}
		if (i == 1)
		{	K[i*N+1] = p*.5;
		}
		if (i == 2)
		{	K[i*N+2] = p*Al_*pow(lx_e,2);
		}
		if (i == 3)
		{	K[i*N+3] = p*.5;
		}
		if (i == 4)
		{	K[i*N+4] = p*.5;
		}
		if (i == 5)
		{	K[i*N+5] = p*Al_*pow(lx_e,2);
		}
	}
}

void buildMe(double *K, double A_, double rho_, double lx_e, double Al_)
// Brought no time in the mass matrix
{	const double p =  (rho_*A_*lx_e);
	int N = 6;
	for (int i = 0; i < N; ++i)
	{	if (i == 0)
		{	K[i*N+0] = p*.5;
		}
		if (i == 1)
		{	K[i*N+1] = p*.5;
		}
		if (i == 2)
		{	K[i*N+2] = p*Al_*pow(lx_e,2);
		}
		if (i == 3)
		{	K[i*N+3] = p*.5;
		}
		if (i == 4)
		{	K[i*N+4] = p*.5;
		}
		if (i == 5)
		{	K[i*N+5] = p*Al_*pow(lx_e,2);
		}
	}
}

void buildKe(double *K, double lx_e, double A_, double E_, double I_)
{	double EI(E_*I_), a(A_*E_/lx_e), b(12*EI/pow(lx_e, 3)), 
			c(6*EI/pow(lx_e, 2)), d(4*EI/lx_e), e(2*EI/lx_e);
	int N = 6;
	for (int i = 0; i < N; ++i)
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

// Build elemental force vector
void buildFe(double *K, double lx_e, double qx_, double qy_)
{	K[0] = K[3] = lx_e*qx_/2;
	K[1] = K[4] = lx_e*qy_/2;
	K[2] = (qy_*lx_e*lx_e)/12;
	K[5] = -(qy_*lx_e*lx_e)/12;
}

// Build Kg - al*M
void buildKeff(double *K, double *K_eff, double *M, double co1_, int Nvar_, int buf_)
{	for (int i = 0; i < Nvar_; ++i)
	{	int pnt = (4+buf_) + i*(9+buf_);
		K_eff[pnt] += (co1_*M[i]);
		K[pnt] += (co1_*M[i]);
	}
}

// Assigns an array to a value
void assignArr(double *K, double V, int N)
{	fill(K, K+N, V);
}

void updateKglb(double *K, int Nvar_, int buf_)
{	if (MPI::mpi_rank==(MPI::mpi_size-1))
	{	int pnt = (Nvar_ - 3)*(9+buf_);
		int pnt2 = (Nvar_ - 6)*(9+buf_);
		for (int i = 0; i < 3*(9+buf_); ++i)
		{	K[pnt2+i] =  K[pnt+i];
		}
		for (int i = 0; i < 3*(9+buf_); ++i)
		{	K[pnt+i] =  0;
		}
		for (int i = Nvar_-3; i < Nvar_; ++i)
		{	K[i*(9+buf_) + 12] = 1;
		}
	}
}

// Updates Force vector for each timestep
void updateVars(double *F, double lx_e, double qx_, double qy_,
	int Nx_g, int Nvar_, int step, int nite_)
{	double coef = double(step)/nite_;
	double *fe = new double[6]();
	buildFe(fe, lx_e, coef*qx_, coef*qy_);
	buildFg(F, fe, Nx_g, Nvar_);
	delete [] fe;
}

// Updates Force vector for each timestep in parallel
void updateParVars(double *F, double lx_e, double qx_, double qy_,
	int Nx_g, int Nvar_, int step, int nite_)
{	double coef = double(step)/nite_;
	double *fe = new double[6]();
	buildFe(fe, lx_e, coef*qx_, coef*qy_);
	buildFgBand(F, fe, coef, Nx_g, Nvar_);
	delete [] fe;
}