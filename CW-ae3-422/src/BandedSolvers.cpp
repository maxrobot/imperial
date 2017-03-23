#include <fstream>
#include <string.h>
#include "Solvers.hpp"
#include "Common.hpp"
#include "CommonMPI.hpp"
#include "BuildFunction.hpp"
#include "Output.hpp"

using namespace std;

void solveSparseExplicit(double *K, double *M, double *F, double lx_e,
	double qx_, double qy_, int Nvar_, int Nx_g, int nite_, int nout_,
	int buf_, std::string test)
{	double *U 		= new double[Nvar_]();
	double *Un_g	= new double[Nvar_]();
	double *Mt_		= new double[Nvar_]();
	double *MKU_o	= new double[Nvar_]();

	int kl = 4;
	int ku = 4;
    int info = 0;
    int *ipiv = new int[Nvar_];

    // Make the initial matrix...
    for (int i = 0; i < Nvar_; ++i)
    {	int pnt = 4 + (i*9);
    	K[pnt] -= 2*M[i];
    }    	

	// =================== Create S Matrix ====================//
	// Start marching through time...
	for (int i = 0; i < nite_; ++i)
	{	F77NAME(dcopy)(Nvar_, M, 1, Mt_, 1); // Save Un-1
		assignArr(MKU_o, 0, Nvar_);
		
		// Calculate MK_o*U{n}
		F77NAME(dgbmv)('n', Nvar_, Nvar_, 4, 4, 1, K, 9, U, 1, 0, MKU_o, 1);

		// Update Variables
		assignArr(F, 0., Nvar_);
		updateVars(F, lx_e, qx_, qy_, Nx_g, Nvar_, i, nite_);

		// Calculate updated M*U{n+1} = S
		for (int i = 0; i < Nvar_; ++i)
		{	double sum = F[i] - MKU_o[i] - Mt_[i]*Un_g[i];
			Un_g[i] = U[i]; // Update the old variable...
			U[i] = sum;
		}

	    F77NAME(dgbsv)(Nvar_, 0, 0, 1, Mt_, 1, ipiv, U, Nvar_, info);
	}
	writeVec(U, Nx_g, 1, test);
}

void solveParSparseExplicit(double *K, double *M, double *F, double lx_e,
	double qx_, double qy_, int Nvar_g, int Nvar_, int Nghost_, int Sghost_,
	int Nx_g, int Nx_, int nite_, int nout_, int buf_, std::string test)
{	// Make the initial matrix, and free unused k before continuing...
	double *MK_o	= new double[Nghost_ * 9]();
	parMatSum(K, M, MK_o, Nvar_, Nghost_);
    delete[] K;

    // Make the rest...
	double *U 		= new double[Nghost_]();
	double *MKU_o	= new double[Nghost_]();
	double *Un_g	= new double[Nvar_]();
	double *Minv_	= new double[Nvar_]();
	double d1, d2;
    for (int i = 0; i < Nvar_; ++i)
    {	Minv_[i] = 1/M[i];
    }

	// =================== Create S Matrix ====================//
	// Start marching through time...
	for (int i = 0; i < nite_; ++i)
	{	
		// Calculate MK_o*U{n}
		F77NAME(dgbmv)('n', Nghost_, Nghost_, 4, 4, 1, MK_o, 9, U, 1, 0, MKU_o, 1);
		MPI::exchangeVecConts(MKU_o, Nghost_, Sghost_);

		// Update Variables
		assignArr(F, 0., Nvar_);
		updateParVars(F, lx_e, qx_, qy_, Nx_g, Nvar_, i, nite_);

		for (int i = 0; i < Nvar_; ++i)
		{	F[i] -= M[i]*Un_g[i];
		}
		parVecSub(F, MKU_o, Nvar_);

		// Update old U into Un_g
		parVecCopy(Un_g, U, Nvar_);

		// Calculate updated M*U{n+1} = S
		parMatSolve(Minv_, F, U, Nvar_);
		MPI::copyVecConts(U, Nghost_);
	}
	parVecCopy(Un_g, U, Nvar_);
	writeParVec(Un_g, Nx_g, Nvar_g, Nvar_, 1, test);
}

void solveSparseImplicit(double *K, double *M, double *F, double lx_e, double qx_,
	double qy_, double dt_, int Nvar_, int Nx_g, int nite_, int nout_, int buf_,
	std::string test)
{	// ===================== Build Tables =====================//
	double *K_eff 		= new double[Nvar_ * (9+buf_)]();
	double *U			= new double[Nvar_]();
	double *Ud			= new double[Nvar_]();
	double *Udd			= new double[Nvar_]();
	double *tmp			= new double[Nvar_]();
	double sum;

	// Constanst and coefficients
	const double beta_(.25);
	const double gamma_(.5);
	const double co1_(1/(beta_*(dt_*dt_)));
	const double co2_(1/(beta_*dt_));
	const double co3_((1/(2*beta_))-1);
	const double co4_(dt_*(1-gamma_));
	const double co5_(dt_*gamma_);

	F77NAME(dcopy)(Nvar_*(9+buf_), K, 1, K_eff, 1);
	// Build Keff
	for (int i = 0; i < Nvar_; ++i)
	{	int pnt = (4+buf_) + i*(9+buf_);
		K_eff[pnt] += (co1_*M[i]);
		K[pnt] += (co1_*M[i]);
	}

    int info = 0;
    int *ipiv = new int[Nvar_];

	for (int i = 0; i < nite_; ++i)
	{	// Create Dynamic Force
		updateVars(F, lx_e, qx_, qy_, Nx_g, Nvar_, i, nite_);

		// Sum U, Ud and, Udd with coefficients, the add forces to sum...
		for (int i = 0; i < Nvar_; ++i)
		{	sum = ((co1_*U[i]) + (co2_*Ud[i]) + (co3_*Udd[i]));
			F[i] += M[i]*sum;
		}

		// Solve Keff*U_{n+1} = S
	    F77NAME(dgbsv)(Nvar_, 4, 4, 1, K, 9+buf_, ipiv, F, Nvar_, info);

		// Update K to contain only the K_eff as desgv overwrites...
		F77NAME(dcopy)(Nvar_*(9+buf_), K_eff, 1, K, 1);

		// Now update Udd
		F77NAME(dcopy)(Nvar_, Udd, 1, tmp, 1);   
		for (int i = 0; i < Nvar_; ++i)
		{	sum = (co2_*Ud[i]) + (co3_*Udd[i]);
			Udd[i] = (co1_*(F[i]-U[i])) - sum;
		}

		// Now update Ud
		for (int i = 0; i < Nvar_; ++i)
		{	sum = (co4_*tmp[i]) + (co5_*Udd[i]);
			Ud[i] = Ud[i] + sum;
		}

		// Now update U
		F77NAME(dcopy)(Nvar_, F, 1, U, 1);
	}
	writeVec(U, Nx_g, 1, test);
}

void solveParSparseImplicit(double *K, double *M, double *F, double lx_e, double qx_,
	double qy_, double dt_, int Nvar_g, int Nvar_, int Nghost_, int Nx_g, int Nx_,
	int nite_, int nout_, int buf_, std::string test)
{	// ===================== Build Tables =====================//
	double *K_eff 		= new double[Nvar_ * (9+buf_)]();
	double *U			= new double[Nvar_]();
	double *Ud			= new double[Nvar_]();
	double *Udd			= new double[Nvar_]();
	double *tmp			= new double[Nvar_]();
	double sum;

	// Constanst and coefficients
	const double beta_(.25);
	const double gamma_(.5);
	const double co1_(1/(beta_*(dt_*dt_)));
	const double co2_(1/(beta_*dt_));
	const double co3_((1/(2*beta_))-1);
	const double co4_(dt_*(1-gamma_));
	const double co5_(dt_*gamma_);

    const int kl   = 4;             // Number of lower diagonals
    const int ku   = 4;             // Number of upper diagonals
    const int lda  = 1 + 2*kl + 2*ku; // Leading dimension (num of rows)
    const int ldb  = Nvar_;
    const int nrhs = 1;
    const int ja   = 1;             // Offset index in global array (col)
    const int ib   = 1;             // Offset index in global array (row)
    const int lwork= (Nvar_+ku)*(kl+ku)+6*(kl+ku)*(kl+2*ku)
                        + max(nrhs*(Nvar_+2*kl+4*ku), 1);

    // Allocate factorisation matrices...
    double* wk   = new double[lwork]();   // Local workspace
    int*    ipiv = new int[Nghost_]();         // Local pivoting array
    int info = 0;

	// Build Keff
	updateKglb(K, Nvar_, buf_);
	F77NAME(dcopy)(Nvar_*(9+buf_), K, 1, K_eff, 1);
	buildKeff(K, K_eff, M, co1_, Nvar_, buf_);
    
    // Create descriptors for matrix and RHS vector storage
    int desca[7];
    desca[0] = 501;         // Type is a banded matrix 1-by-P
    desca[1] = MPI::ctx;    // Context
    desca[2] = Nghost_;     // Problem size
    desca[3] = Nvar_;       // Blocking
    desca[4] = 0;           // Process row/column
    desca[5] = lda;         // Local size
    desca[6] = 0;           // Reserved

    int descb[7];
    descb[0] = 502;         // Type is a banded matrix P-by-1 (RHS)
    descb[1] = MPI::ctx;    // Context
    descb[2] = Nghost_;     // Problem size
    descb[3] = Nvar_;       // Blocking
    descb[4] = 0;           // Process row/column
    descb[5] = Nvar_;       // Local size
    descb[6] = 0;           // Reserved
    
	for (int i = 0; i < nite_; ++i)
	{	// Create Dynamic Force
		updateParVars(F, lx_e, qx_, qy_, Nx_g, Nvar_, i, nite_);

		// Update RHS
		for (int i = 0; i < Nvar_; ++i)
		{	sum = ( (co1_*U[i]) + (co2_*Ud[i]) + (co3_*Udd[i]) );
			F[i] += M[i]*sum;
		}

	    // Solver for parallel system matrix vector (RHS vector replaced by solution)
	    F77NAME(pdgbsv)(Nghost_, kl, ku, nrhs, K, ja, desca, ipiv, F, ib, descb, wk, lwork, &info);

		// Update K to contain only the K_eff as dgbsv overwrites...
		F77NAME(dcopy)(Nvar_*(9+buf_), K_eff, 1, K, 1);

		// Now update Udd
		F77NAME(dcopy)(Nvar_, Udd, 1, tmp, 1);
		for (int i = 0; i < Nvar_; ++i)
		{	sum = (co2_*Ud[i]) + (co3_*Udd[i]);
			Udd[i] = (co1_*(F[i]-U[i])) - sum;
		}

		// Now update Ud
		for (int i = 0; i < Nvar_; ++i)
		{	sum = (co4_*tmp[i]) + (co5_*Udd[i]);
			Ud[i] = Ud[i] + sum;
		}

		// Now update U
		F77NAME(dcopy)(Nvar_, F, 1, U, 1);
	}
	writeParVec(U, Nx_g, Nvar_*MPI::mpi_size, Nvar_, 1, test);
}