#include "Solvers.hpp"
#include "Common.hpp"
#include "CommonMPI.hpp"
#include "BuildFunction.hpp"

using namespace std;

void parMatSum(double *K, double *M, double *MK, int Nvar_, int Nghost_)
{	if (MPI::mpi_rank==0)
	{	// Copy K into MK
		for (int i = 0; i < Nvar_*9; ++i)
		{	MK[i] = K[i];
		}
		for (int i = 0; i < Nvar_; ++i)
	    {	int pnt = 4 + (i*9);
	    	MK[pnt] -= 2*M[i];
	    }    	
	}
	else if (MPI::mpi_rank!=0)
	{	// Copy K into MK
		for (int i = 0; i < Nvar_*9; ++i)
		{	MK[i+9] = K[i];
		}
		for (int i = 0; i < Nvar_; ++i)
	    {	int pnt = 4 + (i+1)*9 ;
	    	MK[pnt] -= 2*M[i];
	    }    	
	}
}

void parVecCopy(double *M, double *N, int Nvar_)
{	if (MPI::mpi_rank==0)
	{	for (int i = 0; i < Nvar_; ++i)
		{	M[i] = N[i];
		}
	}
	if (MPI::mpi_rank==(MPI::mpi_size-1))
	{	for (int i = 0; i < Nvar_; ++i)
		{	M[i] = N[i+1];
		}
	}
	if (MPI::mpi_rank>0 && MPI::mpi_rank>(MPI::mpi_size-1))
	{	for (int i = 0; i < Nvar_; ++i)
		{	M[i] = N[i+1];
		}
	}
}

void parVecSum(double *M, double *N, double al, int Nvar_)
{	if (MPI::mpi_rank==0)
	{	for (int i = 0; i < Nvar_; ++i)
		{	M[i] += al*N[i];
		}
	}
	if (MPI::mpi_rank==(MPI::mpi_size-1))
	{	for (int i = 0; i < Nvar_; ++i)
		{	M[i] += al*N[i+1];
		}
	}
	if (MPI::mpi_rank>0 && MPI::mpi_rank>(MPI::mpi_size-1))
	{	for (int i = 0; i < Nvar_; ++i)
		{	M[i] += al*N[i+1];
		}
	}
}

void parMatSolve(double *M, double *N, double *S, int Nvar_)
{	if (MPI::mpi_rank==0)
	{	for (int i = 0; i < Nvar_; ++i)
		{	S[i] = M[i]*N[i];
		}
	}
	if (MPI::mpi_rank==(MPI::mpi_size-1))
	{	for (int i = 0; i < Nvar_; ++i)
		{	S[i+1] += M[i]*N[i];
		}
	}
	if (MPI::mpi_rank>0 && MPI::mpi_rank>(MPI::mpi_size-1))
	{	for (int i = 0; i < Nvar_; ++i)
		{	S[i+1] += M[i]*N[i];
		}
	}
}

void solveSparseExplicit(double *K, double *M, double *F, double lx_e,
	double qx_, double qy_, int Nvar_, int Nx_g, int nite_, int nout_,
	int buf_, std::string test)
{	
	double *MK_o	= new double[Nvar_ * 9]();
	double *U 		= new double[Nvar_]();
	double *Un_g	= new double[Nvar_]();
	double *S		= new double[Nvar_]();
	double *Mt_		= new double[Nvar_]();
	double *MKU_o	= new double[Nvar_]();
	double *MU_o	= new double[Nvar_]();

	int kl = 4;
	int ku = 4;
    int info = 0;
    int *ipiv = new int[Nvar_];

    // Make the initial matrix...
	F77NAME(dcopy)(Nvar_*(9+buf_), K, 1, MK_o, 1);

    for (int i = 0; i < Nvar_; ++i)
    {	int pnt = 4 + (i*9);
    	MK_o[pnt] -= 2*M[i];
    }    	

	// =================== Create S Matrix ====================//
	// Start marching through time...
	// for (int i = 0; i <= nite_; ++i)
	for (int i = 0; i < 3; ++i)
	{	F77NAME(dcopy)(Nvar_, M, 1, Mt_, 1); // Save Un-1
		assignArr(MKU_o, 0., Nvar_);
		
		// Calculate MK_o*U{n}
		F77NAME(dgbmv)('n', Nvar_, Nvar_, 4, 4, 1, MK_o, 9, U, 1, 0, MKU_o, 1);
		// Update Variables
		assignArr(F, 0., Nvar_);
		updateVars(F, lx_e, qx_, qy_, Nx_g, Nvar_, i, nite_);

		// Calculate M*U{n-1}
		for (int i = 0; i < Nvar_; ++i)
		{	MU_o[i] = Mt_[i]*Un_g[i];
		}

		for (int i = 0; i < Nvar_; ++i)
		{	double sum = F[i] - MKU_o[i] - MU_o[i];
			S[i] = sum;
		}

		// Calculate updated M*U{n+1} = S
	    F77NAME(dgbsv)(Nvar_, 0, 0, 1, Mt_, 1, ipiv, S, Nvar_, info);
		showVec(S, Nvar_);
		// showVec(F, Nvar_);
		// showVec(Un_g, Nvar_);
		// showMat(MK_o,9, Nvar_);
		// Now update vars for repeat
		F77NAME(dcopy)(Nvar_, U, 1, Un_g, 1); // Save Un-1
		F77NAME(dcopy)(Nvar_, S, 1, U, 1); // Update Un
	}
	writeVec(U, Nx_g, 1, test);
}

void solveParSparseExplicit(double *K, double *M, double *F, double lx_e,
	double qx_, double qy_, int Nvar_g, int Nvar_, int Nghost_, int Nx_g,
	int Nx_, int nite_, int nout_, int buf_, std::string test)
{	
	double *MK_o	= new double[Nghost_ * 9]();
	double *U 		= new double[Nghost_]();
	double *MKU_o	= new double[Nghost_]();
	double *Un_g	= new double[Nvar_]();
	double *S		= new double[Nvar_]();
	double *Mt_		= new double[Nvar_]();
	double *Minv_	= new double[Nvar_]();
	double *MU_o	= new double[Nvar_]();
	double d1, d2;

	int kl = 0;
	int ku = 0;
    int info = 0;
    const int lda = 1 + 2*kl + 2*ku;
    const int ldb  = Nvar_;
    int nrhs = 1;
    const int lwork= (Nvar_+ku)*(kl+ku)+6*(kl+ku)*(kl+2*ku)
                    + max(nrhs*(Nvar_+2*kl+4*ku), 1);
    double* wk   = new double[lwork]();   // Local workspace
    int *ipiv = new int[Nvar_];

    // Make the initial matrix...
	parMatSum(K, M, MK_o, Nvar_, Nghost_);
    for (int i = 0; i < Nvar_; ++i)
    {	Minv_[i] = 1/M[i];
    }    	


	// =================== Create S Matrix ====================//
	// Start marching through time...
	for (int i = 0; i <= nite_; ++i)
	// for (int i = 0; i < 3; ++i)
	{	F77NAME(dcopy)(Nvar_, M, 1, Mt_, 1); // Save Un-1
		assignArr(MKU_o, 0., Nghost_);

		
		// Calculate MK_o*U{n}
		F77NAME(dgbmv)('n', Nghost_, Nghost_, 4, 4, 1, MK_o, 9, U, 1, 0, MKU_o, 1);
		MPI::exchangeVecConts(MKU_o, Nghost_);
		MPI::exchangeVecConts(U, Nghost_);
		// for (int i = 0; i < MPI::mpi_size; ++i)

		// Update Variables
		assignArr(F, 0., Nvar_);
		updateParVars(F, lx_e, qx_, qy_, Nx_g, Nvar_, i, nite_);

		// Calculate M*U{n-1}

		for (int i = 0; i < Nvar_; ++i)
		{	double sum = F[i] - MU_o[i];
			S[i] = sum;
		}
		parVecSum(S, MKU_o, -1., Nvar_);
		for (int i = 0; i < Nvar_; ++i)
		{	MU_o[i] = Mt_[i]*Un_g[i];
		}

		// Calculate updated M*U{n+1} = S
		parMatSolve(Minv_, S, U, Nvar_);

		// for (int i = 0; i < MPI::mpi_size; ++i)
		// {	if (MPI::mpi_rank==0)
		// 	{	
		// 		// showMat(MK_o, 9, Nghost_);
		// 		// showVec(MKU_o, Nvar_);
		// 		// showVec(F, Nvar_);
		// 		// showVec(Un_g, Nvar_);
		// 	}
		// }
		// Now update vars for repeat
		// F77NAME(dcopy)(Nvar_, U, 1, Un_g, 1); // Save Un-1
		parVecCopy(Un_g, U, Nvar_);

		F77NAME(dcopy)(Nvar_, S, 1, U, 1); // Update Un
	}
	writeVec(U, Nx_g, 1, test);
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
	double *tmp2		= new double[Nvar_]();
	double *S			= new double[Nvar_]();

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

	// for (int i = 0; i < 2; ++i)
	for (int i = 0; i < nite_; ++i)
	{	// Create Dynamic Force
		assignArr(F, 0., Nvar_);
		assignArr(S, 0., Nvar_);
		assignArr(tmp, 0., Nvar_);
		assignArr(tmp2, 0., Nvar_);
		updateVars(F, lx_e, qx_, qy_, Nx_g, Nvar_, i, nite_);


		// Sum U, Ud and, Udd with coefficients
		for (int i = 0; i < Nvar_; ++i)
		{	double sum = (co1_*U[i]) + (co2_*Ud[i]) + (co3_*Udd[i]);
			tmp2[i] = sum;
		}

		// Multiple mass by sum
		for (int i = 0; i < Nvar_; ++i)
		{	S[i] = M[i]*tmp2[i];
		}

		// Add forces to sum...
		for (int i = 0; i < Nvar_; ++i)
		{	double sum = S[i] + F[i];
			S[i] = sum;
		}

		// Solve Keff*U_{n+1} = S
	    F77NAME(dgbsv)(Nvar_, 4, 4, 1, K, 9+buf_, ipiv, S, Nvar_, info);

		// Update K to contain only the K_eff as desgv overwrites...
		F77NAME(dcopy)(Nvar_*(9+buf_), K_eff, 1, K, 1);

		// Now update Udd
		F77NAME(dcopy)(Nvar_, Udd, 1, tmp, 1);   
		for (int i = 0; i < Nvar_; ++i)
		{	double sum = (co1_*(S[i]-U[i])) - (co2_*Ud[i]) -
				(co3_*Udd[i]);
			Udd[i] = sum;
		}

		// Now update Ud
		for (int i = 0; i < Nvar_; ++i)
		{	double sum = U[i] + (co4_*tmp[i]) + (co5_*Udd[i]);
			Ud[i] = sum;
		}

		// Now update U
		F77NAME(dcopy)(Nvar_, S, 1, U, 1);
	}
	writeVec(U, Nx_g, 1, test);
}