#include "Solvers.hpp"
#include "Common.hpp"
#include "CommonMPI.hpp"
#include "BuildFunction.hpp"
#include "Output.hpp"

using namespace std;

// Static Solver Launch
void runSolver(double lx_e, double A_,
	double E_, double I_, double qx_, double qy_, int Nvar_,
	int Nx_g,std::string scheme_)
{	printInfo(Nx_g, 0, Nvar_, 0, scheme_);
	// ================ Initialise Local Vars. ================//
	const int buf_(4);	  	// Buffer
	if (MPI::mpi_size>1)
	{	printMessage("ERROR: static is single processor only!");
	}
	// Matrices
	double *F_g	= new double[Nvar_]();
	double *K_g	= new double[Nvar_ * (9+buf_)]();
	double *K_e	= new double[6*6]();
	double *F_e	= new double[6]();

	// ============== Create Elemental K Matrix ===============//
	buildKe(K_e, lx_e, A_, E_, I_);
	buildFe(F_e, lx_e, qx_, qy_);
	
	buildKgBand(K_g, K_e, Nvar_, Nx_g, buf_);
	buildFg(F_g, F_e, Nx_g, Nvar_);

	// ==================== Free Memory =======================//
	delete[] K_e;
	delete[] F_e;

	// ==================== Run Solver ========================//
	solveStatic(K_g, F_g, Nvar_, 9+buf_, Nx_g, "task_1");
}

// Dynamic Explicit Launch
void runSolver(double dt_, double lx_e, double A_,
	double E_, double I_, double rho_, double qx_, double qy_,
	int Nvar_, int Nvar_e, int Nghost_, int Sghost_, int Nx_g,
	int Nx_, int nite_, int nout_, const int buf_, string sparse_,
	string scheme_)
{	printInfo(Nx_g, Nx_, Nvar_, Nvar_e, scheme_);
	// ================ Initialise Local Vars. ================//
	const double Al_(1./24);  	// Constant Alpha

	double *F_g	= new double[Nvar_e]();
	double *K_e	= new double[6*6]();
	double *M_e	= new double[6*6]();

	// ============== Create Elemental K Matrix ===============//
	if (buf_ == 4 || buf_ == 8)
	{	buildMe(M_e, A_, rho_, lx_e, Al_);
	}
	if (buf_ == 0)	
	{	buildMe(M_e, A_, rho_, lx_e, Al_, dt_);
	}
	
	buildKe(K_e, lx_e, A_, E_, I_);

	if (sparse_=="none")
	{	// Matrices
		double *K_g = new double[Nvar_e*Nvar_e]();
		double *M_g	= new double[Nvar_e*Nvar_e]();

	// ===================== Build Tables =====================//
		buildKg(K_g, K_e, Nvar_e, Nx_g);
		buildKg(M_g, M_e, Nvar_e, Nx_g);

	// ==================== Free Memory =======================//
		delete[] K_e;
		delete[] M_e;

		solveExplicit(K_g, M_g, F_g, lx_e, qx_, qy_, Nvar_e,
			Nx_g, nite_, nout_, buf_, "task_2");
	}
	else if (sparse_=="sparse")
	{	// Matrices
		double *K_g = new double[Nvar_e*(9+buf_)]();
		double *M_g	= new double[Nvar_e]();

		if (MPI::mpi_size==1)
		{	// ===================== Build Tables =====================//
			buildKgBand(K_g, K_e, Nvar_e, Nx_g, buf_);
			buildMgBand(M_g, M_e, Nvar_e, Nx_g, buf_);

			// ==================== Free Memory =======================//
			delete[] K_e;
			delete[] M_e;

			// ==================== Run Solver ========================//
			if (buf_==4)
			{	solveSparseImplicit(K_g, M_g, F_g, lx_e, qx_, qy_,
					dt_, Nvar_e, Nx_g, nite_, nout_, buf_, "task_3");
			}
			else
			{	solveSparseExplicit(K_g, M_g, F_g, lx_e, qx_, qy_,
					Nvar_e, Nx_g, nite_, nout_, buf_,"task_2");
			}
		}

		else if (MPI::mpi_size>1)
		{	// ===================== Build Tables =====================//
			buildKgBandPar(K_g, K_e, Nvar_e, Nx_g, buf_);
			buildMgBandPar(M_g, M_e, Nvar_e, Nx_g, buf_);
			
			// ==================== Free Memory =======================//
			delete[] K_e;
			delete[] M_e;

			// ==================== Run Solver ========================//
			if (buf_==8)
			{	solveParSparseImplicit(K_g, M_g, F_g, lx_e, qx_, qy_,
					dt_, Nvar_, Nvar_e, Nghost_, Nx_g, Nx_, nite_,
					nout_, buf_, "task_5");
			}
			else
			{	solveParSparseExplicit(K_g, M_g, F_g, lx_e, qx_, qy_,
					Nvar_, Nvar_e, Nghost_, Sghost_, Nx_g, Nx_, nite_,
					nout_, buf_,"task_4");
			}
		}
	}
}

void solveStatic(double *K, double *F, int Nvar_, int ldab, int Nx_, std::string test)
{	const int nrhs = 1;
    int info = 0;
    int *ipiv = new int[Nvar_];
    int kl = 4;
    int ku = 4;
    int ldb = Nvar_;

    // Use blas to solve system...
    F77NAME(dgbsv)(Nvar_, kl, ku, nrhs, K, ldab, ipiv, F, ldb, info);
    writeVec(F, Nx_, test);
	delete[] ipiv;
}

void solveExplicit(double *K, double *M, double *F, double lx_e,
	double qx_, double qy_, int Nvar_, int Nx_g, int nite_, int nout_,
	int buf_, std::string test)
{	double *U 		= new double[Nvar_]();
	double *MK_o	= new double[Nvar_ * Nvar_]();
	double *Un_g	= new double[Nvar_]();
	double *S		= new double[Nvar_]();
	double *MKU_o	= new double[Nvar_]();
	double *MU_o	= new double[Nvar_]();

    int info = 0;
    int *ipiv = new int[Nvar_];

    for (int i = 0; i < Nvar_*Nvar_; ++i)
    {	MK_o[i] = K[i] - (2*M[i]);
    }    	

	// =================== Create S Matrix ====================//
	// Start marching through time...
	for (int i = 0; i < nite_; ++i)
	{	// Calculate MK_o*U{n}
		F77NAME(dgemv)('n', Nvar_, Nvar_, 1, MK_o, Nvar_, U, 1, 0, MKU_o, 1);
		
		assignArr(F, 0., Nvar_);
		updateVars(F, lx_e, qx_, qy_, Nx_g, Nvar_, i, nite_);

		// Calculate M*U{n-1}
		for (int i = 0; i < Nvar_; ++i)
		{	int pnt = i*Nvar_ + i;
			MU_o[i] = M[pnt]*Un_g[i];
		}

		for (int i = 0; i < Nvar_; ++i)
		{	S[i] = F[i] - MKU_o[i] - MU_o[i];
		}

		// Calculate updated M*U{n+1} = S
		F77NAME(dgesv)(Nvar_, 1, M, Nvar_, ipiv, S, Nvar_, info);
		// Now update vars for repeat
		F77NAME(dcopy)(Nvar_, U, 1, Un_g, 1); // Save Un-1
		F77NAME(dcopy)(Nvar_, S, 1, U, 1); // Update Un
	}
	writeVec(U, Nx_g, 1, test);
}

void solveImplicit(double *K, double *M, double *F, double lx_e, double qx_,
	double qy_, double dt_, int Nvar_, int Nx_g, int nite_, int nout_,
    std::string test)
{	// ===================== Build Tables =====================//
	double *K_eff 		= new double[Nvar_ * Nvar_]();
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


	// Build Keff
	for (int i = 0; i < Nvar_*Nvar_; ++i)
	{	double sum = (co1_*M[i]) + K[i];
		K_eff[i] = sum;
		K[i] = sum;
	}

    int info = 0;
    int *ipiv = new int[Nvar_];

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
		F77NAME(dgemv)('n', Nvar_, Nvar_, 1, M, Nvar_, tmp2, 1, 0, S, 1);
		// Add forces to sum...
		for (int i = 0; i < Nvar_; ++i)
		{	double sum = S[i] + F[i];
			S[i] = sum;
		}

		// Solve Keff*U_{n+1} = S
		F77NAME(dgesv)(Nvar_, 1, K, Nvar_, ipiv, S, Nvar_, info);

		// Update K to contain only the K_eff as desgv overwrites...
		for (int i = 0; i < Nvar_*Nvar_; ++i)
		{	K[i] = K_eff[i];
		}

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

void parMatSum(double *K, double *M, double *MK, int Nvar_, int Nghost_)
{	int Sghost_ = 3;
	int cnt = 4;
	if (MPI::mpi_rank==0)
	{	// Copy K into MK
		for (int i = 0; i < Nvar_*9; ++i)
		{	MK[i] = K[i];
		}
		for (int i = 0; i < Nvar_; ++i)
	    {	int pnt = cnt + (i*9);
	    	MK[pnt] -= 2*M[i];
	    }    	
	}
	else if (MPI::mpi_rank!=0)
	{	// Copy K into MK
		for (int i = 0; i < Nvar_*9; ++i)
		{	MK[i+(Sghost_*9)] = K[i];
		}
		for (int i = 0; i < Nvar_; ++i)
	    {	int pnt = cnt + (i+Sghost_)*9 ;
	    	MK[pnt] -= 2*M[i];
	    }    	
	}
}

void parVecCopy(double *M, double *N, int Nvar_)
{	int Sghost_ = 3;
	if (MPI::mpi_rank==0)
	{	for (int i = 0; i < Nvar_; ++i)
		{	M[i] = N[i];
		}
	}
	if (MPI::mpi_rank==(MPI::mpi_size-1))
	{	for (int i = 0; i < Nvar_; ++i)
		{	M[i] = N[i+Sghost_];
		}
	}
	if (MPI::mpi_rank>0 && MPI::mpi_rank<(MPI::mpi_size-1))
	{	for (int i = 0; i < Nvar_; ++i)
		{	M[i] = N[i+Sghost_];
		}
	}
}

void parVecAdd(double *M, double *N, int Nvar_)
{	int Sghost_ = 3;
	if (MPI::mpi_rank==0)
	{	for (int i = 0; i < Nvar_; ++i)
		{	M[i] += N[i];
		}
	}
	if (MPI::mpi_rank==(MPI::mpi_size-1))
	{	for (int i = 0; i < Nvar_; ++i)
		{	M[i] += N[i+Sghost_];
		}
	}
	if (MPI::mpi_rank>0 && MPI::mpi_rank<(MPI::mpi_size-1))
	{	for (int i = 0; i < Nvar_; ++i)
		{	M[i] += N[i+Sghost_];
		}
	}
}

void parVecSub(double *M, double *N, int Nvar_)
{	int Sghost_ = 3;
	if (MPI::mpi_rank==0)
	{	for (int i = 0; i < Nvar_; ++i)
		{	M[i] -= N[i];
		}
	}
	if (MPI::mpi_rank==(MPI::mpi_size-1))
	{	for (int i = 0; i < Nvar_; ++i)
		{	M[i] -= N[i+Sghost_];
		}
	}
	if (MPI::mpi_rank>0 && MPI::mpi_rank<(MPI::mpi_size-1))
	{	for (int i = 0; i < Nvar_; ++i)
		{	M[i] -= N[i+Sghost_];
		}
	}
}

void parMatSolve(double *M, double *S, double *U, int Nvar_)
{	int Sghost_ = 3;
	if (MPI::mpi_rank==0)
	{	for (int i = 0; i < Nvar_; ++i)
		{	U[i] = M[i]*S[i];
		}
	}
	if (MPI::mpi_rank==(MPI::mpi_size-1))
	{	for (int i = 0; i < Nvar_; ++i)
		{	U[i+Sghost_] = M[i]*S[i];
		}
	}
	if (MPI::mpi_rank>0 && MPI::mpi_rank<(MPI::mpi_size-1))
	{	for (int i = 0; i < Nvar_; ++i)
		{	U[i+Sghost_] = M[i]*S[i];
		}
	}
}