#include "Solvers.hpp"
#include "Common.hpp"
#include "CommonMPI.hpp"
#include "BuildFunction.hpp"

using namespace std;

// Static Solver Launch
void runSolver(double *K_e, double *U_g, double *F_g, double lx_e,
	double A_, double E_, double I_, double qx_, double qy_, int Nvar_,
	int Nx_g)
{	// ================ Initialise Local Vars. ================//
	const int buf_(4);	  	// Buffer
	if (MPI::mpi_size>1)
	{	printMessage("ERROR: static is single processor only!");
	}
	
	// Matrices
	double *K_g	= new double[Nvar_ * (9+buf_)]();
	double *F_e	= new double[6]();

	// ============== Create Elemental K Matrix ===============//
	buildKele(K_e, lx_e, A_, E_, I_);
	buildFele(F_e, lx_e, qx_, qy_);
	
	buildKglbSparse(K_g, K_e, Nvar_, Nx_g, buf_);
	buildFglb(F_g, F_e, Nx_g, Nvar_);

	// =================== Solve System =======================//
	solveStatic(K_g, F_g, Nvar_, 9+buf_, Nx_g, "task1_");
}

// Dynamic Explicit Launch
void runSolver(double *K_e, double *U_g, double *F_g, double dt_, 
	double lx_e, double A_, double E_, double I_, double rho_,
	double qx_,	double qy_,	int Nvar_, int Nvar_e, int Nx_g, int nite_,
	int nout_, const int buf_, string sparse_)
{	// ================ Initialise Local Vars. ================//
	const double Al_(1./24);  	// Constant Alpha

	double *M_e	= new double[6*6]();

	// ============== Create Elemental K Matrix ===============//
	if (buf_ == 4)
	{	buildMele(M_e, A_, rho_, lx_e, Al_);
	}
	if (buf_ == 0)	
	{	buildMele(M_e, A_, rho_, lx_e, Al_, dt_);
	}
	
	buildKele(K_e, lx_e, A_, E_, I_);

	if (sparse_=="none")
	{	// Matrices
		double *K_g = new double[Nvar_*Nvar_]();
		double *M_g	= new double[Nvar_*Nvar_]();

	// ===================== Build Tables =====================//
		buildKglb(K_g, K_e, Nvar_, Nx_g);
		buildKglb(M_g, M_e, Nvar_, Nx_g);

		solveExplicit(K_g, M_g, F_g, U_g, lx_e, qx_, qy_, Nvar_,
			Nx_g, nite_, nout_, buf_,"task2_");
	}
	else if (sparse_=="sparse")
	{	// Matrices
		double *K_ 	= new double[Nvar_*(9+buf_)]();
		double *K_g = new double[Nvar_e*(9+buf_)]();
		double *M_g	= new double[Nvar_]();
		
	// ===================== Build Tables =====================//
		if (MPI::mpi_size==1)
		{	buildSparse(K_g, K_e, Nvar_, Nx_g, buf_);
			buildMglbSparse(M_g, M_e, Nvar_, Nx_g, buf_);
		}
		else if (MPI::mpi_size>1)
		{	buildSparse(K_, K_e, Nvar_, Nx_g, buf_);
			buildBandSparse(K_g, K_e, Nvar_e, Nx_g, buf_);
		}

		if (MPI::mpi_rank==0)
		{	//showMat(K_e, 6);
			showMat(K_, (9+buf_), Nvar_);
			// showMat(K_g, (9+buf_), Nvar_);
			MPI_Barrier;
		}
		showMat(K_g, (9+buf_), Nvar_e);
		// for (int i = 0; i < MPI::mpi_size; ++i)
		// {	if (MPI::mpi_rank==i)
		// 	{	showMat(K_g, (9+buf_), Nvar_);
		// 		MPI_Barrier;
		// 	}
		// }
		
	// ==================== Run Solver ========================//
		if (buf_==4)
		{	solveSparseImplicit(K_g, M_g, F_g, U_g, lx_e, qx_, qy_,
				dt_, Nvar_, Nx_g, nite_, nout_, buf_, "task3_");
		}
		else
		{	solveSparseExplicit(K_g, M_g, F_g, U_g, lx_e, qx_, qy_,
				Nvar_, Nx_g, nite_, nout_, buf_,"task2_");
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

void solveExplicit(double *K, double *M, double *F, double *U, double lx_e,
	double qx_, double qy_, int Nvar_, int Nx_g, int nite_, int nout_,
	int buf_, std::string test)
{	
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
	for (int i = 0; i <= nite_; ++i)
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

void solveImplicit(double *K, double *M, double *F, double *U, double lx_e,
    double qx_, double qy_, double dt_, int Nvar_, int Nx_g, int nite_, int nout_,
    std::string test)
{	// ===================== Build Tables =====================//
	double *K_eff 		= new double[Nvar_ * Nvar_]();
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

	for (int i = 0; i < 2; ++i)
	// for (int i = 0; i < nite_; ++i)
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