#include <iostream>
#include <mpi.h>
#include <time.h>

#include "Common.hpp"
#include "CommonMPI.hpp"
#include "InitFile.hpp"
#include "BuildFunction.hpp"
#include "Solvers.hpp"

using namespace std;
// using namespace MPI;

int main(int argc, char *argv[])
{	// ======== Declaration of all things good and MPI ======== //
	MPI_Init(&argc, &argv);

	MPI::initMpiStuff();
	// ================ Reading of Inputs =====================//

	ifstream param_file;
	param_file.open(argv[1], ifstream::in);
	if (!param_file)
	{	printMessage("Unable to open input file");
		exit (EXIT_FAILURE);
	}   

	// Global Variables #############################################
	string eq_("none");			// Type of eq., static or dynamic
	string scheme_("none");		// Type of integration, explicit or 
	string sparse_("none");		// Type of integration, explicit or 

	int T_(0);	         		// Simulation length (s)
	int nite_(0);         		// Number of time steps
	int nout_(0);         		// Number of output interval
	int Nx_g(0);          		// Number of global elements
	int Nx_(0);          		// Number of local elements
	int Nvar_(0);				// Number of variables global in domain
	int Nvar_e(0);				// Number of variables local in domain
	double lx_g(0);       		// Length of global domain
	double lx_e(0);       		// Length of global domain
	double dt_(0);        		// Time step
	double dx_(0);        		// Mesh size

	double rho_(0);       		// Density
	double I_(0);         		// Second moment area
	double E_(0);         		// Youngs modulus
	double qx_(0);			    // Axial uniform load
	double qy_(0);        		// Traverse uniform load
	double b_(0);         		// Cross-sectional width
	double h_(0);         		// Cross-sectional height
	double A_(0);         		// Cross-sectional area

	// End of Global Variables ######################################
	readParamFile(param_file, &T_, &nite_, &Nx_g, &nout_, &lx_g, &E_,
		&rho_, &b_, &h_, &qx_, &qy_, &eq_, &scheme_, &sparse_);
	initVars(&b_, &h_, &A_, &I_, &E_, &dt_, &Nvar_, &Nvar_e, &Nx_g,
		&Nx_, &T_, &nite_);
	lx_e = lx_g/Nx_g;		// Local element length
	if (MPI::mpi_rank==0)
	{
		cout << Nx_g << "  " << Nx_ << "    " << Nvar_ << "  " << Nvar_e << endl;
	}
	// ===================== Build Tables =====================//
	double *F_g	= new double[Nvar_]();
	double *U_g	= new double[Nvar_]();
	double *K_e	= new double[6*6]();
	
	if (eq_=="static")
	{	runSolver(K_e, U_g, F_g, lx_e, A_, E_, I_, qx_, qy_, Nvar_, Nx_g);
	}
/*
	if (eq_=="dynamic")
	{	if (scheme_=="explicit")
		{	const int buf_(0);
			runSolver(K_e, U_g, F_g, dt_, lx_e, A_, E_, I_, rho_, qx_, qy_, Nvar_,
				Nx_g, nite_, nout_, buf_, sparse_);
		}
		if (scheme_=="implicit")
		{	const int buf_(4);
			runSolver(K_e, U_g, F_g, dt_, lx_e, A_, E_, I_, rho_, qx_, qy_, Nvar_,
				Nx_g, nite_, nout_, buf_, sparse_);
		}
		else
		{	printMessage("Please Choose Integration Scheme. (explicit/implicit)");
			exit(EXIT_FAILURE);
		}
	}
*/

	MPI_Finalize();
	return 0;
}

// =======
// 	MPI::initMpiDomain(&Nx_g, &Nx_);

// 	// ================ Initialise Local Vars. ================//
// 	initVars(&lx_g, &lx_e, &b_, &h_, &A_, &I_, &E_, &dt_, &Nvar_,
// 		&Nx_g, &Nx_, &T_,	&nite_);

// 	if (eq_=="static" && MPI::mpi_size==1)
// 	{	
// 		const double buf_(4);	  	// Buffer
		
// 		// ===================== Build Tables =====================//
// 		// Matrices
// 		double *K_g 		= allocateDbl(Nvar_*(9+buf_));

// 		// Vectors
// 		double *F_g			= allocateDbl(Nvar_);
// 		double *U_g			= allocateDbl(Nvar_);

// 		double *K_e			= allocateDbl(6*6);
// 		double *F_e	 		= allocateDbl(6);


// 		// ============== Create Elemental K Matrix ===============//
// 		buildKele(K_e, lx_e, A_, E_, I_);
// 		buildFele(F_e, lx_e, qx_, qy_);

// 		buildKglbSparse(K_g, K_e, Nvar_, Nx_g, buf_);
// 		buildFglb(F_g, F_e, Nx_g, Nvar_);

// 		// ================== Clean up Space ======================//
// 		delete[] K_e;
// 		delete[] F_e;

// 		// =================== Solve System =======================//
// 		solveStatic(K_g, F_g, Nvar_, 9+buf_, Nx_g, "task1");
// 	}
// 	else if (eq_=="dynamic" && scheme_=="explicit")
// 	{	const double Al_(1./24);  	// Constant Alpha
// 		const double buf_(0);	  	// Buffer

// 		// Vectors
// 		double *F_g			= allocateDbl(Nvar_);
// 		double *F_orig		= allocateDbl((Nx_g-1)*3);
// 		double *U_g			= allocateDbl(Nvar_);

// 		double *M_e			= allocateDbl(6*6);
// 		double *K_e			= allocateDbl(6*6);
// 		double *F_e			= allocateDbl(6);


// 		// ============== Create Elemental K Matrix ===============//
// 		buildMele(M_e, A_, rho_, lx_e, Al_, dt_);
// 		buildKele(K_e, lx_e, A_, E_, I_);
// 		buildFele(F_e, lx_e, qx_, qy_);


// 		if (sparse_=="none")
// 		{	// Matrices
// 			double *K_g 		= allocateDbl(Nvar_*Nvar_);
// 			double *K_ 			= allocateDbl(Nvar_*Nvar_);
// 			double *M_g			= allocateDbl(Nvar_*Nvar_);
		
// 		// ===================== Build Tables =====================//
// 			buildKglb(K_g, K_e, Nvar_, Nx_g);
// 			buildKglb(M_g, M_e, Nvar_, Nx_g);

// 			if (MPI::mpi_size==1)
// 			{	solveExplicit(K_g, M_g, F_g, U_g, lx_e, qx_, qy_, Nvar_,
// 					Nx_g, nite_, nout_, buf_,"task2_");
// 			}
// 			else if (MPI::mpi_size>1)
// 			{	printMessage("Error: Full matrices only possible in serial.");
// 			}
// 		}
// 		else if (sparse_=="sparse")
// 		{	// Matrices
// 			double *K_g 		= allocateDbl(Nvar_*(9+buf_));
// 			double *K_ 			= allocateDbl(Nvar_*Nvar_);
// 			double *M_g			= allocateDbl(Nvar_);

// 		// ===================== Build Tables =====================//
// 			if (MPI::mpi_size==1)
// 			{	// Build Matrices
// 				buildSparse(K_g, K_e, Nvar_, Nx_, buf_);
// 				buildMglbSparse(M_g, M_e, Nvar_, Nx_, buf_);

// 				solveSparseExplicit(K_g, M_g, F_g, U_g, lx_e, qx_, qy_,
// 					Nvar_, Nx_g, nite_, nout_, buf_,"task2_");

// 			}
// 			else if (MPI::mpi_size>1)
// 			{	// Build matrices
// 				buildBandSparse(K_g, K_e, Nvar_, Nx_, buf_);
// 				buildBandFglb(F_g, F_e, Nx_g, Nvar_);
// 				buildMglbSparse(M_g, M_e, Nvar_, Nx_, buf_);
// 			}

// 			// if (MPI::mpi_rank==0)
// 			// {	
// 			// 	buildFglb(F_orig, F_e, Nx_g, (Nx_g-1)*3);
// 			// 	showVec(F_orig, (Nx_g-1)*3);
// 			// }
// 			// for (int i = 0; i < MPI::mpi_size; ++i)
// 			// {	if (MPI::mpi_rank==i)
// 			// 	{	showMat(K_g, 9+buf_, Nvar_);
// 			// 	}
// 			// 	MPI_Barrier(MPI::mpi_comm);
// 			// }
// 			for (int i = 0; i < MPI::mpi_size; ++i)
// 			{	if (MPI::mpi_rank==i)
// 				{	showVec(F_g, Nvar_);
// 					showMat(K_g)
// 				}
// 				MPI_Barrier(MPI::mpi_comm);
// 			}
// 		}		
// 	}

// 	// else if (eq_=="dynamic" && scheme_=="implicit")
// 	// {	const double buf_(4);	  	// Buffer
// 	// 	const double Al_(1./24);  	// Constant Alpha

// 	// 	// Vectors
// 	// 	double *F_g			= allocateDbl(Nvar_);
// 	// 	double *U_g			= allocateDbl(Nvar_);

// 	// 	double *M_e			= allocateDbl(6*6);
// 	// 	double *K_e			= allocateDbl(6*6);
// 	// 	// ============== Create Elemental K Matrix ===============//
// 	// 	buildMele(M_e, A_, rho_, lx_e, Al_);
// 	// 	buildKele(K_e, lx_e, A_, E_, I_);


// 	// 	// ==================== Run Solver ========================//
		
// 	// 	if (sparse_=="none")
// 	// 	{	// Matrices
// 	// 		double *K_g 		= allocateDbl(Nvar_*Nvar_);
// 	// 		double *K_ 			= allocateDbl(Nvar_*Nvar_);
// 	// 		double *M_g			= allocateDbl(Nvar_*Nvar_);

// 	// 	// ===================== Build Tables =====================//
// 	// 		buildKglb(K_g, K_e, Nvar_, Nx_g);
// 	// 		buildKglb(M_g, M_e, Nvar_, Nx_g);

// 	// 	// ==================== Run Solver ========================//
// 	// 		solveImplicit(K_g, M_g, F_g, U_g, lx_e, qx_, qy_, dt_,
// 	// 			Nvar_, Nx_g, nite_, nout_, "task3_");
// 	// 	}
// 	// 	else if (sparse_=="sparse")
// 	// 	{	// Matrices
// 	// 		double *K_g 		= allocateDbl(Nvar_*(9+buf_));
// 	// 		double *K_ 			= allocateDbl(Nvar_*Nvar_);
// 	// 		double *M_g			= allocateDbl(Nvar_);

// 	// 	// ===================== Build Tables =====================//
// 	// 		buildSparse(K_g, K_e, Nvar_, Nx_g, buf_);
// 	// 		buildMglbSparse(M_g, M_e, Nvar_, Nx_g, buf_);

// 	// 	// ==================== Run Solver ========================//
// 	// 		solveSparseImplicit(K_g, M_g, F_g, U_g, lx_e, qx_, qy_,
// 	// 			dt_, Nvar_, Nx_g, nite_, nout_, buf_, "task3_");
// 	// 	}

// 	// }
// 	// if (eq_=="dynamic" && scheme_=="none")
// 	// {	printMessage("Please Choose Integration Scheme. (explicit/implicit)");
// 	// 	exit(EXIT_FAILURE);
// 	// }
// >>>>>>> master
// 	else
// 	{ 	if (MPI::mpi_rank==0)
// 		{	printMessage("Error: please check parameter file!");
// 		}
// 	}
// <<<<<<< HEAD	