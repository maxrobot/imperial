#include <iostream>
#include <mpi.h>
#include <time.h>

#include "Common.hpp"
#include "CommonMPI.hpp"
#include "InitFile.hpp"
#include "BuildFunction.hpp"
#include "Solvers.hpp"

using namespace std;

int main(int argc, char *argv[])
{	// =========== Declaration of all things MPI ==============//
	int retval = MPI_Init(&argc, &argv);


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
	double lx_e = lx_g/Nx_g;		// Local element length
	if (MPI::mpi_rank==0)
	{
		cout << Nvar_ << "  " << Nvar_e << endl;
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
	else
	{
		exit(EXIT_FAILURE);
	}
*/	
	MPI_Finalize();
	return 0;
}