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
{	// ======== Declaration of all things good and MPI ======== //
	MPI_Init(&argc, &argv);

	MPI::initMpiStuff();
	MPI::initMpiDomain();
	// MPI::initCblacsStuff();
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
	string sparse_("none");		// Type of matrices, sparse or 

	int T_(0);	         		// Simulation length (s)
	int nite_(0);         		// Number of time steps
	int nout_(0);         		// Number of output interval
	int Nx_g(0);          		// Number of global elements
	int Nx_(0);          		// Number of local elements
	int Nvar_(0);				// Number of variables global in domain
	int Nvar_e(0);				// Number of variables local in domain
	int Nghost_(0);				// Number of ghost variables local in domain
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
	initVars(&b_, &h_, &A_, &I_, &E_, &dt_, &Nvar_, &Nvar_e, &Nghost_,
		&Nx_g, &Nx_, &T_, &nite_);
	lx_e = lx_g/Nx_g;		// Local element length

	// ===================== Build Tables =====================//
	double *K_e	= new double[6*6]();

	if (eq_=="static")
	{	runSolver(K_e, lx_e, A_, E_, I_, qx_, qy_, Nvar_, Nx_g);
	}
	if (eq_=="dynamic")
	{	if (scheme_=="explicit")
		{	const int buf_(0);
			runSolver(K_e, dt_, lx_e, A_, E_, I_, rho_, qx_, qy_, Nvar_,
				Nvar_e, Nghost_, Nx_g, Nx_, nite_, nout_, buf_, sparse_);
		}
		if (scheme_=="implicit")
		{	if (MPI::mpi_size==1)
			{	const int buf_(4);
				runSolver(K_e, dt_, lx_e, A_, E_, I_, rho_, qx_, qy_, Nvar_,
					Nvar_e, Nghost_, Nx_g, Nx_, nite_, nout_, buf_, sparse_);
			}
			if (MPI::mpi_size>1)
			{	const int buf_(8);
				runSolver(K_e, dt_, lx_e, A_, E_, I_, rho_, qx_, qy_, Nvar_,
					Nvar_e, Nghost_, Nx_g, Nx_, nite_, nout_, buf_, sparse_);
			}
		}
		if (scheme_=="none")
		{	printMessage("Please Choose Integration Scheme. (explicit/implicit)");
		}
	}
	MPI_Finalize();
	return 0;
}