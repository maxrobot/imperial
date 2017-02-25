#include <iostream>
#include <stdio.h>
#include <string.h>

#include "InitFile.hpp"
#include "Common.hpp"
#include "BuildFunction.hpp"
#include "GlobalVars.hpp"
#include "Memory.hpp"

using namespace std;

int main(int argc, char *argv[])
{	// ================ Reading of Inputs =====================//
	ifstream param_file;
	param_file.open(argv[1], ifstream::in);
	if (!param_file)
	{	cout << "Unable to open input file" << endl;
		exit (EXIT_FAILURE);
	}   

	// Global Variables #############################################
	int nite_(0);         		// Number of time steps
	int Nx_g(0);          		// Number of global elements
	int Nvar_(0);				// Number of variables in domain
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

	const double Al_(1./12);  	// Constant Alpha
	// End of Global Variables ######################################

	readParamFile(param_file, nite_, Nx_g, lx_g, E_, rho_, b_, h_, qx_, qy_);
	// ===================== Build Tables =====================//
	double *K_g	 		= allocateDbl(Nvar_*Nvar_);
	double *F_g			= allocateDbl(Nvar_);
	double *U_g			= allocateDbl(Nvar_);

	// ================ Initialise Local Vars. ================//
	double K_e[6*6] = {};
	double F_e[6] = {};
	double lx_e = lx_g/Nx_g;		// Local element length

	initVars(b_, h_, A_, I_, Nvar_, Nx_g);

	// ============== Create Elemental K Matrix ===============//
	buildKele(K_e, lx_e, A_, E_, I_);
	buildKglb(K_g, K_e, Nvar_, Nx_g);
	buildFele(F_e, lx_e, qx_, qy_);
	buildFglb(F_g, F_e, Nx_g);

	showMat(K_g, Nvar_);
	showVec(K_e, Nvar_);

	return 0;
}

	// for (int i = 0; i < Nx_g; ++i)
	// {	int pnt = i*3;
	// 	cout << K_g[pnt] << " " << K_g[pnt+1] << " " << K_g[pnt+2] << " " << endl;
	// 	cout << pnt << " " << pnt+1 << " " << pnt+2 << " " << endl;
	// 	// myfile << K_g[pnt] << " " << K_g[pnt+1] << " " << K_g[pnt+2] << endl;

	// }

    // const int nrhs = 1;
    // int info = 0;
    // int*    ipiv = new int[Nvar_];

	// showVec(F_g, Nvar_);
	// F77NAME(dgesv)(Nvar_, nrhs, K_g, Nvar_, ipiv, F_g, Nvar_, info);
	// cout << endl;
	// showVec(F_g, Nvar_);
