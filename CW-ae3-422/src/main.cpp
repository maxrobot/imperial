#include <iostream>

#include "InitFile.hpp"
#include "Common.hpp"
#include "BuildFunction.hpp"
#include "GlobalVars.hpp"

using namespace std;

// Global Variables #############################################
int nite_(0);         		// Number of time steps
int Nx_g(0);          		// Number of global elements
double lx_g(0);       		// Length of global domain
double dt_(0);        		// Time step
double dx_(0);        		// Mesh size

double rho_(0);       		// Density
double I_(0);         		// Second moment area
double E_(0);         		// Youngs modulus
double b_(0);         		// Cross-sectional width
double h_(0);         		// Cross-sectional height
double A_(0);         		// Cross-sectional area

const double Al_(1./12);  	// Constant Alpha
// End of Global Variables ######################################

int main(int argc, char *argv[])
{	// ================ Reading of Inputs =====================//
	ifstream param_file;
	param_file.open(argv[1], ifstream::in);
	if (!param_file)
	{	cout << "Unable to open input file" << endl;
		exit (EXIT_FAILURE);
	}   

	readParamFile(param_file);
	initVars();
	// ================ Initialise Local Vars. ===============//
	double K_e[6*6] = {};
	double F_e[6] = {};
	double lx_e = lx_g/Nx_g;		// Local element length
	cout << A_ << " " << E_ <<  " " << Nx_g << " " << lx_g << endl;

	// ============= Create Elemental K Matrix ===============//
	showMatrix(K_e, 6);
	buildKele(K_e, lx_e);
	cout << endl;
	showMatrix(K_e, 6);
	// ============= Create Elemental F Matrix ===============//
	// buildFele(F_e, lx_e)


	return 0;
}