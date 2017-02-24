#include <iostream>

#include "InitFile.hpp"
#include "Common.hpp"
#include "BuildFunction.hpp"
#include "GlobalVars.hpp"
#include "Memory.hpp"

using namespace std;

// Define LAPACK Shit ###########################################
#define F77NAME(x) x##_
extern "C" {
    double F77NAME(dnrm2)(const int& N, const double *X, const int& incX);
    void F77NAME(daxpy)(const int& N, const double& alpha, const double *X,
                             const int& incX, double *Y, const int& incY);
    void F77NAME(dgemv)(const char& trans, const int& m, const int& n,
                    const double& alpha, const double* a, const int& lda,
                    const double* x, const int& incx, const double& beta,
                    double* y, const int& incy);
    void F77NAME(dgesv)(const int& n, const int& nrhs, const double * A,
                    const int& lda, int * ipiv, double * B, const int& ldb,
                    int& info);
}

// Global Variables #############################################
int nite_(0);         		// Number of time steps
int Nx_g(0);          		// Number of global elements
int Nvars_(0);				// Number of variables in domain
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
	// ===================== Build Tables =====================//
	double *K_g			= allocateDbl(Nvars_*Nvars_);
	double *F_g			= allocateDbl(Nvars_);
	double *U_g			= allocateDbl(Nvars_);

	// ================ Initialise Local Vars. ================//
	double K_e[6*6] = {};
	double F_e[6] = {};
	double lx_e = lx_g/Nx_g;		// Local element length

	// ============== Create Elemental K Matrix ===============//
	buildKele(K_e, lx_e);
	buildFele(F_e, lx_e);
	buildKglb(K_g, K_e);
	buildFglb(F_g, F_e);
	// double test[9] = {3,1,1,2};
	// double test2[3] = {9,8};
	// double test3[3] = {0,0};

	return 0;
}
