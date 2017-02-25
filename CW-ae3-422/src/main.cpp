#include <iostream>
#include <stdio.h>
#include <string.h>

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

int main(int argc, char *argv[])
{	// ================ Reading of Inputs =====================//
	ifstream param_file;
	param_file.open(argv[1], ifstream::in);
	if (!param_file)
	{	cout << "Unable to open input file" << endl;
		exit (EXIT_FAILURE);
	}   

	readParamFile(param_file);
	// ===================== Build Tables =====================//
	double *K_g	 		= allocateDbl(Nvar_*Nvar_);
	double *F_g			= allocateDbl(Nvar_);
	double *U_g			= allocateDbl(Nvar_);

	// ================ Initialise Local Vars. ================//
	double K_e[6*6] = {};
	double F_e[6] = {};
	double lx_e = lx_g/Nx_g;		// Local element length

	initVars();
	zeroMat(K_g, Nvar_*Nvar_);
	zeroMat(F_g, Nvar_);
	zeroMat(U_g, Nvar_);

	// ============== Create Elemental K Matrix ===============//
	buildKele(K_e, lx_e);
	buildFele(F_e, lx_e);
	buildKglb(K_g, K_e);
	buildFglb(F_g, F_e);


 //    const int nrhs = 1;
 //    int info = 0;
 //    int*    ipiv = new int[Nvar_];





	// showVec(F_g, Nvar_);
	// F77NAME(dgesv)(Nvar_, nrhs, K_g, Nvar_, ipiv, F_g, Nvar_, info);
	// cout << endl;
	// showVec(F_g, Nvar_);

	return 0;
}
