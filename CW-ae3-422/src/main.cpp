#include <iostream>

#include "InitFile.hpp"
#include "Common.hpp"
#include "BuildFunction.hpp"
#include "GlobalVars.hpp"
#include "Solvers.hpp"
#include "Memory.hpp"

#include <cblas.h>

using namespace std;

void solve(double *K, double *F, int Nvar_, int ldab, int Nx_)
{	const int nrhs = 1;
    int info = 0;
    int *ipiv = new int[Nvar_];
    int kl = 4;
    int ku = 4;
    int ldb = Nvar_;

    // Use blas to solve system...
    F77NAME(dgbsv)(Nvar_, kl, ku, nrhs, K, ldab, ipiv, F, ldb, info);
	writeVec(F, Nx_, "output");

	delete[] ipiv;
}

int main(int argc, char *argv[])
{	// ================ Reading of Inputs =====================//
	ifstream param_file;
	param_file.open(argv[1], ifstream::in);
	if (!param_file)
	{	cout << "Unable to open input file" << endl;
		exit (EXIT_FAILURE);
	}   

	// Global Variables #############################################
	string eq_("none");			// Type of eq., static or dynamic

	int T_(0);	         		// Simulation length (s)
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

	const double Al_(1./24);  	// Constant Alpha
	const double buf(4);	  	// Buffer
	// End of Global Variables ######################################
	readParamFile(param_file, &T_, &nite_, &Nx_g, &lx_g, &E_,
		&rho_, &b_, &h_, &qx_, &qy_, &eq_);

	if (eq_=="static")
	{	
		// ================ Initialise Local Vars. ================//
		initVars(&b_, &h_, &A_, &I_, &E_, &Nvar_, &Nx_g);
		double lx_e = lx_g/Nx_g;		// Local element length
		
		// ===================== Build Tables =====================//
		// Matrices
		double *K_g 		= allocateDbl(Nvar_*(9+buf));

		// Vectors
		double *F_g			= allocateDbl(Nvar_);
		double *U_g			= allocateDbl(Nvar_);

		double *K_e			= allocateDbl(6*6);
		double *F_e	 		= allocateDbl(6);


		// ============== Create Elemental K Matrix ===============//
		buildKele(K_e, lx_e, A_, E_, I_);
		buildKglbSparse(K_g, K_e, Nvar_, Nx_g, buf);
		buildFele(F_e, lx_e, qx_, qy_);
		buildFglb(F_g, F_e, Nx_g);
		// showVec(F_g, Nvar_);
		// showMat(K_g,9+buf,Nvar_);


		// =================== Solve System =======================//
		solve(K_g, F_g, Nvar_, 9+buf, Nx_g);
	    // const int nrhs = 1;
	    // int info = 0;
	    // int*    ipiv = new int[Nvar_];
	    // int kl = 4;
	    // int ku = 4;
	    // int ldab = 1 + 2*kl + ku;
	    // int ldb = Nvar_;
		// F77NAME(dgbsv)(Nvar_, kl, ku, nrhs, K_g, ldab, ipiv, F_g, ldb, info);
		// writeVec(F_g, Nx_g, "output");
	}
	else if (eq_=="dynamic")
	{	// ================ Initialise Local Vars. ================//
		initVars(&b_, &h_, &A_, &I_, &E_, &dt_, &Nvar_, &Nx_g, &T_, &nite_);
		double lx_e = lx_g/Nx_g;		// Local element length
		// ===================== Build Tables =====================//
		// Matrices
		double *K_g 		= allocateDbl(Nvar_*Nvar_);
		double *M_g			= allocateDbl(Nvar_*Nvar_);
		double *eye			= allocateDbl(Nvar_*Nvar_);

		// Vectors
		double *F_g			= allocateDbl(Nvar_);
		double *U_g			= allocateDbl(Nvar_);
		double *Un_g		= allocateDbl(Nvar_);
		double *S_g			= allocateDbl(Nvar_);

		double *M_e			= allocateDbl(6*6);
		double *K_e			= allocateDbl(6*6);

		// ============== Create Elemental K Matrix ===============//
		buildEye(eye, Nvar_);
		buildMele(M_e, A_, rho_, lx_e, Al_, dt_);
		buildKele(K_e, lx_e, A_, E_, I_);

		buildKglb(K_g, K_e, Nvar_, Nx_g);
		buildKglb(M_g, M_e, Nvar_, Nx_g);

		double *MK_o	= allocateDbl(Nvar_*Nvar_);
		double *MKU_o		= allocateDbl(Nvar_);
		double *MU_o	= allocateDbl(Nvar_);

	    int info = 0;
	    int *ipiv = new int[Nvar_];

	    for (int i = 0; i < Nvar_*Nvar_; ++i)
		{	MK_o[i] = 2*M_g[i] - K_g[i];
		}
		
		// =================== Create S Matrix ====================//
		// Start marching through time...
		for (int i = 0; i <= nite_; ++i)
		{	// Make C: output = L_c
			assignArr(F_g, 0., Nvar_);
			updateVars(F_g, lx_e, qx_, qy_, Nx_g, i, nite_);

			// Calculate M*U{n-1}
			F77NAME(dgemv)('n', Nvar_, Nvar_, 1, M_g, Nvar_, Un_g, 1, 0, MU_o, 1);

			// Calculate MK_o*U{n}
			F77NAME(dgemv)('n', Nvar_, Nvar_, 1, MK_o, Nvar_, U_g, 1, 0, MKU_o, 1);

			for (int i = 0; i < Nvar_; ++i)
			{	S_g[i] = F_g[i] + MKU_o[i] - MU_o[i];
			}

			// Calculate updated M*U{n+1} = S
			F77NAME(dgesv)(Nvar_, 1, M_g, Nvar_, ipiv, S_g, Nvar_, info);

			// Now update vars for repeat
			F77NAME(dcopy)(Nvar_, U_g, 1, Un_g, 1); // Save Un-1
			F77NAME(dcopy)(Nvar_, S_g, 1, U_g, 1); // Update Un
			if (i%10==0)
			{	writeVec(U_g, Nx_g, i, "output");
			}
		}
	}
	else
	{
		exit(EXIT_FAILURE);
	}

	return 0;
}