#include <iostream>

#include "InitFile.hpp"
#include "Common.hpp"
#include "BuildFunction.hpp"
#include "GlobalVars.hpp"
#include "Memory.hpp"
#include "Solvers.hpp"

#include <cblas.h>

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

	const double Al_(1./12);  	// Constant Alpha
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

		// =================== Solve System =======================//
	    const int nrhs = 1;
	    int info = 0;
	    int*    ipiv = new int[Nvar_];
	    int kl = 4;
	    int ku = 4;
	    int ldab = 1 + 2*kl + ku;
	    int ldb = Nvar_;
		F77NAME(dgbsv)(Nvar_, kl, ku, nrhs, K_g, ldab, ipiv, F_g, ldb, info);
		writeVec(F_g, Nx_g, "output");
	}
	else if (eq_=="dynamic")
	{	// Solve Mu^{n+1} =  
		// ================ Initialise Local Vars. ================//
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
		double *F_e	 		= allocateDbl(6);

		// ============== Create Elemental K Matrix ===============//
		buildEye(eye, Nvar_);
		buildMele(M_e, A_, rho_, lx_e, Al_, dt_);
		buildKele(K_e, lx_e, A_, E_, I_);
		buildFele(F_e, lx_e, qx_, qy_);

		buildKglb(K_g, K_e, Nvar_, Nx_g);
		buildFglb(F_g, F_e, Nx_g);
		buildKglb(M_g, M_e, Nvar_, Nx_g);

		showMat(M_e, 6);	

		/*double *L_c		= allocateDbl(Nvar_);
		double *M_gt	= allocateDbl(Nvar_*Nvar_);
		double *F_gt	= allocateDbl(Nvar_);

	    int info = 0;
	    int *ipiv = new int[Nvar_];

		// =================== Create S Matrix ====================//
		// Start marching through time...
		// for (int i = 0; i < 1; ++i)
		for (int i = 0; i < nite_; ++i)
		{	// Make C: output = L_c
			// showMat	(M_g, Nvar_);
			F77NAME(dgemv)('n', Nvar_, Nvar_, 1, M_g, Nvar_, Un_g, 1, 0, L_c, 1);
			// showMat(M_g, Nvar_);

			// Make B: output = M_gt
			F77NAME(dcopy)(Nvar_*Nvar_, M_g, 1, M_gt, 1);
			// showMat(M_gt, Nvar_);
			F77NAME(dgemm)('n', 'n', Nvar_, Nvar_, Nvar_, dt_*dt_, eye, Nvar_, K_g,
					Nvar_, -2, M_gt, Nvar_);
			// showMat(M_gt, Nvar_);

			// Combine B and C: output = L_c
			// showMat(M_gt, Nvar_);
			// showVec(L_c, Nvar_);
			F77NAME(dgemv)('n', Nvar_, Nvar_, 1, M_gt, Nvar_, U_g, 1, -1, L_c, 1);
			
			// Make L_a: output = L_c
			F77NAME(dcopy)(Nvar_, F_g, 1, F_gt, 1);
			// showVec(F_gt, Nvar_);
			// showVec(L_c, Nvar_);
			F77NAME(dgemv)('n', Nvar_, Nvar_, dt_*dt_, eye, Nvar_, F_gt, 1, -1, L_c, 1);
			// showVec(L_c, Nvar_);

			// Calculate updated U_g
			// showVec(U_g	, Nvar_);
			F77NAME(dgesv)(Nvar_, 1, M_g, Nvar_, ipiv, L_c, Nvar_, info);
			if (i%5==0)
			{	showVec(L_c, Nvar_);
				}

			// Now update vars for repeat
			F77NAME(dcopy)(Nvar_, U_g, 1, Un_g, 1);

		}*/
	}
	else
	{
		exit(EXIT_FAILURE);
	}

	return 0;
}