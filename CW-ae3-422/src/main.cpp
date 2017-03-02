#include <iostream>

#include "InitFile.hpp"
#include "Common.hpp"
#include "BuildFunction.hpp"
#include "GlobalVars.hpp"
#include "Solvers.hpp"
#include "Memory.hpp"

using namespace std;

void solveDynamic(double *K, double *M, double *F, double *U, double lx_e,
	double qx_, double qy_, int Nvar_, int Nx_g, int nite_, int nout_, std::string test)
{	
	double *eye			= allocateDbl(Nvar_*Nvar_);
	double *Un_g		= allocateDbl(Nvar_);
	double *S			= allocateDbl(Nvar_);
	double *MK_o	= allocateDbl(Nvar_*Nvar_);
	double *MKU_o	= allocateDbl(Nvar_);
	double *MU_o	= allocateDbl(Nvar_);

	buildEye(eye, Nvar_);

    int info = 0;
    int *ipiv = new int[Nvar_];

    for (int i = 0; i < Nvar_*Nvar_; ++i)
	{	MK_o[i] = 2*M[i] - K[i];
	}

	// =================== Create S Matrix ====================//
	// Start marching through time...
	for (int i = 0; i <= nite_; ++i)
	{	// Make C: output = L_c
		assignArr(F, 0., Nvar_);
		updateVars(F, lx_e, qx_, qy_, Nx_g, i, nite_);

		// Calculate M*U{n-1}
		F77NAME(dgemv)('n', Nvar_, Nvar_, 1, M, Nvar_, Un_g, 1, 0, MU_o, 1);

		// Calculate MK_o*U{n}
		F77NAME(dgemv)('n', Nvar_, Nvar_, 1, MK_o, Nvar_, U, 1, 0, MKU_o, 1);

		for (int i = 0; i < Nvar_; ++i)
		{	S[i] = F[i] + MKU_o[i] - MU_o[i];
		}

		// Calculate updated M*U{n+1} = S
		F77NAME(dgesv)(Nvar_, 1, M, Nvar_, ipiv, S, Nvar_, info);

		// Now update vars for repeat
		F77NAME(dcopy)(Nvar_, U, 1, Un_g, 1); // Save Un-1
		F77NAME(dcopy)(Nvar_, S, 1, U, 1); // Update Un
		if ((i%nout_)==0)
			writeVec(U, Nx_g, i, test);
		// showVec(U, Nvar_);
		// }
	}
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
	int nout_(0);         		// Number of output interval
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
	readParamFile(param_file, &T_, &nite_, &Nx_g, &nout_, &lx_g, &E_,
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

		// ================== Clean up Space ======================//
		delete[] K_e;
		delete[] F_e;

		// =================== Solve System =======================//
		solveStatic(K_g, F_g, Nvar_, 9+buf, Nx_g, "output");
	}
	else if (eq_=="dynamic")
	{	// ================ Initialise Local Vars. ================//
		initVars(&b_, &h_, &A_, &I_, &E_, &dt_, &Nvar_, &Nx_g, &T_,
			&nite_);
		double lx_e = lx_g/Nx_g;		// Local element length
		// ===================== Build Tables =====================//
		// Matrices
		double *K_g 		= allocateDbl(Nvar_*Nvar_);
		double *M_g			= allocateDbl(Nvar_*Nvar_);

		// Vectors
		double *F_g			= allocateDbl(Nvar_);
		double *U_g			= allocateDbl(Nvar_);


		double *M_e			= allocateDbl(6*6);
		double *K_e			= allocateDbl(6*6);

		// ============== Create Elemental K Matrix ===============//
		buildMele(M_e, A_, rho_, lx_e, Al_, dt_);
		buildKele(K_e, lx_e, A_, E_, I_);

		buildKglb(K_g, K_e, Nvar_, Nx_g);
		buildKglb(M_g, M_e, Nvar_, Nx_g);

		// ================== Clean up Space ======================//
		delete[] M_e;
		delete[] K_e;
		// double *MK_o	= allocateDbl(Nvar_*Nvar_);
		// double *MKU_o	= allocateDbl(Nvar_);
		// double *MU_o	= allocateDbl(Nvar_);

	 //    int info = 0;
	 //    int *ipiv = new int[Nvar_];

	 //    for (int i = 0; i < Nvar_*Nvar_; ++i)
		// {	MK_o[i] = 2*M_g[i] - K_g[i];
		// }
		
		// // =================== Create S Matrix ====================//
		// // Start marching through time...
		// for (int i = 0; i <= nite_; ++i)
		// {	// Make C: output = L_c
		// 	assignArr(F_g, 0., Nvar_);
		// 	updateVars(F_g, lx_e, qx_, qy_, Nx_g, i, nite_);

		// 	// Calculate M*U{n-1}
		// 	F77NAME(dgemv)('n', Nvar_, Nvar_, 1, M_g, Nvar_, Un_g, 1, 0, MU_o, 1);

		// 	// Calculate MK_o*U{n}
		// 	F77NAME(dgemv)('n', Nvar_, Nvar_, 1, MK_o, Nvar_, U_g, 1, 0, MKU_o, 1);

		// 	for (int i = 0; i < Nvar_; ++i)
		// 	{	S_g[i] = F_g[i] + MKU_o[i] - MU_o[i];
		// 	}

		// 	// Calculate updated M*U{n+1} = S
		// 	F77NAME(dgesv)(Nvar_, 1, M_g, Nvar_, ipiv, S_g, Nvar_, info);

		// 	// Now update vars for repeat
		// 	F77NAME(dcopy)(Nvar_, U_g, 1, Un_g, 1); // Save Un-1
		// 	F77NAME(dcopy)(Nvar_, S_g, 1, U_g, 1); // Update Un
		// 	if (i%nout_==0)
		// 	{	writeVec(U_g, Nx_g, i, "output");
		// 	}
		// }
		solveDynamic(K_g, M_g, F_g, U_g, lx_e, qx_, qy_, Nvar_, Nx_g, nite_, nout_, "output");
	}
	else
	{
		exit(EXIT_FAILURE);
	}

	return 0;
}