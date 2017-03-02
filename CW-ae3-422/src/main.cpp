#include <iostream>

#include "Memory.hpp"
#include "Common.hpp"
#include "InitFile.hpp"
#include "BuildFunction.hpp"
#include "Solvers.hpp"

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
	string scheme_("implicit");	// Type of integration, explicit or 
								// implicit

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

	// End of Global Variables ######################################
	readParamFile(param_file, &T_, &nite_, &Nx_g, &nout_, &lx_g, &E_,
		&rho_, &b_, &h_, &qx_, &qy_, &eq_);
	if (eq_=="static")
	{	
		// ================ Initialise Local Vars. ================//
		initVars(&b_, &h_, &A_, &I_, &E_, &dt_, &Nvar_, &Nx_g, &T_,
			&nite_);
		double lx_e = lx_g/Nx_g;		// Local element length
		const double buf(4);	  	// Buffer
		
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
	else if (eq_=="dynamic" && scheme_=="explicit")
	{	// ================ Initialise Local Vars. ================//
		initVars(&b_, &h_, &A_, &I_, &E_, &dt_, &Nvar_, &Nx_g, &T_,
			&nite_);
		double lx_e = lx_g/Nx_g;
		const double Al_(1./24);  	// Constant Alpha

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
		solveDynamic(K_g, M_g, F_g, U_g, lx_e, qx_, qy_, Nvar_,
			Nx_g, nite_, nout_, "output");
	}
	else if (eq_=="dynamic" && scheme_=="implicit")
	{	// ================ Initialise Local Vars. ================//
		initVars(&b_, &h_, &A_, &I_, &E_, &dt_, &Nvar_, &Nx_g, &T_,
			&nite_);
		double lx_e = lx_g/Nx_g;		// Local element length
		const double Al_(1./24);  	// Constant Alpha
		const double Bt_(.25);  	// Constant Beta
		const double Gm_(.5);  	// Constant Gamma

		const double coeff1_(1/(Bt_*dt_*dt_));	// Coefficient 1
		const double coeff2_(1/(Bt_*dt_));		// Coefficient 2
		const double coeff3_((1/(2*Bt_)) -1);	// Coefficient 3
		const double coeff4_(dt_*(1-Gm_));		// Coefficient 4
		const double coeff5_(dt_*Gm_);			// Coefficient 5

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

		double *K_eff 		= allocateDbl(Nvar_*Nvar_);
		double *M_gt		= allocateDbl(Nvar_*Nvar_);
		double *Ud_g		= allocateDbl(Nvar_);
		double *Udd_g		= allocateDbl(Nvar_);
		double *Udd_gt		= allocateDbl(Nvar_);
		double *S_g			= allocateDbl(Nvar_);
		double *S_gt		= allocateDbl(Nvar_);

		int info = 0;
	    int *ipiv = new int[Nvar_];
	    const int coeff_ = 1/(Bt_*dt_*dt_);
		for (int i = 0; i < Nvar_*Nvar_; ++i)
		{	M_gt[i] = coeff_*M_g[i];
		}
		// Build K_eff = M/(Bt*Dt_^2) + K
		for (int i = 0; i < Nvar_*Nvar_; ++i)
		{	K_eff[i] = M_gt[i] + K_g[i];
		}


		// for (int i = 0; i < nite_+1; ++i)
		for (int i = 0; i < 100; ++i)
		{
			assignArr(F_g, 0., Nvar_);
			updateVars(F_g, lx_e, qx_, qy_, Nx_g, i+1, nite_);

			for (int i = 0; i < Nvar_; ++i)
			{	S_g[i] = (coeff1_*U_g[i]) + (coeff2_*Ud_g[i]) + (coeff3_*Udd_g[i]);
			}
			// Now multiply M*S = S
			F77NAME(dgemv)('n', Nvar_, Nvar_, 1, M_g, Nvar_, S_g, 1, 0, S_g, 1);
			// showVec(S_g, Nvar_);

			for (int i = 0; i < Nvar_; ++i)
			{	S_gt[i] = S_g[i] + F_g[i];
			}
			F77NAME(dgesv)(Nvar_, 1, K_eff, Nvar_, ipiv, S_gt, Nvar_, info);
			// Update Values Udd_g and Ud_g
			// showVec(S_gt, Nvar_);
			// showVec(Udd_g, Nvar_);
			F77NAME(dcopy)(Nvar_, Udd_g, 1, Udd_gt, 1); // Save Un-1
			// showVec(Udd_gt, Nvar_);
			assignArr(Udd_g, 0., Nvar_);
			for (int i = 0; i < Nvar_; ++i)
			{	Udd_g[i] = coeff1_*(S_gt[i] - U_g[i]) -
					(coeff2_*Ud_g[i]) - (coeff3_*Udd_g[i]);
			}
			assignArr(Ud_g, 0., Nvar_);
			for (int i = 0; i < Nvar_; ++i)
			{	Ud_g[i] = S_g[i] + (coeff4_*Udd_gt[i]) + (coeff5_*Udd_g[i]);
			}
			// Create sum matrix S
			// cout << 1 << endl;
			// showVec(U_g, Nvar_);
			// cout << 2 << endl;
			// showVec(Ud_g, Nvar_);
			// cout << 3 << endl;
			// showVec(Udd_g, Nvar_);
			F77NAME(dcopy)(Nvar_, S_gt, 1, U_g, 1); // Save Un-1
			if ((i%nout_)==0)
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