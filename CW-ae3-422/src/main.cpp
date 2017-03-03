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
	// string scheme_("explicit");	// Type of integration, explicit or 
	string scheme_("none");	// Type of integration, explicit or 
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
		&rho_, &b_, &h_, &qx_, &qy_, &eq_, &scheme_);
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

		// ===================== Build Tables =====================//
		// Matrices
		double *K_g 		= allocateDbl(Nvar_*Nvar_);
		double *K_eff 		= allocateDbl(Nvar_*Nvar_);
		double *M_g			= allocateDbl(Nvar_*Nvar_);

		// Vectors
		double *F_g			= allocateDbl(Nvar_);
		double *U_g			= allocateDbl(Nvar_);
		double *Ud_g		= allocateDbl(Nvar_);
		double *Udd_g		= allocateDbl(Nvar_);
		double *tmp			= allocateDbl(Nvar_);
		double *tmp2		= allocateDbl(Nvar_);
		double *S_g			= allocateDbl(Nvar_);

		double *M_e			= allocateDbl(6*6);
		double *K_e			= allocateDbl(6*6);

		// Constanst and coefficients
		const double beta_(.25);
		const double gamma_(.5);
		const double co1_(1/(beta_*(dt_*dt_)));
		const double co2_(1/(beta_*dt_));
		const double co3_((1/(2*beta_))-1);
		const double co4_(dt_*(1-gamma_));
		const double co5_(dt_*gamma_);

		// ============== Create Elemental K Matrix ===============//
		buildMele(M_e, A_, rho_, lx_e, Al_);
		buildKele(K_e, lx_e, A_, E_, I_);

		buildKglb(K_g, K_e, Nvar_, Nx_g);
		buildKglb(M_g, M_e, Nvar_, Nx_g);

		// Build Keff
		for (int i = 0; i < Nvar_*Nvar_; ++i)
		{	double sum = (co1_*M_g[i]) + K_g[i];
			K_eff[i] = sum;
			K_g[i] = sum;
		}

	    int info = 0;
	    int *ipiv = new int[Nvar_];

		F77NAME(daxpy)(Nvar_*Nvar_, co1_, M_g, 1, K_g, 1);

		for (int i = 0; i < nite_; ++i)
		{	// Create Dynamic Force
			assignArr(F_g, 0., Nvar_);
			assignArr(S_g, 0., Nvar_);
			assignArr(tmp, 0., Nvar_);
			assignArr(tmp2, 0., Nvar_);
			updateVars(F_g, lx_e, qx_, qy_, Nx_g, i, nite_);


			// Sum U, Ud and, Udd with coefficients
			for (int i = 0; i < Nvar_; ++i)
			{	double sum = (co1_*U_g[i]) + (co2_*Ud_g[i]) + (co3_*Udd_g[i]);
				tmp2[i] = sum;
			}

			// Multiple mass by sum
			F77NAME(dgemv)('n', Nvar_, Nvar_, 1, M_g, Nvar_, tmp2, 1, 0, S_g, 1);

			// Add forces to sum...
			for (int i = 0; i < Nvar_; ++i)
			{	double sum = S_g[i] + F_g[i];
				S_g[i] = sum;
			}

			// Solve Keff*U_{n+1} = S
			F77NAME(dgesv)(Nvar_, 1, K_g, Nvar_, ipiv, S_g, Nvar_, info);

			// Update K_g to contain only the K_eff as desgv overwrites...
			for (int i = 0; i < Nvar_*Nvar_; ++i)
			{	K_g[i] = K_eff[i];
			}

			// Now update Udd
			F77NAME(dcopy)(Nvar_, Udd_g, 1, tmp, 1);   
			for (int i = 0; i < Nvar_; ++i)
			{	double sum = (co1_*(S_g[i]-U_g[i])) - (co2_*Ud_g[i]) -
					(co3_*Udd_g[i]);
				Udd_g[i] = sum;
			}

			// Now update Ud
			for (int i = 0; i < Nvar_; ++i)
			{	double sum = U_g[i] + (co4_*tmp[i]) + (co5_*Udd_g[i]);
				Ud_g[i] = sum;
			}

			// Now update U
			F77NAME(dcopy)(Nvar_, S_g, 1, U_g, 1);
			if (i%nout_==0)
			{	writeVec(U_g, Nx_g, i, "output");
			}
		}

	}
	if (eq_=="dynamic" && scheme_=="none")
	{	cout << "Please Choose Integration Scheme. (explicit/implicit)" << endl;
		exit(EXIT_FAILURE);
	}
	else
	{
		exit(EXIT_FAILURE);
	}

	return 0;
}