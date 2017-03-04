#include "Solvers.hpp"
#include "Common.hpp"
#include "BuildFunction.hpp"

using namespace std;

void solveStatic(double *K, double *F, int Nvar_, int ldab, int Nx_, std::string test)
{	const int nrhs = 1;
    int info = 0;
    int *ipiv = new int[Nvar_];
    int kl = 4;
    int ku = 4;
    int ldb = Nvar_;

	showVec(F, Nvar_);

    // Use blas to solve system...
    F77NAME(dgbsv)(Nvar_, kl, ku, nrhs, K, ldab, ipiv, F, ldb, info);
    writeVec(F, Nx_, test);
	delete[] ipiv;
}

void solveExplicit(double *K, double *M, double *F, double *U, double lx_e,
	double qx_, double qy_, int Nvar_, int Nx_g, int nite_, int nout_,
	int buf_, std::string test)
{	
	double *MK_o	= new double[Nvar_ * Nvar_]();
	double *Un_g	= new double[Nvar_]();
	double *S		= new double[Nvar_]();
	double *MKU_o	= new double[Nvar_]();
	double *MU_o	= new double[Nvar_]();

    int info = 0;
    int *ipiv = new int[Nvar_];

    for (int i = 0; i < Nvar_*Nvar_; ++i)
    {	MK_o[i] = K[i] - (2*M[i]);
    }    	

	// =================== Create S Matrix ====================//
	// Start marching through time...
	for (int i = 0; i <= nite_; ++i)
	{	// Calculate MK_o*U{n}
		F77NAME(dgemv)('n', Nvar_, Nvar_, 1, MK_o, Nvar_, U, 1, 0, MKU_o, 1);
		
		assignArr(F, 0., Nvar_);
		updateVars(F, lx_e, qx_, qy_, Nx_g, Nvar_, i, nite_);

		// Calculate M*U{n-1}
		for (int i = 0; i < Nvar_; ++i)
		{	int pnt = i*Nvar_ + i;
			MU_o[i] = M[pnt]*Un_g[i];
		}

		for (int i = 0; i < Nvar_; ++i)
		{	S[i] = F[i] - MKU_o[i] - MU_o[i];
		}

		// Calculate updated M*U{n+1} = S
		F77NAME(dgesv)(Nvar_, 1, M, Nvar_, ipiv, S, Nvar_, info);
		// Now update vars for repeat
		F77NAME(dcopy)(Nvar_, U, 1, Un_g, 1); // Save Un-1
		F77NAME(dcopy)(Nvar_, S, 1, U, 1); // Update Un
	}
	writeVec(U, Nx_g, 1, test);
}

void solveImplicit(double *K, double *M, double *F, double *U, double lx_e,
    double qx_, double qy_, double dt_, int Nvar_, int Nx_g, int nite_, int nout_,
    std::string test)
{	// ===================== Build Tables =====================//
	double *K_eff 		= new double[Nvar_ * Nvar_]();
	double *Ud			= new double[Nvar_]();
	double *Udd			= new double[Nvar_]();
	double *tmp			= new double[Nvar_]();
	double *tmp2		= new double[Nvar_]();
	double *S			= new double[Nvar_]();

	// Constanst and coefficients
	const double beta_(.25);
	const double gamma_(.5);
	const double co1_(1/(beta_*(dt_*dt_)));
	const double co2_(1/(beta_*dt_));
	const double co3_((1/(2*beta_))-1);
	const double co4_(dt_*(1-gamma_));
	const double co5_(dt_*gamma_);


	// Build Keff
	for (int i = 0; i < Nvar_*Nvar_; ++i)
	{	double sum = (co1_*M[i]) + K[i];
		K_eff[i] = sum;
		K[i] = sum;
	}

    int info = 0;
    int *ipiv = new int[Nvar_];

	for (int i = 0; i < 2; ++i)
	// for (int i = 0; i < nite_; ++i)
	{	// Create Dynamic Force
		assignArr(F, 0., Nvar_);
		assignArr(S, 0., Nvar_);
		assignArr(tmp, 0., Nvar_);
		assignArr(tmp2, 0., Nvar_);
		updateVars(F, lx_e, qx_, qy_, Nx_g, Nvar_, i, nite_);


		// Sum U, Ud and, Udd with coefficients
		for (int i = 0; i < Nvar_; ++i)
		{	double sum = (co1_*U[i]) + (co2_*Ud[i]) + (co3_*Udd[i]);
			tmp2[i] = sum;
		}

		// Multiple mass by sum
		F77NAME(dgemv)('n', Nvar_, Nvar_, 1, M, Nvar_, tmp2, 1, 0, S, 1);
		// Add forces to sum...
		for (int i = 0; i < Nvar_; ++i)
		{	double sum = S[i] + F[i];
			S[i] = sum;
		}

		// Solve Keff*U_{n+1} = S
		F77NAME(dgesv)(Nvar_, 1, K, Nvar_, ipiv, S, Nvar_, info);

		// Update K to contain only the K_eff as desgv overwrites...
		for (int i = 0; i < Nvar_*Nvar_; ++i)
		{	K[i] = K_eff[i];
		}

		// Now update Udd
		F77NAME(dcopy)(Nvar_, Udd, 1, tmp, 1);   
		for (int i = 0; i < Nvar_; ++i)
		{	double sum = (co1_*(S[i]-U[i])) - (co2_*Ud[i]) -
				(co3_*Udd[i]);
			Udd[i] = sum;
		}

		// Now update Ud
		for (int i = 0; i < Nvar_; ++i)
		{	double sum = U[i] + (co4_*tmp[i]) + (co5_*Udd[i]);
			Ud[i] = sum;
		}

		// Now update U
		F77NAME(dcopy)(Nvar_, S, 1, U, 1);
	}
	writeVec(U, Nx_g, 1, test);
}