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

    // Use blas to solve system...
    F77NAME(dgbsv)(Nvar_, kl, ku, nrhs, K, ldab, ipiv, F, ldb, info);
    writeVec(F, Nx_, test);
	delete[] ipiv;
}

void solveDynamic(double *K, double *M, double *F, double *U, double lx_e,
	double qx_, double qy_, int Nvar_, int Nx_g, int nite_, int nout_,
	std::string test)
{	
	double *eye		= new double[Nvar_ * Nvar_]();
	double *MK_o	= new double[Nvar_ * Nvar_]();
	double *Un_g	= new double[Nvar_]();
	double *S		= new double[Nvar_]();
	double *MKU_o	= new double[Nvar_]();
	double *MU_o	= new double[Nvar_]();

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
		{	writeVec(U, Nx_g, i, test);
		}
	}
}