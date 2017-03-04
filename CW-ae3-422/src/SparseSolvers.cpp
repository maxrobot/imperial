#include "Solvers.hpp"
#include "Common.hpp"
#include "BuildFunction.hpp"

using namespace std;


void solveSparseExplicit(double *K, double *M, double *F, double *U, double lx_e,
	double qx_, double qy_, int Nvar_, int Nx_g, int nite_, int nout_,
	int buf_, std::string test)
{	
	double *MK_o	= new double[Nvar_ * 9]();
	double *Un_g	= new double[Nvar_]();
	double *S		= new double[Nvar_]();
	double *Mt_		= new double[Nvar_]();
	double *MKU_o	= new double[Nvar_]();
	double *MU_o	= new double[Nvar_]();

	int kl = 4;
	int ku = 4;
    int info = 0;
    int *ipiv = new int[Nvar_];

    // Make the initial matrix...
	for (int i = 0; i < Nvar_; ++i)
	{	for (int j = 0; j < 9; ++j)
		{	int pnt = (j+4) + (i*13);
			int pnt2 = j + (i*9);
			MK_o[pnt2] = K[pnt];
		}
	}

    for (int i = 0; i < Nvar_; ++i)
    {	int pnt = 4 + (i*9);
    	MK_o[pnt] -= 2*M[i];
    }    	
	F77NAME(dgbmv)('n', Nvar_, 9, 4, 4, 1, MK_o, Nvar_, U, 1, 0, MKU_o, 1);

	showMat(MK_o, 9, Nvar_);
	showVec(U, Nvar_);
	showVec(MKU_o, Nvar_);
	// =================== Create S Matrix ====================//
	// Start marching through time...
	// for (int i = 0; i <= nite_; ++i)
	// for (int i = 0; i <= 5; ++i)
	// {	// Make C: output = L_c
	// 	F77NAME(dcopy)(Nvar_, M, 1, Mt_, 1); // Save Un-1
	// 	assignArr(F, 0., Nvar_);
	// 	updateVars(F, lx_e, qx_, qy_, Nx_g, Nvar_, i, nite_);

	// 	// Calculate M*U{n-1}
	// 	for (int i = 0; i < Nvar_; ++i)
	// 	{	MU_o[i] = Mt_[i]*Un_g[i];
	// 	}

	// 	// Calculate MK_o*U{n}
	// 	// F77NAME(dgbmv)('n', Nvar_, 9, 4, 4, 1, MK_o, Nvar_, U, 1, 0, MKU_o, 1);
	// 	// for (int i = 0; i < Nvar_; ++i)
	// 	// {	double sum = F[i] - MKU_o[i] - MU_o[i];
	// 	// 	S[i] = sum;
	// 	// }
	// 	// showMat(MK_o,9, Nvar_);
	// 	showVec(F, Nvar_);
	// 	showVec(S, Nvar_);

	// 	// Calculate updated M*U{n+1} = S
	// 	// showMat(MK_o, 9, Nvar_);
	//     F77NAME(dgbsv)(Nvar_, 0, 0, 1, Mt_, 1, ipiv, S, Nvar_, info);
	// 	// showVec(S, Nvar_);

	// 	// Now update vars for repeat
	// 	F77NAME(dcopy)(Nvar_, U, 1, Un_g, 1); // Save Un-1
	// 	F77NAME(dcopy)(Nvar_, S, 1, U, 1); // Update Un
	// 	// if ((i%nout_)==0)
	// 	// {	showVec(M, Nvar_);
	// 	// }
	// }
	writeVec(U, Nx_g, 1, test);
}