#include "Solvers.hpp"
#include "Common.hpp"

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