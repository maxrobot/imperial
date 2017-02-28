#include "Solvers.hpp"

void solveStatic(double *A, double* b, int Nvar_)
{	// =================== Solve System =======================//
    const int nrhs = 1;
    int info = 0;
    int *ipiv = new int(Nvar_);
    int kl = 4;
    int ku = 4;
    int ldab = 1 + 2*kl + ku;
    int ldb = Nvar_;
	F77NAME(dgbsv)(Nvar_, kl, ku, nrhs, A, ldab, ipiv, b, ldb, info);
}