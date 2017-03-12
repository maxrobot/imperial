#ifndef SOLVERS_HPP_INCLUDED
#define SOLVERS_HPP_INCLUDED

#include <iostream>

// // Define LAPACK Shit ###########################################
// #define F77NAME(x) x##_
// extern "C" {
//     double F77NAME(dnrm2)(const int &N, const double *X, const int &incX);
//     void F77NAME(daxpy)(const int &N, const double &alpha, const double *X,
//                              const int &incX, double *Y, const int &incY);
//     void F77NAME(dgemv)(const char& trans, const int& m, const int& n,
//                     const double& alpha, const double* a, const int& lda,
//                     const double* x, const int& incx, const double& beta,
//                     double* y, const int& incy);
//     void F77NAME(dgesv)(const int& n, const int& nrhs, const double * A,
//                     const int& lda, int * ipiv, double * B, const int& ldb,
//                     int& info);
//     void F77NAME(dgbsv)(const int &N, const int &kl, const int &ku, const int &nrhs,
//                     const double *A, const int &ldab, const int *ipiv,
//                     const double *b, const int &ldb, const int &info);
// }

void runSolver(double *K_e, double *U_g, double *F_g, double lx_e,
    double A_, double E_, double I_, double qx_, double qy_, int Nvar_,
    int Nx_g);

void runSolver(double *K_e, double *U_g, double *F_g, double dt_, 
    double lx_e, double A_, double E_, double I_, double rho_,
    double qx_, double qy_, int Nvar_, int Nvar_e, int Nx_g, int nite_,
    int nout_, const int buf_, std::string sparse_);

void solveStatic(double *K, double *F, int Nvar_, int ldab, int Nx_,
    std::string test);

void solveExplicit(double *K, double *M, double *F, double *U, double lx_e,
    double qx_, double qy_, int Nvar_, int Nx_g, int nite_, int nout_,
    int buf_, std::string test);

void solveSparseExplicit(double *K, double *M, double *F, double *U, double lx_e,
    double qx_, double qy_, int Nvar_, int Nx_g, int nite_, int nout_,
    int buf_, std::string test);

void solveParallelSparseExplicit(double *K, double *M, double *F, double *U, double lx_e,
    double qx_, double qy_, int Nvar_, int Nx_g, int nite_, int nout_,
    int buf_, std::string test);

void solveImplicit(double *K, double *M, double *F, double *U, double lx_e,
    double qx_, double qy_, double dt_, int Nvar_, int Nx_g, int nite_, int nout_,
    std::string test);

void solveSparseImplicit(double *K, double *M, double *F, double *U, double lx_e,
    double qx_, double qy_, double dt_, int Nvar_, int Nx_g, int nite_, int nout_,
    int buf_, std::string test);

#endif // SOLVERS_HPP_INCLUDED