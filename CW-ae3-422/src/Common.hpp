#ifndef COMMON_HPP_INCLUDED
#define COMMON_HPP_INCLUDED

#include <iostream>
#include <cmath>
#include <fstream>
#include <sstream>
#include <iomanip>


// Define LAPACK Shit ###########################################
#define F77NAME(x) x##_
extern "C" {
    double F77NAME(ddot)(const int &N, const double *X, const int &incx, 
                const double *Y, const int &incy);
    double F77NAME(dnrm2)(const int &N, const double *X, const int &incX);
    void F77NAME(dcopy)(const int &N, double *DX, const int &incx,
                double *DY, const int &incy);   
    void F77NAME(dscal)(const int &N, const double &alpha, const double *X,
                const int &incX);
    void F77NAME(daxpy)(const int &N, const double &alpha, const double *X,
                const int &incX, double *Y, const int &incY);
    void F77NAME(dgemv)(const char& trans, const int& m, const int& n,
                const double& alpha, const double* a, const int& lda,
                const double* x, const int& incx, const double& beta,
                double* y, const int& incy);
    void F77NAME(dgemm)(const char &transa, const char &transb, const int &m,
                const int &n, const int &k, const double &alpha, double *A,
                const int &lda, double *B, const int &ldb, const double &beta,
                double *C, const int &ldc);
    void F77NAME(dgesv)(const int& n, const int& nrhs, const double * A,
                const int& lda, int * ipiv, double * B, const int& ldb,
                int& info);
    void F77NAME(dgbsv)(const int &N, const int &kl, const int &ku, const int &nrhs,
                const double *A, const int &ldab, const int *ipiv,
                const double *b, const int &ldb, const int &info);
}

void showVec(double *M, int N);

void showMat(double *M, int N);
void showMat(double *M, int N, int O);

void writeVec(double *M, int N, std::string test);
void writeVec(double *M, int N, int step, std::string test);

#endif // COMMON_HPP_INCLUDED
