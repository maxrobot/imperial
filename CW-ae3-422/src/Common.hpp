#ifndef COMMON_HPP_INCLUDED
#define COMMON_HPP_INCLUDED

#include <iostream>
#include <cmath>
#include <fstream>
#include <sstream>
#include <iomanip>
// #include <stdio.h>

// #include <string.h>

// Define LAPACK Shit ###########################################
#define F77NAME(x) x##_
extern "C" {
    double F77NAME(dnrm2)(const int& N, const double *X, const int& incX);
    void F77NAME(daxpy)(const int& N, const double& alpha, const double *X,
                             const int& incX, double *Y, const int& incY);
    void F77NAME(dgemv)(const char& trans, const int& m, const int& n,
                    const double& alpha, const double* a, const int& lda,
                    const double* x, const int& incx, const double& beta,
                    double* y, const int& incy);
    void F77NAME(dgesv)(const int& n, const int& nrhs, const double * A,
                    const int& lda, int * ipiv, double * B, const int& ldb,
                    int& info);
}

void showVec(double *M, int N);

void showMat(double *M, int N);

void writeVec(double *M, int N, std::string test);

#endif // COMMON_HPP_INCLUDED
