#ifndef SOLVERS_HPP_INCLUDED
#define SOLVERS_HPP_INCLUDED

#include <iostream>

void runSolver(double *K_e, double lx_e, double A_, double E_, double I_,
    double qx_, double qy_, int Nvar_, int Nx_g);

void runSolver(double *K_e, double dt_, double lx_e, double A_, double E_,
    double I_, double rho_, double qx_, double qy_, int Nvar_, int Nvar_e,
    int Nghost_, int Nx_g, int Nx_, int nite_, int nout_, const int buf_,
    std::string sparse_);

void solveStatic(double *K, double *F, int Nvar_, int ldab, int Nx_,
    std::string test);

void solveExplicit(double *K, double *M, double *F, double lx_e,
    double qx_, double qy_, int Nvar_, int Nx_g, int nite_, int nout_,
    int buf_, std::string test);

void solveSparseExplicit(double *K, double *M, double *F, double lx_e,
    double qx_, double qy_, int Nvar_, int Nx_g, int nite_, int nout_,
    int buf_, std::string test);

void solveParSparseExplicit(double *K, double *M, double *F, double lx_e,
    double qx_, double qy_, int Nvar_, int Nghost_, int Nx_g, int nite_, int nout_,
    int buf_, std::string test);

void solveImplicit(double *K, double *M, double *F, double lx_e, double qx_,
    double qy_, double dt_, int Nvar_, int Nx_g, int nite_, int nout_,
    std::string test);

void solveSparseImplicit(double *K, double *M, double *F, double lx_e, double qx_,
    double qy_, double dt_, int Nvar_, int Nx_g, int nite_, int nout_, int buf_,
    std::string test);

#endif // SOLVERS_HPP_INCLUDED