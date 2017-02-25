#ifndef BUILDFUNCTION_HPP_INCLUDED
#define BUILDFUNCTION_HPP_INCLUDED

#include <iostream>
#include <cmath>

void buildKglb(double *Kg, double *ke, int Nvar_, int Nx_g);
void buildFglb(double *Kg, double *ke, int Nx_g);

void buildKele(double *K, double lx_e, double A_, double E_, double I_);
void buildFele(double *K, double lx_e, double qx_, double qy_);

void zeroMat(double *K, int N);

#endif // BUILDFUNCTION_HPP_INCLUDED
