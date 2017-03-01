#ifndef BUILDFUNCTION_HPP_INCLUDED
#define BUILDFUNCTION_HPP_INCLUDED

#include <iostream>
#include <cmath>

void buildEye(double *K, int 	Nvar_);

void buildKglb(double *Kg, double *ke, int Nvar_, int Nx_g);
void buildKglbSparse(double *Kg, double *ke, int Nvar_, int Nx_g, int buf);
void buildFglb(double *Kg, double *ke, int Nx_g);

void buildMele(double *K, double A_, double rho_, double lx_e, double Al_, double dt_);
void buildKele(double *K, double lx_e, double A_, double E_, double I_);
void buildFele(double *K, double lx_e, double qx_, double qy_);

void assignArr(double *K, double V, int N);

#endif // BUILDFUNCTION_HPP_INCLUDED
