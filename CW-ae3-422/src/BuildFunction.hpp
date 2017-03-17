#ifndef BUILDFUNCTION_HPP_INCLUDED
#define BUILDFUNCTION_HPP_INCLUDED

#include <iostream>
#include <cmath>

void buildEye(double *K, int 	Nvar_);

void buildKglb(double *Kg, double *ke, int Nvar_, int Nx_g);
void buildSparse(double *Kg, double *ke, int Nvar_, int Nx_g, int buf);
void buildSparsePar(double *Kg, double *ke, int Nvar_, int Nx_g, int buf);
void buildKglbSparse(double *Kg, double *ke, int Nvar_, int Nx_g, int buf);
void buildMglbSparse(double *Kg, double *ke, int Nvar_, int Nx_g, int buf);
void buildMglbPar(double *Kg, double *ke, int Nvar_, int Nx_g, int buf);
void buildBandFglb(double *Kg, double *ke, int Nx_g, int Nvar_);
void buildFglb(double *Kg, double *ke, int Nx_g, int Nvar_);


void buildMele(double *K, double A_, double rho_, double lx_e, double Al_,
	double dt_);
void buildMele(double *K, double A_, double rho_, double lx_e, double Al_);
void buildKele(double *K, double lx_e, double A_, double E_, double I_);
void buildFele(double *K, double lx_e, double qx_, double qy_);

void assignArr(double *K, double V, int N);

void updateVars(double *F, double lx_e, double qx_,	double qy_, int Nx_g,
	int Nvar_, int step, int nite_);
void updateParVars(double *F, double lx_e, double qx_,double qy_, int Nx_g,
	int Nvar_, int step, int nite_);

#endif // BUILDFUNCTION_HPP_INCLUDED
