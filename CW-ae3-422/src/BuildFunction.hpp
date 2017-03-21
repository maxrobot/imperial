#ifndef BUILDFUNCTION_HPP_INCLUDED
#define BUILDFUNCTION_HPP_INCLUDED

#include <iostream>
#include <cmath>

void buildEye(double *K, int 	Nvar_);

void buildKg(double *Kg, double *ke, int Nvar_, int Nx_g);
void buildKgBand(double *Kg, double *ke, int Nvar_, int Nx_g, int buf);
void buildKgBandPar(double *Kg, double *ke, int Nvar_, int Nx_g, int buf);
void buildKgBand(double *Kg, double *ke, int Nvar_, int Nx_g, int buf);

void buildMgBand(double *Kg, double *ke, int Nvar_, int Nx_g, int buf);
void buildMgBandPar(double *Kg, double *ke, int Nvar_, int Nx_g, int buf);

void buildFgBand(double *Kg, double *ke, double coeff, int Nx_g, int Nvar_);
void buildFg(double *Kg, double *ke, int Nx_g, int Nvar_);


void buildMe(double *K, double A_, double rho_, double lx_e, double Al_,
	double dt_);
void buildMe(double *K, double A_, double rho_, double lx_e, double Al_);
void buildKe(double *K, double lx_e, double A_, double E_, double I_);
void buildFe(double *K, double lx_e, double qx_, double qy_);
void buildKeff(double *K, double *K_eff, double *M, double co1_, int Nvar_, int buf_);

void assignArr(double *K, double V, int N);
void updateKglb(double *K, int Nvar_, int buf_);

void updateVars(double *F, double lx_e, double qx_,	double qy_, int Nx_g,
	int Nvar_, int step, int nite_);
void updateParVars(double *F, double lx_e, double qx_,double qy_, int Nx_g,
	int Nvar_, int step, int nite_);

#endif // BUILDFUNCTION_HPP_INCLUDED
