#ifndef PARAMFILE_HPP_INCLUDED
#define PARAMFILE_HPP_INCLUDED

#include <iostream>
#include <fstream>
#include <cmath>
#include <algorithm>

using namespace std;

void readParamFile(ifstream &in_run_input_file, int *T_, int *nite_, int *Nx_g,
    int *nout_, double *lx_g, double *E_, double *rho_, double *b_, double *h_,
    double *qx_, double *qy_, string *eq_, string *scheme_, string *sparse_);

inline void StringToLowerCase(string & str)
{ transform(str.begin(), str.end(), str.begin(), ::tolower);
}

inline void StringToUpperCase(string & str)
{ transform(str.begin(), str.end(), str.begin(), ::toupper);
}

bool readLine(string & str, string & option_name, string & option_value);

void analyzeLine(string &keyword, const char *value, int *T_, int *nite_, int *Nx_g,
    int *nout_, double *lx_g, double *E_, double *rho_, double *b_, double *h_,
    double *qx_, double *qy_, string *eq_, string *scheme_, string *sparse_);

string trim(string & str);

string replaceTabsAndReturns(string & str);

void initVars(double *b_, double *h_, double *A_, double *I_,  double *E_,
	double *dt_, int *Nvar_, int *Nvar_e, int *Nghost_, int *Sghost_, int *Nx_g,
	int *Nx_, int *T_, int *nite_);

#endif // PARAMFILE_HPP_INCLUDED
