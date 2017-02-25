#ifndef PARAMFILE_HPP_INCLUDED
#define PARAMFILE_HPP_INCLUDED

#include <iostream>
#include <fstream>
#include <cmath>
#include <algorithm>

using namespace std;

void readParamFile(ifstream &in_run_input_file, int *nite_, int *Nx_g, double *lx_g,
      double *E_, double *rho_, double *b_, double *h_, double *qx_, double *qy_);

/*!
 * \brief utility function for converting strings to uppercase
 * \param[in, out] str - string we want to convert
 */
inline void StringToLowerCase(string & str)
{ transform(str.begin(), str.end(), str.begin(), ::tolower);
}

inline void StringToUpperCase(string & str)
{ transform(str.begin(), str.end(), str.begin(), ::toupper);
}
/*!
 * \brief utility function for converting strings to uppercase
 * \param[in] str - string we want a copy of converted to uppercase
 * \returns a copy of str in uppercase
 */
// inline string StringToLowerCase(const string & str) {
//   string upp_str(str);
//   transform(upp_str.begin(), upp_str.end(), upp_str.begin(), ::tolower);
//   return upp_str;
// }

bool readLine(string & str, string & option_name, string & option_value);

void analyzeLine(string & keyword, const char *value, int *nite_, int *Nx_g, double *lx_g,
      double *E_, double *rho_, double *b_, double *h_, double *qx_, double *qy_);

string trim(string & str);

string replaceTabsAndReturns(string & str);

void initVars(double *b_, double *h_, double *A_, double *I_, int *Nvar_, int *Nx_g);

#endif // PARAMFILE_HPP_INCLUDED
