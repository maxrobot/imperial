#ifndef OUTPUT_HPP_INCLUDED
#define OUTPUT_HPP_INCLUDED

#include <iostream>
#include <cmath>
#include <fstream>
#include <sstream>
#include <iomanip>

void printHeader();

void printInfo(int Nx_g, int Nx_, int Nvar_, int Nvar_e, std::string sparse_);

void writeVec(double *M, int N, std::string test);
void writeVec(double *M, int N, int step, std::string test);
void writeParVec(double *M, int N, int Nvar_g, int Nvar_, int step, std::string test);

void printMessage(std::string message);

#endif // OUTPUT_HPP_INCLUDED