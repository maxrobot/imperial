#ifndef GLOBALVARS_HPP_INCLUDED
#define GLOBALVARS_HPP_INCLUDED

  // #include "commonMPI.hpp"

  // using namespace MPI;
  using namespace std;

  extern int Nx_g;          // Number of global elements
  extern int nite_;         // Number of time steps
  extern double lx_g;       // Length of global domain
  extern double dt_;        // Time step
  extern double dx_;        // Mesh size

  extern double rho_;       // Density
  extern double I_;         // Second moment area
  extern double E_;         // Youngs modulus
  extern double b_;         // Cross-sectional width
  extern double h_;         // Cross-sectional height
  extern double A_;         // Cross-sectional area


#endif // GLOBALVARS_HPP_INCLUDED