#include <iostream>
#include "Output.hpp"

using namespace std;

void printHeader(){
  cout << "       *************************************************************" << endl;
  cout << "       *                                                           *" << endl;
  cout << "       *     __      ______  ____    ____    ___     ____    ____  *" << endl;
  cout << "       *    / /     / ____/ / __ \\  / __ \\  /   |   / __ \\  / __ \\ *" << endl;
  cout << "       *   / /     / __/   / / / / / /_/ / / /| |  / /_/ / / / / / *" << endl;
  cout << "       *  / /___  / /___  / /_/ / / ____/ / ___ | / _, _/ / /_/ /  *" << endl;
  cout << "       * /_____/ /_____/  \\____/ /_/     /_/  |_|/_/ |_| /_____/   *" << endl;
  cout << "       *                                                           *" << endl;
  cout << "       *         \033[0m \033[1mL\033[0mattic\033[0m\033[1mE\033[0m b\033[0m\033[1mO\033[0mltzmann \033[0m \033[1mP\033[0ml\033[0m\033[1mA\033[0mtefo\033[0m\033[1mR\033[0mm \033[0m\033[1mD\033[0mevelopment         *" << endl;
  cout << "       *************************************************************" << endl;
  cout << "       *                   Copyright 2015, CERFACS                 *" << endl;
  cout << "       *************************************************************\033[0m" << endl;
}

void printInfo(){
  cout << "\n";
  cout << "       ============================================================="<< endl;
  cout << "       =                     Simulation parameters                 ="<< endl;
  cout << "       ============================================================="<< endl;
  // cout << "           " << std::left << std::setw(32) << "LBM modelisation: "         << std::setw(7) << LBMmodel_ << endl;
  // cout << "           " << std::left << std::setw(32) << "Regularisation: "         << std::setw(7) << BoolToString(reg_) << endl;
  // cout << "           " << std::left << std::setw(32) << "Reference Temperature: "    << std::setw(7) << Tref_ << endl;
  // cout << "           " << std::left << std::setw(32) << "Kinematic Viscosity: "      << std::setw(7) << Nu_ << endl;
  // cout << "           " << std::left << std::setw(32) << "Thermal Diffusivity: "      << std::setw(7) << kappa_ << endl;
  // cout << "           " << std::left << std::setw(32) << "Volumic force (S.I. units)" << std::setw(7) << gx_ << " " << gy_ << endl;
  // cout << "           " << std::left << std::setw(32) << "Number of iteration: "      << std::setw(7) << nite_ << endl;
  // cout << "           " << std::left << std::setw(32) << "Result output frequency: "  << std::setw(7) << tResults_ << endl;
  // cout << "           " << std::left << std::setw(32) << "Configuration file: "       << std::setw(7) << config_file << endl;
  // cout << "           " << std::left << std::setw(32) << "R gas constant: "           << std::setw(7) << Rgas_ << endl;
  // cout << "           " << std::left << std::setw(32) << "Gamma: "                    << std::setw(7) << Gamma_ << endl;
  // cout << "           " << std::left << std::setw(32) << "Aero distribution dev. order: " << std::setw(7) << aeroOrder_ << endl;
  // cout << "           " << std::left << std::setw(32) << "Convergence criterion: "    << std::setw(7) << convCriterion_ << endl;
  // cout << "           " << std::left << std::setw(32) << "Residuals comp. frequency: "  << std::setw(7) << tResiduals_ << endl;
  // cout << "           " << std::left << std::setw(32) << "Number of MPI Procs: "  << std::setw(7) << mpi_size << endl;
  // cout << "           " << std::left << std::setw(32) << "MPI Decomposition (x, y): "  << std::setw(0) << mpi_dims[0] << "," << mpi_dims[1] << endl;
  cout << "       ============================================================="<< endl;
}