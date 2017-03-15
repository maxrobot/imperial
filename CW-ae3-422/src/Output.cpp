#include <iostream>
#include "Output.hpp"
#include "CommonMPI.hpp"

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

void printInfo(int Nx_g, int Nx_, int Nvar_, int Nvar_e)
{ if (MPI::mpi_rank==0)
  { cout << "\n";
    cout << "       ============================================================="<< endl;
    cout << "       =                     Simulation parameters                 ="<< endl;
    cout << "       ============================================================="<< endl;
    cout << "           " << std::left << std::setw(32) << "Global element numbers: "         << std::setw(7) << Nx_g << endl;
    cout << "           " << std::left << std::setw(32) << "Local element numbers: "         << std::setw(7) << Nx_ << endl;
    cout << "           " << std::left << std::setw(32) << "Global nodes numbers: "         << std::setw(7) << Nvar_ << endl;
    cout << "           " << std::left << std::setw(32) << "Local nodes numbers: "         << std::setw(7) << Nvar_e << endl;
    cout << "           " << std::left << std::setw(32) << "Number of processes: "         << std::setw(7) << MPI::mpi_size << endl;
    cout << "       ============================================================="<< endl;
  }
}