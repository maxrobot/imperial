#include <iostream>
#include "Output.hpp"
#include "CommonMPI.hpp"

using namespace std;

void printHeader()
{ MPI_Barrier(MPI_COMM_WORLD);
  if (MPI::mpi_rank==0)
  { cout << "       *************************************************************" << endl;
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
  MPI_Barrier(MPI_COMM_WORLD);
}

void printInfo(int Nx_g, int Nx_, int Nvar_, int Nvar_e, std::string sparse_)
{ MPI_Barrier(MPI_COMM_WORLD);
  if (MPI::mpi_rank==0)
  { cout << "\n";
    cout << "       ============================================================="<< endl;
    cout << "       =                     Simulation parameters                 ="<< endl;
    cout << "       ============================================================="<< endl;
    cout << "           " << std::left << std::setw(32) << "Simulation type: "         << std::setw(7) << sparse_ << endl;
    cout << "           " << std::left << std::setw(32) << "Number of processes: "         << std::setw(7) << MPI::mpi_size << endl;
    cout << "           " << std::left << std::setw(32) << "Global elements: "         << std::setw(7) << Nx_g << endl;
    cout << "           " << std::left << std::setw(32) << "Local elements: "         << std::setw(7) << Nx_ << endl;
    cout << "           " << std::left << std::setw(32) << "Global nodes: "         << std::setw(7) << Nvar_ << endl;
    cout << "           " << std::left << std::setw(32) << "Local nodes: "         << std::setw(7) << Nvar_e << endl;
    cout << "       ============================================================="<< endl;
  }
 MPI_Barrier(MPI_COMM_WORLD);
}

void writeVec(double *M, int N, std::string test)
{ ofstream myfile;
  myfile.open ("./output/data/" + test + ".txt");
  // myfile.open ("./output/data/" + test + ".txt");
  myfile << setprecision(10) << 0. << setw(20) << 0.  <<  setw(20) << 0.  << endl;
  for (int i = 0; i < N-1; ++i)
  { int pnt = i*3;
    myfile << setprecision(10) << M[pnt] << setw(20) << M[pnt+1] <<  setw(20) << M[pnt+2] << endl;
  }
  myfile << setprecision(10) << 0. << setw(20) << 0.  <<  setw(20) << 0.  << endl;
  myfile.close();
}

void writeVec(double *M, int N, int step, std::string test)
{ ofstream myfile;
  char numstr[21]; 
  sprintf(numstr, "%d", step);
  myfile.open ("./output/data/" + test + numstr + ".txt");
  myfile << setprecision(10) << 0. << setw(20) << 0.  <<  setw(20) << 0.  << endl;
  for (int i = 0; i < N-1; ++i)
  { int pnt = i*3;
    myfile << setprecision(10) << M[pnt] << setw(20) << M[pnt+1] <<  setw(20) << M[pnt+2] << endl;
  }
  myfile << setprecision(10) << 0. << setw(20) << 0.  <<  setw(20) << 0.  << endl;
  myfile.close();
}

void writeParVec(double *M, int N, int Nvar_g, int Nvar_, int step, std::string test)
{ ofstream myfile;
  char numstr[21]; 
  sprintf(numstr, "%d", step);
  double *output = new double[Nvar_g]();
  MPI::gatherData(output, M, Nvar_g, Nvar_);

  if (MPI::mpi_rank==0)
  { //showVec(output, Nvar_g);
    myfile.open ("./output/data/" + test + numstr + ".txt");
    myfile << setprecision(10) << 0. << setw(20) << 0.  <<  setw(20) << 0.  << endl;
    for (int i = 0; i < N-1; ++i)
    { int pnt = i*3;
      myfile << setprecision(10) << output[pnt] << setw(20) << output[pnt+1] <<  setw(20) << output[pnt+2] << endl;
    }
    myfile << setprecision(10) << 0. << setw(20) << 0.  <<  setw(20) << 0.  << endl;
    myfile.close();
  }
}

void printMessage(std::string message)
{ if (MPI::mpi_rank==0)
  { std::cout << message << std::endl;
  }
}