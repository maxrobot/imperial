#include "InitFile.hpp"
#include "Common.hpp"
#include "CommonMPI.hpp"

using namespace std;

// Moves through param file and reads each line
void readParamFile(ifstream &in_run_input_file, int *T_, int *nite_, int *Nx_g,
    int *nout_, double *lx_g, double *E_, double *rho_, double *b_, double *h_,
    double *qx_, double *qy_, string *eq_, string *scheme_, string *sparse_)
{ 
  string text_line;
  string keyword, value;
  // Now read in parameters
  while (getline (in_run_input_file, text_line))
  { if (readLine(text_line, keyword, value)){
      analyzeLine(keyword, value.c_str(), T_, nite_, Nx_g, nout_, lx_g, E_,
        rho_, b_, h_, qx_, qy_, eq_, scheme_, sparse_);
    }
  }
}

// Look at each line and see if the is something we need
void analyzeLine(string &keyword, const char *value, int *T_, int *nite_, int *Nx_g,
    int *nout_, double *lx_g, double *E_, double *rho_, double *b_, double *h_,
    double *qx_, double *qy_, string *eq_, string *scheme_, string *sparse_)
{ 
  if (!keyword.compare("equation")){
    string tmp_value = string(value);
    *eq_ = string(tmp_value.c_str());
  }
  if (!keyword.compare("integration")){
    string tmp_value = string(value);
    *scheme_ = string(tmp_value.c_str());
  }
  if (!keyword.compare("matrix format")){
    string tmp_value = string(value);
    *sparse_ = string(tmp_value.c_str());
  }
  if (!keyword.compare("simulation time (s)"))
  { sscanf(value,"%d",T_);
  }
  if (!keyword.compare("output spacing"))
  { sscanf(value,"%d",nout_);
  }
  if (!keyword.compare("number of iterations"))
  { sscanf(value,"%d",nite_);
  }
  if (!keyword.compare("number of elements"))
  { string tmp_value = string(value);
    sscanf(value,"%d",Nx_g);
    if (*Nx_g % 2)
    { printMessage("Please choose even number of elements!");
      exit(EXIT_FAILURE);
    }
  }
  if (!keyword.compare("beam length (mm)"))
  { string tmp_value = string(value);
    sscanf(value,"%lf",lx_g);
  }
  if (!keyword.compare("youngs modulus (mpa)"))
  { string tmp_value = string(value);
    sscanf(value,"%lf",E_);
  }
  if (!keyword.compare("rho (ton/mm^3)")) sscanf(value,"%lf",rho_);
  if (!keyword.compare("cross-sectional width (mm)"))  sscanf(value,"%lf",b_);
  if (!keyword.compare("cross-sectional height (mm)")) sscanf(value,"%lf",h_);
  if (!keyword.compare("axial uniform load (n/mm)")) sscanf(value,"%lf",qx_);     
  if (!keyword.compare("traverse uniform load (n/mm)")) sscanf(value,"%lf",qy_);
}


bool readLine(string & str, string & keyword, string & value)
{ const string delimiters(" \t\n\r");
  string::size_type pos, last_pos;

  if(str.empty()) return false; // empty str 

  replaceTabsAndReturns(str);
  pos = str.find_first_of("#");
  if (pos == 0)  return false;  // a full comment line
  if (pos != string::npos) str.erase(pos); // remove comment at end if necessary

  pos = str.find("//");
  if (pos != string::npos) str.erase(pos); // remove comment at end if necessary
    
  // find the : sign and split string
  string name_part, value_part;
  pos = str.find(":");
  if (pos == string::npos) return false;

  name_part = str.substr(0, pos);
  value_part = str.substr(pos+1, string::npos);

  // remove left white space or tabs
  last_pos = name_part.find_first_not_of(delimiters, 0);
  pos = name_part.find_first_of(delimiters, last_pos);
  keyword = name_part.substr(last_pos, name_part.length() - last_pos);
  keyword = trim(keyword);
  StringToLowerCase(keyword);

  // remove left white space or tabs
  last_pos = value_part.find_first_not_of(delimiters, 0);
  pos = value_part.find_first_of(delimiters, last_pos);
  value = value_part.substr(last_pos, value_part.length() - last_pos);
  value = trim(value);
  if (keyword.compare("configuration file")) StringToLowerCase(value);

  return true;
}

string trim(string & str)
{ if(str.empty())
      return str;
  size_t firstScan = str.find_first_not_of(' ');
  size_t first     = firstScan == string::npos ? str.length() : firstScan;
  size_t last      = str.find_last_not_of(' ');
  return str.substr(first, last-first+1);
}

string replaceTabsAndReturns(string & str)
{ size_t pos;
  if(str.empty())
      return str;

  while (1){
    pos = str.find_first_of("\t");
    if (pos == string::npos) break;
    str.replace(pos,1," ");      
  }

  while (1){
    pos = str.find_first_of("\r");
    if (pos == string::npos) break;
    str.replace(pos,1," ");      
  }

  return str;
}

void initVars(double *b_, double *h_, double *A_, double *I_, double *E_,
  double *dt_, int *Nvar_, int *Nvar_e, int *Nghost_, int *Sghost_, int *Nx_g,
  int *Nx_, int *T_, int *nite_)
{ // Initialise most value stuffs
  *A_ = *b_ * *h_;                   // Cross-sectional Area Calculation
  *I_ = (*b_ * pow(*h_,3.))/12;      // Second Moments of area calculation
  *E_ = *E_;
  *dt_ = double(*T_)/ *nite_;
  *Nvar_ = (*Nx_g-1)*3;          // Number of variables in global matrices excluding boundaries

  // Initialise domain decomposition using bool to ensure it does so correclty...
  
  int pnt = 0;

  // First check number of mpi processes
  if (MPI::mpi_size%2 && MPI::mpi_size>1)
  { printMessage("Error: odd  number of proc, choose even amount!");
    exit(EXIT_FAILURE);
  }
  bool safe_decomp = false;

  if (MPI::mpi_size>1)
  { while(safe_decomp!=true)
    { 
      *Nvar_ = (*Nx_g-1)*3;          // Number of variables in global matrices excluding boundaries
      *Nvar_e = (*Nvar_ + (MPI::mpi_size-1)*3)/MPI::mpi_size;
      if (*Nvar_e%2)
      { printMessage("Error: Domain decomposed irregularly!");
        printMessage("Trying to fix!");
        safe_decomp = false; // Flick switch to say domain not happy
        Nx_g += 2;
        pnt +=1;
      }
      else
      {
        safe_decomp = true;
      }
      if (pnt==3)
      { printMessage("Fix didn't work, exiting!");
        exit(EXIT_FAILURE);
      }
      *Nvar_ = (*Nx_g-1)*3;          // Number of variables in global matrices excluding boundaries
      *Nvar_e = (*Nvar_ + (MPI::mpi_size-1)*3)/MPI::mpi_size;
      *Nx_ = *Nvar_e/3;
    }
  } 
  else if (MPI::mpi_size==1)
  { *Nvar_e = *Nvar_;
    *Nx_ = *Nx_g;
  }  
  // Initialise ghost cells we use 3 ghost cells as 4th makes no contribution 
  *Sghost_ = 3;
  if (MPI::mpi_size==1)
  { *Nghost_ = *Nvar_e;
  }
  if (MPI::mpi_size>1)
  { if (MPI::mpi_rank==0 || MPI::mpi_rank==(MPI::mpi_size-1))
    { *Nghost_ = *Nvar_e+ *Sghost_;
    }
    if (MPI::mpi_rank>0 && MPI::mpi_rank<(MPI::mpi_size-1))
    { *Nghost_ = *Nvar_e+(2* *Sghost_);
    }
  }  
}