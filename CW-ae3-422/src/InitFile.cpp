#include "InitFile.hpp"
#include "GlobalVars.hpp"

using namespace std;

void readParamFile(ifstream& in_run_input_file)
{ 
  string text_line;
  string keyword, value;

  // Now read in parameters
  while (getline (in_run_input_file, text_line)) {
    if (readLine(text_line, keyword, value)){
      analyzeLine(keyword, value.c_str());
    }
  }
}

void analyzeLine(string & keyword, const char * value){
  if (!keyword.compare("number of iterations"))
  { sscanf(value,"%d",&nite_);
  }
  if (!keyword.compare("number of elements"))
  { string tmp_value = string(value);
    sscanf(value,"%d",&Nx_g);
  }
  if (!keyword.compare("beam length (mm)"))
  { string tmp_value = string(value);
    sscanf(value,"%lf",&lx_g);
  }
  if (!keyword.compare("youngs modulus (mpa)"))
  { string tmp_value = string(value);
    sscanf(value,"%lf",&E_);
  }
  if (!keyword.compare("rho (ton/mm^3)")) sscanf(value,"%lf",&rho_);
  if (!keyword.compare("cross-sectional width (mm)"))  sscanf(value,"%lf",&b_);
  if (!keyword.compare("cross-sectional height (mm)")) sscanf(value,"%lf", &h_);      
  if (!keyword.compare("axial uniform load (n/mm)")) sscanf(value,"%lf", &qx_);     
  if (!keyword.compare("traverse uniform load (n/mm)")) sscanf(value,"%lf", &qy_); 
}


bool readLine(string & str, string & keyword, string & value) {

  const string delimiters(" \t\n\r");
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
{
  if(str.empty())
      return str;
  size_t firstScan = str.find_first_not_of(' ');
  size_t first     = firstScan == string::npos ? str.length() : firstScan;
  size_t last      = str.find_last_not_of(' ');
  return str.substr(first, last-first+1);
}

string replaceTabsAndReturns(string & str)
{
  size_t pos;
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

// Initialise local and global arrays...
void initVars()
{ A_ = b_*h_;                   // Cross-sectional Area Calculation
  I_ = (b_*pow(h_,3.))/12;      // Second Moments of area calculation
  Nvar_ = (Nx_g+1)*3;          // Number of variables in global matrices
}