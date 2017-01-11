#include "TROOT.h"
#include "src/Example.hh"
#include <stdio.h>
#include <stdlib.h>
#include <stdexcept>
#include <fstream>
#include <sstream>
#include <iostream>
#include <string>
#include <math.h>
#include <iomanip>
#include <assert.h>


using namespace std;

int main(int argc, char *argv[]) {

  //NOTE: - Time is in ps
  //      - Angles/phases in degrees

  double tmean           =  0.4101;           //ps^-1
  double DeltaGammaTau   =  1.43e-2;        //Delta(Gamma)*tau
  double DeltaM          =  1.18e+1;        // ps^-1 ?? changed from e+4 to e+10 and u.m. is hs-1
  double Mod_lambda_f    =  1.0;            //Module of lambda_f
  double Phi_lambda_f    = 10.0;            //Phase of lambda_f (in degrees)
  double ResModel_Bias   = 0.0;             //ps. Bias  of resolution model (Gaussian function)
  double ResModel_Width  = 0.001;         //ps. Width of resolution model (Gaussian function)
  double Omega           = 0.0;             //Average mistag probabity (omega_D0 + omega_D0bar)/2
  double DeltaOmega      = 0.0;             //Mistag probability difference (omega_D0 - omega_D0bar)
  int Nevents            = 5000000;           //Number of events to be generated
  int seed               = 8918793;         //Seed
  const char* Option     = "SingleSided";   //Model option: SingleSided or DoubleSided
  const char* OutputFile = "Plots/Example"; //Ouput file name

  if(argc == 14) {
    tmean          = atof(argv[1]);
    DeltaGammaTau  = atof(argv[2]);
    DeltaM         = atof(argv[3]);
    Mod_lambda_f   = atof(argv[4]);
    Phi_lambda_f   = atof(argv[5]);
    ResModel_Bias  = atof(argv[6]);
    ResModel_Width = atof(argv[7]);
    Omega          = atof(argv[8]);
    DeltaOmega     = atof(argv[9]);
    Nevents        = atoi(argv[10]);
    seed           = atoi(argv[11]);
    Option         = argv[12];    
    OutputFile     = argv[13];
  }
  else {
    cout << endl;
    cout << "Wrong number of input parameters" << endl;
    cout << endl;
    assert(false);
  }

  if(TString(Option) != TString("SingleSided") && 
     TString(Option) != TString("DoubleSided")) {
    cout << endl;
    cout << "The Option variable can only have the values:" << endl;
    cout << " - SingleSided," << endl;
    cout << " - DoubleSided." << endl;
    cout << "Check your inputs. Exiting now!!!" << 
    cout << endl;
    assert(false);
  }
  
  
  //int seed = rand() % 100000 + 1;
  cout << endl;
  cout << "=============================================" << endl;
  cout << "The function arguments are:"                   << endl;
  cout << "---------------------------------------------" << endl;
  cout << "tau                  = " << tmean          << " ps"      << endl;
  cout << "Delta(Lambda)*tau    = " << DeltaGammaTau                << endl;
  cout << "DeltaM               = " << DeltaM         << " ps^-1"   << endl;
  cout << "|lambda_f|           = " << Mod_lambda_f                 << endl;
  cout << "Phase(lambda_f)      = " << Phi_lambda_f   << " degrees" << endl;
  cout << "Res Model Bias       = " << ResModel_Bias  << " ps"      << endl;
  cout << "Res Model Width      = " << ResModel_Width << " ps"      << endl;
  cout << "omega                = " << Omega                        << endl;
  cout << "Delta(omega)         = " << DeltaOmega                   << endl;
  cout << "# to be generated    = " << Nevents                      << endl;
  
  cout << "Seed                 = " << seed                         << endl;
  cout << "Option               = " << Option                       << endl;
  cout << "OutputFile           = " << OutputFile                   << endl;
  cout << "=============================================" << endl;
	}
	myfile.close();
  
  Example(tmean,
	  DeltaGammaTau,
	  DeltaM,
	  Mod_lambda_f,
	  Phi_lambda_f,
	  ResModel_Bias,
	  ResModel_Width,
	  Omega,
	  DeltaOmega,
	  Nevents,
	  seed,
	  Option,
	  OutputFile);

  return 0;  

}
