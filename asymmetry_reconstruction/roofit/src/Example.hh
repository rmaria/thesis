//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Fri Feb 18 16:10:55 2011 by ROOT version 5.26/00
// from TTree T/T
// found on file: hits_DAFNE_nomTraject_posiBeam_noVtitl.root
//////////////////////////////////////////////////////////

#ifndef Example_hh
#define Example_hh

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TH1.h>
#include <TH2.h>
#include <TH2D.h>
#include <TString.h>
#include <TLorentzVector.h>
#include <TLorentzRotation.h>
#include <TVector3.h>
#include <TString.h>





//===============================================================
void Example(
		 double &argLfFit,
	     double &argLfFitError,
	     double &minusLog,
	     double &edm,
	     double &stat,
		 double tmean           =  40.1,
	     double DeltaGammaTau   =  1.43e-2,
	     double DeltaM          =  1.18e+4,
	     double Mod_lambda_f    =  1.0,
	     double Phi_lambda_f    = 45.0,
	     double ResModel_Bias   = 0.0,
	     double ResModel_Width  = 1.0e-12,
	     double Omega           = 0.0,
	     double DeltaOmega      = 0.0,
	     int Nevents            = 10000,
	     int seed               = 9383995,
	     const char* Option     = "SingleSided",
	     const char* OutputFile = "Plots/Example"
	     );
//===============================================================

#endif

