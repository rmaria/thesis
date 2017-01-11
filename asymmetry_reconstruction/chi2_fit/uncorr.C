  // --------------------------------------------------------------------
//
// Set of functions to evaluate the measurement precision of CP-violation
//  parameters in the D0-D0bar system with time dependent asymmetry.
// The method is based on a toy Monte-Carlo:
//  - asymmetry returns the value of the asymmetry for a given Delta(t) between
//     the two D0s,
//  - simuAsymmetry returns the evaluation of the lambda_f parameter 
//     for a single toy-experiment with a given number of D0-D0bar pairs, 
//     uncertainty on Delta(t) and mistag probability,
//  - evalAsymmetry returns the lambda_f evaluation distribution over a given 
//     number of toy-experiments.
//
// Inputs are specified either as arguments of function 
//  or in the simuASymmetry function.
//
// Usage in Root:
//   Root> .L simuAsymmetry.C+ // compilation step
// Run a single toy-experiment with 10000 pairs / 20 bins / generated Argument = 0.1 (degrés) / time uncertainty = 0.01e-12, mistag = 5 % / verbose: 
//   Root> simuAsymmetry(10000, 20, 0.2, .01e-12, 0.05, 1) 
//   Root> evalAsymmetry(100, 10000, 20, 3., .01e-12, 0.05) // average over 100 toy-experiments
//
// Created: JB 2012/03/20  
// Modified: IRB May 2012 
//
// --------------------------------------------------------------------

#include "TF1.h"
#include "TH1.h"
#include "TH2.h"
#include "TCanvas.h"
#include "TRandom.h"
#include "TStyle.h"
#include "TMath.h"
#include "RooHist.h"


// Global Input parameters
/*
Double_t x = 0.005; // deltam/gamma
Double_t y = 0.015; // deltag/2/gamma
Double_t width = 1./410.e-15; // (J) hbar(1.e-34 J, 197 326 963 eV)/tau
// le boost intervient quand on veut interpreter la resolution necessaire sur delta_t en delta_z (detecteur)
//Double_t deltag = y*2.*width; // PDG average 1.5
//Double_t deltam = x*width; // PDG average 2.39
Double_t deltag = 2*y*width; // PDG average 0.015*2*gamma
Double_t deltam = x*width; // PDG average 2.39e10

Double_t modLf_gen =  1.;
*/
Double_t tau = 410e-15;;
Double_t width = 1./tau;
Double_t DeltaGammaTau=1.44e-2;      //1.44e-2 pdg 2013
Double_t deltag = DeltaGammaTau/tau;
Double_t deltam = 1.18e-2;        //1.18e-2 pdg 2013
Double_t modLf_gen =  1.; 


// ---------------------------------------------------------------------------
Double_t asymmetry( Double_t *x, Double_t *par) {

  Double_t deltat = x[0];
  
  Double_t deltag = par[0];
  Double_t deltam = par[1];  
  Double_t modLf = par[2];
  Double_t argLf = par[3];
  Double_t mistag = par[4];
  
  Double_t absLf2 = modLf*modLf;
  Double_t reLf = modLf *cos(argLf*TMath::Pi()/180.);
  Double_t imLf = modLf *sin(argLf*TMath::Pi()/180.);

  Double_t hplus  = 1. + exp(deltag*TMath::Abs(deltat));
  Double_t hminus = 1. - exp(deltag*TMath::Abs(deltat));
  Double_t dilution = 1.-2*mistag;

  return 2.*dilution*exp(deltag*TMath::Abs(deltat)/2.)
	*( (absLf2-1.)*cos(deltam*deltat) + 2.*imLf*sin(deltam*deltat) )
	/( (1.+absLf2)*hplus + 2.*hminus*reLf );


}


// ---------------------------------------------------------------------------
void simuAsymmetry( Double_t &argLf, Double_t &errargLf, Int_t nGen=100, Int_t nBins=20, Double_t argLf_gen=0., Double_t rangefit=10, Double_t timeUncertainty= 1., Double_t mistag = 0.05, Bool_t verbose = kFALSE){
 
  Double_t dilution = 1.-2*mistag;

  if( verbose ) {
    printf( "\nPhysical parameters:\n");
//    printf( " x = %f, y = %f\n", x, y);
    printf( " Delta(Gamma) = %e\n", deltag); 
    printf( " Delta(Mass) = %e\n", deltam); 
    printf( " Module(Lambda_f) generated at %f\n", modLf_gen);
	  
    printf("\nParameter to be estimated:\n"); 
    printf( " Arg(Lambda_f) generated at %f\n", argLf_gen); 
	  
    printf( "\nRange of the fit vs. Delta_t in ps = +/- %f\n", rangefit);
    
    printf("\nDetector imperfection:\n");
    printf( " uncertainty on Delta(t) = %e\n", timeUncertainty); 
    printf( " mistag probability = %f\n", mistag); 
    printf( " dilution factor = %f\n", dilution); 
  }
  
  // Range of measured Delta_t

    Double_t range_deltat[2] = {0.,rangefit};

  // Function describing theoretical asymmetry

  TF1 *fgentheo = new TF1( "fgentheo", asymmetry, range_deltat[0], range_deltat[1], 5);
  fgentheo->SetParameters( deltag, deltam, modLf_gen, argLf_gen, 0.);
  

  // Function used to simulate the asymmetry

  TF1 *fgen = new TF1( "fgen", asymmetry, range_deltat[0], range_deltat[1], 5);
  fgen->SetParameters( deltag, deltam, modLf_gen, argLf_gen, mistag);
  
  // Histograms:
  
  TH1F *hd0 = new TH1F("hd0", "number of D0", nBins, range_deltat[0], range_deltat[1]);
  hd0->SetXTitle("#Deltat (ps)");

  TH1F *hd0bar = new TH1F("hd0bar", "number of anti-D0", nBins, range_deltat[0], range_deltat[1]);
  hd0bar->SetXTitle("#Deltat (ps)");

  // Histogram used to store the simulated asymmetry
  //  without the measurement uncertainty
  
  /*
  TH1F *hsim = new TH1F("hsim", "Simulated asymmetry", nBins, range_deltat[0], range_deltat[1]);
  hsim->SetXTitle("#Deltat (ps)");
  hsim->Sumw2();
  */
  TH1F *htotal = new TH1F("htotal", "Total nb of Ds",  nBins, range_deltat[0], range_deltat[1]);
  htotal->SetXTitle("#Deltat (ps)");
  
 // ======================================
 // Simulation step
 // ======================================

  if( verbose ) printf( "\nSimulation of %d Delta(t)\n", nGen);
  Double_t deltat, A, probBar, weight;
    Int_t countD0 = 0;
    Int_t countD0bar = 0;
  for( int i=0; i<nGen; i++) {
    //deltat = gRandom->Uniform(range_deltat[0], range_deltat[1]);
    // Generate one Delta(t)
    deltat = gRandom->Exp( 1./width);
                 
    // Compute asymmetry
    A = fgen->Eval( deltat);
                 //printf("\n A = %.10f\n",A);
    // Decide if we measure a D0 or a D0bar    
    probBar = (1.+A)/2.;
                 //printf("\n probBar = %.16f { probBar = (1+A) / 2 }\n",probBar);
   Double_t random_decide_D0 = gRandom->Uniform();
               // printf("\n random_decide_D0 = %.16f  {if random > probBar => antiD0 if not D0}\n",random_decide_D0);
   
    
    if( random_decide_D0 > probBar ) {
      weight = 1.;                    //should be -1
               // printf("\nWe have an antiD0\n");
     countD0++;
    }
    else {
      weight = -1.;                  //should be 1
              // printf("\nWe have a D0\n");
      countD0bar++;
    }
    
    Double_t uncertainty_add = gRandom->Gaus(0,timeUncertainty);
                // printf("\n uncertainty_add = %.16f\n",uncertainty_add);         
    // Take into account the time measurement error
    deltat += uncertainty_add;
               // printf("\n New deltat = %.16f\n",deltat);                      
               // printf("\n-----------------------------------\n");
               // printf("\n-----------------------------------\n"); 
    // Fill Histogram
    if( weight>0. ) hd0->Fill(deltat, 1.);                     //should be hd0
    else if( weight<=0.) hd0bar->Fill(deltat, 1.);                //should be hd0bar
    
    
    //I think it makes the opposite, asym = gamma - gammaBar / +
    //hsim->Fill( deltat, weight);
    htotal->Fill( deltat, 1.);
  }
    RooHist* hsim = new RooHist(*hd0bar, *hd0);
    hsim->SetLineColor(2);
    hsim->SetMarkerColor(2);
    hsim->SetName("hsim");
    
  
  printf("\nThe total number of D0 = %d and antiD0 = %d\n",countD0,countD0bar);
  //hsim->Divide( htotal);


 // ======================================
 // Fit step
 // ======================================
  if( verbose ) printf( "\nFit results:\n");
  // Function used to fit the simulated asymmetry
//
  TF1 *ffit = new TF1("ffit", asymmetry, range_deltat[0], range_deltat[1], 5);
// ffit->SetParameters( deltag, deltam, modLf_gen, argLf_gen, mistag);
  ffit->SetParameters( deltag, deltam, modLf_gen, 1.5, mistag);
  ffit->FixParameter( 0, deltag);
  ffit->FixParameter( 1, deltam);
  ffit->FixParameter( 2, modLf_gen);
  ffit->FixParameter( 4, mistag);
  ffit->SetParNames( "#Delta#Gamma", "#DeltaM", "Mod(#lambda_f)", "Arg(#lambda_f)", "mistag");

  hsim->Fit( "ffit","OQR");
//  hsim->Fit( "ffit","OQEM");
  argLf = ffit->GetParameter(3);
  errargLf = ffit->GetParError(3);
//
  
  // ======================================
  // Display step
  // ======================================
  gStyle->SetOptFit(1);
  if( verbose ) {
    TCanvas *c = new TCanvas( "c", "Asymmetry simulation", 1);
    c->Divide(2,2);
    
    c->cd(1);
    fgentheo->SetLineColor(2);
    fgentheo->DrawCopy();
/*    fgen->DrawCopy("same");    
    fgen->SetParameters( deltag, deltam, 1., 0.);
    fgen->DrawCopy("same");
    fgen->SetParameters( deltag, deltam, 1., TMath::Pi()/2.);
    fgen->DrawCopy("same");
    fgen->SetParameters( deltag, deltam, 1., TMath::Pi()/4);
    fgen->DrawCopy("same");  
*/      

    c->cd(2); 
  //  hd0->SetFillColor(2);
    hd0->SetLineColor(2);
    hd0->DrawCopy();
    hd0bar->SetLineColor(3);
    hd0bar->DrawCopy("same");
    gPad->SetLogy();
	  
    c->cd(3);
 //   hsim->SetMinimum( fgentheo->Eval( range_deltat[0] ) );
 //   hsim->SetMaximum( fgentheo->Eval( range_deltat[1] ) );
    //hsim->SetMinimum( -0.02  );
    //hsim->SetMaximum( 0.02);
    hsim->SetLineColor(1);
    //hsim->DrawCopy("e");       //23.07.2014
   hsim->SetMinimum(-2.5);
   hsim->SetMaximum(2.5);
   hsim->Draw("AP");
    //ffit->Draw("same");
    //ffit->FixParameter( 2, 2.);
   
   // fgentheo->SetLineColor(2);
   // fgentheo->DrawCopy("same");
    
    
    c->cd(4);
    htotal->DrawCopy();
    gPad->SetLogy();
    
    
  }
  printf( " Arg(Lambda_f) = %f +/- %f,  relative error = %.3f\n", argLf, errargLf, (argLf-argLf_gen)/argLf_gen);
  
  // ======================================
  // End
  // ======================================
//  hsim->Delete();
  htotal->Delete();
  hd0->Delete();
  hd0bar->Delete();

}

// ---------------------------------------------------------------------------
void simuAsymmetry( Int_t nGen=100, Int_t nBins=20, Double_t argLf_gen=0.0, Double_t rangefit=10., Double_t timeUncertainty= 1., Double_t mistag = 0.05, Bool_t verbose = kTRUE){
  //
  // Created: JB 2012/03/20   / IRB 2012/05/10
  
  Double_t argLf, errargLf;
  simuAsymmetry( argLf, errargLf, nGen, nBins, argLf_gen, rangefit, timeUncertainty, mistag, verbose);
   
}

// ---------------------------------------------------------------------------
void evalAsymmetry( Int_t nExp=100, Int_t nGen=100, Int_t nBins=20, Double_t argLf_gen=0.0, Double_t rangefit=10., Double_t timeUncertainty= 1., Double_t mistag = 0.05){

  
  Double_t argLf, errargLf;
 printf("\n nExp = %d and nGen = %d",nExp,nGen);
 printf("\n argLf_gen = %f and rangefit = %f and timeUncertainty = %.15f \n",argLf_gen,rangefit,timeUncertainty);

// Range of measured Delta_t
    Double_t range_deltat[2] = {0.,rangefit};
	
	
// Function describing theoretical asymmetry	
  TF1 *fgentheor = new TF1( "fgentheor", asymmetry, range_deltat[0], range_deltat[1], 5);
  fgentheor->SetParameters( deltag, deltam, modLf_gen, argLf_gen, 0.);	
	
  Char_t title[200];
  sprintf( title, "Distribution for %d D^0", nGen);//, #sigma_#Deltat=%e, mistag=%.3f", nGen, timeUncertainty, mistag);
 TH1F *hArgLf = new TH1F( "harglf", title, 60, 0.,0.);
  //TH1F *hArgLf = new TH1F( "harglf", title, 60, 0., argLf_gen+10.);
  hArgLf->SetXTitle( "Arg(#lambda_f)");
  
  // Loop on experiments
  for( Int_t iExp=0; iExp<nExp; iExp++) {
    simuAsymmetry( argLf, errargLf, nGen, nBins, argLf_gen, rangefit, timeUncertainty, mistag, kFALSE);
    hArgLf->Fill( argLf);
    cout<<"Experiment number "<<iExp<<endl;
  }
  
  // Display the result: 	
	
  printf("\nSimulating %d experiments\n", nExp);
  printf("Arg(Lambda_f): average = %f, bias = %f, absolute uncertainty = %f, relative uncertainty = %f\n", 
		 hArgLf->GetMean(), hArgLf->GetMean()-argLf_gen, hArgLf->GetRMS(), hArgLf->GetRMS()/hArgLf->GetMean() ); 

  // Function describing the result of the fit:	
  TF1 *fgenfit = new TF1( "fgenfit", asymmetry, range_deltat[0], range_deltat[1], 5);
  fgenfit->SetParameters( deltag, deltam, modLf_gen, hArgLf->GetRMS(), mistag);

  TCanvas *cexp = new TCanvas( "cexp", "Asymmetry smulation over experiment", 1);
  cexp->Divide(2,2);
	
  cexp->cd(1);
  hArgLf->DrawCopy();
	
  cexp->cd(2);
  fgenfit->DrawCopy();
  fgentheor->SetLineColor(2);
  fgentheor->DrawCopy("same");	
  	

// ======================================
// End
// ======================================	
 // hArgLf->Delete();

}
