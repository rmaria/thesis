#include "TF1.h"
#include "TH1F.h"
#include "TCanvas.h"
#include "TRandom3.h"
#include "TStyle.h"
#include "TMath.h"
#include "Math/Minimizer.h"
#include "Math/Factory.h"
#include "Math/Functor.h"
#include "TError.h"
#include <iostream>
#include <vector>

// Did not produce the good results :(
// uses the Gamma and Gammabar equations instead of the asymmetry equation, just to compare with the Roofit method



// Global Input parameters

Double_t xx = 0.005; 
Double_t yy = 0.01; 
Double_t width = 1./410.e-15; 
Double_t deltag = yy*2.*width; // PDG average 0.015*2*gamma
Double_t deltam = xx*width; // PDG average 2.39e10
Double_t timeUncertainty= 0.2e-12;
Double_t mistag = 0.0;
Double_t gamma1 = width*(1+yy);               
Double_t modLf_gen =  1.;

vector<double> *measurement_flavour = new vector<double>;
vector<double> *measurement_deltat = new vector<double>; 


// Histograms:

 TCanvas *cannevas = new TCanvas( "cannevas", "asymmetry measurement", 1);

 TH1F *hd0 = new TH1F("hd0", "number of D0", 30,-10.e-12, 10.e-12);
 TH1F *hd0bar = new TH1F("hd0bar", "number of anti-D0", 30,-10.e-12, 10.e-12);
 TH1F *htotal = new TH1F("htotal", "Total nb of Ds", 30,-10.e-12, 10.e-12);
 TH1F *hsim = new TH1F("hsim", "Data asymmetry", 30,-10.e-12, 10.e-12);




// ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

Double_t asymmetry( Double_t *x, Double_t *par) {
 

  Double_t deltat = x[0];
  
  Double_t deltag = par[0];
  Double_t deltam = par[1];  
  Double_t modLf = par[2];
  Double_t argLf = par[3];
  Double_t mistag = par[4];
  
  Double_t absLf2 = modLf*modLf;
  Double_t reLf = modLf *cos(argLf*TMath::DegToRad());
  Double_t imLf = modLf *sin(argLf*TMath::DegToRad());

  Double_t hplus  = 1. + exp(deltag*TMath::Abs(deltat));     // in the article no "abs"????
  Double_t hminus = 1. - exp(deltag*TMath::Abs(deltat));
	
  Double_t dilution = 1.-2*mistag; 

  return 2.*dilution*exp(deltag*TMath::Abs(deltat)/2.)
	*( (absLf2-1.)*cos(deltam*deltat) + 2.*imLf*sin(deltam*deltat) )
	/( (1.+absLf2)*hplus + 2.*hminus*reLf );
}

// ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
double pdfExponential( double *x, double *par )

//    The PDF = (1-w)*Gamma+w*GammaBar
//              where Gamma    = the function for time evolution of a D0 decaying to a final state f (correlated case)
//                    GammaBar = the function for time evolution of a D0bar decaying to a final state f  

{
  Double_t deltat = x[0];

  Double_t deltag = par[0];
  Double_t deltam = par[1];  
  Double_t modLf = par[2];
  Double_t argLf = par[3];
  Double_t mistag = par[4];
 
  Double_t hplus  = 1. + exp(deltag*TMath::Abs(deltat));
  Double_t hminus = 1. - exp(deltag*TMath::Abs(deltat));
  Double_t absLf2 = modLf*modLf;
  Double_t reLf = modLf *cos(argLf*TMath::DegToRad());
  Double_t imLf = modLf *sin(argLf*TMath::DegToRad());	
  
return (1-mistag) * ( exp(-gamma1*TMath::Abs(deltat))*( hplus/2 + reLf*hminus/(1+absLf2)
        + exp(deltag*TMath::Abs(deltat)/2.)*( (1-absLf2)*cos(deltam*deltat)/(1+absLf2) - (2*imLf)*sin(deltam*deltat)/(1+absLf2))) )
        + mistag * ( exp(-gamma1*TMath::Abs(deltat))*(hplus/2 + reLf*hminus/(1+absLf2)
        + exp(deltag*TMath::Abs(deltat)/2.)*( -(1-absLf2)*cos(deltam*deltat)/(1+absLf2) + (2*imLf)*sin(deltam*deltat)/(1+absLf2))));
  
}

// ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
double pdfExponentialBar( double *x, double *par )

//        The Pdf = (1-w)*GammaBar+w*Gamma

{
  Double_t deltat = x[0];
  
  Double_t deltag = par[0];
  Double_t deltam = par[1];  
  Double_t modLf = par[2];
  Double_t argLf = par[3];
  Double_t mistag = par[4]; 

  Double_t hplus  = 1. + exp(deltag*TMath::Abs(deltat));     
  Double_t hminus = 1. - exp(deltag*TMath::Abs(deltat));
  Double_t absLf2 = modLf*modLf;
  Double_t reLf = modLf *cos(argLf*TMath::DegToRad());
  Double_t imLf = modLf *sin(argLf*TMath::DegToRad());	
  
	return mistag * ( exp(-gamma1*TMath::Abs(deltat))*( hplus/2 + reLf*hminus/(1+absLf2)
														   + exp(deltag*TMath::Abs(deltat)/2.)*( (1-absLf2)*cos(deltam*deltat)/(1+absLf2) - (2*imLf)*sin(deltam*deltat)/(1+absLf2))) )
	+ (1-mistag) * ( exp(-gamma1*TMath::Abs(deltat))*(hplus/2 + reLf*hminus/(1+absLf2)
												  + exp(deltag*TMath::Abs(deltat)/2.)*( -(1-absLf2)*cos(deltam*deltat)/(1+absLf2) + (2*imLf)*sin(deltam*deltat)/(1+absLf2))));
	
	
  
}

// ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
double likelyhood(const double *param)

//    Likelyhood function, a sum over -log from all the pdfExponential or pdfExponentialBar, the choose depending 
//        whether we have generated a D0 or a Dobar  

{
  Double_t par[5];
  par[0] = deltag;                  //deltag
  par[1] = deltam;                      //deltam
  par[2] = 1. ;                           //modLf
  par[3] = param[0];                     //THE FREE PARAMETER
  par[4] = mistag;                         //mistag
  Double_t like = 0.;
  for(int i=0; i<(int)measurement_flavour->size(); i++)
     {
      if(measurement_flavour->at(i) > 0.)
         like += -log( pdfExponential( &(measurement_deltat->at(i)), par ) );    
      else if(measurement_flavour->at(i) < 0.)                                                         
         like += -log(pdfExponentialBar( &(measurement_deltat -> at(i)), par));  
     }
 return like;
}

// ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
void generateMeasurements(Int_t nGen, Double_t argLf_gen)
{
	
	
// Generate the measurements
  
  //Double_t dilution = 1.-2*mistag;

	Double_t deltat=0., A=0., probBar=0., weight=1.;
	Int_t countd0=0, countd0bar=0 ;
	
	
  for( int i=0; i<nGen; i++) 
  {
// generate a deltat:
//    deltat = gRandom->Exp( 1./width);
	  deltat =  (-1./width)*log( gRandom->Uniform());
// decide if D0signal has decayed before (then deltat<0) or after (then deltat>0) D0tag:  	  
	if( gRandom->Uniform()>.5) deltat *= -1;
	  
// Compute asymmetry
//	in order to know how many D0 vs. D0bar have to be generated, 
//	taken into account the possible mistake due to mistag:
// (1+A)/2 = D0/(D0+D0bar)	  
	TF1 *fgen = new TF1( "fgen", asymmetry, -10.e-12,10.e-12, 5);	
	fgen->SetParameters( deltag, deltam, modLf_gen, argLf_gen, mistag);	  
    A = fgen->Eval( deltat);
// decide whether it is a D0 ar a D0bar, taken asymmetry into account:  
    probBar = (1.+A)/2.;
    if( gRandom->Uniform() > probBar ) {
      weight = -1.;
	  countd0bar++;
    }
    else {
      weight = 1.;
	  countd0++;
    }
	  
	  
    // Take into account the time measurement error
    deltat += gRandom->Gaus(0,timeUncertainty);
	  
// Store simulated measurements:
	measurement_deltat->push_back(deltat); 
	measurement_flavour->push_back(weight);  

// Fill Histogram
    if( weight>0. ) hd0->Fill(deltat, 1.);
	else if( weight<=0.) hd0bar->Fill(deltat, 1.);
	hsim->Fill( deltat, weight);
	htotal->Fill( deltat, 1.);

  }
	
	hsim->Divide( htotal);

	
	cout << measurement_flavour->size() << " measurements were generated with arg_lambda_f = " << argLf_gen << endl;
	cout << countd0 << " D0 were generated taken mistag into account, and "<< countd0bar<<" anti-D0"<< endl;

}
//------------------------------------------------------------------------------
void drawResult(double argLf_gen, double arg_fit) {
	
	cannevas->Divide(2,2);
	
	cannevas->cd(1);
// Generated D0 and D0bar as a function of time:
	hd0->SetXTitle("#Deltat (ps)");
	hd0->SetLineColor(2);
    hd0->DrawCopy();
	hd0bar->SetXTitle("#Deltat (ps)");
    hd0bar->SetLineColor(3);
    hd0bar->DrawCopy("same");
    gPad->SetLogy();
	
	cannevas->cd(2);
// Generated data asymmetry:	
	hsim->SetXTitle("#Deltat (ps)");
	hsim->Sumw2();
	hsim->SetAxisRange(-0.1,0.1,"Y");
	hsim->DrawCopy("e");
	
	cannevas->cd(3);
// Function describing the theoretical generated asymmetry:	
	TF1 *fgentheo = new TF1( "fgentheo", asymmetry, -10.e-12, 10.e-12, 5);
	fgentheo->SetParameters( deltag, deltam, modLf_gen, argLf_gen, 0.);
	fgentheo->SetLineWidth(.5);
	
	TF1 *fcopy;
		
	fgentheo->SetLineColor(1);    
	fcopy = fgentheo->DrawCopy();
	fcopy->GetHistogram()->SetXTitle("#Delta t");
	fcopy->GetHistogram()->SetYTitle("Asymmetry");
	fcopy->GetHistogram()->SetTitle("Theoretical (generated) asymmetry");

	
	
}


// ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
Int_t NumericalMinimization(Int_t nGen=1000, Double_t argLf_gen=6)
{
   
  const char * minName = "Minuit";
  const char *algoName = "";
	
  ROOT::Math::Minimizer* min = 
  ROOT::Math::Factory::CreateMinimizer(minName, algoName);
  
  // set tolerance , etc...
  min->SetStrategy(2); // for Minuit
  min->SetMaxFunctionCalls(1000000); // for Minuit/Minuit2 
  min->SetMaxIterations(10000);  // for GSL 
  min->SetTolerance(0.0000001);
  min->SetPrintLevel(1);
  
  // create function wrapper for minmizer
  // a IMultiGenFunction type 
  ROOT::Math::Functor flike(&likelyhood,1);    
  double step[1] = {0.000001};	
  // starting point
  
  cout << "Generating" << endl;
  generateMeasurements(nGen,argLf_gen);
  


  min->SetFunction(flike);
  
  // Set the free variables to be minimized!
  min->SetLimitedVariable(0,"Arg_Lf",5., step[0],5.,48.);    
  
  // do the minimization and get indicators on fit quality
  bool convergence = min->Minimize(); 
  
  int convergenceStatus = min->Status();
  double distToMin = min->Edm();
  bool errorValid = min->IsValidError();
  cout << "Fit has converge: " << convergence << " with status: " << convergenceStatus << endl;
  cout << "Error computation is valid: " << errorValid << endl;
  cout << "Expected distance to minimum: " << distToMin << endl;
  
  // Perform better uncertainty computation
  min->Hesse();
    
  // Store result with uncertainties
  const double *estimation = min->X();
  const double *uncertainty = min->Errors();
  cout << "Estimated param: " << estimation[0] << " +/- " 
  << uncertainty[0] << std::endl;
  std::cout << "Minimum of likelyhood: " << min->MinValue() << std::endl;
  
  // Display step:
 
//   if( convergence ) 
 //  {
    drawResult(argLf_gen, estimation[0]);                                            
//   }   
 
  // Final message
  if ( fabs(estimation[0]-argLf_gen)  < .1  )                                         //instead of trueTau
    std::cout << "Minimizer " << minName << " - " << algoName 
    << "   converged to the right minimum" << std::endl;
  else {
    std::cout << "Minimizer " << minName << " - " << algoName 
    << "   failed to converge !!!" << std::endl;
    Error("NumericalMinimization","fail to converge");
  }
  
  return convergenceStatus;
}




