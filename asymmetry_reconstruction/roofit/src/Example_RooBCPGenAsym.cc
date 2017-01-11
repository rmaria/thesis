#include "Example.hh"
#include <TStyle.h>
#include <TCanvas.h>
#include <TROOT.h>
#include <TH1.h>
#include <TH2.h>
#include <TLegend.h>
#include <TTree.h>
#include <TLeaf.h>
#include <TMath.h>
#include <TLorentzVector.h>
#include <TVector2.h>
#include <TVector3.h>
#include <TRotation.h>
#include <TLorentzRotation.h>
#include <TLine.h>
#include <TArrow.h>
#include <TSystem.h>
#include <TF1.h>
#include <TApplication.h>
#include <TRandom.h>
#include <TLatex.h>

//RooFit
#include <RooRealConstant.h>
#include <RooAbsReal.h>
#include <RooRealVar.h>
#include <RooGaussModel.h>
#include <RooBCPGenDecay.h>   //rmaria commented
//#include "/Users/rmaria/work/programs/learning/roofit/robert_bu/Example/mypdf2/inc/myRooBDecay.h"     //rmaria added
//#include "/Users/rmaria/work/programs/learning/roofit/robert_bu/Example/mypdf/inc/myDecay.h"     //rmaria added
#include <RooFormulaVar.h>
#include <RooCategory.h>
#include <RooAddPdf.h>
#include <RooSimultaneous.h>
#include <RooGlobalFunc.h>
#include <RooPlot.h>
#include <RooFitResult.h>
#include "RooAddition.h"
#include <RooArgSet.h>
//#include <RooBDecay.h>   //rmaria commented
//#include <myRooBDecay.h>     //rmaria added
#include <RooCategory.h>
#include <RooCmdArg.h>
#include <RooDataSet.h>
#include <RooFitResult.h>
#include <RooGaussian.h>
#include <RooGaussModel.h>
#include <RooMinuit.h>
#include <RooNLLVar.h>
#include <RooPlot.h>
#include <RooProdPdf.h>
#include "RooProduct.h"
#include <RooRandom.h>
#include <RooRealVar.h>
#include "RooRealSumPdf.h"
#include <RooRealConstant.h>
#include <RooResolutionModel.h>
#include <RooTruthModel.h>
#include "TH1.h"    //rmaria
#include <RooHist.h> //rmaria

//C++, C
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
using namespace RooFit;

//============================================================================
void Example(double tmean,
	     double DeltaGammaTau,
	     double DeltaM,
	     double Mod_lambda_f,
	     double Phi_lambda_f,
	     double ResModel_Bias,
	     double ResModel_Width,
	     double Omega,
	     double DeltaOmega,
	     int Nevents,
	     int seed,
	     const char* Option,
	     const char* OutputFile)
{

  bool IsSingleSided = true;
  if(TString(Option) == TString("DoubleSided")) IsSingleSided = false;

  Double_t ModLf  = Mod_lambda_f;
  Double_t PhiLf  = Phi_lambda_f*TMath::Pi()/180.0;
 // Double_t tau = tmean;
  //Double_t dm     = DeltaM;


  RooRealVar ModlambdaF("ModlambdaF","ModlambdaF",ModLf,0,10.);
  RooRealVar PhilambdaF("PhilambdaF","PhilambdaF",PhiLf,-TMath::Pi(),TMath::Pi());

  //rmaria mod for myDecay 
 RooAbsReal* avgC = new RooFormulaVar("avgC",
										"avgC",
										"(1.0 - pow(ModlambdaF,2))/(1.0 + pow(ModlambdaF,2))",
										RooArgList(ModlambdaF));
 RooAbsReal* avgS = new RooFormulaVar("avgS",
										"avgS",
										"-2.0*ModlambdaF*sin(PhilambdaF)/(1.0 + pow(ModlambdaF,2))",    
										RooArgList(ModlambdaF,PhilambdaF));
 // end rmaria mod for myDecay

  
  
  //rmaria coeff for myDecay to work
 
 
 RooCategory D0flav("D0flav","D0flav flavour");
 D0flav.defineType("D0",    1);
 D0flav.defineType("D0bar",-1);
 
 RooRealVar bias("bias","title",ResModel_Bias,-10.,10.,"ps"); 
 RooRealVar sigma("sigma","sigma",ResModel_Width,1.e-8,10.,"ps");
 RooRealVar t("t","t",0.,5.,"ps"); 
 RooRealVar tau("tau","tau",tmean,0,10.,"ps");
 RooRealVar dm("dm","dm",DeltaM,0.,10.,"ps^{-1}");
 RooRealVar omega("omega","omega",Omega,0.,0.5);
 RooRealVar deltaomega("deltaomega","deltaomega",DeltaOmega,0.,0.5);
 RooRealVar mu("mu","mu",0.,0.,0.5);
 RooGaussModel model("model","Gaussian Resolution Model",t,bias,sigma);
 
 
 RooBCPGenDecay pdfD0("pdfD0", "title", 
			           t,D0flav,
			           tau, dm,
			           omega, 
			           *avgC, *avgS,
			           deltaomega,
                       mu,
			           model, RooBCPGenDecay::SingleSided) ; 
			           
			           

 
  TRandom r(seed);
  
  
  
  int N_D0    = r.Poisson(Nevents);
  
  
  RooDataSet* data_D0 = pdfD0.generate(RooArgSet(t,D0flav),N_D0) ;
  //RooDataSet* data_D0bar = pdfD0bar.generate(t,500000) ;
  
  cout << data_D0->numEntries() << endl;
    
  //  ---- FIT --------
  RooFitResult* res = pdfD0.fitTo(*data_D0,Save());
  res->Print("v");  

  int nBin = 100; 
  
  
 // -- Define histograms for D0 and D0bar
  
  TH1F *histD0 = new TH1F("histD0","histo",nBin, 0., 5.);
  histD0->Sumw2();
  TH1F *histD0bar = new TH1F("histD0bar","histo",nBin, 0., 5.);
  histD0bar->Sumw2();
  
  //RooArgSet* set;
  
  for(int i=0; i<data_D0->numEntries(); i++)
  	{
    
  	const RooArgSet* set = data_D0->get(i);
    Double_t time = ((RooRealVar*)set->find(t.GetName()))->getVal();
    int flav      = ((RooCategory*)set->find(D0flav.GetName()))->getIndex();
    
    //cout << i+1 << "  " << time << "  " << flav << endl;
    
  	if(flav == 1)       histD0->Fill(time);
  	else if(flav == -1) histD0bar->Fill(time);
  	
  	}
  
  
  //---- The asymmetry -----	
  RooHist* asym2 = new RooHist(*histD0bar, *histD0);
  asym2->SetName("asym2");
  
  
  
  TCanvas *c = new TCanvas("c","c",800,600);
  c->Divide(1,2);
  c->cd(1);
  double Maximum = TMath::Max(histD0->GetMaximum(),
                              histD0bar->GetMaximum());
  histD0->SetMaximum(Maximum*(1.0 + 0.10));
  histD0bar->SetMaximum(Maximum*(1.0 + 0.10));
  histD0->SetLineColor(4);
  histD0bar->SetLineColor(2);
  histD0->Draw();
  histD0bar->Draw("same");
  c->cd(2);
  asym2->Draw("AP");
  c->SaveAs("D0.eps");
 
  			           

  
  return;

}
//============================================================================


