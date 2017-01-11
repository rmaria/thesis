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
#include <TLegend.h>


//RooFit
#include <RooRealConstant.h>
#include <RooAbsReal.h>
#include <RooRealVar.h>
#include <RooGaussModel.h>
#include <RooBCPGenDecay.h>   //rmaria commented
//#include "/Users/rmaria/work/programs/learning/roofit/robert_bu/Example/mypdf2/inc/myRooBDecay.h"     //rmaria added
#include "/Users/rmaria/work/programs/roofit/learning_again/test22.06.15/mypdf/inc/myDecay.h"     //rmaria added
#include <RooFormulaVar.h>
#include <RooCategory.h>
#include <RooAddPdf.h>
#include <RooSimultaneous.h>
#include <RooGlobalFunc.h>
#include <RooPlot.h>
#include <RooFitResult.h>
#include <RooAddition.h>     //23.07.2014
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
#include <RooRealSumPdf.h>      //23.07.2014
#include <RooRealConstant.h>
#include <RooResolutionModel.h>
#include <RooTruthModel.h>
#include "TH1.h"    //rmaria
#include <RooHist.h> //rmaria
#include <RooMCStudy.h>

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


//time us rmaria 25.06.15
#include <sys/time.h>
#include <unistd.h>

using namespace std;
using namespace RooFit;



int nrD0, nrD0bar;
Double_t modLf_gen =  1.;

Double_t asymmetry( Double_t *x, Double_t *par) {

  Double_t t = x[0];
  
  Double_t deltag = par[0];
  Double_t deltam = par[1];  
  Double_t modLf = par[2];
  Double_t argLf = par[3];
  Double_t misstag = par[4];
  
  Double_t absLf2 = modLf*modLf;
  Double_t reLf = modLf *cos(argLf*TMath::Pi()/180.);
  Double_t imLf = modLf *sin(argLf*TMath::Pi()/180.);

  Double_t hplus  = 1. + exp(deltag*TMath::Abs(t));
  Double_t hminus = 1. - exp(deltag*TMath::Abs(t));
  Double_t dilution = 1.-2*misstag;

  return 2.*dilution*exp(deltag*TMath::Abs(t)/2.)
	*( (absLf2-1.)*cos(deltam*t) + 2.*imLf*sin(deltam*t) )
	/( (1.+absLf2)*hplus + 2.*hminus*reLf );


}


//============================================================================
void Example(
         double &argLfFit,
	     double &argLfFitError,
	     double &minusLog,
	     double &edm,
	     double &stat,
         double tmean,
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
	     const char* OutputFile
	     )
{


  //seed, rmaria 25.06.15
  //struct timeval end;                                 // usecond does not work
  
  //gettimeofday(&end, NULL);				// usecond does not work
  
  //srand (end.tv_usec); 				//usecond does not work	

  //gettimeofday(&end, NULL);                               //usecond does not work
  //cout << "------> rand band = " << end.tv_usec << endl; // usecond does not work
  

  //print verification information
  
  cout << "---------- Verification informations ------- " << endl;
  cout << "Seed in the Example.cc is = " << seed << endl;
  cout << "tmean verif = " << tmean << endl;
  cout << "DeltaGammaTau = " << DeltaGammaTau << endl;
  cout << "DeltaM verif = " << DeltaM << endl;
  cout << "Mod_lambda_f verif = " << Mod_lambda_f << endl;
  cout << "Phi_lambda_f = " << Phi_lambda_f << endl;
  

  bool IsSingleSided = true;
  if(TString(Option) == TString("DoubleSided")) 
  	{
	cout << "???????????????????" << endl;
	cout << "???????????????????" << endl;
	IsSingleSided = false;
	}
  else
  	{
	cout << "!!!!!!!!!!!!!!!!!!!!" << endl;
	cout << "!!!!!!!!!!!!!!!!!!!!" << endl;
	}		
  Double_t ModLf  = Mod_lambda_f;
  Double_t PhiLf  = Phi_lambda_f*TMath::Pi()/180.0;
  Double_t dg = DeltaGammaTau/tmean;
 // Double_t tau = tmean;
  //Double_t dm     = DeltaM;


  RooRealVar ModlambdaF("ModlambdaF","ModlambdaF",ModLf);   
  RooRealVar PhilambdaF("PhilambdaF","PhilambdaF",PhiLf,-TMath::Pi(),TMath::Pi());

  //rmaria mod for myDecay
  
 RooAbsReal* _avgC = new RooFormulaVar("avgC",
										"avgC",
										"(1.0 - pow(ModlambdaF,2))/(1.0 + pow(ModlambdaF,2))",
										RooArgList(ModlambdaF));
 RooAbsReal* _avgS = new RooFormulaVar("avgS",
										"avgS",
										"2.0*ModlambdaF*sin(PhilambdaF)/(1.0 + pow(ModlambdaF,2))",    
										RooArgList(ModlambdaF,PhilambdaF));

 RooAbsReal* _avgCh = new RooFormulaVar("avgCh",
										"avgCh",
										"1.0",
										RooArgList(ModlambdaF));
 RooAbsReal* _avgSh = new RooFormulaVar("_avgSh",
										"_avgSh",
										"-2.0* ModlambdaF*cos(PhilambdaF)/(1.0 + pow(ModlambdaF,2))",    
										RooArgList(ModlambdaF,PhilambdaF));
 // end rmaria mod for myDecay

 
  
  //rmaria coeff for myDecay to work
 
 
 RooCategory D0flav("D0flav","D0flav flavour");
 D0flav.defineType("D0",1);
 D0flav.defineType("D0bar",-1);
 
 RooRealVar bias("bias","title",ResModel_Bias,"ps"); 
 RooRealVar sigma("sigma","sigma",ResModel_Width,"ps");
 RooRealVar t("t","t",0.,20.,"ps"); 
 RooRealVar tau("tau","tau",tmean,"ps");
 RooRealVar dm("dm","dm",DeltaM,"ps^{-1}");
 RooRealVar dGamma("dGamma","dGamma",dg,"ps^{-1}");
 RooRealVar omega("omega","omega",Omega);
 RooRealVar deltaomega("deltaomega","deltaomega",DeltaOmega);
 RooRealVar mu("mu","mu",0.);
 RooGaussModel model("model","Gaussian Resolution Model",t,bias,sigma);
 
 
  /*
   RooBCPGenDecay pdfD0("pdfD0", "title", 
			           t,D0flav,
			           tau, dm,
			           omega, 
			           *_avgC, *_avgS,
			           deltaomega,
                       mu,
			           model, RooBCPGenDecay::SingleSided) ;
 */
  
 myDecay pdfD0_my("pdfD0_my", "title", 
			           t,D0flav,
			           tau, dm,
			           dGamma,
			           omega, 
			           *_avgC, *_avgS,
			           *_avgCh, *_avgSh,
			           deltaomega,
                       mu,
			           model, myDecay::SingleSided) ;			           

 
  TRandom r(seed);
 
  
  int N_D0    = r.Poisson(Nevents);
 //int N_D0 = Nevents;
 
 
  
  
  RooDataSet* data_D0_my = pdfD0_my.generate(RooArgSet(t,D0flav),N_D0) ;  
  cout << "My Data D0 = " << data_D0_my->numEntries() << endl;
  
  
 // RooDataSet* data_D0 = pdfD0.generate(RooArgSet(t,D0flav),N_D0) ;  
 // cout << "Data D0    = " << data_D0->numEntries() << endl;
   
  
   
   
    
  //----Fit here ----
  RooFitResult* res = pdfD0_my.fitTo(*data_D0_my,Save(),Range(0.,5.));
  res->Print("v");
  
  cout << "<><><><><><><><><><><><>" << endl;
   
  RooRealVar* par1_fitresult = (RooRealVar*)res->floatParsFinal().find("PhilambdaF");
  
  //try other param (even if not free in the fit for the moment)
  // RooRealVar* par2_fitresult = (RooRealVar*)res->floatParsFinal().find("ModlambdaF");
  
  
  //declared here
 // double modLfFit;
  
  argLfFit = par1_fitresult->getVal();
  //modLfFit = par2_fitresult->getVal();
  argLfFitError = par1_fitresult->getError();
   edm = res->edm();
   stat = res->status();
   minusLog = res->minNll();
  const TMatrixDSym& cor = res->correlationMatrix() ;
  const TMatrixDSym& cov = res->covarianceMatrix() ;
  //Double_t correlation(const char* argLfFit, const char* modLfFit) const;
  
  /*
  if(stat == 0)
  	{
  ofstream myfile;
  myfile.open ("no_error_status.txt");
  myfile << "The Fit has converged" << endl;	
  myfile << "edm = " << edm << endl;
  myfile << "status = " << stat << endl;
  myfile << "-log(L) at minimum = " << minusLog << endl ;    
  myfile << "argLf = " << argLfFit << endl;
  myfile << "correlations between parameters"<<endl;
  //myfile << res->correlation("argLfFit","modLfFit") << endl;
  myfile << "<><><><><><><><><><><><>" << endl;
  myfile.close();
 	}
 
 if(stat != 0)
 	{
 	 
 	ofstream myfile;
  	myfile.open ("error_status.txt");
  	myfile << "Error status != OK\n";
  	myfile << "seed = "<<seed<<endl;
  	myfile << "argLf = " << argLfFit << endl;
  	myfile << "error =" << argLfFitError << endl;
 	myfile << "-log(L) at minimum = " << minusLog << endl;
 	myfile << "edm = " << edm << endl;
  	myfile.close();
 	
 	} 
 */
 /*
 //Let's write a TTree
  TTree* tree = new TTree("tree","tree") ;
  Double_t* result_argLF = new Double_t ;
  Double_t* error_argLf = new Double_t ;
  tree->Branch("result",result_argLF,"result/D") ;
  tree->Branch("error",error_argLf,"error/D") ;
  
  *result_argLF = argLfFit ;
  *error_argLf = argLfFitError ;
  
  tree->Fill() ;
  tree->SaveAs("AltNtuple.root");
  */
  
  
 /*
 // ------------------------------------------ OLD METHOD ----------------------------------//
 // -----------------------------------------------------------------------------------------//
 
 Double_t argLf, errargLf;
 Double_t argLf2, errargLf2;
 
 //Variables used
 Double_t misstag = Omega;
 Double_t rangefit = 5.;
 Double_t deltag = dGamma.getVal();
 Double_t deltam = dm.getVal();
 Double_t modLf_gen = ModlambdaF.getVal();
 Double_t argLf_gen = Phi_lambda_f;
 Double_t nGen = N_D0;
 Double_t nBins = 50;
 Double_t timeUncertainty = ResModel_Width;
 Double_t width = 1./tmean;
 
 Double_t dilution = 1.-2*misstag;

  
  
  // Range of measured Delta_t

    Double_t range_deltat[2] = {0., rangefit};

  // Function describing theoretical asymmetry

  TF1 *fgentheo = new TF1( "fgentheo", asymmetry, range_deltat[0], range_deltat[1], 5);
  fgentheo->SetParameters( deltag, deltam, modLf_gen, argLf_gen, 0.);
  

  // Function used to simulate the asymmetry

  TF1 *fgen = new TF1( "fgen", asymmetry, range_deltat[0], range_deltat[1], 5);
  fgen->SetParameters( deltag, deltam, modLf_gen, argLf_gen, misstag);
  
  // Histograms:
  
  TH1F *hd0 = new TH1F("hd0", "number of D0", nBins, range_deltat[0], range_deltat[1]);
  hd0->SetXTitle("#Deltat (ps)");

  TH1F *hd0bar = new TH1F("hd0bar", "number of anti-D0", nBins, range_deltat[0], range_deltat[1]);
  hd0bar->SetXTitle("#Deltat (ps)");

  // Histogram used to store the simulated asymmetry
  //  without the measurement uncertainty


  
  TH1F *hsim = new TH1F("hsim", "Simulated asymmetry", nBins, range_deltat[0], range_deltat[1]);
  hsim->SetXTitle("#Deltat (ps)");
  hsim->SetLineColor(1);
  hsim->Sumw2();
  
  

  
  TH1F *htotal = new TH1F("htotal", "Total nb of Ds",  nBins, range_deltat[0], range_deltat[1]);
  htotal->SetXTitle("#Deltat (ps)");

 // ======================================
 // Simulation step
 // ======================================

  
  Double_t deltat, A, probBar, weight;
    Int_t countD0 = 0;
    Int_t countD0bar = 0;
    nrD0 =0;
    nrD0bar = 0;
  for( int i=0; i<nGen; i++) {
    //deltat = gRandom->Uniform(range_deltat[0], range_deltat[1]);
    // Generate one Delta(t)
    deltat = gRandom->Exp( tau.getVal());
              
    A = fgen->Eval( deltat);
    
    
    probBar = (1.+A)/2.;
               
   Double_t random_decide_D0 = gRandom->Uniform();
             
    
    if( random_decide_D0 > probBar ) {
      weight = 1.; //should be -1
              
     countD0++;
     nrD0++;
    }
    else {
      weight = -1.;         //should be 1
             
      countD0bar++;
      nrD0bar++;
    }
    
    Double_t uncertainty_add = gRandom->Gaus(0,timeUncertainty);
               
    deltat += uncertainty_add;
               
    if( weight>0. ) hd0->Fill(deltat, 1.);            //should be hd0
    else if( weight<=0.) hd0bar->Fill(deltat, 1.);        //should be hd0bar
    hsim->Fill( deltat, -weight);
    htotal->Fill( deltat, 1.);
  }
  printf("\nThe total number of D0 = %d and antiD0 = %d\n",countD0,countD0bar);
  hsim->Divide( htotal);



  RooHist* hsim2 = new RooHist(*hd0bar, *hd0);
  hsim2->SetLineColor(4);
  hsim2->SetMarkerColor(4);
  hsim2->SetName("hsim2");
 

 

  TF1 *ffit = new TF1("ffit", asymmetry, range_deltat[0], range_deltat[1], 5);
// ffit->SetParameters( deltag, deltam, modLf_gen, argLf_gen, misstag);
  ffit->SetParameters( deltag, deltam, modLf_gen, 1.5, misstag);
  ffit->FixParameter( 0, deltag);
  ffit->FixParameter( 1, deltam);
  ffit->FixParameter( 2, modLf_gen);
  ffit->SetParLimits(3,-180.,+180.);   // if the fit fails doesn't help
  ffit->FixParameter( 4, misstag);
  ffit->SetParNames( "#Delta#Gamma", "#DeltaM", "Mod(#lambda_f)", "Arg(#lambda_f)", "misstag");
  ffit->SetLineColor(kGreen);

 TF1 *ffit2 = new TF1("ffit2", asymmetry, range_deltat[0], range_deltat[1], 5);
// ffit->SetParameters( deltag, deltam, modLf_gen, argLf_gen, misstag);
  ffit2->SetParameters( deltag, deltam, modLf_gen, 1.5, misstag);
  ffit2->FixParameter( 0, deltag);
  ffit2->FixParameter( 1, deltam);
  ffit2->FixParameter( 2, modLf_gen);
  ffit2->SetParLimits(3,-180.,+180.);      // if the fit fails doesn't help
  ffit2->FixParameter( 4, misstag);
  ffit2->SetParNames( "#Delta#Gamma", "#DeltaM", "Mod(#lambda_f)", "Arg(#lambda_f)", "misstag");
  ffit2->SetLineColor(kBlue);
 
 
  hsim->Fit( ffit,"R");
//  hsim->Fit( "ffit","OQEM");
  argLf = ffit->GetParameter(3);
  errargLf = ffit->GetParError(3);
   cout << "-------------------------------------------->>" << endl;
   cout << "---->>ArgLf = " << argLf << " +- " << errargLf << endl;
   
   
 hsim2->Fit( ffit2,"R");
 argLf2 = ffit2->GetParameter(3);
 errargLf2 = ffit2->GetParError(3);
   cout << "--------------------------------------------->>" << endl;
   cout << "----->>ArgLf SECOND = " << argLf2 << " +- " << errargLf2 << endl; 
  
  
  
 // -- Define histograms for D0 and D0bar
  
  TH1F *histD0_my = new TH1F("histD0_my","histo",nBins, range_deltat[0],range_deltat[1]);
  histD0_my->SetLineColor(4);
  histD0_my->Sumw2();
  TH1F *histD0bar_my = new TH1F("histD0bar_my","histo",nBins, range_deltat[0],range_deltat[1]);
  histD0bar_my->SetLineColor(4);
  histD0bar_my->Sumw2();
  
  //RooArgSet* set;
  
  int nD0 = 0;
  int nD0bar = 0;
  
  
  /*
  for(int i=0; i<data_D0_my->numEntries(); i++)
  	{
    
  	const RooArgSet* set = data_D0_my->get(i);
    Double_t time = ((RooRealVar*)set->find(t.GetName()))->getVal();
    int flav      = ((RooCategory*)set->find(D0flav.GetName()))->getIndex();
    TString lab      = ((RooCategory*)set->find(D0flav.GetName()))->getLabel();
    //cout << i+1 << "  " << time << "  " << flav << endl;
    
  	if(lab == TString("D0"))
  	     {
  	       nD0++;
  	       histD0_my->Fill(time);
  	     }
  	else if(lab == TString("D0bar")) 
  	      {
  	      nD0bar++;
  	      histD0bar_my->Fill(time);
  	      }
  	
  	}
  
  cout << "ROOFIT nb D0 = " << nD0 << endl;
  cout << "ROOFIT nb D0bar = " << nD0bar << endl;
  
  
  
  //---- The asymmetry -----	
  RooHist* asym2_my = new RooHist(*histD0bar_my, *histD0_my);
  asym2_my->SetLineColor(4);
  asym2_my->SetMarkerColor(4);
  asym2_my->SetName("asym2_my");
  
  */
  
  /*
  TH1F *histD0 = new TH1F("histD0","histo",nBin, 0., 5.);
  histD0->SetLineColor(2);
  histD0->Sumw2();
  TH1F *histD0bar = new TH1F("histD0bar","histo",nBin, 0., 5.);
  histD0bar->SetLineColor(2);
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
  asym2->SetLineColor(2);
  asym2->SetMarkerColor(2);
  asym2->SetName("asym2");
  
  
  //old asymmetry function
  Double_t reLf = ModLf *cos(PhiLf);
  Double_t imLf = ModLf *sin(PhiLf);
  Double_t absLf2 = ModLf * ModLf;
  //Double_t dilution = 1.;
  
  
  
  TF1* func2 = new TF1("func2","2.*[0]*exp( [1] * x / 2. )*( ( [2]-1. )*cos([3] * x) + 2.*sin([4]*TMath::Pi()/180.)*sin([3] * x) )/( (1.+[2])* (1. + exp([1] * x) ) + 2.* (1. - exp([1] * x) ) * cos([4]*TMath::Pi()/180.) )  ",0.,5.);
  func2->SetLineColor(kBlue);
  func2->SetParameter(0,dilution);
  func2->SetParameter(1,dGamma.getVal());
  func2->SetParameter(2,absLf2);
  func2->SetParameter(3,dm.getVal());
  func2->SetParameter(4,Phi_lambda_f);
  
  
  TF1* func3 = new TF1("func3","( sin( [0]*TMath::Pi()/180 )*sin([1]*x) )/( cosh([2]*x/2) - cos( [0]*TMath::Pi()/180) * sinh([2]*x/2) )");
  func3->SetLineColor(2);
  func3->SetParameter(0,Phi_lambda_f);
  func3->SetParameter(1,dm.getVal());
  func3->SetParameter(2,dGamma.getVal());
  
  
  //TF1* func2 = new TF1("func2", asymmetry, range_t[0], range_t[1], 5);
  //func2->SetParameters( dGamma.getVal(), dm.getVal(), ModLf, Phi_lambda_f, 0);
  
  
  
  
  TCanvas *c = new TCanvas("c","c",800,600);
  c->Divide(2,2);
  c->cd(1);
  
  hd0->SetLineColor(2);
  hd0->Draw();
  histD0_my->SetLineColor(4);
  histD0_my->Draw("same");
 // histD0->Draw("same");
 
 TLegend *legend = new TLegend(0.5,0.7,0.70,0.9);
 legend->SetHeader("D0 comparation");
 legend->AddEntry("hd0","red Old Histo","l");
 legend->AddEntry("histD0_my","blue Roofit Histo","l");
 legend->Draw();
  
  c->cd(2);
  hd0bar->SetLineColor(2);
  hd0bar->Draw();
  histD0bar_my->SetLineColor(4);
  histD0bar_my->Draw("same");
  TLegend *legend1 = new TLegend(0.5,0.7,0.70,0.9);
  legend1->SetHeader("antiD0 comparation");
  legend1->AddEntry("hd0bar","red Old Histo","l");
  legend1->AddEntry("histD0bar_my","blue Roofit Histo","l");
  legend1->Draw();
  
  c->cd(3);
  
  asym2_my->SetMinimum(-2.5);
  asym2_my->SetMaximum(2.5);
  asym2_my->Draw("AP");
  func2->SetLineColor(1);
  func2->Draw("same");
  func3->Draw("same");
  TLegend *legend2 = new TLegend(0.3,0.7,0.7,0.9);
  legend1->SetHeader("Asym");
  legend2->AddEntry("asym2_my","Roofit data","l");
  legend2->AddEntry("func2","Th. Asym","l");
  legend2->SetEntrySeparation(0.0025);
  //legend2->SetTextFont(82);
  //legend2->SetTextSize(0.05);
  legend2->Draw();
  
  c->cd(4);
  hsim->SetMinimum(-2.5);
  hsim->SetMaximum(2.5);
  hsim->Draw("AP");
  TLegend *legend3 = new TLegend(0.3,0.7,0.70,0.9);
  legend3->SetHeader("Data and #chi^2 fit");
  legend3->AddEntry("hsim","Old method","l");
  
  legend3->Draw();
  
  //asym2->SetLineColor(4);
  //asym2->Draw("P");
  c->SaveAs("D0.root");          
  
  TCanvas *cc = new TCanvas("cc","cc",800,600);
  cc->Divide(2,2);
  cc->cd(1);
  asym2_my->SetLineColor(2);
  asym2_my->Draw("AP");
  hsim->SetLineColor(4);
  hsim->Draw("e same");
  cc->SaveAs("comp_asym.root");
  */
  
  delete data_D0_my;
  
  //return;

}
//============================================================================

