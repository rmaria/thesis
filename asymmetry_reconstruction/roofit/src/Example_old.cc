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
//#include <RooBDecay.h>   //rmaria commented
//#include "/Users/rmaria/work/programs/learning/roofit/robert_bu/Example/mypdf2/inc/myRooBDecay.h"     //rmaria added
#include "/Users/rmaria/work/programs/learning/roofit/robert_bu/Example/mypdf/inc/myDecay.h"     //rmaria added
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
	     double GenAsym,
	     int seed,
	     const char* Option,
	     const char* OutputFile)
{

  bool IsSingleSided = true;
  if(TString(Option) == TString("DoubleSided")) IsSingleSided = false;

  myDecay::DecayType Opt_num = myDecay::SingleSided;

  //The Time variable:
  double R_t[2];
  double N_tau = 10.0;
  if(IsSingleSided) {
    Opt_num = myDecay::SingleSided;
    R_t[0] = 0.0;
    //R_t[1] = tmean*N_tau;
    R_t[1] = 5.;
  }
  else {
    Opt_num = myDecay::DoubleSided;
    R_t[0] = -tmean*N_tau;
    R_t[1] =  tmean*N_tau;
  }
  RooRealVar t("t","t",R_t[0],R_t[1],"ps");
  if(!IsSingleSided) t.SetTitle("#Deltat");

  //Resolution Model parameters:
  Double_t TheBias = ResModel_Bias;
  Double_t TheWidth = ResModel_Width;
  RooRealVar bias("bias",  "Resolution model Bias",TheBias,-10.0,10.0,"ps");
  RooRealVar sigma("sigma","Resolution model #sigma",TheWidth,0.0001,100.0,"ps");
  RooGaussModel ResModel("ResModel","Gaussian Resolution Model",t,bias,sigma);

  //Truth PDF parameters:
  Double_t Tau = tmean;
  Double_t dgamma = DeltaGammaTau/tmean;
  Double_t DM     = DeltaM;
  Double_t ModLf  = Mod_lambda_f;
  Double_t PhiLf  = Phi_lambda_f*TMath::Pi()/180.0;

  RooRealVar tau("tau","Average Decay time",Tau,0.0,500.0,"ps");
  RooRealVar DGamma("DGamma","#Delta#Gamma",dgamma,-10.0,10.0,"ps^{-1}");     //it was -100.0,100.0
  RooRealVar dm("dm","#DeltaM",DM,0.0,10.0,"ps^{-1}");                       //it was 0.0, 100.0 
  RooRealVar ModlambdaF("ModlambdaF","#lambda_{f}",ModLf,0.0,10.0,"");     //it was 0.0, 1000.0
  RooRealVar PhilambdaF("PhilambdaF","#lambda_{f}",PhiLf,-TMath::Pi(),TMath::Pi(),"rad"); //it was -180.0,180.0

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

  //D0->f pdf:
  //Coefficient of Sinh:
  RooAbsReal* f1_D0 = new RooFormulaVar("f1_D0",
										"f1_D0",
										"-2.0*ModlambdaF*cos(PhilambdaF)/(1.0 + pow(ModlambdaF,2))",   
										RooArgList(PhilambdaF,ModlambdaF));
  //Coefficient of Cos:
  RooAbsReal* f2_D0 = new RooFormulaVar("f2_D0",
										"f2_D0",
										"(1.0 - pow(ModlambdaF,2))/(1.0 + pow(ModlambdaF,2))",
										RooArgList(ModlambdaF));
  //Coefficient of Sin:
  RooAbsReal* f3_D0 = new RooFormulaVar("f3_D0",
										"f3_D0",
										"-2.0*ModlambdaF*sin(PhilambdaF)/(1.0 + pow(ModlambdaF,2))",    
										RooArgList(ModlambdaF,PhilambdaF));
  
  
  
  myDecay D0Pdf("D0Pdf","#Gamma(D^{0} #rightarrow f) pdf",
		  t,tag,
		  tau, dm,avgMistag,a,b,delMistag,mu,model,Opt_num);
		  //RooRealConstant::value(1),
		  //RooRealConstant::value(0),
		  //RooRealConstant::value((1.0 - pow(ModLf,2))/(1.0 + pow(ModLf,2))),
		  //RooRealConstant::value(-2.0*ModLf*sin(PhiLf)/(1.0 + pow(ModLf,2))),
		  //*f1_D0,
		  //*f2_D0,
		  //*f3_D0,
		  //dm,
		  //ResModel,
		  //Opt_num);

  //D0bar->f pdf:
  //Coefficient of Sinh:
  RooAbsReal* f1_D0bar = new RooFormulaVar("f1_D0bar",
					   "f1_D0bar",
					   "-2.0*ModlambdaF*cos(PhilambdaF)/(1.0 + pow(ModlambdaF,2))",  
					   RooArgList(PhilambdaF,ModlambdaF));
  //Coefficient of Cos:
  RooAbsReal* f2_D0bar = new RooFormulaVar("f2_D0bar",
					   "f2_D0bar",
					   "-(1.0 - pow(ModlambdaF,2))/(1.0 + pow(ModlambdaF,2))",
					   RooArgList(ModlambdaF));
  //Coefficient of Sin:
  RooAbsReal* f3_D0bar = new RooFormulaVar("f3_D0bar",
					   					   "f3_D0bar",
					   					   "2.0*ModlambdaF*sin(PhilambdaF)/(1.0 + pow(ModlambdaF,2))",
/*					   					   RooArgList(ModlambdaF,PhilambdaF));
  myDecay D0barPdf("D0barPdf","#Gamma(#bar{D}^{0} #rightarrow f) pdf",
		     t,
		     tau,DGamma,
		     RooRealConstant::value(1),
		     //RooRealConstant::value(0),
		     //RooRealConstant::value(-(1.0 - pow(ModLf,2))/(1.0 + pow(ModLf,2))),
		     //RooRealConstant::value(+2.0*ModLf*sin(PhiLf)/(1.0 + pow(ModLf,2))),
		     *f1_D0bar,
		     *f2_D0bar,
		     *f3_D0bar,
		     dm,
		     ResModel,
		     Opt_num);


  //Inclusion of mistagging effects:
  Double_t TheOmega      = Omega;
  Double_t TheDeltaOmega = DeltaOmega;
  RooRealVar omega("omega","Average mistag probability",TheOmega,0.0,1.0,"");
  RooRealVar Domega("Domega","mistag probability difference",TheDeltaOmega,-1.0,1.0,"");
  RooAbsReal* omega_D0 = new RooFormulaVar("omega_D0",
					   "omega_D0",
					   "omega + 0.5*Domega",
					   RooArgList(omega,Domega));
  RooAbsReal* omega_D0bar = new RooFormulaVar("omega_D0bar",
					      "omega_D0bar",
					      "omega - 0.5*Domega",
					      RooArgList(omega,Domega));
 /*
  //D0-tagged -> f pdf:
  RooAddPdf D0TaggedPdf("D0TaggedPdf",
			"D^{0}(tagged) #rightarrow f) pdf",
			RooArgList(D0barPdf,D0Pdf),
			RooArgList(*omega_D0));

  //D0-tagged -> f pdf:
  RooAddPdf D0barTaggedPdf("D0barTaggedPdf",
			   "#bar{D}^{0}(tagged) #rightarrow f) pdf",
			   RooArgList(D0Pdf,D0barPdf),
			   RooArgList(*omega_D0bar));

  //D0 flavor:
  RooCategory D0flav("D0flav","D0 flavour");
  D0flav.defineType("D0",    1);
  D0flav.defineType("D0bar",-1);

  RooSimultaneous FullDDecayPdf("FullDDecayPdf",
				"Full D^{0}/#bar{D}^{0} #rightarrow f decay pdf",
				D0flav);
  FullDDecayPdf.addPdf(D0TaggedPdf,   "D0");
  FullDDecayPdf.addPdf(D0barTaggedPdf,"D0bar");
  FullDDecayPdf.Print("v");


  TRandom r(seed);
  double N_D0_gen    = 0.5*Nevents*(1.0 + GenAsym);
  double N_D0bar_gen = 0.5*Nevents*(1.0 - GenAsym);

  cout << endl;
  cout << "Generating " << Nevents << " (N_D0 + N_D0bar) events, with asymmetry " 
       << GenAsym << " (N_D0 - N_D0bar)/(N_D0 + N_D0bar)" << endl;
  cout << "Generating " << N_D0_gen    << " D0    events" << endl;
  cout << "Generating " << N_D0bar_gen << " D0bar events" << endl;
  cout << endl;

 // int N_D0    = r.Poisson(N_D0_gen);
 // int N_D0bar = r.Poisson(N_D0bar_gen);

  cout << endl;
  cout << "N D0    = " << N_D0 << endl;
  cout << "N D0bar = " << N_D0bar << endl;
  cout << "N       = " << N_D0 + N_D0bar << endl;
  cout << "A       = " << (N_D0-N_D0bar)/(N_D0+N_D0bar) << endl;
  cout << endl;
*/

  //Generation of D0 events:
  //RooDataSet* data_D0    = D0TaggedPdf.generate(RooArgSet(t),N_D0);
  
  
  //Int_t Nevents = 5000;
 
  //TRandom r(seed);
  
  
  
  //int N_D0    = r.Poisson(Nevents);
  
  
 // RooDataSet* dataD0 = D0Pdf.generate(RooArgSet(t,tag),5000) ;
  
  
  // RooDataSet* data_D0    = D0Pdf.generate(RooArgSet(t,tag),500000);
  //Generation of D0 events:
 // RooDataSet* data_D0bar = D0barTaggedPdf.generate(RooArgSet(t),N_D0bar);
  //RooDataSet* data_D0bar = D0barPdf.generate(RooArgSet(t),N_D0bar);


  /*
  //Combine Samples:
  RooDataSet* combData = new RooDataSet("combData","combined D0 and D0bar data",
					t,
					Index(D0flav),
					Import("D0",   *data_D0),
					Import("D0bar",*data_D0bar));
  combData->Print("v");
*/

  //defining the variables which are constant in Fit:
  /*
  bias.setConstant(true);
  sigma.setConstant(true);
  tau.setConstant(true);
  DGamma.setConstant(true);
  dm.setConstant(true);
  ModlambdaF.setConstant(true);
  PhilambdaF.setConstant(false);
  //omega.setConstant(true);
 // Domega.setConstant(true);
  //RooFitResult* res = FullDDecayPdf.fitTo(*combData,Save());
  //res->Print("v");

  */
  /*
  //rmaria begin
  
  int nBin = 1000; 
  
  TH1F *funcD0 = new TH1F("funcD0","histo",nBin, R_t[0], R_t[1]);
  TH1F *funcD0bar = new TH1F("funcD0bar","histo",nBin,R_t[0], R_t[1]);
  TH1F *Funcasym = new TH1F("Funcasym","histo",nBin, R_t[0], R_t[1]);
  
  for(int i=0; i< nBin; i++)
  	{
  	Double_t tval = funcD0->GetBinCenter(i+1);
  	t.setVal(tval);
  	Double_t funcval = D0TaggedPdf.getVal(RooArgSet(t));
  	funcD0->SetBinContent(i+1,funcval);
  	
    funcval = D0barTaggedPdf.getVal(RooArgSet(t));
  	funcD0bar->SetBinContent(i+1,funcval);
  	}
  	
  Double_t integral1 = funcD0->Integral("width");
  Double_t integral2 = funcD0bar->Integral("width");
  if(integral1 > 0.)
  	{
  	funcD0->Scale(1./integral1);
  	}	
  if(integral2 > 0.)
  	{
  	funcD0bar->Scale(1./integral2);
  	}
  
  //funcD0->Scale(1.0/funcD0->Integral("width"));
  //funcD0bar->Scale(1.0/funcD0bar->Integral("width"));
  
  for(int i=0; i< nBin; i++)
  	{
  	Double_t valD0 = funcD0->GetBinContent(i+1);
  	Double_t valD0bar = funcD0bar->GetBinContent(i+1);
  	
  	Double_t difference = valD0bar - valD0;
  	Double_t sum = valD0bar + valD0;
  	Double_t asym_bin = difference / sum;
  	Funcasym->SetBinContent(i+1,asym_bin);
  	
  	}
  
  //Plot the 2 pdf inside a histogram (?)
  Double_t xD0 = 0;
  Double_t xD0bar = 0.;
  
  //take the normalized values
  //Double_t xnormValD0 =D0Pdf.getVal(xD0);
  //Double_t xnormValD0bar =D0barPdf.getVal(xD0bar);
  //-----------
  
  
  TH1F *histD0 = new TH1F("histD0","histo",nBin, 0., 5.);
  histD0->Sumw2();
  TH1F *histD0bar = new TH1F("histD0bar","histo",nBin, 0., 5.);
  histD0bar->Sumw2();
  
  //RooArgSet* set;
  
  for(int i=0; i<data_D0->numEntries(); i++)
  	{
    
  	const RooArgSet* set = data_D0->get(i);
    Double_t time = ((RooRealVar*)set->find(t.GetName()))->getVal();
  	histD0->Fill(time);
  	
  	}
  for(int i=0; i<data_D0bar->numEntries(); i++)
  	{
  	const RooArgSet* set = data_D0bar->get(i);
    Double_t time = ((RooRealVar*)set->find(t.GetName()))->getVal();
  	histD0bar->Fill(time);
  	
  	}
  
  RooHist asym2(*histD0bar, *histD0);
  asym2.SetName("asym2");
  TFile f("D0.root","RECREATE") ;
  
  histD0->Write();
  histD0bar->Write();
  asym2.Write();
  funcD0->Write();
  funcD0bar->Write();
  Funcasym->Write();
  /*  Is it necessary to write a ttree?
  TTree *mytree = new TTree("T","title");
  Double_t D0time[Nevents] ;
  Double_t antiD0time[Nevents] ;
  Double_t D0antiD0time[Nevents] ;
  
  TBranch *gaga1 = mytree->Branch("D0",&D0time,"D0time/D");
  TBranch *gaga2 = mytree->Branch("antiD0",&antiD0time,"antiD0/D");
  TBranch *gaga3 = mytree->Branch("D0antiD0",&D0antiD0time,"D0antiD0/D");
  
  */
  
  /*
  RooPlot *timeFrame_bubu = t.frame(t.getMin(),t.getMax(),100);
  timeFrame_bubu->SetTitle("Decay BUBU time distribution, projection D^{0}");
  D0TaggedPdf.plotOn(timeFrame_bubu,LineColor(kBlue));
  timeFrame_bubu->Draw();
  
  timeFrame_bubu->Write() ;
  
  
  f.Close() ;
  
  
  //rmaria end
*/
/*

  TString EPSName = TString(OutputFile) + TString(".root");   //rmaria it was .eps
  TString EPSNameO = EPSName + TString("[");
  TString EPSNameC = EPSName + TString("]");

  TCanvas* c1 = new TCanvas("c1","c1");
  c1->SetFillColor(10);
  c1->SetFrameFillColor(10);
  c1->SetTickx(1);
  c1->SetTicky(1);
  c1->SetLeftMargin(0.15);
  c1->SetBottomMargin(0.15);
  c1->SetRightMargin(0.10);

  c1->Print(EPSNameO.Data());
#if 0
  c1->Clear();
  RooPlot *timeFrame = t.frame(t.getMin(),t.getMax(),100);
  timeFrame->SetTitle("Decay time distribution");
  combData->plotOn(timeFrame,LineColor(kBlack),MarkerColor(kBlack));
  FullDDecayPdf.plotOn(timeFrame,LineColor(kRed));
  timeFrame->Draw();
  c1->Print(EPSName.Data());

#endif
  c1->Clear();
  c1->SetLogy(1);
  RooPlot *timeFrame_1 = t.frame(t.getMin(),t.getMax(),100);
  timeFrame_1->SetTitle("Decay time distribution, projection D^{0}");
  //combData->plotOn(timeFrame_1,Cut("D0flav==D0flav::D0"),LineColor(kBlue),MarkerColor(kBlue));
  //FullDDecayPdf.plotOn(timeFrame_1,Slice(D0flav,"D0"),ProjWData(D0flav,*combData),LineColor(kGreen+2));
  D0TaggedPdf.plotOn(timeFrame_1,LineColor(kBlue));
  timeFrame_1->Draw();
  c1->Print(EPSName.Data());
  c1->SetLogy(0);

  c1->Clear();
  c1->SetLogy(1);
  RooPlot *timeFrame_2 = t.frame(t.getMin(),t.getMax(),100);
  timeFrame_2->SetTitle("Decay time distribution, projection #bar{D}^{0}");
  //combData->plotOn(timeFrame_2,Cut("D0flav==D0flav::D0bar"),LineColor(kRed),MarkerColor(kRed));
  //FullDDecayPdf.plotOn(timeFrame_2,Slice(D0flav,"D0bar"),ProjWData(D0flav,*combData),LineColor(kOrange-2));
  D0barTaggedPdf.plotOn(timeFrame_2,LineColor(kRed));
  timeFrame_2->Draw();
  c1->Print(EPSName.Data());
  c1->SetLogy(0);

  c1->Clear();
  c1->SetLogy(1);
  RooPlot *timeFrame_3 = t.frame(t.getMin(),t.getMax(),100);
  timeFrame_3->SetTitle("Decay time distribution, projection D^{0}");
  D0TaggedPdf.plotOn(timeFrame_3,LineColor(kBlue));
  D0barTaggedPdf.plotOn(timeFrame_3,LineColor(kRed));
  timeFrame_3->Draw();
  c1->Print(EPSName.Data());
  c1->SetLogy(0);

  c1->Print(EPSNameC.Data());
  
*/
  return;

}
//============================================================================


