#include "/home/aperez/BaBarAnalyzes/bTosGamma/AnalysisCode/ThePi0EtaMassPdf/src/RooDoubleCBShape.hh"
#include "/home/aperez/BaBarAnalyzes/bTosGamma/AnalysisCode/ThePi0EtaMassPdf/src/RooBifurcatedCBShape.hh"
#include "/home/aperez/BaBarAnalyzes/bTosGamma/AnalysisCode/ThePi0EtaMassPdf/src/RooPolyAndDecay.hh"
#include "/home/aperez/BaBarAnalyzes/bTosGamma/AnalysisCode/ThePi0EtaMassPdf/src/RooPolyAndExpDecay.hh"
#include "/home/aperez/BaBarAnalyzes/bTosGamma/AnalysisCode/ThePi0EtaMassPdf/src/RooPolyAndExpDecayV2.hh"
#include "/home/aperez/BaBarAnalyzes/bTosGamma/AnalysisCode/ThePi0EtaMassPdf/src/RooPolyAndExpDecayV3.hh"
#include "BXsGamma_Analysis.hh"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TROOT.h>
#include <TH1.h>
#include <TLegend.h>
#include <TTree.h>
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
#include <TTree.h>
#include <TGraphErrors.h>
#include <TFitResultPtr.h>
#include <TFitResult.h>
#include <TMatrixDSym.h>
#include <TMatrixDSymEigen.h>
#include <RooRealVar.h>
#include <RooDataSet.h>
#include <RooCBShape.h>
#include <RooArgusBG.h>
#include <RooAddPdf.h>
#include <RooAbsPdf.h>
#include <RooFitResult.h>
#include <RooPlot.h>
#include <RooLinkedList.h>
#include <RooGaussian.h>
#include <RooPolynomial.h>
#include <RooProdPdf.h>
#include <RooBifurGauss.h>
#include <RooHist.h>
#include <RooChebychev.h>
#include <RooGlobalFunc.h>
#include <RooConstVar.h>
#include <RooAbsReal.h>

//TMVA
#include "TMVA/Factory.h"
#include "TMVA/Tools.h"
#include "TMVA/Reader.h"
#include "TMVA/MethodCuts.h"

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
void BXsGamma_Analysis::Loop(const char* outputFile)
{

  double epsilon = 1.0e-8;

  TLatex* latex = new TLatex();
  latex->SetTextAlign(12);
  latex->SetTextSize(0.03);
  latex->SetTextColor(1);

  //Bkg-control region:
  double R_EgBRF_Bkg[2];
  R_EgBRF_Bkg[0] = 1.3;
  R_EgBRF_Bkg[1] = 1.8;

  //Bkg-control region2:
  double R_EgBRF_Bkg2[2];
  R_EgBRF_Bkg2[0] = 2.7;
  R_EgBRF_Bkg2[1] = 3.0;

  //Signal Region
  double R_EgBRF_Sig[2];
  R_EgBRF_Sig[0] = 1.8;
  R_EgBRF_Sig[1] = 2.7;

  TStyle* style = new TStyle();
  style->SetPalette(1);

  RooFitResult* resFit_mES_all = NULL;
  bool DoGlobalFit = true;
  if(_FitResFile != TString("XXX")) {
    DoGlobalFit = false;
    TFile f_res(_FitResFile.Data());
    resFit_mES_all = (RooFitResult*)f_res.Get("resFit_mES_all");
    f_res.Close();
    resFit_mES_all->Print("v");
  }

  TH1D* h_EgBRF_AllBkgPrediction_tot_corr                     = NULL;
  TH2D* h_EgBRF_AllBkgPrediction_Cov_AllStat                  = NULL;
  TH2D* h_EgBRF_AllBkgPrediction_Cov_AllHEEfficCorr           = NULL;
  TH2D* h_EgBRF_AllBkgPrediction_Cov_NonPi0EtaBkgBrecoCorr    = NULL;
  TH2D* h_EgBRF_AllBkgPrediction_Cov_Pi0EtaBkgMCYield         = NULL;
  TH2D* h_EgBRF_AllBkgPrediction_Cov_Pi0EtaBkgRpi0Reta        = NULL;
  TH2D* h_EgBRF_AllBkgPrediction_Cov_Pi0EtaBkgRpi0CRetaC      = NULL;
  TH2D* h_EgBRF_AllBkgPrediction_Cov_Pi0EtaBkgmESFitYields    = NULL;
  TH2D* h_EgBRF_AllBkgPrediction_Cov_Pi0EtaBkgmESFixX1        = NULL;
  TH2D* h_EgBRF_AllBkgPrediction_Cov_Pi0EtaBkgmESFixX2        = NULL;
  TH2D* h_EgBRF_AllBkgPrediction_Cov_Pi0EtaBkgmESFixX3        = NULL;
  TH2D* h_EgBRF_AllBkgPrediction_Cov_Pi0EtaBkgmESFix_alphaP   = NULL;
  TH2D* h_EgBRF_AllBkgPrediction_Cov_Pi0EtaBkgmESFix_BBComb   = NULL;
  TH2D* h_EgBRF_AllBkgPrediction_Cov_Pi0EtaBkgmESFix_slope    = NULL;
  TH2D* h_EgBRF_AllBkgPrediction_Cov_Pi0EtaBkgmESFitBias      = NULL;
  TH2D* h_EgBRF_AllBkgPrediction_Cov_Pi0EtaBkgMggFitValSyst   = NULL;
  TH2D* h_EgBRF_AllBkgPrediction_Cov_Pi0EtaBkgLE_Base         = NULL;
  TH2D* h_EgBRF_AllBkgPrediction_Cov_Pi0EtaBkgULE             = NULL;
  TH2D* h_EgBRF_AllBkgPrediction_Cov_OmegaBkgAlpha            = NULL;
  TH2D* h_EgBRF_AllBkgPrediction_Cov_etaprimeBkgAlpha         = NULL;
  TH2D* h_EgBRF_AllBkgPrediction_Cov_BkgSLAlpha               = NULL;
  TH2D* h_EgBRF_AllBkgPrediction_Cov_BremsBkgMaterial         = NULL;
  TH2D* h_EgBRF_AllBkgPrediction_Cov_ElecBkgTrackIneffic      = NULL;
  TH2D* h_EgBRF_AllBkgPrediction_Cov_AntiNeutronsBkgAlpha     = NULL;
  TH2D* h_EgBRF_AllBkgPrediction_Cov_BDT                      = NULL;
  TH2D* h_EgBRF_AllBkgPrediction_Cov_Pi0EtaVeto_MassLineShape = NULL;
  TH2D* h_EgBRF_AllBkgPrediction_Cov_Pi0EtaVeto_2ndGammaMult  = NULL;

  TH2D* h_mESFitYield_vs_EgBRF_Cov_Stat            = NULL;
  TH2D* h_mESFitYield_vs_EgBRF_Cov_mESFixX1        = NULL;
  TH2D* h_mESFitYield_vs_EgBRF_Cov_mESFixX2        = NULL;
  TH2D* h_mESFitYield_vs_EgBRF_Cov_mESFixX3        = NULL;
  TH2D* h_mESFitYield_vs_EgBRF_Cov_mESFix_alphaP   = NULL;
  TH2D* h_mESFitYield_vs_EgBRF_Cov_mESFix_BBComb   = NULL;
  TH2D* h_mESFitYield_vs_EgBRF_Cov_mESFix_slope    = NULL;
  TH2D* h_mESFitYield_vs_EgBRF_Cov_mESFitBias      = NULL;

  TH2D* h_SignalYield_vs_EgBRF_Cov_Stat                     = NULL;
  TH2D* h_SignalYield_vs_EgBRF_Cov_BBkgStat                 = NULL;
  TH2D* h_SignalYield_vs_EgBRF_Cov_AllHEEfficCorr           = NULL;
  TH2D* h_SignalYield_vs_EgBRF_Cov_BrecoCorr                = NULL;
  TH2D* h_SignalYield_vs_EgBRF_Cov_Pi0EtaBkgMCYield         = NULL;
  TH2D* h_SignalYield_vs_EgBRF_Cov_Pi0EtaBkgRpi0Reta        = NULL;
  TH2D* h_SignalYield_vs_EgBRF_Cov_Pi0EtaBkgRpi0CRetaC      = NULL;
  TH2D* h_SignalYield_vs_EgBRF_Cov_Pi0EtaBkgmESFitYields    = NULL;
  TH2D* h_SignalYield_vs_EgBRF_Cov_mESFixX1                 = NULL;
  TH2D* h_SignalYield_vs_EgBRF_Cov_mESFixX2                 = NULL;
  TH2D* h_SignalYield_vs_EgBRF_Cov_mESFixX3                 = NULL;
  TH2D* h_SignalYield_vs_EgBRF_Cov_mESFix_alphaP            = NULL;
  TH2D* h_SignalYield_vs_EgBRF_Cov_mESFix_BBComb            = NULL;
  TH2D* h_SignalYield_vs_EgBRF_Cov_mESFix_slope             = NULL;
  TH2D* h_SignalYield_vs_EgBRF_Cov_mESFitBias               = NULL;
  TH2D* h_SignalYield_vs_EgBRF_Cov_Pi0EtaBkgMggFitValSyst   = NULL;
  TH2D* h_SignalYield_vs_EgBRF_Cov_Pi0EtaBkgLE_Base         = NULL;
  TH2D* h_SignalYield_vs_EgBRF_Cov_Pi0EtaBkgULE             = NULL;
  TH2D* h_SignalYield_vs_EgBRF_Cov_OmegaBkgAlpha            = NULL;
  TH2D* h_SignalYield_vs_EgBRF_Cov_etaprimeBkgAlpha         = NULL;
  TH2D* h_SignalYield_vs_EgBRF_Cov_BkgSLAlpha               = NULL;
  TH2D* h_SignalYield_vs_EgBRF_Cov_BremsBkgMaterial         = NULL;
  TH2D* h_SignalYield_vs_EgBRF_Cov_ElecBkgTrackIneffic      = NULL;
  TH2D* h_SignalYield_vs_EgBRF_Cov_AntiNeutronsBkgAlpha     = NULL;
  TH2D* h_SignalYield_vs_EgBRF_Cov_BDT                      = NULL;
  TH2D* h_SignalYield_vs_EgBRF_Cov_Pi0EtaVeto_MassLineShape = NULL;
  TH2D* h_SignalYield_vs_EgBRF_Cov_Pi0EtaVeto_2ndGammaMult  = NULL;

  cout << endl;
  cout << "Get BB-background predictions from file: " << _BkgPredFile.Data() << endl;
  cout << endl;
  TFile f_Bkg(_BkgPredFile.Data());
  //Bkg-prediction Cov matrices:
  h_EgBRF_AllBkgPrediction_tot_corr = (TH1D*)f_Bkg.Get("h_EgBRF_AllBkgPrediction_tot_corr");
  h_EgBRF_AllBkgPrediction_tot_corr->SetDirectory(0);

  h_EgBRF_AllBkgPrediction_Cov_AllStat = (TH2D*)f_Bkg.Get("h_EgBRF_AllBkgPrediction_Cov_AllStat");
  h_EgBRF_AllBkgPrediction_Cov_AllStat->SetDirectory(0);

  h_EgBRF_AllBkgPrediction_Cov_AllHEEfficCorr = (TH2D*)f_Bkg.Get("h_EgBRF_AllBkgPrediction_Cov_AllHEEfficCorr");
  h_EgBRF_AllBkgPrediction_Cov_AllHEEfficCorr->SetDirectory(0);

  h_EgBRF_AllBkgPrediction_Cov_NonPi0EtaBkgBrecoCorr = (TH2D*)f_Bkg.Get("h_EgBRF_AllBkgPrediction_Cov_NonPi0EtaBkgBrecoCorr");
  h_EgBRF_AllBkgPrediction_Cov_NonPi0EtaBkgBrecoCorr->SetDirectory(0);

  h_EgBRF_AllBkgPrediction_Cov_Pi0EtaBkgMCYield = (TH2D*)f_Bkg.Get("h_EgBRF_AllBkgPrediction_Cov_Pi0EtaBkgMCYield");
  h_EgBRF_AllBkgPrediction_Cov_Pi0EtaBkgMCYield->SetDirectory(0);

  h_EgBRF_AllBkgPrediction_Cov_Pi0EtaBkgRpi0Reta = (TH2D*)f_Bkg.Get("h_EgBRF_AllBkgPrediction_Cov_Pi0EtaBkgRpi0Reta");
  h_EgBRF_AllBkgPrediction_Cov_Pi0EtaBkgRpi0Reta->SetDirectory(0);

  h_EgBRF_AllBkgPrediction_Cov_Pi0EtaBkgRpi0CRetaC = (TH2D*)f_Bkg.Get("h_EgBRF_AllBkgPrediction_Cov_Pi0EtaBkgRpi0CRetaC");
  h_EgBRF_AllBkgPrediction_Cov_Pi0EtaBkgRpi0CRetaC->SetDirectory(0);

  h_EgBRF_AllBkgPrediction_Cov_Pi0EtaBkgmESFitYields = (TH2D*)f_Bkg.Get("h_EgBRF_AllBkgPrediction_Cov_Pi0EtaBkgmESFitYields");
  h_EgBRF_AllBkgPrediction_Cov_Pi0EtaBkgmESFitYields->SetDirectory(0);

  h_EgBRF_AllBkgPrediction_Cov_Pi0EtaBkgmESFixX1 = (TH2D*)f_Bkg.Get("h_EgBRF_AllBkgPrediction_Cov_Pi0EtaBkgmESFixX1");
  h_EgBRF_AllBkgPrediction_Cov_Pi0EtaBkgmESFixX1->SetDirectory(0);

  h_EgBRF_AllBkgPrediction_Cov_Pi0EtaBkgmESFixX2 = (TH2D*)f_Bkg.Get("h_EgBRF_AllBkgPrediction_Cov_Pi0EtaBkgmESFixX2");
  h_EgBRF_AllBkgPrediction_Cov_Pi0EtaBkgmESFixX2->SetDirectory(0);

  h_EgBRF_AllBkgPrediction_Cov_Pi0EtaBkgmESFixX3 = (TH2D*)f_Bkg.Get("h_EgBRF_AllBkgPrediction_Cov_Pi0EtaBkgmESFixX3");
  h_EgBRF_AllBkgPrediction_Cov_Pi0EtaBkgmESFixX3->SetDirectory(0);

  h_EgBRF_AllBkgPrediction_Cov_Pi0EtaBkgmESFix_alphaP = (TH2D*)f_Bkg.Get("h_EgBRF_AllBkgPrediction_Cov_Pi0EtaBkgmESFix_alphaP");
  h_EgBRF_AllBkgPrediction_Cov_Pi0EtaBkgmESFix_alphaP->SetDirectory(0);

  h_EgBRF_AllBkgPrediction_Cov_Pi0EtaBkgmESFix_BBComb = (TH2D*)f_Bkg.Get("h_EgBRF_AllBkgPrediction_Cov_Pi0EtaBkgmESFix_BBComb");
  h_EgBRF_AllBkgPrediction_Cov_Pi0EtaBkgmESFix_BBComb->SetDirectory(0);

  h_EgBRF_AllBkgPrediction_Cov_Pi0EtaBkgmESFix_slope = (TH2D*)f_Bkg.Get("h_EgBRF_AllBkgPrediction_Cov_Pi0EtaBkgmESFix_slope");
  h_EgBRF_AllBkgPrediction_Cov_Pi0EtaBkgmESFix_slope->SetDirectory(0);

  h_EgBRF_AllBkgPrediction_Cov_Pi0EtaBkgmESFitBias = (TH2D*)f_Bkg.Get("h_EgBRF_AllBkgPrediction_Cov_Pi0EtaBkgmESFitBias");
  h_EgBRF_AllBkgPrediction_Cov_Pi0EtaBkgmESFitBias->SetDirectory(0);

  h_EgBRF_AllBkgPrediction_Cov_Pi0EtaBkgMggFitValSyst = (TH2D*)f_Bkg.Get("h_EgBRF_AllBkgPrediction_Cov_Pi0EtaBkgMggFitValSyst");
  h_EgBRF_AllBkgPrediction_Cov_Pi0EtaBkgMggFitValSyst->SetDirectory(0);

  h_EgBRF_AllBkgPrediction_Cov_Pi0EtaBkgLE_Base = (TH2D*)f_Bkg.Get("h_EgBRF_AllBkgPrediction_Cov_Pi0EtaBkgLE_Base");
  h_EgBRF_AllBkgPrediction_Cov_Pi0EtaBkgLE_Base->SetDirectory(0);

  h_EgBRF_AllBkgPrediction_Cov_Pi0EtaBkgULE = (TH2D*)f_Bkg.Get("h_EgBRF_AllBkgPrediction_Cov_Pi0EtaBkgULE");
  h_EgBRF_AllBkgPrediction_Cov_Pi0EtaBkgULE->SetDirectory(0);

  h_EgBRF_AllBkgPrediction_Cov_OmegaBkgAlpha = (TH2D*)f_Bkg.Get("h_EgBRF_AllBkgPrediction_Cov_OmegaBkgAlpha");
  h_EgBRF_AllBkgPrediction_Cov_OmegaBkgAlpha->SetDirectory(0);

  h_EgBRF_AllBkgPrediction_Cov_etaprimeBkgAlpha = (TH2D*)f_Bkg.Get("h_EgBRF_AllBkgPrediction_Cov_etaprimeBkgAlpha");
  h_EgBRF_AllBkgPrediction_Cov_etaprimeBkgAlpha->SetDirectory(0);

  h_EgBRF_AllBkgPrediction_Cov_BkgSLAlpha = (TH2D*)f_Bkg.Get("h_EgBRF_AllBkgPrediction_Cov_BkgSLAlpha");
  h_EgBRF_AllBkgPrediction_Cov_BkgSLAlpha->SetDirectory(0);

  h_EgBRF_AllBkgPrediction_Cov_BremsBkgMaterial = (TH2D*)f_Bkg.Get("h_EgBRF_AllBkgPrediction_Cov_BremsBkgMaterial");
  h_EgBRF_AllBkgPrediction_Cov_BremsBkgMaterial->SetDirectory(0);

  h_EgBRF_AllBkgPrediction_Cov_ElecBkgTrackIneffic = (TH2D*)f_Bkg.Get("h_EgBRF_AllBkgPrediction_Cov_ElecBkgTrackIneffic");
  h_EgBRF_AllBkgPrediction_Cov_ElecBkgTrackIneffic->SetDirectory(0);

  h_EgBRF_AllBkgPrediction_Cov_AntiNeutronsBkgAlpha = (TH2D*)f_Bkg.Get("h_EgBRF_AllBkgPrediction_Cov_AntiNeutronsBkgAlpha");
  h_EgBRF_AllBkgPrediction_Cov_AntiNeutronsBkgAlpha->SetDirectory(0);

  h_EgBRF_AllBkgPrediction_Cov_BDT = (TH2D*)f_Bkg.Get("h_EgBRF_AllBkgPrediction_Cov_BDT");
  h_EgBRF_AllBkgPrediction_Cov_BDT->SetDirectory(0);

  h_EgBRF_AllBkgPrediction_Cov_Pi0EtaVeto_MassLineShape = (TH2D*)f_Bkg.Get("h_EgBRF_AllBkgPrediction_Cov_Pi0EtaVeto_MassLineShape");
  h_EgBRF_AllBkgPrediction_Cov_Pi0EtaVeto_MassLineShape->SetDirectory(0);

  h_EgBRF_AllBkgPrediction_Cov_Pi0EtaVeto_2ndGammaMult = (TH2D*)f_Bkg.Get("h_EgBRF_AllBkgPrediction_Cov_Pi0EtaVeto_2ndGammaMult");
  h_EgBRF_AllBkgPrediction_Cov_Pi0EtaVeto_2ndGammaMult->SetDirectory(0);

  //mES fit Peaking Yield Cov matrices:
  h_mESFitYield_vs_EgBRF_Cov_Stat = (TH2D*)f_Bkg.Get("h_EgBRF_AllBkgPrediction_Cov_BDT");
  h_mESFitYield_vs_EgBRF_Cov_Stat->SetName("h_mESFitYield_vs_EgBRF_Cov_Stat");
  h_mESFitYield_vs_EgBRF_Cov_Stat->SetTitle("Cov(N_{i},N_{j}) for m_{ES} peaking yield from Stat");
  h_mESFitYield_vs_EgBRF_Cov_Stat->SetDirectory(0);

  h_mESFitYield_vs_EgBRF_Cov_mESFixX1 = (TH2D*)f_Bkg.Get("h_EgBRF_AllBkgPrediction_Cov_BDT");
  h_mESFitYield_vs_EgBRF_Cov_mESFixX1->SetName("h_mESFitYield_vs_EgBRF_Cov_mESFixX1");
  h_mESFitYield_vs_EgBRF_Cov_mESFixX1->SetTitle("Cov(N_{i},N_{j}) for m_{ES} peaking yield from Fix-par X_{1}");
  h_mESFitYield_vs_EgBRF_Cov_mESFixX1->SetDirectory(0);

  h_mESFitYield_vs_EgBRF_Cov_mESFixX2 = (TH2D*)f_Bkg.Get("h_EgBRF_AllBkgPrediction_Cov_BDT");
  h_mESFitYield_vs_EgBRF_Cov_mESFixX2->SetName("h_mESFitYield_vs_EgBRF_Cov_mESFixX2");
  h_mESFitYield_vs_EgBRF_Cov_mESFixX2->SetTitle("Cov(N_{i},N_{j}) for m_{ES} peaking yield from Fix-par X_{2}");
  h_mESFitYield_vs_EgBRF_Cov_mESFixX2->SetDirectory(0);

  h_mESFitYield_vs_EgBRF_Cov_mESFixX3 = (TH2D*)f_Bkg.Get("h_EgBRF_AllBkgPrediction_Cov_BDT");
  h_mESFitYield_vs_EgBRF_Cov_mESFixX3->SetName("h_mESFitYield_vs_EgBRF_Cov_mESFixX3");
  h_mESFitYield_vs_EgBRF_Cov_mESFixX3->SetTitle("Cov(N_{i},N_{j}) for m_{ES} peaking yield from Fix-par X_{3}");
  h_mESFitYield_vs_EgBRF_Cov_mESFixX3->SetDirectory(0);

  h_mESFitYield_vs_EgBRF_Cov_mESFix_alphaP = (TH2D*)f_Bkg.Get("h_EgBRF_AllBkgPrediction_Cov_BDT");
  h_mESFitYield_vs_EgBRF_Cov_mESFix_alphaP->SetName("h_mESFitYield_vs_EgBRF_Cov_mESFix_alphaP");
  h_mESFitYield_vs_EgBRF_Cov_mESFix_alphaP->SetTitle("Cov(N_{i},N_{j}) for m_{ES} peaking yield from Fix-par #alpha^{P}");
  h_mESFitYield_vs_EgBRF_Cov_mESFix_alphaP->SetDirectory(0);

  h_mESFitYield_vs_EgBRF_Cov_mESFix_BBComb = (TH2D*)f_Bkg.Get("h_EgBRF_AllBkgPrediction_Cov_BDT");
  h_mESFitYield_vs_EgBRF_Cov_mESFix_BBComb->SetName("h_mESFitYield_vs_EgBRF_Cov_mESFix_BBComb");
  h_mESFitYield_vs_EgBRF_Cov_mESFix_BBComb->SetTitle("Cov(N_{i},N_{j}) for m_{ES} peaking yield from Fix-par BB-comb");
  h_mESFitYield_vs_EgBRF_Cov_mESFix_BBComb->SetDirectory(0);

  h_mESFitYield_vs_EgBRF_Cov_mESFix_slope = (TH2D*)f_Bkg.Get("h_EgBRF_AllBkgPrediction_Cov_BDT");
  h_mESFitYield_vs_EgBRF_Cov_mESFix_slope->SetName("h_mESFitYield_vs_EgBRF_Cov_mESFix_slope");
  h_mESFitYield_vs_EgBRF_Cov_mESFix_slope->SetTitle("Cov(N_{i},N_{j}) for m_{ES} peaking yield from Fix-par Argus-slope");
  h_mESFitYield_vs_EgBRF_Cov_mESFix_slope->SetDirectory(0);

  h_mESFitYield_vs_EgBRF_Cov_mESFitBias = (TH2D*)f_Bkg.Get("h_EgBRF_AllBkgPrediction_Cov_BDT");
  h_mESFitYield_vs_EgBRF_Cov_mESFitBias->SetName("h_mESFitYield_vs_EgBRF_Cov_mESFitBias");
  h_mESFitYield_vs_EgBRF_Cov_mESFitBias->SetTitle("Cov(N_{i},N_{j}) for m_{ES} peaking yield from Fit-bias");
  h_mESFitYield_vs_EgBRF_Cov_mESFitBias->SetDirectory(0);

  //Signal Yield Cov matrices:
  h_SignalYield_vs_EgBRF_Cov_Stat = (TH2D*)f_Bkg.Get("h_EgBRF_AllBkgPrediction_Cov_BDT");
  h_SignalYield_vs_EgBRF_Cov_Stat->SetName("h_SignalYield_vs_EgBRF_Cov_Stat");
  h_SignalYield_vs_EgBRF_Cov_Stat->SetTitle("Cov(S_{i},S_{j}) Data stat error");
  h_SignalYield_vs_EgBRF_Cov_Stat->SetDirectory(0);

  h_SignalYield_vs_EgBRF_Cov_BBkgStat = (TH2D*)f_Bkg.Get("h_EgBRF_AllBkgPrediction_Cov_BDT");
  h_SignalYield_vs_EgBRF_Cov_BBkgStat->SetName("h_SignalYield_vs_EgBRF_Cov_BBkgStat");
  h_SignalYield_vs_EgBRF_Cov_BBkgStat->SetTitle("Cov(S_{i},S_{j}) BBkg stat error");
  h_SignalYield_vs_EgBRF_Cov_BBkgStat->SetDirectory(0);

  h_SignalYield_vs_EgBRF_Cov_AllHEEfficCorr = (TH2D*)f_Bkg.Get("h_EgBRF_AllBkgPrediction_Cov_BDT");
  h_SignalYield_vs_EgBRF_Cov_AllHEEfficCorr->SetName("h_SignalYield_vs_EgBRF_Cov_AllHEEfficCorr");
  h_SignalYield_vs_EgBRF_Cov_AllHEEfficCorr->SetTitle("Cov(S_{i},S_{j}) BBkg EH-#gamma Effic syst");
  h_SignalYield_vs_EgBRF_Cov_AllHEEfficCorr->SetDirectory(0);

  h_SignalYield_vs_EgBRF_Cov_BrecoCorr = (TH2D*)f_Bkg.Get("h_EgBRF_AllBkgPrediction_Cov_BDT");
  h_SignalYield_vs_EgBRF_Cov_BrecoCorr->SetName("h_SignalYield_vs_EgBRF_Cov_BrecoCorr");
  h_SignalYield_vs_EgBRF_Cov_BrecoCorr->SetTitle("Cov(S_{i},S_{j}) BBkg Breco corr. syst");
  h_SignalYield_vs_EgBRF_Cov_BrecoCorr->SetDirectory(0);

  h_SignalYield_vs_EgBRF_Cov_Pi0EtaBkgMCYield = (TH2D*)f_Bkg.Get("h_EgBRF_AllBkgPrediction_Cov_BDT");
  h_SignalYield_vs_EgBRF_Cov_Pi0EtaBkgMCYield->SetName("h_SignalYield_vs_EgBRF_Cov_Pi0EtaBkgMCYield");
  h_SignalYield_vs_EgBRF_Cov_Pi0EtaBkgMCYield->SetTitle("Cov(S_{i},S_{j}) BBkg #pi^{0}/#eta analysis MC-yield error");
  h_SignalYield_vs_EgBRF_Cov_Pi0EtaBkgMCYield->SetDirectory(0);

  h_SignalYield_vs_EgBRF_Cov_Pi0EtaBkgRpi0Reta = (TH2D*)f_Bkg.Get("h_EgBRF_AllBkgPrediction_Cov_BDT");
  h_SignalYield_vs_EgBRF_Cov_Pi0EtaBkgRpi0Reta->SetName("h_SignalYield_vs_EgBRF_Cov_Pi0EtaBkgRpi0Reta");
  h_SignalYield_vs_EgBRF_Cov_Pi0EtaBkgRpi0Reta->SetTitle("Cov(S_{i},S_{j}) BBkg #pi^{0}/#eta analysis R(#pi^{0})/R(#eta) error");
  h_SignalYield_vs_EgBRF_Cov_Pi0EtaBkgRpi0Reta->SetDirectory(0);

  h_SignalYield_vs_EgBRF_Cov_Pi0EtaBkgRpi0CRetaC = (TH2D*)f_Bkg.Get("h_EgBRF_AllBkgPrediction_Cov_BDT");
  h_SignalYield_vs_EgBRF_Cov_Pi0EtaBkgRpi0CRetaC->SetName("h_SignalYield_vs_EgBRF_Cov_Pi0EtaBkgRpi0CRetaC");
  h_SignalYield_vs_EgBRF_Cov_Pi0EtaBkgRpi0CRetaC->SetTitle("Cov(S_{i},S_{j}) BBkg #pi^{0}/#eta analysis R^{C}(#pi^{0})/R^{C}(#eta) error");
  h_SignalYield_vs_EgBRF_Cov_Pi0EtaBkgRpi0CRetaC->SetDirectory(0);

  h_SignalYield_vs_EgBRF_Cov_Pi0EtaBkgmESFitYields = (TH2D*)f_Bkg.Get("h_EgBRF_AllBkgPrediction_Cov_BDT");
  h_SignalYield_vs_EgBRF_Cov_Pi0EtaBkgmESFitYields->SetName("h_SignalYield_vs_EgBRF_Cov_Pi0EtaBkgmESFitYields");
  h_SignalYield_vs_EgBRF_Cov_Pi0EtaBkgmESFitYields->SetTitle("Cov(S_{i},S_{j}) BBkg #pi^{0}/#eta analysis m_{ES}-fit yield errors");
  h_SignalYield_vs_EgBRF_Cov_Pi0EtaBkgmESFitYields->SetDirectory(0);

  h_SignalYield_vs_EgBRF_Cov_mESFixX1 = (TH2D*)f_Bkg.Get("h_EgBRF_AllBkgPrediction_Cov_BDT");
  h_SignalYield_vs_EgBRF_Cov_mESFixX1->SetName("h_SignalYield_vs_EgBRF_Cov_mESFixX1");
  h_SignalYield_vs_EgBRF_Cov_mESFixX1->SetTitle("Cov(S_{i},S_{j}) m_{ES}-fix X_{1} error");
  h_SignalYield_vs_EgBRF_Cov_mESFixX1->SetDirectory(0);

  h_SignalYield_vs_EgBRF_Cov_mESFixX2 = (TH2D*)f_Bkg.Get("h_EgBRF_AllBkgPrediction_Cov_BDT");
  h_SignalYield_vs_EgBRF_Cov_mESFixX2->SetName("h_SignalYield_vs_EgBRF_Cov_mESFixX2");
  h_SignalYield_vs_EgBRF_Cov_mESFixX2->SetTitle("Cov(S_{i},S_{j}) m_{ES}-fix X_{2} error");
  h_SignalYield_vs_EgBRF_Cov_mESFixX2->SetDirectory(0);

  h_SignalYield_vs_EgBRF_Cov_mESFixX3 = (TH2D*)f_Bkg.Get("h_EgBRF_AllBkgPrediction_Cov_BDT");
  h_SignalYield_vs_EgBRF_Cov_mESFixX3->SetName("h_SignalYield_vs_EgBRF_Cov_mESFixX3");
  h_SignalYield_vs_EgBRF_Cov_mESFixX3->SetTitle("Cov(S_{i},S_{j}) m_{ES}-fix X_{3} error");
  h_SignalYield_vs_EgBRF_Cov_mESFixX3->SetDirectory(0);

  h_SignalYield_vs_EgBRF_Cov_mESFix_alphaP = (TH2D*)f_Bkg.Get("h_EgBRF_AllBkgPrediction_Cov_BDT");
  h_SignalYield_vs_EgBRF_Cov_mESFix_alphaP->SetName("h_SignalYield_vs_EgBRF_Cov_mESFix_alphaP");
  h_SignalYield_vs_EgBRF_Cov_mESFix_alphaP->SetTitle("Cov(S_{i},S_{j}) m_{ES}-fix #alpha^{P} error");
  h_SignalYield_vs_EgBRF_Cov_mESFix_alphaP->SetDirectory(0);

  h_SignalYield_vs_EgBRF_Cov_mESFix_BBComb = (TH2D*)f_Bkg.Get("h_EgBRF_AllBkgPrediction_Cov_BDT");
  h_SignalYield_vs_EgBRF_Cov_mESFix_BBComb->SetName("h_SignalYield_vs_EgBRF_Cov_mESFix_BBComb");
  h_SignalYield_vs_EgBRF_Cov_mESFix_BBComb->SetTitle("Cov(S_{i},S_{j}) m_{ES}-fix BB-comb error");
  h_SignalYield_vs_EgBRF_Cov_mESFix_BBComb->SetDirectory(0);

  h_SignalYield_vs_EgBRF_Cov_mESFix_slope = (TH2D*)f_Bkg.Get("h_EgBRF_AllBkgPrediction_Cov_BDT");
  h_SignalYield_vs_EgBRF_Cov_mESFix_slope->SetName("h_SignalYield_vs_EgBRF_Cov_mESFix_slope");
  h_SignalYield_vs_EgBRF_Cov_mESFix_slope->SetTitle("Cov(S_{i},S_{j}) m_{ES}-fix Argus-slope error");
  h_SignalYield_vs_EgBRF_Cov_mESFix_slope->SetDirectory(0);

  h_SignalYield_vs_EgBRF_Cov_mESFitBias = (TH2D*)f_Bkg.Get("h_EgBRF_AllBkgPrediction_Cov_BDT");
  h_SignalYield_vs_EgBRF_Cov_mESFitBias->SetName("h_SignalYield_vs_EgBRF_Cov_mESFitBias");
  h_SignalYield_vs_EgBRF_Cov_mESFitBias->SetTitle("Cov(S_{i},S_{j}) m_{ES}-fit Bias syst.");
  h_SignalYield_vs_EgBRF_Cov_mESFitBias->SetDirectory(0);

  h_SignalYield_vs_EgBRF_Cov_Pi0EtaBkgMggFitValSyst = (TH2D*)f_Bkg.Get("h_EgBRF_AllBkgPrediction_Cov_BDT");
  h_SignalYield_vs_EgBRF_Cov_Pi0EtaBkgMggFitValSyst->SetName("h_SignalYield_vs_EgBRF_Cov_Pi0EtaBkgMggFitValSyst");
  h_SignalYield_vs_EgBRF_Cov_Pi0EtaBkgMggFitValSyst->SetTitle("Cov(S_{i},S_{j}) BBkg #pi^{0}/#eta analysis, M_{#gamma#gamma} Val. syst.");
  h_SignalYield_vs_EgBRF_Cov_Pi0EtaBkgMggFitValSyst->SetDirectory(0);

  h_SignalYield_vs_EgBRF_Cov_Pi0EtaBkgLE_Base = (TH2D*)f_Bkg.Get("h_EgBRF_AllBkgPrediction_Cov_BDT");
  h_SignalYield_vs_EgBRF_Cov_Pi0EtaBkgLE_Base->SetName("h_SignalYield_vs_EgBRF_Cov_Pi0EtaBkgLE_Base");
  h_SignalYield_vs_EgBRF_Cov_Pi0EtaBkgLE_Base->SetTitle("Cov(S_{i},S_{j}) BBkg #epsilon^{Base}_{LE} error");
  h_SignalYield_vs_EgBRF_Cov_Pi0EtaBkgLE_Base->SetDirectory(0);

  h_SignalYield_vs_EgBRF_Cov_Pi0EtaBkgULE = (TH2D*)f_Bkg.Get("h_EgBRF_AllBkgPrediction_Cov_BDT");
  h_SignalYield_vs_EgBRF_Cov_Pi0EtaBkgULE->SetName("h_SignalYield_vs_EgBRF_Cov_Pi0EtaBkgULE");
  h_SignalYield_vs_EgBRF_Cov_Pi0EtaBkgULE->SetTitle("Cov(S_{i},S_{j}) BBkg #epsilon_{ULE} error");
  h_SignalYield_vs_EgBRF_Cov_Pi0EtaBkgULE->SetDirectory(0);

  h_SignalYield_vs_EgBRF_Cov_OmegaBkgAlpha = (TH2D*)f_Bkg.Get("h_EgBRF_AllBkgPrediction_Cov_BDT");
  h_SignalYield_vs_EgBRF_Cov_OmegaBkgAlpha->SetName("h_SignalYield_vs_EgBRF_Cov_OmegaBkgAlpha");
  h_SignalYield_vs_EgBRF_Cov_OmegaBkgAlpha->SetTitle("Cov(S_{i},S_{j}) BBkg #alpha(#omega) errors");
  h_SignalYield_vs_EgBRF_Cov_OmegaBkgAlpha->SetDirectory(0);

  h_SignalYield_vs_EgBRF_Cov_etaprimeBkgAlpha = (TH2D*)f_Bkg.Get("h_EgBRF_AllBkgPrediction_Cov_BDT");
  h_SignalYield_vs_EgBRF_Cov_etaprimeBkgAlpha->SetName("h_SignalYield_vs_EgBRF_Cov_etaprimeBkgAlpha");
  h_SignalYield_vs_EgBRF_Cov_etaprimeBkgAlpha->SetTitle("Cov(S_{i},S_{j}) BBkg #alpha(#eta') errors");
  h_SignalYield_vs_EgBRF_Cov_etaprimeBkgAlpha->SetDirectory(0);

  h_SignalYield_vs_EgBRF_Cov_BkgSLAlpha = (TH2D*)f_Bkg.Get("h_EgBRF_AllBkgPrediction_Cov_BDT");
  h_SignalYield_vs_EgBRF_Cov_BkgSLAlpha->SetName("h_SignalYield_vs_EgBRF_Cov_BkgSLAlpha");
  h_SignalYield_vs_EgBRF_Cov_BkgSLAlpha->SetTitle("Cov(S_{i},S_{j}) BBkg #alpha(SL) errors");
  h_SignalYield_vs_EgBRF_Cov_BkgSLAlpha->SetDirectory(0);

  h_SignalYield_vs_EgBRF_Cov_BremsBkgMaterial = (TH2D*)f_Bkg.Get("h_EgBRF_AllBkgPrediction_Cov_BDT");
  h_SignalYield_vs_EgBRF_Cov_BremsBkgMaterial->SetName("h_SignalYield_vs_EgBRF_Cov_BremsBkgMaterial");
  h_SignalYield_vs_EgBRF_Cov_BremsBkgMaterial->SetTitle("Cov(S_{i},S_{j}) BBkg Bremss. Material syst.");
  h_SignalYield_vs_EgBRF_Cov_BremsBkgMaterial->SetDirectory(0);

  h_SignalYield_vs_EgBRF_Cov_ElecBkgTrackIneffic = (TH2D*)f_Bkg.Get("h_EgBRF_AllBkgPrediction_Cov_BDT");
  h_SignalYield_vs_EgBRF_Cov_ElecBkgTrackIneffic->SetName("h_SignalYield_vs_EgBRF_Cov_ElecBkgTrackIneffic");
  h_SignalYield_vs_EgBRF_Cov_ElecBkgTrackIneffic->SetTitle("Cov(S_{i},S_{j}) BBkg e^{#pm} Track-ineffic. syst.");
  h_SignalYield_vs_EgBRF_Cov_ElecBkgTrackIneffic->SetDirectory(0);

  h_SignalYield_vs_EgBRF_Cov_AntiNeutronsBkgAlpha = (TH2D*)f_Bkg.Get("h_EgBRF_AllBkgPrediction_Cov_BDT");
  h_SignalYield_vs_EgBRF_Cov_AntiNeutronsBkgAlpha->SetName("h_SignalYield_vs_EgBRF_Cov_AntiNeutronsBkgAlpha");
  h_SignalYield_vs_EgBRF_Cov_AntiNeutronsBkgAlpha->SetTitle("Cov(S_{i},S_{j}) BBkg #alpha(#bar{n^{0}}) errors");
  h_SignalYield_vs_EgBRF_Cov_AntiNeutronsBkgAlpha->SetDirectory(0);

  h_SignalYield_vs_EgBRF_Cov_BDT = (TH2D*)f_Bkg.Get("h_EgBRF_AllBkgPrediction_Cov_BDT");
  h_SignalYield_vs_EgBRF_Cov_BDT->SetName("h_SignalYield_vs_EgBRF_Cov_BDT");
  h_SignalYield_vs_EgBRF_Cov_BDT->SetTitle("Cov(S_{i},S_{j}) BBkg BDT syst");
  h_SignalYield_vs_EgBRF_Cov_BDT->SetDirectory(0);

  h_SignalYield_vs_EgBRF_Cov_Pi0EtaVeto_MassLineShape = (TH2D*)f_Bkg.Get("h_EgBRF_AllBkgPrediction_Cov_BDT");
  h_SignalYield_vs_EgBRF_Cov_Pi0EtaVeto_MassLineShape->SetName("h_SignalYield_vs_EgBRF_Cov_Pi0EtaVeto_MassLineShape");
  h_SignalYield_vs_EgBRF_Cov_Pi0EtaVeto_MassLineShape->SetTitle("Cov(S_{i},S_{j}) BBkg #pi^{0}/#eta veto, mass-line-shape syst");
  h_SignalYield_vs_EgBRF_Cov_Pi0EtaVeto_MassLineShape->SetDirectory(0);

  h_SignalYield_vs_EgBRF_Cov_Pi0EtaVeto_2ndGammaMult = (TH2D*)f_Bkg.Get("h_EgBRF_AllBkgPrediction_Cov_BDT");
  h_SignalYield_vs_EgBRF_Cov_Pi0EtaVeto_2ndGammaMult->SetName("h_SignalYield_vs_EgBRF_Cov_Pi0EtaVeto_2ndGammaMult");
  h_SignalYield_vs_EgBRF_Cov_Pi0EtaVeto_2ndGammaMult->SetTitle("Cov(S_{i},S_{j}) BBkg #pi^{0}/#eta veto, 2nd-#gamma Mult. syst");
  h_SignalYield_vs_EgBRF_Cov_Pi0EtaVeto_2ndGammaMult->SetDirectory(0);
  f_Bkg.Close();
  h_EgBRF_AllBkgPrediction_tot_corr->SetLineColor(2);
  h_EgBRF_AllBkgPrediction_tot_corr->SetMarkerColor(2);
  h_EgBRF_AllBkgPrediction_tot_corr->SetStats(false);

  for(int i=0;i<h_mESFitYield_vs_EgBRF_Cov_mESFixX1->GetXaxis()->GetNbins();i++) {
    for(int j=0;j<h_mESFitYield_vs_EgBRF_Cov_mESFixX1->GetYaxis()->GetNbins();j++) {
      h_mESFitYield_vs_EgBRF_Cov_Stat->SetBinContent(i+1,j+1,0.0);
      h_mESFitYield_vs_EgBRF_Cov_mESFixX1->SetBinContent(i+1,j+1,0.0);
      h_mESFitYield_vs_EgBRF_Cov_mESFixX2->SetBinContent(i+1,j+1,0.0);
      h_mESFitYield_vs_EgBRF_Cov_mESFixX3->SetBinContent(i+1,j+1,0.0);
      h_mESFitYield_vs_EgBRF_Cov_mESFix_alphaP->SetBinContent(i+1,j+1,0.0);
      h_mESFitYield_vs_EgBRF_Cov_mESFix_BBComb->SetBinContent(i+1,j+1,0.0);
      h_mESFitYield_vs_EgBRF_Cov_mESFix_slope->SetBinContent(i+1,j+1,0.0);
      h_mESFitYield_vs_EgBRF_Cov_mESFitBias->SetBinContent(i+1,j+1,0.0);

      h_SignalYield_vs_EgBRF_Cov_Stat->SetBinContent(i+1,j+1,0.0);
      h_SignalYield_vs_EgBRF_Cov_BBkgStat->SetBinContent(i+1,j+1,0.0);
      h_SignalYield_vs_EgBRF_Cov_AllHEEfficCorr->SetBinContent(i+1,j+1,0.0);
      h_SignalYield_vs_EgBRF_Cov_BrecoCorr->SetBinContent(i+1,j+1,0.0);
      h_SignalYield_vs_EgBRF_Cov_Pi0EtaBkgMCYield->SetBinContent(i+1,j+1,0.0);
      h_SignalYield_vs_EgBRF_Cov_Pi0EtaBkgRpi0Reta->SetBinContent(i+1,j+1,0.0);
      h_SignalYield_vs_EgBRF_Cov_Pi0EtaBkgRpi0CRetaC->SetBinContent(i+1,j+1,0.0);
      h_SignalYield_vs_EgBRF_Cov_Pi0EtaBkgmESFitYields->SetBinContent(i+1,j+1,0.0);
      h_SignalYield_vs_EgBRF_Cov_mESFixX1->SetBinContent(i+1,j+1,0.0);
      h_SignalYield_vs_EgBRF_Cov_mESFixX2->SetBinContent(i+1,j+1,0.0);
      h_SignalYield_vs_EgBRF_Cov_mESFixX3->SetBinContent(i+1,j+1,0.0);
      h_SignalYield_vs_EgBRF_Cov_mESFix_alphaP->SetBinContent(i+1,j+1,0.0);
      h_SignalYield_vs_EgBRF_Cov_mESFix_BBComb->SetBinContent(i+1,j+1,0.0);
      h_SignalYield_vs_EgBRF_Cov_mESFix_slope->SetBinContent(i+1,j+1,0.0);
      h_SignalYield_vs_EgBRF_Cov_mESFitBias->SetBinContent(i+1,j+1,0.0);
      h_SignalYield_vs_EgBRF_Cov_Pi0EtaBkgMggFitValSyst->SetBinContent(i+1,j+1,0.0);
      h_SignalYield_vs_EgBRF_Cov_Pi0EtaBkgLE_Base->SetBinContent(i+1,j+1,0.0);
      h_SignalYield_vs_EgBRF_Cov_Pi0EtaBkgULE->SetBinContent(i+1,j+1,0.0);
      h_SignalYield_vs_EgBRF_Cov_OmegaBkgAlpha->SetBinContent(i+1,j+1,0.0);
      h_SignalYield_vs_EgBRF_Cov_etaprimeBkgAlpha->SetBinContent(i+1,j+1,0.0);
      h_SignalYield_vs_EgBRF_Cov_BkgSLAlpha->SetBinContent(i+1,j+1,0.0);
      h_SignalYield_vs_EgBRF_Cov_BremsBkgMaterial->SetBinContent(i+1,j+1,0.0);
      h_SignalYield_vs_EgBRF_Cov_ElecBkgTrackIneffic->SetBinContent(i+1,j+1,0.0);
      h_SignalYield_vs_EgBRF_Cov_AntiNeutronsBkgAlpha->SetBinContent(i+1,j+1,0.0);
      h_SignalYield_vs_EgBRF_Cov_BDT->SetBinContent(i+1,j+1,0.0);
      h_SignalYield_vs_EgBRF_Cov_Pi0EtaVeto_MassLineShape->SetBinContent(i+1,j+1,0.0);
      h_SignalYield_vs_EgBRF_Cov_Pi0EtaVeto_2ndGammaMult->SetBinContent(i+1,j+1,0.0);
    }
  }

  TString HistName, HistTitle, command;
  char ytitle[100];

  int Ndivitions      = 605;
  double TickLength   = 0.05;
  double TitleOffSet  = 0.5;
  double TheSize  = 0.07;
  double TheSize2 = 0.10;
  double TheBottomMargin_Dn = 0.35;
  double Fraction_Pad = 0.40;

  if(fChain == 0) return;

  int _PrintFreq = 1000;
  Long64_t EventsProcessed = 0;
  Long64_t nbytes = 0, nb = 0;
  Long64_t _nentries = fChain->GetEntries();

  ///////////////////////
  //Tag-side variables://
  ///////////////////////
  TheSize = 0.05;

  double mES_threshold_comb = 5.29;
  double R_mES[2];
  R_mES[0] = 5.22;
  R_mES[1] = 5.29;
  double BinSize_mES = 0.001;
  int N_mES = int((R_mES[1]-R_mES[0])/BinSize_mES);
  TH1D h_mES("h_mES",
	     "B_{tag} m_{ES}",
	     N_mES,R_mES[0],R_mES[1]);
  h_mES.SetXTitle("m_{ES} (GeV/c^{2})");
  h_mES.GetXaxis()->CenterTitle(true);
  sprintf(ytitle,"Events / (%.3f GeV/c^{2})",h_mES.GetBinWidth(0));
  h_mES.SetYTitle(ytitle);
  h_mES.GetYaxis()->CenterTitle(true);
  h_mES.SetLineColor(4);
  h_mES.SetLineWidth(2);
  h_mES.SetMinimum(0.0);
  h_mES.GetXaxis()->SetTitleSize(TheSize);
  h_mES.GetXaxis()->SetLabelSize(TheSize);
  h_mES.GetYaxis()->SetTitleSize(TheSize);
  h_mES.GetYaxis()->SetLabelSize(TheSize);

  TH1D h_mES_BrecoTM("h_mES_BrecoTM",
		     "B_{tag} m_{ES} (Breco-TM)",
		     N_mES,R_mES[0],R_mES[1]);
  h_mES_BrecoTM.SetXTitle("m_{ES} (GeV/c^{2})");
  h_mES_BrecoTM.GetXaxis()->CenterTitle(true);
  sprintf(ytitle,"Events / (%.3f GeV/c^{2})",h_mES_BrecoTM.GetBinWidth(0));
  h_mES_BrecoTM.SetYTitle(ytitle);
  h_mES_BrecoTM.GetYaxis()->CenterTitle(true);
  h_mES_BrecoTM.SetLineColor(2);
  h_mES_BrecoTM.SetLineWidth(2);
  h_mES_BrecoTM.SetMinimum(0.0);
  h_mES_BrecoTM.GetXaxis()->SetTitleSize(TheSize);
  h_mES_BrecoTM.GetXaxis()->SetLabelSize(TheSize);
  h_mES_BrecoTM.GetYaxis()->SetTitleSize(TheSize);
  h_mES_BrecoTM.GetYaxis()->SetLabelSize(TheSize);


  ///////////////////////
  //Signal-side photon://
  ///////////////////////
  double Stupid_min_SigPhE = 0.1;

  double R_SigPhoton_EBrest_Fit[2];
  R_SigPhoton_EBrest_Fit[0] = h_EgBRF_AllBkgPrediction_tot_corr->GetXaxis()->GetXmin();
  R_SigPhoton_EBrest_Fit[1] = h_EgBRF_AllBkgPrediction_tot_corr->GetXaxis()->GetXmax();
  int Nbins_SigPhoton_EBrest_Fit = h_EgBRF_AllBkgPrediction_tot_corr->GetXaxis()->GetNbins();

  //Reconstructed photon energy:
  TH1D h_SigPhoton_EBrest_Fit("h_SigPhoton_EBrest_Fit",
			      "E*(#gamma)_{BRF} distribution",
			      Nbins_SigPhoton_EBrest_Fit,
			      R_SigPhoton_EBrest_Fit[0],
			      R_SigPhoton_EBrest_Fit[1]);
  h_SigPhoton_EBrest_Fit.SetXTitle("E*(#gamma)_{BRF} (GeV)");
  h_SigPhoton_EBrest_Fit.GetXaxis()->CenterTitle(true);
  sprintf(ytitle,"Events / (%.1f MeV)",1.0e+3*h_SigPhoton_EBrest_Fit.GetBinWidth(0));
  h_SigPhoton_EBrest_Fit.SetYTitle(ytitle);
  h_SigPhoton_EBrest_Fit.GetYaxis()->CenterTitle(true);
  h_SigPhoton_EBrest_Fit.SetLineColor(4);
  h_SigPhoton_EBrest_Fit.SetLineWidth(2);
  h_SigPhoton_EBrest_Fit.SetMinimum(Stupid_min_SigPhE);
  h_SigPhoton_EBrest_Fit.GetXaxis()->SetTitleSize(TheSize);
  h_SigPhoton_EBrest_Fit.GetXaxis()->SetLabelSize(TheSize);
  h_SigPhoton_EBrest_Fit.GetYaxis()->SetTitleSize(TheSize);
  h_SigPhoton_EBrest_Fit.GetYaxis()->SetLabelSize(TheSize);

  TH1D h_SigPhoton_EBrest_Fit_BrecoTM("h_SigPhoton_EBrest_Fit_BrecoTM",
				      "E*(#gamma)_{BRF} distribution (Breco-TM)",
				      Nbins_SigPhoton_EBrest_Fit,
				      R_SigPhoton_EBrest_Fit[0],
				      R_SigPhoton_EBrest_Fit[1]);
  h_SigPhoton_EBrest_Fit_BrecoTM.SetXTitle("E*(#gamma)_{BRF} (GeV)");
  h_SigPhoton_EBrest_Fit_BrecoTM.GetXaxis()->CenterTitle(true);
  sprintf(ytitle,"Events / (%.1f MeV)",1.0e+3*h_SigPhoton_EBrest_Fit_BrecoTM.GetBinWidth(0));
  h_SigPhoton_EBrest_Fit_BrecoTM.SetYTitle(ytitle);
  h_SigPhoton_EBrest_Fit_BrecoTM.GetYaxis()->CenterTitle(true);
  h_SigPhoton_EBrest_Fit_BrecoTM.SetLineColor(2);
  h_SigPhoton_EBrest_Fit_BrecoTM.SetLineWidth(2);
  h_SigPhoton_EBrest_Fit_BrecoTM.SetMinimum(Stupid_min_SigPhE);
  h_SigPhoton_EBrest_Fit_BrecoTM.GetXaxis()->SetTitleSize(TheSize);
  h_SigPhoton_EBrest_Fit_BrecoTM.GetXaxis()->SetLabelSize(TheSize);
  h_SigPhoton_EBrest_Fit_BrecoTM.GetYaxis()->SetTitleSize(TheSize);
  h_SigPhoton_EBrest_Fit_BrecoTM.GetYaxis()->SetLabelSize(TheSize);

  TH1D h_SigPhoton_EBrest_Fit_BrecoTM_SigGammaTM("h_SigPhoton_EBrest_Fit_BrecoTM_SigGammaTM",
						 "E*(#gamma)_{BRF} distribution (Breco-TM, signal-#gamma-TM)",
						 Nbins_SigPhoton_EBrest_Fit,
						 R_SigPhoton_EBrest_Fit[0],
						 R_SigPhoton_EBrest_Fit[1]);
  h_SigPhoton_EBrest_Fit_BrecoTM_SigGammaTM.SetXTitle("E*(#gamma)_{BRF} (GeV)");
  h_SigPhoton_EBrest_Fit_BrecoTM_SigGammaTM.GetXaxis()->CenterTitle(true);
  sprintf(ytitle,"Events / (%.1f MeV)",1.0e+3*h_SigPhoton_EBrest_Fit_BrecoTM_SigGammaTM.GetBinWidth(0));
  h_SigPhoton_EBrest_Fit_BrecoTM_SigGammaTM.SetYTitle(ytitle);
  h_SigPhoton_EBrest_Fit_BrecoTM_SigGammaTM.GetYaxis()->CenterTitle(true);
  h_SigPhoton_EBrest_Fit_BrecoTM_SigGammaTM.SetLineColor(kGreen+2);
  h_SigPhoton_EBrest_Fit_BrecoTM_SigGammaTM.SetLineWidth(2);
  h_SigPhoton_EBrest_Fit_BrecoTM_SigGammaTM.SetMinimum(Stupid_min_SigPhE);
  h_SigPhoton_EBrest_Fit_BrecoTM_SigGammaTM.GetXaxis()->SetTitleSize(TheSize);
  h_SigPhoton_EBrest_Fit_BrecoTM_SigGammaTM.GetXaxis()->SetLabelSize(TheSize);
  h_SigPhoton_EBrest_Fit_BrecoTM_SigGammaTM.GetYaxis()->SetTitleSize(TheSize);
  h_SigPhoton_EBrest_Fit_BrecoTM_SigGammaTM.GetYaxis()->SetLabelSize(TheSize);

  TH1D h_SigPhoton_EBrest_True_BrecoTM_SigGammaTM("h_SigPhoton_EBrest_True_BrecoTM_SigGammaTM",
						  "E*^{True}(#gamma)_{BRF} distribution (Breco-TM, signal-#gamma-TM)",
						  Nbins_SigPhoton_EBrest_Fit,
						  R_SigPhoton_EBrest_Fit[0],
						  R_SigPhoton_EBrest_Fit[1]);
  h_SigPhoton_EBrest_True_BrecoTM_SigGammaTM.SetXTitle("E*(#gamma)_{BRF} (GeV)");
  h_SigPhoton_EBrest_True_BrecoTM_SigGammaTM.GetXaxis()->CenterTitle(true);
  sprintf(ytitle,"Events / (%.1f MeV)",1.0e+3*h_SigPhoton_EBrest_True_BrecoTM_SigGammaTM.GetBinWidth(0));
  h_SigPhoton_EBrest_True_BrecoTM_SigGammaTM.SetYTitle(ytitle);
  h_SigPhoton_EBrest_True_BrecoTM_SigGammaTM.GetYaxis()->CenterTitle(true);
  h_SigPhoton_EBrest_True_BrecoTM_SigGammaTM.SetLineColor(6);
  h_SigPhoton_EBrest_True_BrecoTM_SigGammaTM.SetLineWidth(2);
  h_SigPhoton_EBrest_True_BrecoTM_SigGammaTM.SetMinimum(Stupid_min_SigPhE);
  h_SigPhoton_EBrest_True_BrecoTM_SigGammaTM.GetXaxis()->SetTitleSize(TheSize);
  h_SigPhoton_EBrest_True_BrecoTM_SigGammaTM.GetXaxis()->SetLabelSize(TheSize);
  h_SigPhoton_EBrest_True_BrecoTM_SigGammaTM.GetYaxis()->SetTitleSize(TheSize);
  h_SigPhoton_EBrest_True_BrecoTM_SigGammaTM.GetYaxis()->SetLabelSize(TheSize);

  TH1D h_SigPhoton_EBrest_mESFit("h_SigPhoton_EBrest_mESFit",
				 "N_{P} yield from m_{ES} fit vs E*(#gamma)_{BRF}",
				 Nbins_SigPhoton_EBrest_Fit,
				 R_SigPhoton_EBrest_Fit[0],
				 R_SigPhoton_EBrest_Fit[1]);
  h_SigPhoton_EBrest_mESFit.SetXTitle("E*(#gamma)_{BRF} (GeV)");
  h_SigPhoton_EBrest_mESFit.GetXaxis()->CenterTitle(true);
  sprintf(ytitle,"Events / (%.1f MeV)",1.0e+3*h_SigPhoton_EBrest_mESFit.GetBinWidth(0));
  h_SigPhoton_EBrest_mESFit.SetYTitle(ytitle);
  h_SigPhoton_EBrest_mESFit.GetYaxis()->CenterTitle(true);
  h_SigPhoton_EBrest_mESFit.SetLineColor(4);
  h_SigPhoton_EBrest_mESFit.SetMarkerColor(4);
  h_SigPhoton_EBrest_mESFit.SetMarkerStyle(20);
  h_SigPhoton_EBrest_mESFit.SetLineWidth(2);
  h_SigPhoton_EBrest_mESFit.SetMinimum(Stupid_min_SigPhE);
  h_SigPhoton_EBrest_mESFit.GetXaxis()->SetTitleSize(TheSize);
  h_SigPhoton_EBrest_mESFit.GetXaxis()->SetLabelSize(TheSize);
  h_SigPhoton_EBrest_mESFit.GetYaxis()->SetTitleSize(TheSize);
  h_SigPhoton_EBrest_mESFit.GetYaxis()->SetLabelSize(TheSize);
  h_SigPhoton_EBrest_mESFit.SetStats(false);

  TH1D h_SigPhoton_EBrest_mESFit_v2("h_SigPhoton_EBrest_mESFit_v2",
				    "N_{P} yield from m_{ES} fit vs E*(#gamma)_{BRF}",
				    Nbins_SigPhoton_EBrest_Fit,
				    R_SigPhoton_EBrest_Fit[0],
				    R_SigPhoton_EBrest_Fit[1]);
  h_SigPhoton_EBrest_mESFit_v2.SetXTitle("E*(#gamma)_{BRF} (GeV)");
  h_SigPhoton_EBrest_mESFit_v2.GetXaxis()->CenterTitle(true);
  sprintf(ytitle,"Events / (%.1f MeV)",1.0e+3*h_SigPhoton_EBrest_mESFit.GetBinWidth(0));
  h_SigPhoton_EBrest_mESFit_v2.SetYTitle(ytitle);
  h_SigPhoton_EBrest_mESFit_v2.GetYaxis()->CenterTitle(true);
  h_SigPhoton_EBrest_mESFit_v2.SetLineColor(4);
  h_SigPhoton_EBrest_mESFit_v2.SetMarkerColor(4);
  h_SigPhoton_EBrest_mESFit_v2.SetMarkerStyle(20);
  h_SigPhoton_EBrest_mESFit_v2.SetLineWidth(2);
  h_SigPhoton_EBrest_mESFit_v2.SetMinimum(Stupid_min_SigPhE);
  h_SigPhoton_EBrest_mESFit_v2.GetXaxis()->SetTitleSize(TheSize);
  h_SigPhoton_EBrest_mESFit_v2.GetXaxis()->SetLabelSize(TheSize);
  h_SigPhoton_EBrest_mESFit_v2.GetYaxis()->SetTitleSize(TheSize);
  h_SigPhoton_EBrest_mESFit_v2.GetYaxis()->SetLabelSize(TheSize);
  h_SigPhoton_EBrest_mESFit_v2.SetStats(false);

  TH1D h_SigPhoton_EBrest_mESFit_REFSig("h_SigPhoton_EBrest_mESFit_REFSig",
					"Peaking yields of m_{ES} fit vs E*(#gamma)_{BRF}",
					100,R_EgBRF_Sig[0],R_EgBRF_Sig[1]);
  h_SigPhoton_EBrest_mESFit_REFSig.SetXTitle("E*(#gamma)_{BRF} (GeV)");
  h_SigPhoton_EBrest_mESFit_REFSig.GetXaxis()->CenterTitle(true);
  sprintf(ytitle,"Events / (%.1f MeV)",1.0e+3*h_SigPhoton_EBrest_mESFit.GetBinWidth(0));
  h_SigPhoton_EBrest_mESFit_REFSig.SetYTitle(ytitle);
  h_SigPhoton_EBrest_mESFit_REFSig.GetYaxis()->CenterTitle(true);
  h_SigPhoton_EBrest_mESFit_REFSig.SetLineColor(1);
  h_SigPhoton_EBrest_mESFit_REFSig.SetLineWidth(2);
  h_SigPhoton_EBrest_mESFit_REFSig.SetMinimum(Stupid_min_SigPhE);
  h_SigPhoton_EBrest_mESFit_REFSig.GetXaxis()->SetTitleSize(TheSize);
  h_SigPhoton_EBrest_mESFit_REFSig.GetXaxis()->SetLabelSize(TheSize);
  h_SigPhoton_EBrest_mESFit_REFSig.GetYaxis()->SetTitleSize(TheSize);
  h_SigPhoton_EBrest_mESFit_REFSig.GetYaxis()->SetLabelSize(TheSize);
  h_SigPhoton_EBrest_mESFit_REFSig.SetStats(false);

  TH1D h_SigPhoton_EBrest_mESFit_REFBkg("h_SigPhoton_EBrest_mESFit_REFBkg",
					"Peaking yields of m_{ES} fit vs E*(#gamma)_{BRF}",
					100,R_EgBRF_Bkg[0],R_EgBRF_Bkg[1]);
  h_SigPhoton_EBrest_mESFit_REFBkg.SetXTitle("E*(#gamma)_{BRF} (GeV)");
  h_SigPhoton_EBrest_mESFit_REFBkg.GetXaxis()->CenterTitle(true);
  sprintf(ytitle,"Events / (%.1f MeV)",1.0e+3*h_SigPhoton_EBrest_mESFit.GetBinWidth(0));
  h_SigPhoton_EBrest_mESFit_REFBkg.SetYTitle(ytitle);
  h_SigPhoton_EBrest_mESFit_REFBkg.GetYaxis()->CenterTitle(true);
  h_SigPhoton_EBrest_mESFit_REFBkg.SetLineColor(1);
  h_SigPhoton_EBrest_mESFit_REFBkg.SetLineWidth(2);
  h_SigPhoton_EBrest_mESFit_REFBkg.SetMinimum(Stupid_min_SigPhE);
  h_SigPhoton_EBrest_mESFit_REFBkg.GetXaxis()->SetTitleSize(TheSize);
  h_SigPhoton_EBrest_mESFit_REFBkg.GetXaxis()->SetLabelSize(TheSize);
  h_SigPhoton_EBrest_mESFit_REFBkg.GetYaxis()->SetTitleSize(TheSize);
  h_SigPhoton_EBrest_mESFit_REFBkg.GetYaxis()->SetLabelSize(TheSize);
  h_SigPhoton_EBrest_mESFit_REFBkg.SetStats(false);

  TH1D h_SigPhoton_EBrest_mESFit_REFBkg2("h_SigPhoton_EBrest_mESFit_REFBkg2",
					 "Peaking yields of m_{ES} fit vs E*(#gamma)_{BRF}",
					 100,R_EgBRF_Bkg2[0],R_EgBRF_Bkg2[1]);
  h_SigPhoton_EBrest_mESFit_REFBkg2.SetXTitle("E*(#gamma)_{BRF} (GeV)");
  h_SigPhoton_EBrest_mESFit_REFBkg2.GetXaxis()->CenterTitle(true);
  sprintf(ytitle,"Events / (%.1f MeV)",1.0e+3*h_SigPhoton_EBrest_mESFit.GetBinWidth(0));
  h_SigPhoton_EBrest_mESFit_REFBkg2.SetYTitle(ytitle);
  h_SigPhoton_EBrest_mESFit_REFBkg2.GetYaxis()->CenterTitle(true);
  h_SigPhoton_EBrest_mESFit_REFBkg2.SetLineColor(1);
  h_SigPhoton_EBrest_mESFit_REFBkg2.SetLineWidth(2);
  h_SigPhoton_EBrest_mESFit_REFBkg2.SetMinimum(Stupid_min_SigPhE);
  h_SigPhoton_EBrest_mESFit_REFBkg2.GetXaxis()->SetTitleSize(TheSize);
  h_SigPhoton_EBrest_mESFit_REFBkg2.GetXaxis()->SetLabelSize(TheSize);
  h_SigPhoton_EBrest_mESFit_REFBkg2.GetYaxis()->SetTitleSize(TheSize);
  h_SigPhoton_EBrest_mESFit_REFBkg2.GetYaxis()->SetLabelSize(TheSize);
  h_SigPhoton_EBrest_mESFit_REFBkg2.SetStats(false);

  TH1D h_SigPhoton_EBrest_SignalYield("h_SigPhoton_EBrest_SignalYield",
				      "Signal Yield (Bkg-subtraction) vs E*(#gamma)_{BRF}",
				      Nbins_SigPhoton_EBrest_Fit,
				      R_SigPhoton_EBrest_Fit[0],
				      R_SigPhoton_EBrest_Fit[1]);
  h_SigPhoton_EBrest_SignalYield.SetXTitle("E*(#gamma)_{BRF} (GeV)");
  h_SigPhoton_EBrest_SignalYield.GetXaxis()->CenterTitle(true);
  sprintf(ytitle,"Events / (%.1f MeV)",1.0e+3*h_SigPhoton_EBrest_SignalYield.GetBinWidth(0));
  h_SigPhoton_EBrest_SignalYield.SetYTitle(ytitle);
  h_SigPhoton_EBrest_SignalYield.GetYaxis()->CenterTitle(true);
  h_SigPhoton_EBrest_SignalYield.SetLineColor(1);
  h_SigPhoton_EBrest_SignalYield.SetMarkerColor(1);
  h_SigPhoton_EBrest_SignalYield.SetMarkerStyle(20);
  h_SigPhoton_EBrest_SignalYield.SetLineWidth(2);
  h_SigPhoton_EBrest_SignalYield.SetMinimum(Stupid_min_SigPhE);
  h_SigPhoton_EBrest_SignalYield.GetXaxis()->SetTitleSize(TheSize);
  h_SigPhoton_EBrest_SignalYield.GetXaxis()->SetLabelSize(TheSize);
  h_SigPhoton_EBrest_SignalYield.GetYaxis()->SetTitleSize(TheSize);
  h_SigPhoton_EBrest_SignalYield.GetYaxis()->SetLabelSize(TheSize);
  h_SigPhoton_EBrest_SignalYield.SetStats(false);

  TH1D h_SigPhoton_EBrest_SignalYield_StatOnly("h_SigPhoton_EBrest_SignalYield_StatOnly",
					       "Signal Yield (Bkg-subtraction) vs E*(#gamma)_{BRF}",
					       Nbins_SigPhoton_EBrest_Fit,
					       R_SigPhoton_EBrest_Fit[0],
					       R_SigPhoton_EBrest_Fit[1]);
  h_SigPhoton_EBrest_SignalYield_StatOnly.SetXTitle("E*(#gamma)_{BRF} (GeV)");
  h_SigPhoton_EBrest_SignalYield_StatOnly.GetXaxis()->CenterTitle(true);
  sprintf(ytitle,"Events / (%.1f MeV)",1.0e+3*h_SigPhoton_EBrest_SignalYield_StatOnly.GetBinWidth(0));
  h_SigPhoton_EBrest_SignalYield_StatOnly.SetYTitle(ytitle);
  h_SigPhoton_EBrest_SignalYield_StatOnly.GetYaxis()->CenterTitle(true);
  h_SigPhoton_EBrest_SignalYield_StatOnly.SetLineColor(4);
  h_SigPhoton_EBrest_SignalYield_StatOnly.SetMarkerColor(4);
  h_SigPhoton_EBrest_SignalYield_StatOnly.SetMarkerStyle(20);
  h_SigPhoton_EBrest_SignalYield_StatOnly.SetLineWidth(2);
  h_SigPhoton_EBrest_SignalYield_StatOnly.SetMinimum(Stupid_min_SigPhE);
  h_SigPhoton_EBrest_SignalYield_StatOnly.GetXaxis()->SetTitleSize(TheSize);
  h_SigPhoton_EBrest_SignalYield_StatOnly.GetXaxis()->SetLabelSize(TheSize);
  h_SigPhoton_EBrest_SignalYield_StatOnly.GetYaxis()->SetTitleSize(TheSize);
  h_SigPhoton_EBrest_SignalYield_StatOnly.GetYaxis()->SetLabelSize(TheSize);
  h_SigPhoton_EBrest_SignalYield_StatOnly.SetStats(false);

  TH1D h_SigPhoton_EBrest_SignalYield_TrueSubTract("h_SigPhoton_EBrest_SignalYield_TrueSubTract",
						   "Signal Yield (Bkg-subtraction) vs E*(#gamma)_{BRF}",
						   Nbins_SigPhoton_EBrest_Fit,
						   R_SigPhoton_EBrest_Fit[0],
						   R_SigPhoton_EBrest_Fit[1]);
  h_SigPhoton_EBrest_SignalYield_TrueSubTract.SetXTitle("E*(#gamma)_{BRF} (GeV)");
  h_SigPhoton_EBrest_SignalYield_TrueSubTract.GetXaxis()->CenterTitle(true);
  sprintf(ytitle,"Events / (%.1f MeV)",1.0e+3*h_SigPhoton_EBrest_SignalYield_TrueSubTract.GetBinWidth(0));
  h_SigPhoton_EBrest_SignalYield_TrueSubTract.SetYTitle(ytitle);
  h_SigPhoton_EBrest_SignalYield_TrueSubTract.GetYaxis()->CenterTitle(true);
  h_SigPhoton_EBrest_SignalYield_TrueSubTract.SetLineColor(kOrange);
  h_SigPhoton_EBrest_SignalYield_TrueSubTract.SetLineWidth(2);
  h_SigPhoton_EBrest_SignalYield_TrueSubTract.SetMinimum(Stupid_min_SigPhE);
  h_SigPhoton_EBrest_SignalYield_TrueSubTract.GetXaxis()->SetTitleSize(TheSize);
  h_SigPhoton_EBrest_SignalYield_TrueSubTract.GetXaxis()->SetLabelSize(TheSize);
  h_SigPhoton_EBrest_SignalYield_TrueSubTract.GetYaxis()->SetTitleSize(TheSize);
  h_SigPhoton_EBrest_SignalYield_TrueSubTract.GetYaxis()->SetLabelSize(TheSize);
  h_SigPhoton_EBrest_SignalYield_TrueSubTract.SetStats(false);

  TH1D* h_mES_EgBRF[Nbins_SigPhoton_EBrest_Fit];
  for(int i=0;i<Nbins_SigPhoton_EBrest_Fit;i++) {
    double c = h_SigPhoton_EBrest_Fit.GetBinCenter(i+1);
    double w = h_SigPhoton_EBrest_Fit.GetBinWidth(i+1)*0.5;
    sprintf(ytitle,"%.1f",c-w);
    TString Range = TString(" E*(#gamma)_{BRF} in (") + TString(ytitle) + TString(",");
    sprintf(ytitle,"%.1f",c+w);
    Range        += TString(ytitle) + TString(") GeV");

    HistName  = TString("h_mES_EgBRF_") + long(i+1);
    HistTitle = TString("B_{tag} mES for ") + Range;
    h_mES_EgBRF[i] = new TH1D(HistName.Data(),
			      HistTitle.Data(),
			      N_mES,R_mES[0],R_mES[1]);
    h_mES_EgBRF[i]->SetXTitle("m_{ES} (GeV/c^{2})");
    h_mES_EgBRF[i]->GetXaxis()->CenterTitle(true);
    sprintf(ytitle,"Events / (%.3f GeV/c^{2})",h_mES_EgBRF[i]->GetBinWidth(0));
    h_mES_EgBRF[i]->SetYTitle(ytitle);
    h_mES_EgBRF[i]->GetYaxis()->CenterTitle(true);
    h_mES_EgBRF[i]->SetLineColor(4);
    h_mES_EgBRF[i]->SetLineStyle(1);
    h_mES_EgBRF[i]->SetLineWidth(2);
    h_mES_EgBRF[i]->SetMinimum(0.0);
    h_mES_EgBRF[i]->GetXaxis()->SetTitleSize(TheSize);
    h_mES_EgBRF[i]->GetXaxis()->SetLabelSize(TheSize);
    h_mES_EgBRF[i]->GetYaxis()->SetTitleSize(TheSize);
    h_mES_EgBRF[i]->GetYaxis()->SetLabelSize(TheSize);
    h_mES_EgBRF[i]->SetStats(false);
  }

  cout << endl;
  cout << "The number of entries is " << _nentries << endl;
  cout << endl;

  //Energy dependent mES fits:
  TString Name_mES = TString(outputFile) + TString("_mES.out");
  std::ofstream file_mES(Name_mES.Data());

  EventsProcessed = 0;
  nbytes = nb = 0;
  _PrintFreq = 5000;
  for (Long64_t jentry=0; jentry<_nentries;jentry++) {
    Long64_t ientry = LoadTree(jentry);
    if(ientry < 0) break;
    nb = fChain->GetEntry(jentry);
    nbytes += nb;
    
    if(!((jentry+1)%_PrintFreq)) cout << jentry+1 << " Events processed!!!" << endl;
    EventsProcessed = jentry+1;
    
    if(EventsProcessed == _nentries) {
      cout << endl;
      cout << EventsProcessed << " Events processed!!!" << endl;
      cout << endl;
    }

    double The_MES = BTagB_mES[0];
    if(_IsMC == 0 || _IsMC == 2) {
      The_MES *= 1.004;
    }

    if(The_MES < R_mES[0]) continue;
    //if(nPi0EtaCands == 0) continue;

    bool BTM   = false;
    if(_IsMC > 1) BTM   = BCand_TM_v3(jentry);

    int TheBTM = -999;
    if(BTM) TheBTM = 1;
    else    TheBTM = 0;

    bool   GammaTM     = false;
    double Egamma_star = -999.0;
    int    idx_gamma   = -999;
    if(_IsMC > 1) {
      int idx = -999;
      int idx1 = 3;
      int idx2 = 4;
      if(dauLen[idx1] >= 2 &&
	 (abs(mcLund[dauIdx[idx1]]) ==   313 ||
	  abs(mcLund[dauIdx[idx1]]) ==   323 ||
	  abs(mcLund[dauIdx[idx1]]) == 30343 ||
	  abs(mcLund[dauIdx[idx1]]) == 30353
	  ) &&
	 mcLund[dauIdx[idx1] + 1] == 22) idx = idx1;
      else if(dauLen[idx2] >= 2 &&
	      (abs(mcLund[dauIdx[idx2]]) ==   313 ||
	       abs(mcLund[dauIdx[idx2]]) ==   323 ||
	       abs(mcLund[dauIdx[idx2]]) == 30343 ||
	       abs(mcLund[dauIdx[idx2]]) == 30353
	       ) &&
	      mcLund[dauIdx[idx2] + 1] == 22) idx = idx2;
      
      //if(idx == -999) cout << jentry+1 << "  Didn't found signal B" << endl;

      if(idx != -999) {
	double gamma_energy_max = -1.0e+20;
	for(int k=0;k<dauLen[idx];k++) {
	  if(mcLund[dauIdx[idx]+k] != 22) continue;
	  if(gamma_energy_max < mcenergy[dauIdx[idx]+k]) {
	    gamma_energy_max = mcenergy[dauIdx[idx]+k];
	    idx_gamma = dauIdx[idx]+k;
	  }
	}
      }
      if(idx_gamma != -999) GammaTM = true;
      if(GammaTM) {
	//Bsig:
	double Bsig_E_CM  = mcenergyCM[idx];
	double Bsig_Px_CM = mcp3CM[idx]*TMath::Cos(mcphiCM[idx])*sqrt(1.0 - pow(mccosthCM[idx],2));
	double Bsig_Py_CM = mcp3CM[idx]*TMath::Sin(mcphiCM[idx])*sqrt(1.0 - pow(mccosthCM[idx],2));
	double Bsig_Pz_CM = mcp3CM[idx]*mccosthCM[idx];
	TVector3 beta(Bsig_Px_CM,Bsig_Py_CM,Bsig_Pz_CM);
	beta = -(1.0/Bsig_E_CM)*beta;
	
	//Gamma:
	double Gamma_E_CM  = mcenergyCM[idx_gamma];
	double Gamma_Px_CM = mcp3CM[idx_gamma]*TMath::Cos(mcphiCM[idx_gamma])*sqrt(1.0 - pow(mccosthCM[idx_gamma],2));
	double Gamma_Py_CM = mcp3CM[idx_gamma]*TMath::Sin(mcphiCM[idx_gamma])*sqrt(1.0 - pow(mccosthCM[idx_gamma],2));
	double Gamma_Pz_CM = mcp3CM[idx_gamma]*mccosthCM[idx_gamma];
	TLorentzVector p4_Gamma_Star(Gamma_Px_CM,Gamma_Py_CM,Gamma_Pz_CM,Gamma_E_CM);
	p4_Gamma_Star.Boost(beta);
	Egamma_star = p4_Gamma_Star.E();
      }
    }


    int TheGammaTM = -999;
    if(GammaTM) TheGammaTM = 1;
    else        TheGammaTM = 0;

    sprintf(ytitle,"%.10f",The_MES);
    command  = TString(ytitle) + TString("    ");
    sprintf(ytitle,"%.10f",SigPhotonenergyStar[0]);
    command += TString(ytitle) + TString("    ");
    sprintf(ytitle,"%d",TheBTM);
    command += TString(ytitle) + TString("    ");
    sprintf(ytitle,"%d",TheGammaTM);
    command += TString(ytitle);
    file_mES << command.Data() << std::endl;

    h_mES.Fill(The_MES);
    h_SigPhoton_EBrest_Fit.Fill(SigPhotonenergyStar[0]);
    if(BTM) {
      h_mES_BrecoTM.Fill(The_MES);
      h_SigPhoton_EBrest_Fit_BrecoTM.Fill(SigPhotonenergyStar[0]);
      if(GammaTM) {
	h_SigPhoton_EBrest_Fit_BrecoTM_SigGammaTM.Fill(SigPhotonenergyStar[0]);
	h_SigPhoton_EBrest_True_BrecoTM_SigGammaTM.Fill(Egamma_star);
      }
    }

    for(int kkk=0;kkk<Nbins_SigPhoton_EBrest_Fit;kkk++) {
      double c = h_SigPhoton_EBrest_Fit.GetBinCenter(kkk+1);
      double w = h_SigPhoton_EBrest_Fit.GetBinWidth(kkk+1)*0.5;
      if(SigPhotonenergyStar[0] >= c-w && 
	 SigPhotonenergyStar[0] <  c+w) {
	h_mES_EgBRF[kkk]->Fill(The_MES);
	break;
      }
    }

  }
  cout << endl;
  cout << endl;

  file_mES.close();

  double porcent = 0.10;
  double Maximum;
  double Minimum;
  TString EPSName;
  EPSName = TString(outputFile) + TString(".eps");   //rmaria it was .eps
  TString EPSNameO = EPSName + TString("[");
  TString EPSNameC = EPSName + TString("]");

  cout << endl;
  cout << "Saving plotted histograms inside " << EPSName.Data() << " file" << endl;
  cout << endl;

  TCanvas* c1 = new TCanvas("c1","c1");
  c1->SetFillColor(10);
  c1->SetFrameFillColor(10);
  c1->SetTickx(1);
  c1->SetTicky(1);
  c1->SetLeftMargin(0.15);
  c1->SetBottomMargin(0.15);
  c1->SetRightMargin(0.20);
  
  c1->Print(EPSNameO.Data());

  c1->Clear();
  h_mES.Draw();
  if(_IsMC > 1) h_mES_BrecoTM.Draw("same");
  c1->Print(EPSName.Data());

  c1->Clear();
  h_SigPhoton_EBrest_Fit.Draw();
  if(_IsMC > 1) {
    h_SigPhoton_EBrest_Fit_BrecoTM.Draw("same");
    h_SigPhoton_EBrest_Fit_BrecoTM_SigGammaTM.Draw("same");
    h_SigPhoton_EBrest_True_BrecoTM_SigGammaTM.Draw("same");
  }
  c1->Print(EPSName.Data());

  if(_IsMC > 1) {
    c1->Clear();
    h_SigPhoton_EBrest_Fit_BrecoTM.Draw();
    h_SigPhoton_EBrest_Fit_BrecoTM_SigGammaTM.Draw("same");
    h_SigPhoton_EBrest_True_BrecoTM_SigGammaTM.Draw("same");
    c1->Print(EPSName.Data());

    c1->Clear();
    Maximum = TMath::Max(h_SigPhoton_EBrest_Fit_BrecoTM_SigGammaTM.GetMaximum(),
			 h_SigPhoton_EBrest_True_BrecoTM_SigGammaTM.GetMaximum());
    porcent = 0.10;
    h_SigPhoton_EBrest_Fit_BrecoTM_SigGammaTM.SetMaximum(Maximum*(1.0+porcent));
    h_SigPhoton_EBrest_True_BrecoTM_SigGammaTM.SetMaximum(Maximum*(1.0+porcent));
    h_SigPhoton_EBrest_Fit_BrecoTM_SigGammaTM.Draw();
    h_SigPhoton_EBrest_True_BrecoTM_SigGammaTM.Draw("same");
    c1->Print(EPSName.Data());
  }


  int Nbins_mES_fit = 50;

  double mES_mean_peak     = 5.27997e+00;
  double R_mES_mean_peak[2];
  R_mES_mean_peak[0] = 5.26;
  R_mES_mean_peak[1] = 5.30;
  double mES_sigma_peakR     = 2.22514e-03;
  double R_mES_sigma_peakR[2];
  R_mES_sigma_peakR[0] = 0.0001;
  R_mES_sigma_peakR[1] = 0.010;
  double mES_sigma_peakL     = 3.10365e-03;
  double R_mES_sigma_peakL[2];
  R_mES_sigma_peakL[0] = 0.0001;
  R_mES_sigma_peakL[1] = 0.010;
  double mES_alpha_peakL     = 1.78835;
  double R_mES_alpha_peakL[2];
  R_mES_alpha_peakL[0] =  0.01;
  R_mES_alpha_peakL[1] = 10.00;
  double mES_order_peakR = 10.0;
  double mES_order_peakL = 1.45;
  double R_mES_order_peakL[2];
  R_mES_order_peakL[0] = 0.5;
  R_mES_order_peakL[1] = 20.0;

  double mES_rho_comb = 0.889488;
  double R_mES_rho_comb[2];
  R_mES_rho_comb[0] = 0.0;
  R_mES_rho_comb[1] = 1.0;

  double mES_alpha_comb = 0.267391;
  double R_mES_alpha_comb[2];
  R_mES_alpha_comb[0] = -3.0;
  R_mES_alpha_comb[1] =  1.0;

  double mES_Slope_comb   = -1.84988e+01;
  double R_mES_Slope_comb[2];
  R_mES_Slope_comb[0] = -200.0;
  R_mES_Slope_comb[1] =    0.0;
  
  double R_B_TM[2];
  R_B_TM[0] = -1.0;
  R_B_TM[1] =  2.0;

  double R_G_TM[2];
  R_G_TM[0] = -1.0;
  R_G_TM[1] =  2.0;

  double Argus_frac = 0.80;
  double R_Argus_frac[2];
  R_Argus_frac[0] = 0.0;
  R_Argus_frac[1] = 1.0;

  double R_mES_peak[2];
  R_mES_peak[0] = 5.260;
  R_mES_peak[1] = 5.290;
  
  double R_mES_fit[2];
  R_mES_fit[0] = 5.22;
  R_mES_fit[1] = mES_threshold_comb;
  
  RooRealVar mES("mES","m_{ES}",R_mES_fit[0],R_mES_fit[1],"GeV/c^{2}");
  RooRealVar Egamma_BRF("Egamma_BRF","Egamma_BRF",R_SigPhoton_EBrest_Fit[0],R_SigPhoton_EBrest_Fit[1],"GeV");
  RooRealVar B_TM("B_TM","B_TM",R_B_TM[0],R_B_TM[1]);
  RooRealVar G_TM("G_TM","G_TM",R_G_TM[0],R_G_TM[1]);

  RooArgSet dataVars_mES;
  dataVars_mES.add(RooArgSet(mES,Egamma_BRF,B_TM,G_TM));

  RooRealVar mESmeanPeak("mESmeanPeak",    "mESmeanPeak",  mES_mean_peak,  R_mES_mean_peak[0],  R_mES_mean_peak[1]);
  RooRealVar mESsigmaPeakL("mESsigmaPeakL","mESsigmaPeakL",mES_sigma_peakL,R_mES_sigma_peakL[0],R_mES_sigma_peakL[1]);
  RooRealVar mESsigmaPeakR("mESsigmaPeakR","mESsigmaPeakR",mES_sigma_peakR,R_mES_sigma_peakR[0],R_mES_sigma_peakR[1]);
  RooRealVar mESalphaPeakL("mESalphaPeakL","mESalphaPeakL",mES_alpha_peakL,R_mES_alpha_peakL[0],R_mES_alpha_peakL[1]);
  RooRealVar mESorderPeakL("mESorderPeakL","mESorderPeakL",mES_order_peakL,R_mES_order_peakL[0],R_mES_order_peakL[1]);
  RooRealVar mESorderPeakR("mESorderPeakR","mESorderPeakR",mES_order_peakR);
  RooDoubleCBShape mESpdfPeak("mESpdfPeak", "mESpdfPeak",
			      mES,
			      mESmeanPeak,
			      mESsigmaPeakL,
			      mESsigmaPeakR,
			      mESalphaPeakL,
			      mESalphaPeakL,
			      mESorderPeakL,
			      mESorderPeakR);

  mESmeanPeak.setConstant(false);
  mESsigmaPeakL.setConstant(false);
  mESsigmaPeakR.setConstant(false);
  mESalphaPeakL.setConstant(true);
  mESorderPeakL.setConstant(true);
  mESorderPeakR.setConstant(true);


  RooRealVar mESrhoComb("mESrhoComb","mESrhoComb",mES_rho_comb,R_mES_rho_comb[0],R_mES_rho_comb[1]);
  RooRealVar mESalphaComb("mESalphaComb","mESalphaComb",mES_alpha_comb,R_mES_alpha_comb[0],R_mES_alpha_comb[1]);
  RooRealVar mESm0Comb("mESm0Comb","mESm0Comb",5.27,5.25,5.28);
  RooRealVar mESmtComb("mESmtComb","mESm0Comb",mES_threshold_comb);
  RooPolyAndExpDecayV3 mESpdfComb1("mESpdfComb1", "mESpdfComb1",
				   mES,
				   mESrhoComb,
				   mESalphaComb,
				   mESm0Comb,
				   mESmtComb);
  mESrhoComb.setConstant(false);
  mESalphaComb.setConstant(false);
  mESm0Comb.setConstant(true);
  mESmtComb.setConstant(true);

  RooRealVar mESslopeComb("mESslopeComb","mESslopeComb",mES_Slope_comb,R_mES_Slope_comb[0],R_mES_Slope_comb[1]);
  RooRealVar mESthresComb("mESthresComb","mESthresComb",5.28992);
  RooArgusBG mESpdfComb2("mESpdfComb2",    "mESpdfComb2",
			 mES,
			 mESthresComb,
			 mESslopeComb);
  mESslopeComb.setConstant(false);
  mESthresComb.setConstant(true);

  RooRealVar mESFracComb("mESFracComb","mESFracComb",Argus_frac,R_Argus_frac[0],R_Argus_frac[1]);
  RooAddPdf mESpdfComb("mESpdfComb",
		       "mESpdfComb",
		       RooArgList(mESpdfComb2,mESpdfComb1),
		       RooArgList(mESFracComb));
  mESFracComb.setConstant(false);


  RooRealVar nmESPeak("nmESPeak","nmESPeak",0.0,5000);
  RooRealVar nmESComb("nmESComb","nmESComb",0.0,5000);
  RooAddPdf PdfTotal_mES("PdfTotal_mES",
                         "PdfTotal_mES",
                         RooArgList(mESpdfPeak,mESpdfComb),
                         RooArgList(nmESPeak,nmESComb));

  TH1D h_mES_clone(h_mES);
  h_mES_clone.SetYTitle("Pulls");
  h_mES_clone.SetStats(false);
  h_mES_clone.SetMaximum( 6.5);
  h_mES_clone.SetMinimum(-6.5);
  h_mES_clone.GetXaxis()->SetTitleSize(TheSize2);
  h_mES_clone.GetXaxis()->SetLabelSize(TheSize2);
  h_mES_clone.GetYaxis()->SetTitleSize(TheSize2);
  h_mES_clone.GetYaxis()->SetLabelSize(TheSize2);
  h_mES_clone.GetYaxis()->SetTitleOffset(TitleOffSet);
  h_mES_clone.GetXaxis()->SetTickLength(TickLength);

  TLine* l_mES_pull_0 = new TLine(h_mES_clone.GetXaxis()->GetXmin(),0.0,
				  h_mES_clone.GetXaxis()->GetXmax(),0.0);
  l_mES_pull_0->SetLineColor(1);
  l_mES_pull_0->SetLineWidth(1);
  l_mES_pull_0->SetLineStyle(1);
  TLine* l_mES_pull_1 = new TLine(h_mES_clone.GetXaxis()->GetXmin(),2.0,
				  h_mES_clone.GetXaxis()->GetXmax(),2.0);
  l_mES_pull_1->SetLineColor(2);
  l_mES_pull_1->SetLineWidth(2);
  l_mES_pull_1->SetLineStyle(2);
  TLine* l_mES_pull_2 = new TLine(h_mES_clone.GetXaxis()->GetXmin(),-2.0,
				  h_mES_clone.GetXaxis()->GetXmax(),-2.0);
  l_mES_pull_2->SetLineColor(2);
  l_mES_pull_2->SetLineWidth(2);
  l_mES_pull_2->SetLineStyle(2);

  double NumEntries;
  double x1,y1,x2,y2,a,b;
  int idx_1,idx_2;
  double Npeak_estimated,Ncomb_estimated,Alpha_tmp;

  double mES_mean_peak_global;
  double mES_sigma_peakL_global;
  double mES_sigma_peakR_global;
  double mES_alpha_peakL_global;
  double mES_order_peakL_global;
  double mES_order_peakR_global;
  double mES_rho_comb_global;
  double mES_alpha_comb_global;
  double mES_slope_comb_global;
  double mES_threshold_comb_global;
  double mES_Frac_comb_global;

  RooDataSet *data_BrecoTM     = NULL;
  RooDataSet *data_Breco_nonTM = NULL;

  if(DoGlobalFit) {
    Egamma_BRF.setMin(R_SigPhoton_EBrest_Fit[0]);
    Egamma_BRF.setMax(R_SigPhoton_EBrest_Fit[1]);
    if(_IsMC > 1) {
      B_TM.setMin(0.5);
      B_TM.setMax(1.5);
      data_BrecoTM = RooDataSet::read(Name_mES.Data(),dataVars_mES,"Q");
      
      B_TM.setMin(-0.5);
      B_TM.setMax(0.5);
      data_Breco_nonTM = RooDataSet::read(Name_mES.Data(),dataVars_mES,"Q");
    }
    
    B_TM.setMin(-0.5);
    B_TM.setMax(1.5);
    RooDataSet *data_mES_Breco = RooDataSet::read(Name_mES.Data(),dataVars_mES,"Q");
    NumEntries = data_mES_Breco->numEntries();
    
    //Estimate of the peaking and combinatoric yields
    idx_1 = h_mES.FindBin(R_mES_peak[0]);
    idx_2 = h_mES.FindBin(R_mES_peak[1]);
    x1 = h_mES.GetBinCenter(idx_1);
    y1 = h_mES.GetBinContent(idx_1);
    x2 = h_mES.GetBinCenter(idx_2);
    y2 = h_mES.GetBinContent(idx_2);
    a  = (y2-y1)/(x2-x1);
    b  = y2 - a*x2;
    
    Npeak_estimated = 0;
    for(int lll=0;lll<h_mES.GetXaxis()->GetNbins();lll++) {
      if(lll+1 < idx_1 || lll+1 > idx_2) continue;
      double x    = h_mES.GetBinCenter(lll+1);
      double comb = a*x + b;
      double tot  = h_mES.GetBinContent(lll+1);
      
      Npeak_estimated += tot - comb;
    }
    Ncomb_estimated = h_mES.Integral() - Npeak_estimated;
    
    Alpha_tmp = Npeak_estimated/(Npeak_estimated+Ncomb_estimated);
    Npeak_estimated = NumEntries*Alpha_tmp;
    Ncomb_estimated = NumEntries*(1.0 - Alpha_tmp);
    
    cout << endl;
    cout << "Initial estimate of Peaking/Combinatoric yeild:" << endl;
    cout << "Total events       = " << h_mES.GetEntries() << endl;
    cout << "Peaking yield      = " << Npeak_estimated << endl;
    cout << "Combinatoric yield = " << Ncomb_estimated << endl;
    cout << endl;
    
    //nmESPeak.setMin(-NumEntries);
    nmESPeak.setMin(0.0);
    nmESPeak.setMax(5.0*NumEntries);
    nmESPeak.setVal(Npeak_estimated);
    nmESComb.setMin(0.0);
    nmESComb.setMax(5.0*NumEntries);
    nmESComb.setVal(Ncomb_estimated);

    mESFracComb.setVal(Argus_frac);
    mESmeanPeak.setVal(mES_mean_peak);
    mESsigmaPeakL.setVal(mES_sigma_peakL);
    mESsigmaPeakR.setVal(mES_sigma_peakR);
    mESalphaPeakL.setVal(mES_alpha_peakL);
    mESrhoComb.setVal(mES_rho_comb);
    mESalphaComb.setVal(mES_alpha_comb);
    mESslopeComb.setVal(mES_Slope_comb);
    mESmeanPeak.setConstant(false);
    mESsigmaPeakL.setConstant(false);
    mESsigmaPeakR.setConstant(false);
    mESalphaPeakL.setConstant(true);
    mESrhoComb.setConstant(true);
    mESalphaComb.setConstant(true);
    mESslopeComb.setConstant(true);
    mESthresComb.setConstant(true);
    mESFracComb.setConstant(false);

    nmESPeak.setConstant(false);
    nmESComb.setConstant(false);
    resFit_mES_all = PdfTotal_mES.fitTo(*data_mES_Breco,Save());
    resFit_mES_all->SetName("resFit_mES_all");

    mES_mean_peak_global      = mESmeanPeak.getVal();
    mES_sigma_peakL_global    = mESsigmaPeakL.getVal();
    mES_sigma_peakR_global    = mESsigmaPeakR.getVal();
    mES_alpha_peakL_global    = mESalphaPeakL.getVal();
    mES_rho_comb_global       = mESrhoComb.getVal();
    mES_alpha_comb_global     = mESalphaComb.getVal();
    mES_slope_comb_global     = mESslopeComb.getVal();
    mES_threshold_comb_global = mESthresComb.getVal();
    mES_Frac_comb_global      = mESFracComb.getVal();
    mES_order_peakL_global    = mES_order_peakL;
    mES_order_peakR_global    = mES_order_peakR;

    RooPlot *mESframe_all_Breco = mES.frame(R_mES[0],R_mES[1],
					    Nbins_mES_fit);
    mESframe_all_Breco->SetTitle("m_{ES} fit for all B-reco events");
    data_mES_Breco->plotOn(mESframe_all_Breco);
    PdfTotal_mES.plotOn(mESframe_all_Breco,LineColor(kRed));
    RooHist* hpull_mES_all_Breco  = mESframe_all_Breco->pullHist();
    PdfTotal_mES.plotOn(mESframe_all_Breco,Components(mESpdfComb),LineStyle(kDashed),LineColor(kGreen+2));
    PdfTotal_mES.plotOn(mESframe_all_Breco,Components(mESpdfPeak),LineStyle(kDashed),LineColor(kBlue));
    if(_IsMC > 1) {
      data_BrecoTM->plotOn(mESframe_all_Breco,LineColor(kBlue),MarkerColor(kBlue));
      data_Breco_nonTM->plotOn(mESframe_all_Breco,LineColor(kGreen+2),MarkerColor(kGreen+2));
    }
    
    {
      c1->Clear();
      TPad all_mES_pad_Up("all_mES_pad_Up","m_{ES} Pad all UP",0.0,Fraction_Pad,1.0,1.0);
      all_mES_pad_Up.SetFillColor(10);
      all_mES_pad_Up.SetFrameFillColor(10);
      all_mES_pad_Up.SetTickx(1);
      all_mES_pad_Up.SetTicky(1);
      all_mES_pad_Up.SetBottomMargin(0.0);
      all_mES_pad_Up.SetLeftMargin(0.15);
      all_mES_pad_Up.Draw();
      all_mES_pad_Up.cd();
      mESframe_all_Breco->GetXaxis()->CenterTitle(true);
      mESframe_all_Breco->GetYaxis()->CenterTitle(true);
      mESframe_all_Breco->GetXaxis()->SetTitleSize(TheSize);
      mESframe_all_Breco->GetXaxis()->SetLabelSize(TheSize);
      mESframe_all_Breco->GetYaxis()->SetTitleSize(TheSize);
      mESframe_all_Breco->GetYaxis()->SetLabelSize(TheSize);
      mESframe_all_Breco->Draw();
      c1->cd();
      TPad all_mES_pad_Dn("all_mES_pad_Dn","m_{ES} Pad all DOWN",0.0,0.0,1.0,Fraction_Pad);
      all_mES_pad_Dn.SetFillColor(10);
      all_mES_pad_Dn.SetFrameFillColor(10);
      all_mES_pad_Dn.SetTickx(1);
      all_mES_pad_Dn.SetTicky(1);
      all_mES_pad_Dn.SetTopMargin(0.0);
      all_mES_pad_Dn.SetBottomMargin(TheBottomMargin_Dn);
      all_mES_pad_Dn.SetLeftMargin(0.15);
      all_mES_pad_Dn.Draw();
      all_mES_pad_Dn.cd();
      h_mES_clone.Draw("AXIS");
      l_mES_pull_0->Draw();
      l_mES_pull_1->Draw();
      l_mES_pull_2->Draw();
      hpull_mES_all_Breco->Draw("P");
      c1->Print(EPSName.Data());
    }

    if(_IsMC > 0) {
      RooPlot *mESframe_all_Breco_SigTest = mES.frame(R_mES[0],R_mES[1],
						      Nbins_mES_fit);
      mESframe_all_Breco_SigTest->SetTitle("m_{ES} fit for all B-reco events (Signal events)");
      data_BrecoTM->plotOn(mESframe_all_Breco_SigTest,LineColor(kBlue),MarkerColor(kBlue));
      mESpdfPeak.plotOn(mESframe_all_Breco,LineColor(kBlue),LineStyle(kDashed));
#if 0
      //PdfTotal_mES.plotOn(mESframe_all_Breco_SigTest,Components(mESpdfPeak),LineStyle(kDashed),LineColor(kBlue));
      RooHist* hpull_mES_all_Breco_SigTest = mESframe_all_Breco_SigTest->pullHist();
      hpull_mES_all_Breco_SigTest->SetLineColor(4);
      hpull_mES_all_Breco_SigTest->SetMarkerColor(4);
#endif 
      {
	c1->Clear();
	TPad all_mES_pad_Up("all_mES_pad_Up","m_{ES} Pad all UP",0.0,Fraction_Pad,1.0,1.0);
	all_mES_pad_Up.SetFillColor(10);
	all_mES_pad_Up.SetFrameFillColor(10);
	all_mES_pad_Up.SetTickx(1);
	all_mES_pad_Up.SetTicky(1);
	all_mES_pad_Up.SetBottomMargin(0.0);
	all_mES_pad_Up.SetLeftMargin(0.15);
	all_mES_pad_Up.Draw();
	all_mES_pad_Up.cd();
	mESframe_all_Breco_SigTest->GetXaxis()->CenterTitle(true);
	mESframe_all_Breco_SigTest->GetYaxis()->CenterTitle(true);
	mESframe_all_Breco_SigTest->GetXaxis()->SetTitleSize(TheSize);
	mESframe_all_Breco_SigTest->GetXaxis()->SetLabelSize(TheSize);
	mESframe_all_Breco_SigTest->GetYaxis()->SetTitleSize(TheSize);
	mESframe_all_Breco_SigTest->GetYaxis()->SetLabelSize(TheSize);
	mESframe_all_Breco_SigTest->Draw();
	c1->cd();
	TPad all_mES_pad_Dn("all_mES_pad_Dn","m_{ES} Pad all DOWN",0.0,0.0,1.0,Fraction_Pad);
	all_mES_pad_Dn.SetFillColor(10);
	all_mES_pad_Dn.SetFrameFillColor(10);
	all_mES_pad_Dn.SetTickx(1);
	all_mES_pad_Dn.SetTicky(1);
	all_mES_pad_Dn.SetTopMargin(0.0);
	all_mES_pad_Dn.SetBottomMargin(TheBottomMargin_Dn);
	all_mES_pad_Dn.SetLeftMargin(0.15);
	all_mES_pad_Dn.Draw();
	all_mES_pad_Dn.cd();
	h_mES_clone.Draw("AXIS");
	l_mES_pull_0->Draw();
	l_mES_pull_1->Draw();
	l_mES_pull_2->Draw();
	//hpull_mES_all_Breco_SigTest->Draw("P");
	c1->Print(EPSName.Data());
      }
    }

  }
  else {
    mES_mean_peak_global      = ((RooRealVar*)(resFit_mES_all->floatParsFinal()).find("mESmeanPeak"))->getVal();
    mES_sigma_peakL_global    = ((RooRealVar*)(resFit_mES_all->floatParsFinal()).find("mESsigmaPeakL"))->getVal();
    mES_sigma_peakR_global    = ((RooRealVar*)(resFit_mES_all->floatParsFinal()).find("mESsigmaPeakR"))->getVal();
    mES_alpha_peakL_global    = ((RooRealVar*)(resFit_mES_all->constPars()).find("mESalphaPeakL"))->getVal();
    mES_order_peakL_global    = (((RooRealVar*)(resFit_mES_all->constPars()).find("mESorderPeakL")))->getVal();
    mES_order_peakR_global    = (((RooRealVar*)(resFit_mES_all->constPars()).find("mESorderPeakR")))->getVal();

    mES_rho_comb_global       = ((RooRealVar*)(resFit_mES_all->constPars()).find("mESrhoComb"))->getVal();
    mES_alpha_comb_global     = (((RooRealVar*)(resFit_mES_all->constPars()).find("mESalphaComb")))->getVal();

    mES_slope_comb_global     = (((RooRealVar*)(resFit_mES_all->constPars()).find("mESslopeComb")))->getVal();
    mES_threshold_comb_global = (((RooRealVar*)(resFit_mES_all->constPars()).find("mESthresComb")))->getVal();

    mES_Frac_comb_global      = (((RooRealVar*)(resFit_mES_all->floatParsFinal()).find("mESFracComb")))->getVal();

    Egamma_BRF.setMin(R_SigPhoton_EBrest_Fit[0]);
    Egamma_BRF.setMax(R_SigPhoton_EBrest_Fit[1]);
    if(_IsMC > 1) {
      B_TM.setMin(0.5);
      B_TM.setMax(1.5);
      data_BrecoTM = RooDataSet::read(Name_mES.Data(),dataVars_mES,"Q");
      
      B_TM.setMin(-0.5);
      B_TM.setMax(0.5);
      data_Breco_nonTM = RooDataSet::read(Name_mES.Data(),dataVars_mES,"Q");
    }
    
    B_TM.setMin(-0.5);
    B_TM.setMax(1.5);
    RooDataSet *data_mES_Breco = RooDataSet::read(Name_mES.Data(),dataVars_mES,"Q");
    NumEntries = data_mES_Breco->numEntries();
    
    //Estimate of the peaking and combinatoric yields
    idx_1 = h_mES.FindBin(R_mES_peak[0]);
    idx_2 = h_mES.FindBin(R_mES_peak[1]);
    x1 = h_mES.GetBinCenter(idx_1);
    y1 = h_mES.GetBinContent(idx_1);
    x2 = h_mES.GetBinCenter(idx_2);
    y2 = h_mES.GetBinContent(idx_2);
    a  = (y2-y1)/(x2-x1);
    b  = y2 - a*x2;
    
    Npeak_estimated = 0;
    for(int lll=0;lll<h_mES.GetXaxis()->GetNbins();lll++) {
      if(lll+1 < idx_1 || lll+1 > idx_2) continue;
      double x    = h_mES.GetBinCenter(lll+1);
      double comb = a*x + b;
      double tot  = h_mES.GetBinContent(lll+1);
      
      Npeak_estimated += tot - comb;
    }
    Ncomb_estimated = h_mES.Integral() - Npeak_estimated;
    
    Alpha_tmp = Npeak_estimated/(Npeak_estimated+Ncomb_estimated);
    Npeak_estimated = NumEntries*Alpha_tmp;
    Ncomb_estimated = NumEntries*(1.0 - Alpha_tmp);
    
    cout << endl;
    cout << "Initial estimate of Peaking/Combinatoric yeild:" << endl;
    cout << "Total events       = " << h_mES.GetEntries() << endl;
    cout << "Peaking yield      = " << Npeak_estimated << endl;
    cout << "Combinatoric yield = " << Ncomb_estimated << endl;
    cout << endl;
    
    //nmESPeak.setMin(-NumEntries);
    nmESPeak.setMin(0.0);
    nmESPeak.setMax(5.0*NumEntries);
    nmESPeak.setVal(Npeak_estimated);
    nmESComb.setMin(0.0);
    nmESComb.setMax(5.0*NumEntries);
    nmESComb.setVal(Ncomb_estimated);

    mESFracComb.setVal(mES_Frac_comb_global);
    mESmeanPeak.setVal(mES_mean_peak_global);
    mESsigmaPeakL.setVal(mES_sigma_peakL_global);
    mESsigmaPeakR.setVal(mES_sigma_peakR_global);
    mESalphaPeakL.setVal(mES_alpha_peakL_global);
    mESrhoComb.setVal(mES_rho_comb_global);
    mESalphaComb.setVal(mES_alpha_comb_global);
    mESslopeComb.setVal(mES_slope_comb_global);
    mESmeanPeak.setConstant(true);
    mESsigmaPeakL.setConstant(true);
    mESsigmaPeakR.setConstant(true);
    mESalphaPeakL.setConstant(true);
    mESrhoComb.setConstant(true);
    mESalphaComb.setConstant(true);
    mESslopeComb.setConstant(true);
    mESthresComb.setConstant(true);
    mESFracComb.setConstant(false);

    nmESPeak.setConstant(false);
    nmESComb.setConstant(false);
    PdfTotal_mES.fitTo(*data_mES_Breco);

    RooPlot *mESframe_all_Breco = mES.frame(R_mES[0],R_mES[1],
					    Nbins_mES_fit);
    mESframe_all_Breco->SetTitle("m_{ES} fit for all B-reco events");
    data_mES_Breco->plotOn(mESframe_all_Breco);
    PdfTotal_mES.plotOn(mESframe_all_Breco,LineColor(kRed));
    RooHist* hpull_mES_all_Breco  = mESframe_all_Breco->pullHist();
    PdfTotal_mES.plotOn(mESframe_all_Breco,Components(mESpdfComb),LineStyle(kDashed),LineColor(kGreen+2));
    PdfTotal_mES.plotOn(mESframe_all_Breco,Components(mESpdfPeak),LineStyle(kDashed),LineColor(kBlue));
    if(_IsMC > 1) {
      data_BrecoTM->plotOn(mESframe_all_Breco,LineColor(kBlue),MarkerColor(kBlue));
      data_Breco_nonTM->plotOn(mESframe_all_Breco,LineColor(kGreen+2),MarkerColor(kGreen+2));
    }
    
    {
      c1->Clear();
      TPad all_mES_pad_Up("all_mES_pad_Up","m_{ES} Pad all UP",0.0,Fraction_Pad,1.0,1.0);
      all_mES_pad_Up.SetFillColor(10);
      all_mES_pad_Up.SetFrameFillColor(10);
      all_mES_pad_Up.SetTickx(1);
      all_mES_pad_Up.SetTicky(1);
      all_mES_pad_Up.SetBottomMargin(0.0);
      all_mES_pad_Up.SetLeftMargin(0.15);
      all_mES_pad_Up.Draw();
      all_mES_pad_Up.cd();
      mESframe_all_Breco->GetXaxis()->CenterTitle(true);
      mESframe_all_Breco->GetYaxis()->CenterTitle(true);
      mESframe_all_Breco->GetXaxis()->SetTitleSize(TheSize);
      mESframe_all_Breco->GetXaxis()->SetLabelSize(TheSize);
      mESframe_all_Breco->GetYaxis()->SetTitleSize(TheSize);
      mESframe_all_Breco->GetYaxis()->SetLabelSize(TheSize);
      mESframe_all_Breco->Draw();
      c1->cd();
      TPad all_mES_pad_Dn("all_mES_pad_Dn","m_{ES} Pad all DOWN",0.0,0.0,1.0,Fraction_Pad);
      all_mES_pad_Dn.SetFillColor(10);
      all_mES_pad_Dn.SetFrameFillColor(10);
      all_mES_pad_Dn.SetTickx(1);
      all_mES_pad_Dn.SetTicky(1);
      all_mES_pad_Dn.SetTopMargin(0.0);
      all_mES_pad_Dn.SetBottomMargin(TheBottomMargin_Dn);
      all_mES_pad_Dn.SetLeftMargin(0.15);
      all_mES_pad_Dn.Draw();
      all_mES_pad_Dn.cd();
      h_mES_clone.Draw("AXIS");
      l_mES_pull_0->Draw();
      l_mES_pull_1->Draw();
      l_mES_pull_2->Draw();
      hpull_mES_all_Breco->Draw("P");
      c1->Print(EPSName.Data());
    }

  }

  const int Nvariables = 3;
  TString Variable_Names[Nvariables];
  Variable_Names[0] = TString("mESmeanPeak");
  Variable_Names[1] = TString("mESsigmaPeakL");
  Variable_Names[2] = TString("mESsigmaPeakR");

  double Cvals[Nvariables];
  double Evals[Nvariables];
  TMatrixTSym<double> _CovMatrix;
  _CovMatrix.ResizeTo(TMatrixTSym<double>(3));
  for(int i=0;i<Nvariables;i++) {
    RooRealVar* The_var = ((RooRealVar*)(resFit_mES_all->floatParsFinal()).find(Variable_Names[i].Data()));
    Cvals[i] = The_var->getVal();
    Evals[i] = The_var->getError();
    for(int j=0;j<Nvariables;j++) {
      _CovMatrix(i,j) = resFit_mES_all->correlation(Variable_Names[i].Data(),
						    Variable_Names[j].Data());
    }
  }
  cout << "Corr-matrix:" << endl;
  _CovMatrix.Print();
  for(int i=0;i<Nvariables;i++) {
    cout << Variable_Names[i].Data() << " = " 
	 << Cvals[i] << " +/- " << Evals[i] << endl;
    for(int j=0;j<Nvariables;j++) {
      _CovMatrix(i,j) = _CovMatrix(i,j)*Evals[i]*Evals[j];
    }
  }
  cout << "Cov-matrix:" << endl;
  _CovMatrix.Print();

  TMatrixTSym<double> _AMatrix(_CovMatrix);
  _AMatrix.Invert();

  TMatrixDSymEigen _AMatrix_eigen(_AMatrix);
  TVectorD _Eigen_values      = _AMatrix_eigen.GetEigenValues();
  TMatrixD _Eigen_vectors     = _AMatrix_eigen.GetEigenVectors();
  TMatrixD _Eigen_vectors_Inv(_Eigen_vectors);
  _Eigen_vectors_Inv.Invert();
  cout << "Eigen-values:" << endl;
  _Eigen_values.Print();
  cout << "Eigen-vectors:" << endl;
  _Eigen_vectors.Print();
  cout << "Eigen-vectors(Inv):" << endl;
  _Eigen_vectors_Inv.Print();

  TMatrixTSym<double> _Tmp_unit;
  _Tmp_unit.ResizeTo(TMatrixTSym<double>(3));
  for(int i=0;i<Nvariables;i++) {
    for(int j=0;j<Nvariables;j++) {
      _Tmp_unit(i,j) = 0.0;
      for(int k=0;k<Nvariables;k++) {
	_Tmp_unit(i,j) += _Eigen_vectors_Inv(i,k)*_Eigen_vectors(k,j);
      }
    }
  }
  cout << "Unit-matrix:" << endl;
  _Tmp_unit.Print();
  
  TMatrixTSym<double> _Diag;
  _Diag.ResizeTo(TMatrixTSym<double>(3));
  for(int i=0;i<Nvariables;i++) {
    for(int j=0;j<Nvariables;j++) {
      _Diag(i,j) = 0.0;
      for(int k=0;k<Nvariables;k++) {
	for(int l=0;l<Nvariables;l++) {
	  _Diag(i,j) += _Eigen_vectors_Inv(i,k)*_AMatrix(k,l)*_Eigen_vectors(l,j);
	}
      }
    }
  }
  cout << "Diag-matrix:" << endl;
  _Diag.Print();

  double Cvals_p[Nvariables];
  double Evals_p[Nvariables];
  for(int i=0;i<Nvariables;i++) {
    Cvals_p[i] = 0.0;
    for(int k=0;k<Nvariables;k++) {
      Cvals_p[i] += _Eigen_vectors_Inv(i,k)*Cvals[k];
    }
    Evals_p[i] = 1.0/sqrt(_Diag(i,i));
  }
  for(int i=0;i<Nvariables;i++) {
    cout << Variable_Names[i] << "_prime = " 
	 << Cvals_p[i] << " +/- " << Evals_p[i]
	 << endl;
  }
  cout << endl;

  double Cvals_p_tmp[Nvariables];
  double Cvals_tmp[Nvariables];
  std::vector<double> _mESmeanPeak_values_p;
  std::vector<double> _mESsigmaPeakL_values_p;
  std::vector<double> _mESsigmaPeakR_values_p;
  std::vector<double> _mESmeanPeak_values_m;
  std::vector<double> _mESsigmaPeakL_values_m;
  std::vector<double> _mESsigmaPeakR_values_m;
  _mESmeanPeak_values_p.clear();
  _mESsigmaPeakL_values_p.clear();
  _mESsigmaPeakR_values_p.clear();
  _mESmeanPeak_values_m.clear();
  _mESsigmaPeakL_values_m.clear();
  _mESsigmaPeakR_values_m.clear();
  for(int i=0;i<Nvariables;i++) {
    for(int j=0;j<Nvariables;j++) {
      Cvals_p_tmp[j] = Cvals_p[j];
    }
    Cvals_p_tmp[i] = Cvals_p[i] + Evals_p[i];
    for(int j=0;j<Nvariables;j++) {
      Cvals_tmp[j] = 0.0;
      for(int k=0;k<Nvariables;k++) {
	Cvals_tmp[j] += _Eigen_vectors(j,k)*Cvals_p_tmp[k];
      }
    }
    _mESmeanPeak_values_p.push_back(Cvals_tmp[0]);
    _mESsigmaPeakL_values_p.push_back(Cvals_tmp[1]);
    _mESsigmaPeakR_values_p.push_back(Cvals_tmp[2]);

    for(int j=0;j<Nvariables;j++) {
      Cvals_p_tmp[j] = Cvals_p[j];
    }
    Cvals_p_tmp[i] = Cvals_p[i] - Evals_p[i];
    for(int j=0;j<Nvariables;j++) {
      Cvals_tmp[j] = 0.0;
      for(int k=0;k<Nvariables;k++) {
	Cvals_tmp[j] += _Eigen_vectors(j,k)*Cvals_p_tmp[k];
      }
    }
    _mESmeanPeak_values_m.push_back(Cvals_tmp[0]);
    _mESsigmaPeakL_values_m.push_back(Cvals_tmp[1]);
    _mESsigmaPeakR_values_m.push_back(Cvals_tmp[2]);
  }
  
  cout << Variable_Names[0].Data() << " = " << Cvals[0] << endl;
  for(int i=0;i<Nvariables;i++) {
    cout << "("
	 << _mESmeanPeak_values_m[i] << ","
	 << _mESmeanPeak_values_p[i] << ")"
	 << endl;
  }
  cout << endl;
  cout << Variable_Names[1].Data() << " = " << Cvals[1] << endl;
  for(int i=0;i<Nvariables;i++) {
    cout << "("
	 << _mESsigmaPeakL_values_m[i] << ","
	 << _mESsigmaPeakL_values_p[i] << ")"
	 << endl;
  }
  cout << endl;
  cout << Variable_Names[2].Data() << " = " << Cvals[2] << endl;
  for(int i=0;i<Nvariables;i++) {
    cout << "("
	 << _mESsigmaPeakR_values_m[i] << ","
	 << _mESsigmaPeakR_values_p[i] << ")"
	 << endl;
  }
  cout << endl;

  double mES_alpha_PeakL_err = 0.0601726;
  double mES_rho_Comb_err    = 0.0100277;
  double mES_alpha_Comb_err  = 0.0138802;
  //Slope value = -18.5451 +/- 0.569265(stat) +/- 1.8441(EgBRF-var) +/- 2.55114(syst-Data-MC-diff)
  double mES_slope_Comb_err  = sqrt(pow(0.569265,2) + pow(1.8441,2) + pow(2.55114,2));
  double mES_FitBias         = 3.0/100.0;

  const int N_par_config = 7;
  TString ParValNames_p[N_par_config];
  TString ParValNames_m[N_par_config];
  std::vector<double> _mES_mean_peak_configs_p;
  std::vector<double> _mES_sigma_peakR_configs_p;
  std::vector<double> _mES_sigma_peakL_configs_p;
  std::vector<double> _mES_alpha_peakL_configs_p;
  std::vector<double> _mES_rho_comb_configs_p;
  std::vector<double> _mES_alpha_comb_configs_p;
  std::vector<double> _mES_Slope_comb_configs_p;
  _mES_mean_peak_configs_p.clear();
  _mES_sigma_peakR_configs_p.clear();
  _mES_sigma_peakL_configs_p.clear();
  _mES_alpha_peakL_configs_p.clear();
  _mES_rho_comb_configs_p.clear();
  _mES_alpha_comb_configs_p.clear();
  _mES_Slope_comb_configs_p.clear();
  std::vector<double> _mES_mean_peak_configs_m;
  std::vector<double> _mES_sigma_peakR_configs_m;
  std::vector<double> _mES_sigma_peakL_configs_m;
  std::vector<double> _mES_alpha_peakL_configs_m;
  std::vector<double> _mES_rho_comb_configs_m;
  std::vector<double> _mES_alpha_comb_configs_m;
  std::vector<double> _mES_Slope_comb_configs_m;
  _mES_mean_peak_configs_m.clear();
  _mES_sigma_peakR_configs_m.clear();
  _mES_sigma_peakL_configs_m.clear();
  _mES_alpha_peakL_configs_m.clear();
  _mES_rho_comb_configs_m.clear();
  _mES_alpha_comb_configs_m.clear();
  _mES_Slope_comb_configs_m.clear();

  ParValNames_p[0] = TString("m_{ES} (#mu,#sigma_{L,R}) X_{1} + #sigma(X_{1})");
  _mES_mean_peak_configs_p.push_back(_mESmeanPeak_values_p[0]);
  _mES_sigma_peakL_configs_p.push_back(_mESsigmaPeakL_values_p[0]);
  _mES_sigma_peakR_configs_p.push_back(_mESsigmaPeakR_values_p[0]);
  _mES_alpha_peakL_configs_p.push_back(mES_alpha_peakL);
  _mES_rho_comb_configs_p.push_back(mES_rho_comb);
  _mES_alpha_comb_configs_p.push_back(mES_alpha_comb);
  _mES_Slope_comb_configs_p.push_back(mES_Slope_comb);
  ParValNames_m[0] = TString("m_{ES} (#mu,#sigma_{L,R}) X_{1} - #sigma(X_{1})");
  _mES_mean_peak_configs_m.push_back(_mESmeanPeak_values_m[0]);
  _mES_sigma_peakR_configs_m.push_back(_mESsigmaPeakR_values_m[0]);
  _mES_sigma_peakL_configs_m.push_back(_mESsigmaPeakL_values_m[0]);
  _mES_alpha_peakL_configs_m.push_back(mES_alpha_peakL);
  _mES_rho_comb_configs_m.push_back(mES_rho_comb);
  _mES_alpha_comb_configs_m.push_back(mES_alpha_comb);
  _mES_Slope_comb_configs_m.push_back(mES_Slope_comb);

  ParValNames_p[1] = TString("m_{ES} (#mu,#sigma_{L,R}) X_{2} + #sigma(X_{2})");
  _mES_mean_peak_configs_p.push_back(_mESmeanPeak_values_p[1]);
  _mES_sigma_peakL_configs_p.push_back(_mESsigmaPeakL_values_p[1]);
  _mES_sigma_peakR_configs_p.push_back(_mESsigmaPeakR_values_p[1]);
  _mES_alpha_peakL_configs_p.push_back(mES_alpha_peakL);
  _mES_rho_comb_configs_p.push_back(mES_rho_comb);
  _mES_alpha_comb_configs_p.push_back(mES_alpha_comb);
  _mES_Slope_comb_configs_p.push_back(mES_Slope_comb);
  ParValNames_m[1] = TString("m_{ES} (#mu,#sigma_{L,R}) X_{2} - #sigma(X_{2})");
  _mES_mean_peak_configs_m.push_back(_mESmeanPeak_values_m[1]);
  _mES_sigma_peakR_configs_m.push_back(_mESsigmaPeakR_values_m[1]);
  _mES_sigma_peakL_configs_m.push_back(_mESsigmaPeakL_values_m[1]);
  _mES_alpha_peakL_configs_m.push_back(mES_alpha_peakL);
  _mES_rho_comb_configs_m.push_back(mES_rho_comb);
  _mES_alpha_comb_configs_m.push_back(mES_alpha_comb);
  _mES_Slope_comb_configs_m.push_back(mES_Slope_comb);

  ParValNames_p[2] = TString("m_{ES} (#mu,#sigma_{L,R}) X_{3} + #sigma(X_{3})");
  _mES_mean_peak_configs_p.push_back(_mESmeanPeak_values_p[2]);
  _mES_sigma_peakL_configs_p.push_back(_mESsigmaPeakL_values_p[2]);
  _mES_sigma_peakR_configs_p.push_back(_mESsigmaPeakR_values_p[2]);
  _mES_alpha_peakL_configs_p.push_back(mES_alpha_peakL);
  _mES_rho_comb_configs_p.push_back(mES_rho_comb);
  _mES_alpha_comb_configs_p.push_back(mES_alpha_comb);
  _mES_Slope_comb_configs_p.push_back(mES_Slope_comb);
  ParValNames_m[2] = TString("m_{ES} (#mu,#sigma_{L,R}) X_{3} - #sigma(X_{3})");
  _mES_mean_peak_configs_m.push_back(_mESmeanPeak_values_m[2]);
  _mES_sigma_peakR_configs_m.push_back(_mESsigmaPeakR_values_m[2]);
  _mES_sigma_peakL_configs_m.push_back(_mESsigmaPeakL_values_m[2]);
  _mES_alpha_peakL_configs_m.push_back(mES_alpha_peakL);
  _mES_rho_comb_configs_m.push_back(mES_rho_comb);
  _mES_alpha_comb_configs_m.push_back(mES_alpha_comb);
  _mES_Slope_comb_configs_m.push_back(mES_Slope_comb);

  ParValNames_p[3] = TString("m_{ES} #alpha^{P} + #sigma(#alpha^{P})");
  _mES_mean_peak_configs_p.push_back(Cvals[0]);
  _mES_sigma_peakL_configs_p.push_back(Cvals[1]);
  _mES_sigma_peakR_configs_p.push_back(Cvals[2]);
  _mES_alpha_peakL_configs_p.push_back(mES_alpha_peakL_global + mES_alpha_PeakL_err);
  _mES_rho_comb_configs_p.push_back(mES_rho_comb_global);
  _mES_alpha_comb_configs_p.push_back(mES_alpha_comb_global);
  _mES_Slope_comb_configs_p.push_back(mES_slope_comb_global);
  ParValNames_m[3] = TString("m_{ES} #alpha^{P} - #sigma(#alpha^{P})");
  _mES_mean_peak_configs_m.push_back(Cvals[0]);
  _mES_sigma_peakL_configs_m.push_back(Cvals[1]);
  _mES_sigma_peakR_configs_m.push_back(Cvals[2]);
  _mES_alpha_peakL_configs_m.push_back(mES_alpha_peakL_global - mES_alpha_PeakL_err);
  _mES_rho_comb_configs_m.push_back(mES_rho_comb_global);
  _mES_alpha_comb_configs_m.push_back(mES_alpha_comb_global);
  _mES_Slope_comb_configs_m.push_back(mES_slope_comb_global);

  ParValNames_p[4] = TString("m_{ES} #rho^{C} + #sigma(#rho^{C})");
  _mES_mean_peak_configs_p.push_back(Cvals[0]);
  _mES_sigma_peakL_configs_p.push_back(Cvals[1]);
  _mES_sigma_peakR_configs_p.push_back(Cvals[2]);
  _mES_alpha_peakL_configs_p.push_back(mES_alpha_peakL_global);
  _mES_rho_comb_configs_p.push_back(mES_rho_comb_global + mES_rho_Comb_err);
  _mES_alpha_comb_configs_p.push_back(mES_alpha_comb_global);
  _mES_Slope_comb_configs_p.push_back(mES_slope_comb_global);
  ParValNames_m[4] = TString("m_{ES} #rho^{C} - #sigma(#rho^{C})");
  _mES_mean_peak_configs_m.push_back(Cvals[0]);
  _mES_sigma_peakL_configs_m.push_back(Cvals[1]);
  _mES_sigma_peakR_configs_m.push_back(Cvals[2]);
  _mES_alpha_peakL_configs_m.push_back(mES_alpha_peakL_global);
  _mES_rho_comb_configs_m.push_back(mES_rho_comb_global - mES_rho_Comb_err);
  _mES_alpha_comb_configs_m.push_back(mES_alpha_comb_global);
  _mES_Slope_comb_configs_m.push_back(mES_slope_comb_global);

  ParValNames_p[5] = TString("m_{ES} #alpha^{C} + #sigma(#alpha^{C})");
  _mES_mean_peak_configs_p.push_back(Cvals[0]);
  _mES_sigma_peakL_configs_p.push_back(Cvals[1]);
  _mES_sigma_peakR_configs_p.push_back(Cvals[2]);
  _mES_alpha_peakL_configs_p.push_back(mES_alpha_peakL_global);
  _mES_rho_comb_configs_p.push_back(mES_rho_comb_global);
  _mES_alpha_comb_configs_p.push_back(mES_alpha_comb_global + mES_alpha_Comb_err);
  _mES_Slope_comb_configs_p.push_back(mES_slope_comb_global);
  ParValNames_m[5] = TString("m_{ES} #alpha^{C} - #sigma(#alpha^{C})");
  _mES_mean_peak_configs_m.push_back(Cvals[0]);
  _mES_sigma_peakL_configs_m.push_back(Cvals[1]);
  _mES_sigma_peakR_configs_m.push_back(Cvals[2]);
  _mES_alpha_peakL_configs_m.push_back(mES_alpha_peakL_global);
  _mES_rho_comb_configs_m.push_back(mES_rho_comb_global);
  _mES_alpha_comb_configs_m.push_back(mES_alpha_comb_global - mES_alpha_Comb_err);
  _mES_Slope_comb_configs_m.push_back(mES_slope_comb_global);

  ParValNames_p[6] = TString("m_{ES} #psi^{C} + #sigma(#psi^{C})");
  _mES_mean_peak_configs_p.push_back(Cvals[0]);
  _mES_sigma_peakL_configs_p.push_back(Cvals[1]);
  _mES_sigma_peakR_configs_p.push_back(Cvals[2]);
  _mES_alpha_peakL_configs_p.push_back(mES_alpha_peakL_global);
  _mES_rho_comb_configs_p.push_back(mES_rho_comb_global);
  _mES_alpha_comb_configs_p.push_back(mES_alpha_comb_global);
  _mES_Slope_comb_configs_p.push_back(mES_slope_comb_global + mES_slope_Comb_err);
  ParValNames_m[6] = TString("m_{ES} #psi^{C} - #sigma(#psi^{C})");
  _mES_mean_peak_configs_m.push_back(Cvals[0]);
  _mES_sigma_peakL_configs_m.push_back(Cvals[1]);
  _mES_sigma_peakR_configs_m.push_back(Cvals[2]);
  _mES_alpha_peakL_configs_m.push_back(mES_alpha_peakL_global);
  _mES_rho_comb_configs_m.push_back(mES_rho_comb_global);
  _mES_alpha_comb_configs_m.push_back(mES_alpha_comb_global);
  _mES_Slope_comb_configs_m.push_back(mES_slope_comb_global - mES_slope_Comb_err);


  TH1D* h_SigPhoton_EBrest_mESFit_Syst_p[N_par_config];
  TH1D* h_SigPhoton_EBrest_mESFit_Syst_m[N_par_config];
  TH1D* h_SigPhoton_EBrest_mESFit_Syst_Delta[N_par_config];
  for(int i=0;i<N_par_config;i++) {
    HistName  = TString("h_SigPhoton_EBrest_mESFit_Syst_p_") + long(i+1);
    HistTitle = TString("N_{P} yield from m_{ES} fit vs E*(#gamma)_{BRF}, ") + ParValNames_p[i];
    h_SigPhoton_EBrest_mESFit_Syst_p[i] = new TH1D(HistName.Data(),
						   HistTitle.Data(),
						   Nbins_SigPhoton_EBrest_Fit,
						   R_SigPhoton_EBrest_Fit[0],
						   R_SigPhoton_EBrest_Fit[1]);
    h_SigPhoton_EBrest_mESFit_Syst_p[i]->SetXTitle("E*(#gamma)_{BRF} (GeV)");
    h_SigPhoton_EBrest_mESFit_Syst_p[i]->GetXaxis()->CenterTitle(true);
    sprintf(ytitle,"Events / (%.1f MeV)",1.0e+3*h_SigPhoton_EBrest_mESFit.GetBinWidth(0));
    h_SigPhoton_EBrest_mESFit_Syst_p[i]->SetYTitle(ytitle);
    h_SigPhoton_EBrest_mESFit_Syst_p[i]->GetYaxis()->CenterTitle(true);
    h_SigPhoton_EBrest_mESFit_Syst_p[i]->SetLineColor(2);
    h_SigPhoton_EBrest_mESFit_Syst_p[i]->SetLineWidth(2);
    h_SigPhoton_EBrest_mESFit_Syst_p[i]->SetMinimum(Stupid_min_SigPhE);
    h_SigPhoton_EBrest_mESFit_Syst_p[i]->GetXaxis()->SetTitleSize(TheSize);
    h_SigPhoton_EBrest_mESFit_Syst_p[i]->GetXaxis()->SetLabelSize(TheSize);
    h_SigPhoton_EBrest_mESFit_Syst_p[i]->GetYaxis()->SetTitleSize(TheSize);
    h_SigPhoton_EBrest_mESFit_Syst_p[i]->GetYaxis()->SetLabelSize(TheSize);
    h_SigPhoton_EBrest_mESFit_Syst_p[i]->SetStats(false);

    HistName  = TString("h_SigPhoton_EBrest_mESFit_Syst_m_") + long(i+1);
    HistTitle = TString("N_{P} yield from m_{ES} fit vs E*(#gamma)_{BRF}, ") + ParValNames_m[i];
    h_SigPhoton_EBrest_mESFit_Syst_m[i] = new TH1D(HistName.Data(),
						   HistTitle.Data(),
						   Nbins_SigPhoton_EBrest_Fit,
						   R_SigPhoton_EBrest_Fit[0],
						   R_SigPhoton_EBrest_Fit[1]);
    h_SigPhoton_EBrest_mESFit_Syst_m[i]->SetXTitle("E*(#gamma)_{BRF} (GeV)");
    h_SigPhoton_EBrest_mESFit_Syst_m[i]->GetXaxis()->CenterTitle(true);
    sprintf(ytitle,"Events / (%.1f MeV)",1.0e+3*h_SigPhoton_EBrest_mESFit.GetBinWidth(0));
    h_SigPhoton_EBrest_mESFit_Syst_m[i]->SetYTitle(ytitle);
    h_SigPhoton_EBrest_mESFit_Syst_m[i]->GetYaxis()->CenterTitle(true);
    h_SigPhoton_EBrest_mESFit_Syst_m[i]->SetLineColor(kGreen+2);
    h_SigPhoton_EBrest_mESFit_Syst_m[i]->SetLineWidth(2);
    h_SigPhoton_EBrest_mESFit_Syst_m[i]->SetMinimum(Stupid_min_SigPhE);
    h_SigPhoton_EBrest_mESFit_Syst_m[i]->GetXaxis()->SetTitleSize(TheSize);
    h_SigPhoton_EBrest_mESFit_Syst_m[i]->GetXaxis()->SetLabelSize(TheSize);
    h_SigPhoton_EBrest_mESFit_Syst_m[i]->GetYaxis()->SetTitleSize(TheSize);
    h_SigPhoton_EBrest_mESFit_Syst_m[i]->GetYaxis()->SetLabelSize(TheSize);
    h_SigPhoton_EBrest_mESFit_Syst_m[i]->SetStats(false);

    HistName  = TString("h_SigPhoton_EBrest_mESFit_Syst_Delta_") + long(i+1);
    HistTitle = TString("#DeltaN_{P} yield from m_{ES} fit vs E*(#gamma)_{BRF}, ") + ParValNames_p[i];
    h_SigPhoton_EBrest_mESFit_Syst_Delta[i] = new TH1D(HistName.Data(),
						       HistTitle.Data(),
						       Nbins_SigPhoton_EBrest_Fit,
						       R_SigPhoton_EBrest_Fit[0],
						       R_SigPhoton_EBrest_Fit[1]);
    h_SigPhoton_EBrest_mESFit_Syst_Delta[i]->SetXTitle("E*(#gamma)_{BRF} (GeV)");
    h_SigPhoton_EBrest_mESFit_Syst_Delta[i]->GetXaxis()->CenterTitle(true);
    sprintf(ytitle,"Events / (%.1f MeV)",1.0e+3*h_SigPhoton_EBrest_mESFit.GetBinWidth(0));
    h_SigPhoton_EBrest_mESFit_Syst_Delta[i]->SetYTitle(ytitle);
    h_SigPhoton_EBrest_mESFit_Syst_Delta[i]->GetYaxis()->CenterTitle(true);
    h_SigPhoton_EBrest_mESFit_Syst_Delta[i]->SetLineColor(4);
    h_SigPhoton_EBrest_mESFit_Syst_Delta[i]->SetLineWidth(2);
    h_SigPhoton_EBrest_mESFit_Syst_Delta[i]->GetXaxis()->SetTitleSize(TheSize);
    h_SigPhoton_EBrest_mESFit_Syst_Delta[i]->GetXaxis()->SetLabelSize(TheSize);
    h_SigPhoton_EBrest_mESFit_Syst_Delta[i]->GetYaxis()->SetTitleSize(TheSize);
    h_SigPhoton_EBrest_mESFit_Syst_Delta[i]->GetYaxis()->SetLabelSize(TheSize);
    h_SigPhoton_EBrest_mESFit_Syst_Delta[i]->SetStats(false);
  }

  //mES different EgBRF bins:
  RooDataSet *data_mES_BrecoTM_EgBRF[Nbins_SigPhoton_EBrest_Fit];
  RooDataSet *data_mES_Breco_nonTM_EgBRF[Nbins_SigPhoton_EBrest_Fit];
  RooDataSet *data_mES_EgBRF[Nbins_SigPhoton_EBrest_Fit];
  RooPlot *Mggframe_mES_EgBRF[Nbins_SigPhoton_EBrest_Fit];
  RooHist *hpull_mES_EgBRF[Nbins_SigPhoton_EBrest_Fit];
  for(int i=0;i<Nbins_SigPhoton_EBrest_Fit;i++) {
    double c = h_SigPhoton_EBrest_Fit.GetBinCenter(i+1);
    double w = h_SigPhoton_EBrest_Fit.GetBinWidth(i+1)*0.5;
    sprintf(ytitle,"%.1f",c-w);
    TString Range = TString("(") + TString(ytitle) + TString(",");
    sprintf(ytitle,"%.1f",c+w);
    Range        += TString(ytitle) + TString(") GeV");
    
    Egamma_BRF.setMin(c-w);
    Egamma_BRF.setMax(c+w);

    if(_IsMC > 1) {
      B_TM.setMin(0.5);
      B_TM.setMax(1.5);
      data_mES_BrecoTM_EgBRF[i] = RooDataSet::read(Name_mES.Data(),dataVars_mES,"Q");
      
      B_TM.setMin(-0.5);
      B_TM.setMax(0.5);
      data_mES_Breco_nonTM_EgBRF[i] = RooDataSet::read(Name_mES.Data(),dataVars_mES,"Q");
    }

    B_TM.setMin(-0.5);
    B_TM.setMax(1.5);
    data_mES_EgBRF[i] = RooDataSet::read(Name_mES.Data(),dataVars_mES,"Q");
    NumEntries = data_mES_EgBRF[i]->numEntries();

    //Estimate of the peaking and combinatoric yields
    idx_1 = h_mES_EgBRF[i]->FindBin(R_mES_peak[0]);
    idx_2 = h_mES_EgBRF[i]->FindBin(R_mES_peak[1]);
    x1 = h_mES_EgBRF[i]->GetBinCenter(idx_1);
    y1 = h_mES_EgBRF[i]->GetBinContent(idx_1);
    x2 = h_mES_EgBRF[i]->GetBinCenter(idx_2);
    y2 = h_mES_EgBRF[i]->GetBinContent(idx_2);
    a  = (y2-y1)/(x2-x1);
    b  = y2 - a*x2;

    Npeak_estimated = 0;
    for(int lll=0;lll<h_mES_EgBRF[i]->GetXaxis()->GetNbins();lll++) {
      if(lll+1 < idx_1 || lll+1 > idx_2) continue;
      double x    = h_mES_EgBRF[i]->GetBinCenter(lll+1);
      double comb = a*x + b;
      double tot  = h_mES_EgBRF[i]->GetBinContent(lll+1);
      
      Npeak_estimated += tot - comb;
    }
    Ncomb_estimated = h_mES_EgBRF[i]->Integral() - Npeak_estimated;
  
    Alpha_tmp = Npeak_estimated/(Npeak_estimated+Ncomb_estimated);
    Npeak_estimated = NumEntries*Alpha_tmp;
    Ncomb_estimated = NumEntries*(1.0 - Alpha_tmp);
    
    cout << endl;
    cout << "Initial estimate of Peaking/Combinatoric yeild:" << endl;
    cout << "Total events       = " << h_mES_EgBRF[i]->GetEntries() << endl;
    cout << "Peaking yield      = " << Npeak_estimated << endl;
    cout << "Combinatoric yield = " << Ncomb_estimated << endl;
    cout << endl;

    nmESPeak.setMin(0.0);
    nmESPeak.setMax(5.0*NumEntries);
    nmESPeak.setVal(Npeak_estimated);
    nmESComb.setMin(0.0);
    nmESComb.setMax(5.0*NumEntries);
    nmESComb.setVal(Ncomb_estimated);

    mESFracComb.setVal(Argus_frac);
    mESmeanPeak.setVal(mES_mean_peak_global);
    mESsigmaPeakL.setVal(mES_sigma_peakL_global);
    mESsigmaPeakR.setVal(mES_sigma_peakR_global);
    mESalphaPeakL.setVal(mES_alpha_peakL_global);
    mESrhoComb.setVal(mES_rho_comb_global);
    mESalphaComb.setVal(mES_alpha_comb_global);
    mESslopeComb.setVal(mES_slope_comb_global);
    mESmeanPeak.setConstant(true);
    mESsigmaPeakL.setConstant(true);
    mESsigmaPeakR.setConstant(true);
    mESalphaPeakL.setConstant(true);
    mESrhoComb.setConstant(true);
    mESalphaComb.setConstant(true);
    mESslopeComb.setConstant(true);
    mESthresComb.setConstant(true);
    mESFracComb.setConstant(false);
    nmESPeak.setConstant(false);
    nmESComb.setConstant(false);
    PdfTotal_mES.fitTo(*data_mES_EgBRF[i]);

    Mggframe_mES_EgBRF[i] = mES.frame(R_mES[0],R_mES[1],
				      Nbins_mES_fit);
    HistTitle = TString("m_{ES} fit for E*(#gamma)_{BRF} in ") + Range;
    Mggframe_mES_EgBRF[i]->SetTitle(HistTitle.Data());
    data_mES_EgBRF[i]->plotOn(Mggframe_mES_EgBRF[i]);
    PdfTotal_mES.plotOn(Mggframe_mES_EgBRF[i],LineColor(kRed));
    hpull_mES_EgBRF[i] = Mggframe_mES_EgBRF[i]->pullHist();
    PdfTotal_mES.plotOn(Mggframe_mES_EgBRF[i],Components(mESpdfComb),LineStyle(kDashed),LineColor(kGreen+2));
    PdfTotal_mES.plotOn(Mggframe_mES_EgBRF[i],Components(mESpdfPeak),LineStyle(kDashed),LineColor(kBlue));
    if(_IsMC > 1) {
      data_mES_BrecoTM_EgBRF[i]->plotOn(Mggframe_mES_EgBRF[i],LineColor(kBlue),MarkerColor(kBlue));
      data_mES_Breco_nonTM_EgBRF[i]->plotOn(Mggframe_mES_EgBRF[i],LineColor(kGreen+2),MarkerColor(kGreen+2));
    }

    h_SigPhoton_EBrest_mESFit.SetBinContent(i+1,nmESPeak.getVal());
    h_SigPhoton_EBrest_mESFit.SetBinError(i+1,nmESPeak.getError());
    h_SigPhoton_EBrest_mESFit_v2.SetBinContent(i+1,nmESPeak.getVal());

    {
      c1->Clear();
      TPad all_mES_pad_Up("all_mES_pad_Up","m_{ES} Pad all UP",0.0,Fraction_Pad,1.0,1.0);
      all_mES_pad_Up.SetFillColor(10);
      all_mES_pad_Up.SetFrameFillColor(10);
      all_mES_pad_Up.SetTickx(1);
      all_mES_pad_Up.SetTicky(1);
      all_mES_pad_Up.SetBottomMargin(0.0);
      all_mES_pad_Up.SetLeftMargin(0.15);
      all_mES_pad_Up.Draw();
      all_mES_pad_Up.cd();
      Mggframe_mES_EgBRF[i]->GetXaxis()->CenterTitle(true);
      Mggframe_mES_EgBRF[i]->GetYaxis()->CenterTitle(true);
      Mggframe_mES_EgBRF[i]->GetXaxis()->SetTitleSize(TheSize);
      Mggframe_mES_EgBRF[i]->GetXaxis()->SetLabelSize(TheSize);
      Mggframe_mES_EgBRF[i]->GetYaxis()->SetTitleSize(TheSize);
      Mggframe_mES_EgBRF[i]->GetYaxis()->SetLabelSize(TheSize);
      Mggframe_mES_EgBRF[i]->Draw();
      c1->cd();
      TPad all_mES_pad_Dn("all_mES_pad_Dn","m_{ES} Pad all DOWN",0.0,0.0,1.0,Fraction_Pad);
      all_mES_pad_Dn.SetFillColor(10);
      all_mES_pad_Dn.SetFrameFillColor(10);
      all_mES_pad_Dn.SetTickx(1);
      all_mES_pad_Dn.SetTicky(1);
      all_mES_pad_Dn.SetTopMargin(0.0);
      all_mES_pad_Dn.SetBottomMargin(TheBottomMargin_Dn);
      all_mES_pad_Dn.SetLeftMargin(0.15);
      all_mES_pad_Dn.Draw();
      all_mES_pad_Dn.cd();
      h_mES_clone.Draw("AXIS");
      l_mES_pull_0->Draw();
      l_mES_pull_1->Draw();
      l_mES_pull_2->Draw();
      hpull_mES_EgBRF[i]->Draw("P");
      c1->Print(EPSName.Data());
    }

    for(int kkk=0;kkk<N_par_config;kkk++) {
      mESmeanPeak.setVal(_mES_mean_peak_configs_p[kkk]);
      mESsigmaPeakL.setVal(_mES_sigma_peakL_configs_p[kkk]);
      mESsigmaPeakR.setVal(_mES_sigma_peakR_configs_p[kkk]);
      mESalphaPeakL.setVal(_mES_alpha_peakL_configs_p[kkk]);
      mESrhoComb.setVal(_mES_rho_comb_configs_p[kkk]);
      mESalphaComb.setVal(_mES_alpha_comb_configs_p[kkk]);
      mESslopeComb.setVal(_mES_Slope_comb_configs_p[kkk]);
      mESmeanPeak.setConstant(true);
      mESsigmaPeakL.setConstant(true);
      mESsigmaPeakR.setConstant(true);
      mESalphaPeakL.setConstant(true);
      mESrhoComb.setConstant(true);
      mESalphaComb.setConstant(true);
      mESslopeComb.setConstant(true);
      mESthresComb.setConstant(true);
      mESFracComb.setConstant(false);
      nmESPeak.setConstant(false);
      nmESComb.setConstant(false);
      PdfTotal_mES.fitTo(*data_mES_EgBRF[i]);
      h_SigPhoton_EBrest_mESFit_Syst_p[kkk]->SetBinContent(i+1,nmESPeak.getVal());

      mESmeanPeak.setVal(_mES_mean_peak_configs_m[kkk]);
      mESsigmaPeakL.setVal(_mES_sigma_peakL_configs_m[kkk]);
      mESsigmaPeakR.setVal(_mES_sigma_peakR_configs_m[kkk]);
      mESalphaPeakL.setVal(_mES_alpha_peakL_configs_m[kkk]);
      mESrhoComb.setVal(_mES_rho_comb_configs_m[kkk]);
      mESalphaComb.setVal(_mES_alpha_comb_configs_m[kkk]);
      mESslopeComb.setVal(_mES_Slope_comb_configs_m[kkk]);
      mESmeanPeak.setConstant(true);
      mESsigmaPeakL.setConstant(true);
      mESsigmaPeakR.setConstant(true);
      mESalphaPeakL.setConstant(true);
      mESrhoComb.setConstant(true);
      mESalphaComb.setConstant(true);
      mESslopeComb.setConstant(true);
      mESthresComb.setConstant(true);
      mESFracComb.setConstant(false);
      nmESPeak.setConstant(false);
      nmESComb.setConstant(false);
      PdfTotal_mES.fitTo(*data_mES_EgBRF[i]);
      h_SigPhoton_EBrest_mESFit_Syst_m[kkk]->SetBinContent(i+1,nmESPeak.getVal());
    }

  }

  for(int kkk=0;kkk<N_par_config;kkk++) {
    for(int lll=0;lll<h_SigPhoton_EBrest_mESFit_v2.GetXaxis()->GetNbins();lll++) {
      double N_nom = h_SigPhoton_EBrest_mESFit_v2.GetBinContent(lll+1);
      double N_p   = h_SigPhoton_EBrest_mESFit_Syst_p[kkk]->GetBinContent(lll+1);
      double N_m   = h_SigPhoton_EBrest_mESFit_Syst_m[kkk]->GetBinContent(lll+1);
      
      double Delta_p   = N_p   - N_nom;
      double Delta_m   = N_nom - N_m;
      double Delta_ave = 0.5*(Delta_p + Delta_m);
      h_SigPhoton_EBrest_mESFit_Syst_Delta[kkk]->SetBinContent(lll+1,Delta_ave);
    }
  }

  TLegend* leg_2 = new TLegend(0.7,0.7,0.9,0.9);
  leg_2->SetFillColor(10);
  for(int kkk=0;kkk<N_par_config;kkk++) {
    c1->Clear();
    c1->Divide(1,2);
    c1->cd(1);
    gPad->SetFillColor(10);
    gPad->SetFrameFillColor(10);
    gPad->SetTickx(1);
    gPad->SetTicky(1);
    gPad->SetLeftMargin(0.15);
    gPad->SetBottomMargin(0.15);
    gPad->SetRightMargin(0.15);
    Maximum = TMath::Max(h_SigPhoton_EBrest_mESFit_Syst_p[kkk]->GetMaximum(),
			 h_SigPhoton_EBrest_mESFit_Syst_m[kkk]->GetMaximum());
    porcent = 0.10;
    h_SigPhoton_EBrest_mESFit_Syst_p[kkk]->SetMaximum(Maximum*(1.0+porcent));
    h_SigPhoton_EBrest_mESFit_Syst_p[kkk]->Draw();
    h_SigPhoton_EBrest_mESFit_v2.Draw("same");
    h_SigPhoton_EBrest_mESFit_Syst_m[kkk]->Draw("same");
    if(kkk==0) {
      leg_2->AddEntry(h_SigPhoton_EBrest_mESFit_Syst_p[kkk],"Nom+#sigma","l");
      leg_2->AddEntry(&h_SigPhoton_EBrest_mESFit_v2,        "Nom",       "l");
      leg_2->AddEntry(h_SigPhoton_EBrest_mESFit_Syst_m[kkk],"Nom-#sigma","l");
    }
    leg_2->Draw("same");
    c1->cd(2);
    gPad->SetFillColor(10);
    gPad->SetFrameFillColor(10);
    gPad->SetTickx(1);
    gPad->SetTicky(1);
    gPad->SetLeftMargin(0.15);
    gPad->SetBottomMargin(0.15);
    gPad->SetRightMargin(0.15);
    h_SigPhoton_EBrest_mESFit_Syst_Delta[kkk]->Draw();
    c1->Print(EPSName.Data());
  }

  for(int i=0;i<h_mESFitYield_vs_EgBRF_Cov_mESFixX1->GetXaxis()->GetNbins();i++) {
    double Del_mESFitX1_i      = h_SigPhoton_EBrest_mESFit_Syst_Delta[0]->GetBinContent(i+1);
    double Del_mESFitX2_i      = h_SigPhoton_EBrest_mESFit_Syst_Delta[1]->GetBinContent(i+1);
    double Del_mESFitX3_i      = h_SigPhoton_EBrest_mESFit_Syst_Delta[2]->GetBinContent(i+1);
    double Del_mESFit_alphaP_i = h_SigPhoton_EBrest_mESFit_Syst_Delta[3]->GetBinContent(i+1);
    double Del_mESFit_BBComb_i = 0.5*(h_SigPhoton_EBrest_mESFit_Syst_Delta[4]->GetBinContent(i+1) + 
				      h_SigPhoton_EBrest_mESFit_Syst_Delta[5]->GetBinContent(i+1));
    double Del_mESFit_slope_i  = h_SigPhoton_EBrest_mESFit_Syst_Delta[6]->GetBinContent(i+1);
    double N_i                 = h_SigPhoton_EBrest_mESFit.GetBinContent(i+1);
    double ErrStat_i           = h_SigPhoton_EBrest_mESFit.GetBinError(i+1);
    for(int j=0;j<h_mESFitYield_vs_EgBRF_Cov_mESFixX1->GetYaxis()->GetNbins();j++) {
      double Del_mESFitX1_j      = h_SigPhoton_EBrest_mESFit_Syst_Delta[0]->GetBinContent(j+1);
      double Del_mESFitX2_j      = h_SigPhoton_EBrest_mESFit_Syst_Delta[1]->GetBinContent(j+1);
      double Del_mESFitX3_j      = h_SigPhoton_EBrest_mESFit_Syst_Delta[2]->GetBinContent(j+1);
      double Del_mESFit_alphaP_j = h_SigPhoton_EBrest_mESFit_Syst_Delta[3]->GetBinContent(j+1);
      double Del_mESFit_BBComb_j = 0.5*(h_SigPhoton_EBrest_mESFit_Syst_Delta[4]->GetBinContent(j+1) + 
					h_SigPhoton_EBrest_mESFit_Syst_Delta[5]->GetBinContent(j+1));
      double Del_mESFit_slope_j  = h_SigPhoton_EBrest_mESFit_Syst_Delta[6]->GetBinContent(j+1);
      double N_j                 = h_SigPhoton_EBrest_mESFit.GetBinContent(j+1);

      double cov = 0.0;

      //Stat:
      cov = 0.0;
      if(i==j) {
	cov = pow(ErrStat_i,2);
      }
      h_mESFitYield_vs_EgBRF_Cov_Stat->SetBinContent(i+1,j+1,cov);

      //mESFixX1:
      cov  = 0.0;
      cov += Del_mESFitX1_i*Del_mESFitX1_j;
      h_mESFitYield_vs_EgBRF_Cov_mESFixX1->SetBinContent(i+1,j+1,cov);

      //mESFixX2:
      cov  = 0.0;
      cov += Del_mESFitX2_i*Del_mESFitX2_j;
      h_mESFitYield_vs_EgBRF_Cov_mESFixX2->SetBinContent(i+1,j+1,cov);

      //mESFixX3:
      cov  = 0.0;
      cov += Del_mESFitX3_i*Del_mESFitX3_j;
      h_mESFitYield_vs_EgBRF_Cov_mESFixX3->SetBinContent(i+1,j+1,cov);

      //mESFix alphaP:
      cov  = 0.0;
      cov += Del_mESFit_alphaP_i*Del_mESFit_alphaP_j;
      h_mESFitYield_vs_EgBRF_Cov_mESFix_alphaP->SetBinContent(i+1,j+1,cov);

      //mESFix BBComb:
      cov  = 0.0;
      cov += Del_mESFit_BBComb_i*Del_mESFit_BBComb_j;
      h_mESFitYield_vs_EgBRF_Cov_mESFix_BBComb->SetBinContent(i+1,j+1,cov);

      //mESFix slope:
      cov  = 0.0;
      cov += Del_mESFit_slope_i*Del_mESFit_slope_j;
      h_mESFitYield_vs_EgBRF_Cov_mESFix_slope->SetBinContent(i+1,j+1,cov);

      //mES Fit Bias:
      cov  = 0.0;
      cov += N_i*N_j*pow(mES_FitBias,2);
      h_mESFitYield_vs_EgBRF_Cov_mESFitBias->SetBinContent(i+1,j+1,cov);
    }
  }

  c1->Clear();
  h_mESFitYield_vs_EgBRF_Cov_Stat->Draw("colz");
  c1->Print(EPSName.Data());

  c1->Clear();
  h_mESFitYield_vs_EgBRF_Cov_mESFixX1->Draw("colz");
  c1->Print(EPSName.Data());

  c1->Clear();
  h_mESFitYield_vs_EgBRF_Cov_mESFixX2->Draw("colz");
  c1->Print(EPSName.Data());

  c1->Clear();
  h_mESFitYield_vs_EgBRF_Cov_mESFixX3->Draw("colz");
  c1->Print(EPSName.Data());

  c1->Clear();
  h_mESFitYield_vs_EgBRF_Cov_mESFix_alphaP->Draw("colz");
  c1->Print(EPSName.Data());

  c1->Clear();
  h_mESFitYield_vs_EgBRF_Cov_mESFix_BBComb->Draw("colz");
  c1->Print(EPSName.Data());

  c1->Clear();
  h_mESFitYield_vs_EgBRF_Cov_mESFix_slope->Draw("colz");
  c1->Print(EPSName.Data());

  c1->Clear();
  h_mESFitYield_vs_EgBRF_Cov_mESFitBias->Draw("colz");
  c1->Print(EPSName.Data());

  int Nbins_ppp = h_mESFitYield_vs_EgBRF_Cov_Stat->GetXaxis()->GetNbins();
  double Min_diag;
  double Min_corr,Max_corr;
  double Test_sym;

  //Fit-yields Stat error:
  Min_diag = Min_corr = 1.0e+20;
  Max_corr = -1.0e+20;
  Test_sym = -1.0e+20;
  for(int i=0;i<Nbins_ppp;i++) {
    double sigma2 = h_mESFitYield_vs_EgBRF_Cov_Stat->GetBinContent(i+1,i+1);
    if(Min_diag > sigma2) Min_diag = sigma2;
    if(sqrt(sigma2) < epsilon) continue;
    for(int j=0;j<Nbins_ppp;j++) {
      if(i==j) continue;
      double corr  = h_mESFitYield_vs_EgBRF_Cov_Stat->GetBinContent(i+1,j+1);
      corr        /= sqrt(h_mESFitYield_vs_EgBRF_Cov_Stat->GetBinContent(i+1,i+1));
      corr        /= sqrt(h_mESFitYield_vs_EgBRF_Cov_Stat->GetBinContent(j+1,j+1));
      if(Min_corr > corr) Min_corr = corr;
      if(Max_corr < corr) Max_corr = corr;

      double diff  = h_mESFitYield_vs_EgBRF_Cov_Stat->GetBinContent(i+1,j+1);
      diff        -= h_mESFitYield_vs_EgBRF_Cov_Stat->GetBinContent(j+1,i+1);
      if(Test_sym < TMath::Abs(diff)) Test_sym = TMath::Abs(diff);
    }
  }
  cout << endl;
  cout << "mES-Fitted-yields. Stat Error:" << endl;
  cout << "Min diag   = " << Min_diag << endl;
  cout << "Range corr = (" << Min_corr << "," << Max_corr << ")" << endl;
  cout << "Test sym   = " << Test_sym << endl;
  cout << endl;

  //Fit-yields mESFixX1 error:
  Min_diag = Min_corr = 1.0e+20;
  Max_corr = -1.0e+20;
  Test_sym = -1.0e+20;
  for(int i=0;i<Nbins_ppp;i++) {
    double sigma2 = h_mESFitYield_vs_EgBRF_Cov_mESFixX1->GetBinContent(i+1,i+1);
    if(Min_diag > sigma2) Min_diag = sigma2;
    if(sqrt(sigma2) < epsilon) continue;
    for(int j=0;j<Nbins_ppp;j++) {
      if(i==j) continue;
      double corr  = h_mESFitYield_vs_EgBRF_Cov_mESFixX1->GetBinContent(i+1,j+1);
      corr        /= sqrt(h_mESFitYield_vs_EgBRF_Cov_mESFixX1->GetBinContent(i+1,i+1));
      corr        /= sqrt(h_mESFitYield_vs_EgBRF_Cov_mESFixX1->GetBinContent(j+1,j+1));
      if(Min_corr > corr) Min_corr = corr;
      if(Max_corr < corr) Max_corr = corr;

      double diff  = h_mESFitYield_vs_EgBRF_Cov_mESFixX1->GetBinContent(i+1,j+1);
      diff        -= h_mESFitYield_vs_EgBRF_Cov_mESFixX1->GetBinContent(j+1,i+1);
      if(Test_sym < TMath::Abs(diff)) Test_sym = TMath::Abs(diff);
    }
  }
  cout << endl;
  cout << "mES-Fitted-yields. mES-Fix-X1 Error:" << endl;
  cout << "Min diag   = " << Min_diag << endl;
  cout << "Range corr = (" << Min_corr << "," << Max_corr << ")" << endl;
  cout << "Test sym   = " << Test_sym << endl;
  cout << endl;

  //Fit-yields mESFixX2 error:
  Min_diag = Min_corr = 1.0e+20;
  Max_corr = -1.0e+20;
  Test_sym = -1.0e+20;
  for(int i=0;i<Nbins_ppp;i++) {
    double sigma2 = h_mESFitYield_vs_EgBRF_Cov_mESFixX2->GetBinContent(i+1,i+1);
    if(Min_diag > sigma2) Min_diag = sigma2;
    if(sqrt(sigma2) < epsilon) continue;
    for(int j=0;j<Nbins_ppp;j++) {
      if(i==j) continue;
      double corr  = h_mESFitYield_vs_EgBRF_Cov_mESFixX2->GetBinContent(i+1,j+1);
      corr        /= sqrt(h_mESFitYield_vs_EgBRF_Cov_mESFixX2->GetBinContent(i+1,i+1));
      corr        /= sqrt(h_mESFitYield_vs_EgBRF_Cov_mESFixX2->GetBinContent(j+1,j+1));
      if(Min_corr > corr) Min_corr = corr;
      if(Max_corr < corr) Max_corr = corr;

      double diff  = h_mESFitYield_vs_EgBRF_Cov_mESFixX2->GetBinContent(i+1,j+1);
      diff        -= h_mESFitYield_vs_EgBRF_Cov_mESFixX2->GetBinContent(j+1,i+1);
      if(Test_sym < TMath::Abs(diff)) Test_sym = TMath::Abs(diff);
    }
  }
  cout << endl;
  cout << "mES-Fitted-yields. mES-Fix-X2 Error:" << endl;
  cout << "Min diag   = " << Min_diag << endl;
  cout << "Range corr = (" << Min_corr << "," << Max_corr << ")" << endl;
  cout << "Test sym   = " << Test_sym << endl;
  cout << endl;

  //Fit-yields mESFixX3 error:
  Min_diag = Min_corr = 1.0e+20;
  Max_corr = -1.0e+20;
  Test_sym = -1.0e+20;
  for(int i=0;i<Nbins_ppp;i++) {
    double sigma2 = h_mESFitYield_vs_EgBRF_Cov_mESFixX3->GetBinContent(i+1,i+1);
    if(Min_diag > sigma2) Min_diag = sigma2;
    if(sqrt(sigma2) < epsilon) continue;
    for(int j=0;j<Nbins_ppp;j++) {
      if(i==j) continue;
      double corr  = h_mESFitYield_vs_EgBRF_Cov_mESFixX3->GetBinContent(i+1,j+1);
      corr        /= sqrt(h_mESFitYield_vs_EgBRF_Cov_mESFixX3->GetBinContent(i+1,i+1));
      corr        /= sqrt(h_mESFitYield_vs_EgBRF_Cov_mESFixX3->GetBinContent(j+1,j+1));
      if(Min_corr > corr) Min_corr = corr;
      if(Max_corr < corr) Max_corr = corr;

      double diff  = h_mESFitYield_vs_EgBRF_Cov_mESFixX3->GetBinContent(i+1,j+1);
      diff        -= h_mESFitYield_vs_EgBRF_Cov_mESFixX3->GetBinContent(j+1,i+1);
      if(Test_sym < TMath::Abs(diff)) Test_sym = TMath::Abs(diff);
    }
  }
  cout << endl;
  cout << "mES-Fitted-yields. mES-Fix-X3 Error:" << endl;
  cout << "Min diag   = " << Min_diag << endl;
  cout << "Range corr = (" << Min_corr << "," << Max_corr << ")" << endl;
  cout << "Test sym   = " << Test_sym << endl;
  cout << endl;

  //Fit-yields mESFix_alphaP error:
  Min_diag = Min_corr = 1.0e+20;
  Max_corr = -1.0e+20;
  Test_sym = -1.0e+20;
  for(int i=0;i<Nbins_ppp;i++) {
    double sigma2 = h_mESFitYield_vs_EgBRF_Cov_mESFix_alphaP->GetBinContent(i+1,i+1);
    if(Min_diag > sigma2) Min_diag = sigma2;
    if(sqrt(sigma2) < epsilon) continue;
    for(int j=0;j<Nbins_ppp;j++) {
      if(i==j) continue;
      double corr  = h_mESFitYield_vs_EgBRF_Cov_mESFix_alphaP->GetBinContent(i+1,j+1);
      corr        /= sqrt(h_mESFitYield_vs_EgBRF_Cov_mESFix_alphaP->GetBinContent(i+1,i+1));
      corr        /= sqrt(h_mESFitYield_vs_EgBRF_Cov_mESFix_alphaP->GetBinContent(j+1,j+1));
      if(Min_corr > corr) Min_corr = corr;
      if(Max_corr < corr) Max_corr = corr;

      double diff  = h_mESFitYield_vs_EgBRF_Cov_mESFix_alphaP->GetBinContent(i+1,j+1);
      diff        -= h_mESFitYield_vs_EgBRF_Cov_mESFix_alphaP->GetBinContent(j+1,i+1);
      if(Test_sym < TMath::Abs(diff)) Test_sym = TMath::Abs(diff);
    }
  }
  cout << endl;
  cout << "mES-Fitted-yields. mES-Fix-alphaP Error:" << endl;
  cout << "Min diag   = " << Min_diag << endl;
  cout << "Range corr = (" << Min_corr << "," << Max_corr << ")" << endl;
  cout << "Test sym   = " << Test_sym << endl;
  cout << endl;

  //Fit-yields mESFix_BBComb error:
  Min_diag = Min_corr = 1.0e+20;
  Max_corr = -1.0e+20;
  Test_sym = -1.0e+20;
  for(int i=0;i<Nbins_ppp;i++) {
    double sigma2 = h_mESFitYield_vs_EgBRF_Cov_mESFix_BBComb->GetBinContent(i+1,i+1);
    if(Min_diag > sigma2) Min_diag = sigma2;
    if(sqrt(sigma2) < epsilon) continue;
    for(int j=0;j<Nbins_ppp;j++) {
      if(i==j) continue;
      double corr  = h_mESFitYield_vs_EgBRF_Cov_mESFix_BBComb->GetBinContent(i+1,j+1);
      corr        /= sqrt(h_mESFitYield_vs_EgBRF_Cov_mESFix_BBComb->GetBinContent(i+1,i+1));
      corr        /= sqrt(h_mESFitYield_vs_EgBRF_Cov_mESFix_BBComb->GetBinContent(j+1,j+1));
      if(Min_corr > corr) Min_corr = corr;
      if(Max_corr < corr) Max_corr = corr;

      double diff  = h_mESFitYield_vs_EgBRF_Cov_mESFix_BBComb->GetBinContent(i+1,j+1);
      diff        -= h_mESFitYield_vs_EgBRF_Cov_mESFix_BBComb->GetBinContent(j+1,i+1);
      if(Test_sym < TMath::Abs(diff)) Test_sym = TMath::Abs(diff);
    }
  }
  cout << endl;
  cout << "mES-Fitted-yields. mES-Fix-BBComb Error:" << endl;
  cout << "Min diag   = " << Min_diag << endl;
  cout << "Range corr = (" << Min_corr << "," << Max_corr << ")" << endl;
  cout << "Test sym   = " << Test_sym << endl;
  cout << endl;

  //Fit-yields mESFix_slope error:
  Min_diag = Min_corr = 1.0e+20;
  Max_corr = -1.0e+20;
  Test_sym = -1.0e+20;
  for(int i=0;i<Nbins_ppp;i++) {
    double sigma2 = h_mESFitYield_vs_EgBRF_Cov_mESFix_slope->GetBinContent(i+1,i+1);
    if(Min_diag > sigma2) Min_diag = sigma2;
    if(sqrt(sigma2) < epsilon) continue;
    for(int j=0;j<Nbins_ppp;j++) {
      if(i==j) continue;
      double corr  = h_mESFitYield_vs_EgBRF_Cov_mESFix_slope->GetBinContent(i+1,j+1);
      corr        /= sqrt(h_mESFitYield_vs_EgBRF_Cov_mESFix_slope->GetBinContent(i+1,i+1));
      corr        /= sqrt(h_mESFitYield_vs_EgBRF_Cov_mESFix_slope->GetBinContent(j+1,j+1));
      if(Min_corr > corr) Min_corr = corr;
      if(Max_corr < corr) Max_corr = corr;

      double diff  = h_mESFitYield_vs_EgBRF_Cov_mESFix_slope->GetBinContent(i+1,j+1);
      diff        -= h_mESFitYield_vs_EgBRF_Cov_mESFix_slope->GetBinContent(j+1,i+1);
      if(Test_sym < TMath::Abs(diff)) Test_sym = TMath::Abs(diff);
    }
  }
  cout << endl;
  cout << "mES-Fitted-yields. mES-Fix-slope Error:" << endl;
  cout << "Min diag   = " << Min_diag << endl;
  cout << "Range corr = (" << Min_corr << "," << Max_corr << ")" << endl;
  cout << "Test sym   = " << Test_sym << endl;
  cout << endl;

  //Fit-yields mESFitBias error:
  Min_diag = Min_corr = 1.0e+20;
  Max_corr = -1.0e+20;
  Test_sym = -1.0e+20;
  for(int i=0;i<Nbins_ppp;i++) {
    double sigma2 = h_mESFitYield_vs_EgBRF_Cov_mESFitBias->GetBinContent(i+1,i+1);
    if(Min_diag > sigma2) Min_diag = sigma2;
    if(sqrt(sigma2) < epsilon) continue;
    for(int j=0;j<Nbins_ppp;j++) {
      if(i==j) continue;
      double corr  = h_mESFitYield_vs_EgBRF_Cov_mESFitBias->GetBinContent(i+1,j+1);
      corr        /= sqrt(h_mESFitYield_vs_EgBRF_Cov_mESFitBias->GetBinContent(i+1,i+1));
      corr        /= sqrt(h_mESFitYield_vs_EgBRF_Cov_mESFitBias->GetBinContent(j+1,j+1));
      if(Min_corr > corr) Min_corr = corr;
      if(Max_corr < corr) Max_corr = corr;

      double diff  = h_mESFitYield_vs_EgBRF_Cov_mESFitBias->GetBinContent(i+1,j+1);
      diff        -= h_mESFitYield_vs_EgBRF_Cov_mESFitBias->GetBinContent(j+1,i+1);
      if(Test_sym < TMath::Abs(diff)) Test_sym = TMath::Abs(diff);
    }
  }
  cout << endl;
  cout << "mES-Fitted-yields. mES-Fit-bias Error:" << endl;
  cout << "Min diag   = " << Min_diag << endl;
  cout << "Range corr = (" << Min_corr << "," << Max_corr << ")" << endl;
  cout << "Test sym   = " << Test_sym << endl;
  cout << endl;

  for(int i=0;i<h_SigPhoton_EBrest_mESFit.GetXaxis()->GetNbins();i++) {
    double e = 0.0;
    e += h_mESFitYield_vs_EgBRF_Cov_Stat->GetBinContent(i+1,i+1);
    e += h_mESFitYield_vs_EgBRF_Cov_mESFixX1->GetBinContent(i+1,i+1);
    e += h_mESFitYield_vs_EgBRF_Cov_mESFixX2->GetBinContent(i+1,i+1);
    e += h_mESFitYield_vs_EgBRF_Cov_mESFixX3->GetBinContent(i+1,i+1);
    e += h_mESFitYield_vs_EgBRF_Cov_mESFix_alphaP->GetBinContent(i+1,i+1);
    e += h_mESFitYield_vs_EgBRF_Cov_mESFix_BBComb->GetBinContent(i+1,i+1);
    e += h_mESFitYield_vs_EgBRF_Cov_mESFix_slope->GetBinContent(i+1,i+1);
    e += h_mESFitYield_vs_EgBRF_Cov_mESFitBias->GetBinContent(i+1,i+1);
    e  = sqrt(e);
    h_SigPhoton_EBrest_mESFit.SetBinError(i+1,e);
  }


  int Medium = 1 + Nbins_SigPhoton_EBrest_Fit/2;
  cout << endl;
  cout << endl;
  cout << "=====================================================" << endl;
  cout << "N^{P} yields results with syst. errors from mES fits:" << endl;
  cout << "===========================================================================================================================" 
       << endl;
  cout << "Type                ";
  for(int i=0;i<Medium;i++) {
    double c = h_SigPhoton_EBrest_mESFit.GetBinCenter(i+1);
    double w = h_SigPhoton_EBrest_mESFit.GetBinWidth(i+1)*0.5;
    double R1,R2;
    R1 = c - w;
    R2 = c + w;
    sprintf(ytitle,"%.1f",R1);
    cout << "(" << ytitle << ",";
    sprintf(ytitle,"%.1f",R2);
    cout << ytitle << ")   ";
  }
  cout << endl;
  cout << "---------------------------------------------------------------------------------------------------------------------------"
       << endl;
  cout << "Yield               ";
  for(int i=0;i<Medium;i++) {
    sprintf(ytitle,"%.2f",h_SigPhoton_EBrest_mESFit.GetBinContent(i+1));
    cout << ytitle
	 << "        ";
  }
  cout << endl;
  cout << "Err-MCStat          ";
  for(int i=0;i<Medium;i++) {
    sprintf(ytitle,"%.2f",sqrt(h_mESFitYield_vs_EgBRF_Cov_Stat->GetBinContent(i+1,i+1)));
    cout << ytitle
	 << "        ";
  }
  cout << endl;
  cout << "Err-mESFixX1        ";
  for(int i=0;i<Medium;i++) {
    sprintf(ytitle,"%.2f",sqrt(h_mESFitYield_vs_EgBRF_Cov_mESFixX1->GetBinContent(i+1,i+1)));
    cout << ytitle
	 << "        ";
  }
  cout << endl;
  cout << "Err-mESFixX2        ";
  for(int i=0;i<Medium;i++) {
    sprintf(ytitle,"%.2f",sqrt(h_mESFitYield_vs_EgBRF_Cov_mESFixX2->GetBinContent(i+1,i+1)));
    cout << ytitle
	 << "        ";
  }
  cout << endl;
  cout << "Err-mESFixX3        ";
  for(int i=0;i<Medium;i++) {
    sprintf(ytitle,"%.2f",sqrt(h_mESFitYield_vs_EgBRF_Cov_mESFixX3->GetBinContent(i+1,i+1)));
    cout << ytitle
	 << "        ";
  }
  cout << endl;
  cout << "Err-mESFix_alphaP    ";
  for(int i=0;i<Medium;i++) {
    sprintf(ytitle,"%.2f",sqrt(h_mESFitYield_vs_EgBRF_Cov_mESFix_alphaP->GetBinContent(i+1,i+1)));
    cout << ytitle
	 << "        ";
  }
  cout << endl;
  cout << "Err-mESFix_BBComb    ";
  for(int i=0;i<Medium;i++) {
    sprintf(ytitle,"%.2f",sqrt(h_mESFitYield_vs_EgBRF_Cov_mESFix_BBComb->GetBinContent(i+1,i+1)));
    cout << ytitle
	 << "        ";
  }
  cout << endl;
  cout << "Err-mESFix_slope     ";
  for(int i=0;i<Medium;i++) {
    sprintf(ytitle,"%.2f",sqrt(h_mESFitYield_vs_EgBRF_Cov_mESFix_slope->GetBinContent(i+1,i+1)));
    cout << ytitle
	 << "        ";
  }
  cout << endl;
  cout << "Err-mESFitBias       ";
  for(int i=0;i<Medium;i++) {
    sprintf(ytitle,"%.2f",sqrt(h_mESFitYield_vs_EgBRF_Cov_mESFitBias->GetBinContent(i+1,i+1)));
    cout << ytitle
	 << "        ";
  }
  cout << endl;
  cout << "---------------------------------------------------------------------------------------------------------------------------"
       << endl;
  cout << "Err-tot             ";
  for(int i=0;i<Medium;i++) {
    sprintf(ytitle,"%.2f",h_SigPhoton_EBrest_mESFit.GetBinError(i+1));
    cout << ytitle
	 << "        ";
  }
  cout << endl;
  cout << "===========================================================================================================================" 
       << endl;
  cout << "Type                ";
  for(int i=Medium;i<Nbins_SigPhoton_EBrest_Fit;i++) {
    double c = h_SigPhoton_EBrest_mESFit.GetBinCenter(i+1);
    double w = h_SigPhoton_EBrest_mESFit.GetBinWidth(i+1)*0.5;
    double R1,R2;
    R1 = c - w;
    R2 = c + w;
    sprintf(ytitle,"%.1f",R1);
    cout << "(" << ytitle << ",";
    sprintf(ytitle,"%.1f",R2);
    cout << ytitle << ")   ";
  }
  cout << endl;
  cout << "---------------------------------------------------------------------------------------------------------------------------"
       << endl;
  cout << "Yield               ";
  for(int i=Medium;i<Nbins_SigPhoton_EBrest_Fit;i++) {
    sprintf(ytitle,"%.2f",h_SigPhoton_EBrest_mESFit.GetBinContent(i+1));
    cout << ytitle
	 << "        ";
  }
  cout << endl;
  cout << "Err-MCStat          ";
  for(int i=Medium;i<Nbins_SigPhoton_EBrest_Fit;i++) {
    sprintf(ytitle,"%.2f",sqrt(h_mESFitYield_vs_EgBRF_Cov_Stat->GetBinContent(i+1,i+1)));
    cout << ytitle
	 << "        ";
  }
  cout << endl;
  cout << "Err-mESFixX1        ";
  for(int i=Medium;i<Nbins_SigPhoton_EBrest_Fit;i++) {
    sprintf(ytitle,"%.2f",sqrt(h_mESFitYield_vs_EgBRF_Cov_mESFixX1->GetBinContent(i+1,i+1)));
    cout << ytitle
	 << "        ";
  }
  cout << endl;
  cout << "Err-mESFixX2        ";
  for(int i=Medium;i<Nbins_SigPhoton_EBrest_Fit;i++) {
    sprintf(ytitle,"%.2f",sqrt(h_mESFitYield_vs_EgBRF_Cov_mESFixX2->GetBinContent(i+1,i+1)));
    cout << ytitle
	 << "        ";
  }
  cout << endl;
  cout << "Err-mESFixX3        ";
  for(int i=Medium;i<Nbins_SigPhoton_EBrest_Fit;i++) {
    sprintf(ytitle,"%.2f",sqrt(h_mESFitYield_vs_EgBRF_Cov_mESFixX3->GetBinContent(i+1,i+1)));
    cout << ytitle
	 << "        ";
  }
  cout << endl;
  cout << "Err-mESFix_alphaP    ";
  for(int i=Medium;i<Nbins_SigPhoton_EBrest_Fit;i++) {
    sprintf(ytitle,"%.2f",sqrt(h_mESFitYield_vs_EgBRF_Cov_mESFix_alphaP->GetBinContent(i+1,i+1)));
    cout << ytitle
	 << "        ";
  }
  cout << endl;
  cout << "Err-mESFix_BBComb    ";
  for(int i=Medium;i<Nbins_SigPhoton_EBrest_Fit;i++) {
    sprintf(ytitle,"%.2f",sqrt(h_mESFitYield_vs_EgBRF_Cov_mESFix_BBComb->GetBinContent(i+1,i+1)));
    cout << ytitle
	 << "        ";
  }
  cout << endl;
  cout << "Err-mESFix_slope     ";
  for(int i=Medium;i<Nbins_SigPhoton_EBrest_Fit;i++) {
    sprintf(ytitle,"%.2f",sqrt(h_mESFitYield_vs_EgBRF_Cov_mESFix_slope->GetBinContent(i+1,i+1)));
    cout << ytitle
	 << "        ";
  }
  cout << endl;
  cout << "Err-mESFitBias       ";
  for(int i=Medium;i<Nbins_SigPhoton_EBrest_Fit;i++) {
    sprintf(ytitle,"%.2f",sqrt(h_mESFitYield_vs_EgBRF_Cov_mESFitBias->GetBinContent(i+1,i+1)));
    cout << ytitle
	 << "        ";
  }
  cout << endl;
  cout << "---------------------------------------------------------------------------------------------------------------------------"
       << endl;
  cout << "Err-tot             ";
  for(int i=Medium;i<Nbins_SigPhoton_EBrest_Fit;i++) {
    sprintf(ytitle,"%.2f",h_SigPhoton_EBrest_mESFit.GetBinError(i+1));
    cout << ytitle
	 << "        ";
  }
  cout << endl;
  cout << "===========================================================================================================================" 
       << endl;
  cout << endl;

  h_SigPhoton_EBrest_Fit_BrecoTM.SetStats(false);
  h_SigPhoton_EBrest_Fit_BrecoTM.SetLineColor(kGreen+2);
  c1->Clear();
  Maximum = -1.0e+20;
  Minimum =  1.0e+20;
  for(int i=0;i<h_SigPhoton_EBrest_mESFit.GetXaxis()->GetNbins();i++) {
    double v,e;
    v = h_SigPhoton_EBrest_mESFit.GetBinContent(i+1);
    e = h_SigPhoton_EBrest_mESFit.GetBinError(i+1);
    if(Maximum < v+e) Maximum = v+e;
    if(Minimum > v-e) Minimum = v-e;

    v = h_EgBRF_AllBkgPrediction_tot_corr->GetBinContent(i+1);
    e = h_EgBRF_AllBkgPrediction_tot_corr->GetBinError(i+1);
    if(Maximum < v+e) Maximum = v+e;
    if(Minimum > v-e) Minimum = v-e;
  }
  porcent = 0.10;
  h_SigPhoton_EBrest_mESFit.SetMaximum(Maximum + porcent*(Maximum-Minimum));
  h_SigPhoton_EBrest_mESFit.SetMinimum(Minimum - porcent*(Maximum-Minimum));
  h_EgBRF_AllBkgPrediction_tot_corr->SetMaximum(Maximum + porcent*(Maximum-Minimum));
  h_EgBRF_AllBkgPrediction_tot_corr->SetMinimum(Minimum - porcent*(Maximum-Minimum));
  h_SigPhoton_EBrest_mESFit.Draw();
  h_EgBRF_AllBkgPrediction_tot_corr->Draw("same");
  TLine* l_1 = new TLine(h_SigPhoton_EBrest_mESFit.GetXaxis()->GetXmin(),0.0,
			 h_SigPhoton_EBrest_mESFit.GetXaxis()->GetXmax(),0.0);
  l_1->SetLineColor(1);
  l_1->SetLineWidth(2);
  l_1->SetLineStyle(2);
  if(h_SigPhoton_EBrest_mESFit.GetMinimum() < 0.0) l_1->Draw();
  TLine* l_sig1 = new TLine(R_EgBRF_Sig[0],h_SigPhoton_EBrest_mESFit.GetMaximum(),
			    R_EgBRF_Sig[0],0.0);
  l_sig1->SetLineColor(2);
  l_sig1->SetLineWidth(2);
  l_sig1->SetLineStyle(2);
  l_sig1->Draw();
  TLine* l_sig2 = new TLine(R_EgBRF_Sig[1],h_SigPhoton_EBrest_mESFit.GetMaximum(),
			    R_EgBRF_Sig[1],0.0);
  l_sig2->SetLineColor(2);
  l_sig2->SetLineWidth(2);
  l_sig2->SetLineStyle(2);
  l_sig2->Draw();

  TLegend* leg1 = new TLegend(0.7,0.7,0.9,0.9);
  leg1->SetFillColor(10);
  leg1->AddEntry(&h_SigPhoton_EBrest_mESFit,       "m_{ES} fit",    "l");
  leg1->AddEntry(h_EgBRF_AllBkgPrediction_tot_corr,"Bkg prediction","l");
  leg1->Draw("same");
  c1->Print(EPSName.Data());

  TLegend* leg12 = new TLegend(0.7,0.7,0.9,0.9);
  leg12->SetFillColor(10);
  if(_IsMC > 1) {
    c1->Clear();
    h_SigPhoton_EBrest_mESFit.Draw();
    h_EgBRF_AllBkgPrediction_tot_corr->Draw("same");
    h_SigPhoton_EBrest_Fit_BrecoTM.Draw("same");
    if(h_SigPhoton_EBrest_mESFit.GetMinimum() < 0.0) l_1->Draw();
    l_sig1->Draw();
    l_sig2->Draw();
    leg12->AddEntry(&h_SigPhoton_EBrest_mESFit,       "m_{ES} fit",    "l");
    leg12->AddEntry(h_EgBRF_AllBkgPrediction_tot_corr,"Bkg prediction","l");
    leg12->AddEntry(&h_SigPhoton_EBrest_Fit_BrecoTM,  "Breco-TM",      "l");
    leg12->Draw("same");
    c1->Print(EPSName.Data());
  }

  c1->Clear();
  Maximum = -1.0e+20;
  Minimum =  1.0e+20;
  for(int i=0;i<h_SigPhoton_EBrest_mESFit.GetXaxis()->GetNbins();i++) {
    double v,e,c,w;

    c = h_SigPhoton_EBrest_mESFit.GetBinCenter(i+1);
    w = h_SigPhoton_EBrest_mESFit.GetBinWidth(i+1)*0.5;

    if(!(c-w + 1.0e-6 > R_EgBRF_Sig[0] && 
	 c+w - 1.0e-6 < R_EgBRF_Sig[1])) continue;

    v = h_SigPhoton_EBrest_mESFit.GetBinContent(i+1);
    e = h_SigPhoton_EBrest_mESFit.GetBinError(i+1);
    if(Maximum < v+e) Maximum = v+e;
    if(Minimum > v-e) Minimum = v-e;

    v = h_EgBRF_AllBkgPrediction_tot_corr->GetBinContent(i+1);
    e = h_EgBRF_AllBkgPrediction_tot_corr->GetBinError(i+1);
    if(Maximum < v+e) Maximum = v+e;
    if(Minimum > v-e) Minimum = v-e;
  }
  porcent = 0.10;
  h_SigPhoton_EBrest_mESFit_REFSig.SetMaximum(Maximum + porcent*(Maximum-Minimum));
  h_SigPhoton_EBrest_mESFit_REFSig.SetMinimum(Minimum - porcent*(Maximum-Minimum));
  h_SigPhoton_EBrest_mESFit_REFSig.Draw();
  h_SigPhoton_EBrest_mESFit.Draw("same");
  h_EgBRF_AllBkgPrediction_tot_corr->Draw("same");
  TLine* l_2 = new TLine(h_SigPhoton_EBrest_mESFit_REFSig.GetXaxis()->GetXmin(),0.0,
			 h_SigPhoton_EBrest_mESFit_REFSig.GetXaxis()->GetXmax(),0.0);
  l_2->SetLineColor(1);
  l_2->SetLineWidth(2);
  l_2->SetLineStyle(2);
  if(h_SigPhoton_EBrest_mESFit_REFSig.GetMinimum() < 0.0) l_2->Draw();
  leg1->Draw("same");
  c1->Print(EPSName.Data());

  if(_IsMC > 1) {
    c1->Clear();
    h_SigPhoton_EBrest_mESFit_REFSig.Draw();
    h_SigPhoton_EBrest_mESFit.Draw("same");
    h_EgBRF_AllBkgPrediction_tot_corr->Draw("same");
    h_SigPhoton_EBrest_Fit_BrecoTM.Draw("same");
    if(h_SigPhoton_EBrest_mESFit_REFSig.GetMinimum() < 0.0) l_2->Draw();
    leg12->Draw("same");
    c1->Print(EPSName.Data());
  }

  cout << endl;
  cout << "mESFit-X1:" << endl;
  for(int i=0;i<h_SigPhoton_EBrest_mESFit.GetXaxis()->GetNbins();i++) {
    for(int j=0;j<h_SigPhoton_EBrest_mESFit.GetXaxis()->GetNbins();j++) {
      if(i>j) continue;
      double den;
      double corr_N  = h_mESFitYield_vs_EgBRF_Cov_mESFixX1->GetBinContent(i+1,j+1);
      den            = sqrt(h_mESFitYield_vs_EgBRF_Cov_mESFixX1->GetBinContent(i+1,i+1));
      den           *= sqrt(h_mESFitYield_vs_EgBRF_Cov_mESFixX1->GetBinContent(j+1,j+1));
      if(den > epsilon) corr_N /= den;
      else              corr_N  = -999.0;

      double corr_B  = h_EgBRF_AllBkgPrediction_Cov_Pi0EtaBkgmESFixX1->GetBinContent(i+1,j+1);
      den            = sqrt(h_EgBRF_AllBkgPrediction_Cov_Pi0EtaBkgmESFixX1->GetBinContent(i+1,i+1));
      den           *= sqrt(h_EgBRF_AllBkgPrediction_Cov_Pi0EtaBkgmESFixX1->GetBinContent(j+1,j+1));
      if(den > epsilon) corr_B /= den;
      else              corr_B  = -999.0;

      cout << i+1 << "  " << j+1 << "  "
	   << corr_N << "  "
	   << corr_B << "  "
	   << endl;
    }
  }
  cout << endl;

  cout << endl;
  cout << "mESFix-X2:" << endl;
  for(int i=0;i<h_SigPhoton_EBrest_mESFit.GetXaxis()->GetNbins();i++) {
    for(int j=0;j<h_SigPhoton_EBrest_mESFit.GetXaxis()->GetNbins();j++) {
      if(i>j) continue;
      double den;
      double corr_N  = h_mESFitYield_vs_EgBRF_Cov_mESFixX2->GetBinContent(i+1,j+1);
      den            = sqrt(h_mESFitYield_vs_EgBRF_Cov_mESFixX2->GetBinContent(i+1,i+1));
      den           *= sqrt(h_mESFitYield_vs_EgBRF_Cov_mESFixX2->GetBinContent(j+1,j+1));
      if(den > epsilon) corr_N /= den;
      else              corr_N  = -999.0;

      double corr_B  = h_EgBRF_AllBkgPrediction_Cov_Pi0EtaBkgmESFixX2->GetBinContent(i+1,j+1);
      den            = sqrt(h_EgBRF_AllBkgPrediction_Cov_Pi0EtaBkgmESFixX2->GetBinContent(i+1,i+1));
      den           *= sqrt(h_EgBRF_AllBkgPrediction_Cov_Pi0EtaBkgmESFixX2->GetBinContent(j+1,j+1));
      if(den > epsilon) corr_B /= den;
      else              corr_B  = -999.0;

      cout << i+1 << "  " << j+1 << "  "
	   << corr_N << "  "
	   << corr_B << "  "
	   << endl;
    }
  }
  cout << endl;

  cout << endl;
  cout << "mESFix-X3:" << endl;
  for(int i=0;i<h_SigPhoton_EBrest_mESFit.GetXaxis()->GetNbins();i++) {
    for(int j=0;j<h_SigPhoton_EBrest_mESFit.GetXaxis()->GetNbins();j++) {
      if(i>j) continue;
      double den;
      double corr_N  = h_mESFitYield_vs_EgBRF_Cov_mESFixX3->GetBinContent(i+1,j+1);
      den            = sqrt(h_mESFitYield_vs_EgBRF_Cov_mESFixX3->GetBinContent(i+1,i+1));
      den           *= sqrt(h_mESFitYield_vs_EgBRF_Cov_mESFixX3->GetBinContent(j+1,j+1));
      if(den > epsilon) corr_N /= den;
      else              corr_N  = -999.0;

      double corr_B  = h_EgBRF_AllBkgPrediction_Cov_Pi0EtaBkgmESFixX3->GetBinContent(i+1,j+1);
      den            = sqrt(h_EgBRF_AllBkgPrediction_Cov_Pi0EtaBkgmESFixX3->GetBinContent(i+1,i+1));
      den           *= sqrt(h_EgBRF_AllBkgPrediction_Cov_Pi0EtaBkgmESFixX3->GetBinContent(j+1,j+1));
      if(den > epsilon) corr_B /= den;
      else              corr_B  = -999.0;

      cout << i+1 << "  " << j+1 << "  "
	   << corr_N << "  "
	   << corr_B << "  "
	   << endl;
    }
  }
  cout << endl;

  cout << endl;
  cout << "mESFix-alphaP:" << endl;
  for(int i=0;i<h_SigPhoton_EBrest_mESFit.GetXaxis()->GetNbins();i++) {
    for(int j=0;j<h_SigPhoton_EBrest_mESFit.GetXaxis()->GetNbins();j++) {
      if(i>j) continue;
      double den;
      double corr_N  = h_mESFitYield_vs_EgBRF_Cov_mESFix_alphaP->GetBinContent(i+1,j+1);
      den            = sqrt(h_mESFitYield_vs_EgBRF_Cov_mESFix_alphaP->GetBinContent(i+1,i+1));
      den           *= sqrt(h_mESFitYield_vs_EgBRF_Cov_mESFix_alphaP->GetBinContent(j+1,j+1));
      if(den > epsilon) corr_N /= den;
      else              corr_N  = -999.0;

      double corr_B  = h_EgBRF_AllBkgPrediction_Cov_Pi0EtaBkgmESFix_alphaP->GetBinContent(i+1,j+1);
      den            = sqrt(h_EgBRF_AllBkgPrediction_Cov_Pi0EtaBkgmESFix_alphaP->GetBinContent(i+1,i+1));
      den           *= sqrt(h_EgBRF_AllBkgPrediction_Cov_Pi0EtaBkgmESFix_alphaP->GetBinContent(j+1,j+1));
      if(den > epsilon) corr_B /= den;
      else              corr_B  = -999.0;

      cout << i+1 << "  " << j+1 << "  "
	   << corr_N << "  "
	   << corr_B << "  "
	   << endl;
    }
  }
  cout << endl;

  cout << endl;
  cout << "mESFix-BBComb:" << endl;
  for(int i=0;i<h_SigPhoton_EBrest_mESFit.GetXaxis()->GetNbins();i++) {
    for(int j=0;j<h_SigPhoton_EBrest_mESFit.GetXaxis()->GetNbins();j++) {
      if(i>j) continue;
      double den;
      double corr_N  = h_mESFitYield_vs_EgBRF_Cov_mESFix_BBComb->GetBinContent(i+1,j+1);
      den            = sqrt(h_mESFitYield_vs_EgBRF_Cov_mESFix_BBComb->GetBinContent(i+1,i+1));
      den           *= sqrt(h_mESFitYield_vs_EgBRF_Cov_mESFix_BBComb->GetBinContent(j+1,j+1));
      if(den > epsilon) corr_N /= den;
      else              corr_N  = -999.0;

      double corr_B  = h_EgBRF_AllBkgPrediction_Cov_Pi0EtaBkgmESFix_BBComb->GetBinContent(i+1,j+1);
      den         = sqrt(h_EgBRF_AllBkgPrediction_Cov_Pi0EtaBkgmESFix_BBComb->GetBinContent(i+1,i+1));
      den           *= sqrt(h_EgBRF_AllBkgPrediction_Cov_Pi0EtaBkgmESFix_BBComb->GetBinContent(j+1,j+1));
      if(den > epsilon) corr_B /= den;
      else              corr_B  = -999.0;

      cout << i+1 << "  " << j+1 << "  "
	   << corr_N << "  "
	   << corr_B << "  "
	   << endl;
    }
  }
  cout << endl;

  cout << endl;
  cout << "mESFix-slope:" << endl;
  for(int i=0;i<h_SigPhoton_EBrest_mESFit.GetXaxis()->GetNbins();i++) {
    for(int j=0;j<h_SigPhoton_EBrest_mESFit.GetXaxis()->GetNbins();j++) {
      if(i>j) continue;
      double den;
      double corr_N  = h_mESFitYield_vs_EgBRF_Cov_mESFix_slope->GetBinContent(i+1,j+1);
      den            = sqrt(h_mESFitYield_vs_EgBRF_Cov_mESFix_slope->GetBinContent(i+1,i+1));
      den           *= sqrt(h_mESFitYield_vs_EgBRF_Cov_mESFix_slope->GetBinContent(j+1,j+1));
      if(den > epsilon) corr_N /= den;
      else              corr_N  = -999.0;

      double corr_B  = h_EgBRF_AllBkgPrediction_Cov_Pi0EtaBkgmESFix_slope->GetBinContent(i+1,j+1);
      den            = sqrt(h_EgBRF_AllBkgPrediction_Cov_Pi0EtaBkgmESFix_slope->GetBinContent(i+1,i+1));
      den           *= sqrt(h_EgBRF_AllBkgPrediction_Cov_Pi0EtaBkgmESFix_slope->GetBinContent(j+1,j+1));
      if(den > epsilon) corr_B /= den;
      else              corr_B  = -999.0;

      cout << i+1 << "  " << j+1 << "  "
	   << corr_N << "  "
	   << corr_B << "  "
	   << endl;
    }
  }
  cout << endl;

  cout << endl;
  cout << "mESFit-Bias:" << endl;
  for(int i=0;i<h_SigPhoton_EBrest_mESFit.GetXaxis()->GetNbins();i++) {
    for(int j=0;j<h_SigPhoton_EBrest_mESFit.GetXaxis()->GetNbins();j++) {
      if(i>j) continue;
      double den;
      double corr_N  = h_mESFitYield_vs_EgBRF_Cov_mESFitBias->GetBinContent(i+1,j+1);
      den            = sqrt(h_mESFitYield_vs_EgBRF_Cov_mESFitBias->GetBinContent(i+1,i+1));
      den           *= sqrt(h_mESFitYield_vs_EgBRF_Cov_mESFitBias->GetBinContent(j+1,j+1));
      if(den > epsilon) corr_N /= den;
      else              corr_N  = -999.0;

      double corr_B  = h_EgBRF_AllBkgPrediction_Cov_Pi0EtaBkgmESFitBias->GetBinContent(i+1,j+1);
      den            = sqrt(h_EgBRF_AllBkgPrediction_Cov_Pi0EtaBkgmESFitBias->GetBinContent(i+1,i+1));
      den           *= sqrt(h_EgBRF_AllBkgPrediction_Cov_Pi0EtaBkgmESFitBias->GetBinContent(j+1,j+1));
      if(den > epsilon) corr_B /= den;
      else              corr_B  = -999.0;

      cout << i+1 << "  " << j+1 << "  "
	   << corr_N << "  "
	   << corr_B << "  "
	   << endl;
    }
  }
  cout << endl;

  for(int i=0;i<h_SigPhoton_EBrest_mESFit.GetXaxis()->GetNbins();i++) {
    for(int j=0;j<h_SigPhoton_EBrest_mESFit.GetXaxis()->GetNbins();j++) {
      double cov = 0.0;
      double tmp;
      double corr;
      double corr_N,corr_B;
      double den;
      double sigma_Ni,sigma_Nj,sigma_Bi,sigma_Bj;

      //Stat:
      cov  = 0.0;
      cov += h_mESFitYield_vs_EgBRF_Cov_Stat->GetBinContent(i+1,j+1);
      h_SignalYield_vs_EgBRF_Cov_Stat->SetBinContent(i+1,j+1,cov);

      //Stat BBkg:
      cov  = 0.0;
      cov += h_EgBRF_AllBkgPrediction_Cov_AllStat->GetBinContent(i+1,j+1);
      h_SignalYield_vs_EgBRF_Cov_BBkgStat->SetBinContent(i+1,j+1,cov);

      //mESFixX1:
      sigma_Ni = sqrt(h_mESFitYield_vs_EgBRF_Cov_mESFixX1->GetBinContent(i+1,i+1));
      sigma_Nj = sqrt(h_mESFitYield_vs_EgBRF_Cov_mESFixX1->GetBinContent(j+1,j+1));
      sigma_Bi = sqrt(h_EgBRF_AllBkgPrediction_Cov_Pi0EtaBkgmESFixX1->GetBinContent(i+1,i+1));
      sigma_Bj = sqrt(h_EgBRF_AllBkgPrediction_Cov_Pi0EtaBkgmESFixX1->GetBinContent(j+1,j+1));
      corr_N  = h_mESFitYield_vs_EgBRF_Cov_mESFixX1->GetBinContent(i+1,j+1);
      den     = sigma_Ni;
      den    *= sigma_Nj;
      if(den > epsilon) corr_N /= den;
      else              corr_N  = -999.0;

      corr_B  = h_EgBRF_AllBkgPrediction_Cov_Pi0EtaBkgmESFixX1->GetBinContent(i+1,j+1);
      den     = sigma_Bi;
      den    *= sigma_Bj;
      if(den > epsilon) corr_B /= den;
      else              corr_B  = -999.0;

      if(corr_B != -999.0)      corr = corr_B;
      else if(corr_N != -999.0) corr = corr_N;
      else                      corr = 1.0;

      cov  = 0.0;
      cov += corr*(sigma_Ni-sigma_Bi)*(sigma_Nj-sigma_Bj);
      h_SignalYield_vs_EgBRF_Cov_mESFixX1->SetBinContent(i+1,j+1,cov);

      //mESFixX2:
      sigma_Ni = sqrt(h_mESFitYield_vs_EgBRF_Cov_mESFixX2->GetBinContent(i+1,i+1));
      sigma_Nj = sqrt(h_mESFitYield_vs_EgBRF_Cov_mESFixX2->GetBinContent(j+1,j+1));
      sigma_Bi = sqrt(h_EgBRF_AllBkgPrediction_Cov_Pi0EtaBkgmESFixX2->GetBinContent(i+1,i+1));
      sigma_Bj = sqrt(h_EgBRF_AllBkgPrediction_Cov_Pi0EtaBkgmESFixX2->GetBinContent(j+1,j+1));
      corr_N  = h_mESFitYield_vs_EgBRF_Cov_mESFixX2->GetBinContent(i+1,j+1);
      den     = sigma_Ni;
      den    *= sigma_Nj;
      if(den > epsilon) corr_N /= den;
      else              corr_N  = -999.0;

      corr_B  = h_EgBRF_AllBkgPrediction_Cov_Pi0EtaBkgmESFixX2->GetBinContent(i+1,j+1);
      den     = sigma_Bi;
      den    *= sigma_Bj;
      if(den > epsilon) corr_B /= den;
      else              corr_B  = -999.0;

      if(corr_B != -999.0)      corr = corr_B;
      else if(corr_N != -999.0) corr = corr_N;
      else                      corr = 1.0;

      cov  = 0.0;
      cov += corr*(sigma_Ni-sigma_Bi)*(sigma_Nj-sigma_Bj);
      h_SignalYield_vs_EgBRF_Cov_mESFixX2->SetBinContent(i+1,j+1,cov);

      //mESFixX3:
      sigma_Ni = sqrt(h_mESFitYield_vs_EgBRF_Cov_mESFixX3->GetBinContent(i+1,i+1));
      sigma_Nj = sqrt(h_mESFitYield_vs_EgBRF_Cov_mESFixX3->GetBinContent(j+1,j+1));
      sigma_Bi = sqrt(h_EgBRF_AllBkgPrediction_Cov_Pi0EtaBkgmESFixX3->GetBinContent(i+1,i+1));
      sigma_Bj = sqrt(h_EgBRF_AllBkgPrediction_Cov_Pi0EtaBkgmESFixX3->GetBinContent(j+1,j+1));
      corr_N  = h_mESFitYield_vs_EgBRF_Cov_mESFixX3->GetBinContent(i+1,j+1);
      den     = sigma_Ni;
      den    *= sigma_Nj;
      if(den > epsilon) corr_N /= den;
      else              corr_N  = -999.0;

      corr_B  = h_EgBRF_AllBkgPrediction_Cov_Pi0EtaBkgmESFixX3->GetBinContent(i+1,j+1);
      den     = sigma_Bi;
      den    *= sigma_Bj;
      if(den > epsilon) corr_B /= den;
      else              corr_B  = -999.0;

      if(corr_B != -999.0)      corr = corr_B;
      else if(corr_N != -999.0) corr = corr_N;
      else                      corr = 1.0;

      cov  = 0.0;
      cov += corr*(sigma_Ni-sigma_Bi)*(sigma_Nj-sigma_Bj);
      h_SignalYield_vs_EgBRF_Cov_mESFixX3->SetBinContent(i+1,j+1,cov);

      //mESFix_alphaP:
      sigma_Ni = sqrt(h_mESFitYield_vs_EgBRF_Cov_mESFix_alphaP->GetBinContent(i+1,i+1));
      sigma_Nj = sqrt(h_mESFitYield_vs_EgBRF_Cov_mESFix_alphaP->GetBinContent(j+1,j+1));
      sigma_Bi = sqrt(h_EgBRF_AllBkgPrediction_Cov_Pi0EtaBkgmESFix_alphaP->GetBinContent(i+1,i+1));
      sigma_Bj = sqrt(h_EgBRF_AllBkgPrediction_Cov_Pi0EtaBkgmESFix_alphaP->GetBinContent(j+1,j+1));
      corr_N  = h_mESFitYield_vs_EgBRF_Cov_mESFix_alphaP->GetBinContent(i+1,j+1);
      den     = sigma_Ni;
      den    *= sigma_Nj;
      if(den > epsilon) corr_N /= den;
      else              corr_N  = -999.0;

      corr_B  = h_EgBRF_AllBkgPrediction_Cov_Pi0EtaBkgmESFix_alphaP->GetBinContent(i+1,j+1);
      den     = sigma_Bi;
      den    *= sigma_Bj;
      if(den > epsilon) corr_B /= den;
      else              corr_B  = -999.0;

      if(corr_B != -999.0)      corr = corr_B;
      else if(corr_N != -999.0) corr = corr_N;
      else                      corr = 1.0;

      cov  = 0.0;
      cov += corr*(sigma_Ni-sigma_Bi)*(sigma_Nj-sigma_Bj);
      h_SignalYield_vs_EgBRF_Cov_mESFix_alphaP->SetBinContent(i+1,j+1,cov);

      //mESFix_BBComb:
      sigma_Ni = sqrt(h_mESFitYield_vs_EgBRF_Cov_mESFix_BBComb->GetBinContent(i+1,i+1));
      sigma_Nj = sqrt(h_mESFitYield_vs_EgBRF_Cov_mESFix_BBComb->GetBinContent(j+1,j+1));
      sigma_Bi = sqrt(h_EgBRF_AllBkgPrediction_Cov_Pi0EtaBkgmESFix_BBComb->GetBinContent(i+1,i+1));
      sigma_Bj = sqrt(h_EgBRF_AllBkgPrediction_Cov_Pi0EtaBkgmESFix_BBComb->GetBinContent(j+1,j+1));
      corr_N  = h_mESFitYield_vs_EgBRF_Cov_mESFix_BBComb->GetBinContent(i+1,j+1);
      den     = sigma_Ni;
      den    *= sigma_Nj;
      if(den > epsilon) corr_N /= den;
      else              corr_N  = -999.0;

      corr_B  = h_EgBRF_AllBkgPrediction_Cov_Pi0EtaBkgmESFix_BBComb->GetBinContent(i+1,j+1);
      den     = sigma_Bi;
      den    *= sigma_Bj;
      if(den > epsilon) corr_B /= den;
      else              corr_B  = -999.0;

      if(corr_B != -999.0)      corr = corr_B;
      else if(corr_N != -999.0) corr = corr_N;
      else                      corr = 1.0;

      cov  = 0.0;
      cov += corr*(sigma_Ni-sigma_Bi)*(sigma_Nj-sigma_Bj);
      h_SignalYield_vs_EgBRF_Cov_mESFix_BBComb->SetBinContent(i+1,j+1,cov);


      //mESFix_slope:
      sigma_Ni = sqrt(h_mESFitYield_vs_EgBRF_Cov_mESFix_slope->GetBinContent(i+1,i+1));
      sigma_Nj = sqrt(h_mESFitYield_vs_EgBRF_Cov_mESFix_slope->GetBinContent(j+1,j+1));
      sigma_Bi = sqrt(h_EgBRF_AllBkgPrediction_Cov_Pi0EtaBkgmESFix_slope->GetBinContent(i+1,i+1));
      sigma_Bj = sqrt(h_EgBRF_AllBkgPrediction_Cov_Pi0EtaBkgmESFix_slope->GetBinContent(j+1,j+1));
      corr_N  = h_mESFitYield_vs_EgBRF_Cov_mESFix_slope->GetBinContent(i+1,j+1);
      den     = sigma_Ni;
      den    *= sigma_Nj;
      if(den > epsilon) corr_N /= den;
      else              corr_N  = -999.0;

      corr_B  = h_EgBRF_AllBkgPrediction_Cov_Pi0EtaBkgmESFix_slope->GetBinContent(i+1,j+1);
      den     = sigma_Bi;
      den    *= sigma_Bj;
      if(den > epsilon) corr_B /= den;
      else              corr_B  = -999.0;

      if(corr_B != -999.0)      corr = corr_B;
      else if(corr_N != -999.0) corr = corr_N;
      else                      corr = 1.0;

      cov  = 0.0;
      cov += corr*(sigma_Ni-sigma_Bi)*(sigma_Nj-sigma_Bj);
      h_SignalYield_vs_EgBRF_Cov_mESFix_slope->SetBinContent(i+1,j+1,cov);

      //mESFitBias:
      sigma_Ni = sqrt(h_mESFitYield_vs_EgBRF_Cov_mESFitBias->GetBinContent(i+1,i+1));
      sigma_Nj = sqrt(h_mESFitYield_vs_EgBRF_Cov_mESFitBias->GetBinContent(j+1,j+1));
      sigma_Bi = sqrt(h_EgBRF_AllBkgPrediction_Cov_Pi0EtaBkgmESFitBias->GetBinContent(i+1,i+1));
      sigma_Bj = sqrt(h_EgBRF_AllBkgPrediction_Cov_Pi0EtaBkgmESFitBias->GetBinContent(j+1,j+1));
      corr_N  = h_mESFitYield_vs_EgBRF_Cov_mESFitBias->GetBinContent(i+1,j+1);
      den     = sigma_Ni;
      den    *= sigma_Nj;
      if(den > epsilon) corr_N /= den;
      else              corr_N  = -999.0;

      corr_B  = h_EgBRF_AllBkgPrediction_Cov_Pi0EtaBkgmESFitBias->GetBinContent(i+1,j+1);
      den     = sigma_Bi;
      den    *= sigma_Bj;
      if(den > epsilon) corr_B /= den;
      else              corr_B  = -999.0;

      if(corr_B != -999.0)      corr = corr_B;
      else if(corr_N != -999.0) corr = corr_N;
      else                      corr = 1.0;

      cov  = 0.0;
      cov += corr*(sigma_Ni-sigma_Bi)*(sigma_Nj-sigma_Bj);
      h_SignalYield_vs_EgBRF_Cov_mESFitBias->SetBinContent(i+1,j+1,cov);

      //HE-gamma efficiency:
      cov  = 0.0;
      cov += h_EgBRF_AllBkgPrediction_Cov_AllHEEfficCorr->GetBinContent(i+1,j+1);
      h_SignalYield_vs_EgBRF_Cov_AllHEEfficCorr->SetBinContent(i+1,j+1,cov);

      //Breco correction:
      cov  = 0.0;
      cov += h_EgBRF_AllBkgPrediction_Cov_NonPi0EtaBkgBrecoCorr->GetBinContent(i+1,j+1);
      h_SignalYield_vs_EgBRF_Cov_BrecoCorr->SetBinContent(i+1,j+1,cov);

      //pi0/eta bkg, MC-yields:
      cov  = 0.0;
      cov += h_EgBRF_AllBkgPrediction_Cov_Pi0EtaBkgMCYield->GetBinContent(i+1,j+1);
      h_SignalYield_vs_EgBRF_Cov_Pi0EtaBkgMCYield->SetBinContent(i+1,j+1,cov);

      //pi0/eta bkg, Rpi0/Reta:
      cov  = 0.0;
      cov += h_EgBRF_AllBkgPrediction_Cov_Pi0EtaBkgRpi0Reta->GetBinContent(i+1,j+1);
      h_SignalYield_vs_EgBRF_Cov_Pi0EtaBkgRpi0Reta->SetBinContent(i+1,j+1,cov);

      //pi0/eta bkg, Rpi0C/RetaC:
      cov  = 0.0;
      cov += h_EgBRF_AllBkgPrediction_Cov_Pi0EtaBkgRpi0CRetaC->GetBinContent(i+1,j+1);
      h_SignalYield_vs_EgBRF_Cov_Pi0EtaBkgRpi0CRetaC->SetBinContent(i+1,j+1,cov);

      //pi0/eta bkg, mES fit Yields:
      cov  = 0.0;
      cov += h_EgBRF_AllBkgPrediction_Cov_Pi0EtaBkgmESFitYields->GetBinContent(i+1,j+1);
      h_SignalYield_vs_EgBRF_Cov_Pi0EtaBkgmESFitYields->SetBinContent(i+1,j+1,cov);

      //pi0/eta bkg, Mgg Val syst:
      cov  = 0.0;
      cov += h_EgBRF_AllBkgPrediction_Cov_Pi0EtaBkgMggFitValSyst->GetBinContent(i+1,j+1);
      h_SignalYield_vs_EgBRF_Cov_Pi0EtaBkgMggFitValSyst->SetBinContent(i+1,j+1,cov);

      //pi0/eta bkg, LE-gamma efficiency base-line:
      cov  = 0.0;
      cov += h_EgBRF_AllBkgPrediction_Cov_Pi0EtaBkgLE_Base->GetBinContent(i+1,j+1);
      h_SignalYield_vs_EgBRF_Cov_Pi0EtaBkgLE_Base->SetBinContent(i+1,j+1,cov);

      //pi0/eta bkg, ULE-gamma efficiency:
      cov  = 0.0;
      cov += h_EgBRF_AllBkgPrediction_Cov_Pi0EtaBkgULE->GetBinContent(i+1,j+1);
      h_SignalYield_vs_EgBRF_Cov_Pi0EtaBkgULE->SetBinContent(i+1,j+1,cov);

      //Omega corrections:
      cov  = 0.0;
      cov += h_EgBRF_AllBkgPrediction_Cov_OmegaBkgAlpha->GetBinContent(i+1,j+1);
      h_SignalYield_vs_EgBRF_Cov_OmegaBkgAlpha->SetBinContent(i+1,j+1,cov);

      //Eta' corrections:
      cov  = 0.0;
      cov += h_EgBRF_AllBkgPrediction_Cov_etaprimeBkgAlpha->GetBinContent(i+1,j+1);
      h_SignalYield_vs_EgBRF_Cov_etaprimeBkgAlpha->SetBinContent(i+1,j+1,cov);

      //Semi-Leptonic corrections:
      cov  = 0.0;
      cov += h_EgBRF_AllBkgPrediction_Cov_BkgSLAlpha->GetBinContent(i+1,j+1);
      h_SignalYield_vs_EgBRF_Cov_BkgSLAlpha->SetBinContent(i+1,j+1,cov);

      //Bremss Material syst:
      cov  = 0.0;
      cov += h_EgBRF_AllBkgPrediction_Cov_BremsBkgMaterial->GetBinContent(i+1,j+1);
      h_SignalYield_vs_EgBRF_Cov_BremsBkgMaterial->SetBinContent(i+1,j+1,cov);

      //Tracking ineffic. syst:
      cov  = 0.0;
      cov += h_EgBRF_AllBkgPrediction_Cov_ElecBkgTrackIneffic->GetBinContent(i+1,j+1);
      h_SignalYield_vs_EgBRF_Cov_ElecBkgTrackIneffic->SetBinContent(i+1,j+1,cov);

      //anti-n0 corrections:
      cov  = 0.0;
      cov += h_EgBRF_AllBkgPrediction_Cov_AntiNeutronsBkgAlpha->GetBinContent(i+1,j+1);
      h_SignalYield_vs_EgBRF_Cov_AntiNeutronsBkgAlpha->SetBinContent(i+1,j+1,cov);

      //BDT syst:
      cov  = 0.0;
      cov += h_EgBRF_AllBkgPrediction_Cov_BDT->GetBinContent(i+1,j+1);
      h_SignalYield_vs_EgBRF_Cov_BDT->SetBinContent(i+1,j+1,cov);

      //pi0/eta veto, mass line-shape:
      cov  = 0.0;
      cov += h_EgBRF_AllBkgPrediction_Cov_Pi0EtaVeto_MassLineShape->GetBinContent(i+1,j+1);
      h_SignalYield_vs_EgBRF_Cov_Pi0EtaVeto_MassLineShape->SetBinContent(i+1,j+1,cov);

      //pi0/eta veto, 2nd-gamma mult.:
      cov  = 0.0;
      cov += h_EgBRF_AllBkgPrediction_Cov_Pi0EtaVeto_2ndGammaMult->GetBinContent(i+1,j+1);
      h_SignalYield_vs_EgBRF_Cov_Pi0EtaVeto_2ndGammaMult->SetBinContent(i+1,j+1,cov);
    }
  }

  for(int i=0;i<h_SigPhoton_EBrest_mESFit.GetXaxis()->GetNbins();i++) {
    double diff = h_SigPhoton_EBrest_mESFit.GetBinContent(i+1);
    diff -= h_EgBRF_AllBkgPrediction_tot_corr->GetBinContent(i+1);
    h_SigPhoton_EBrest_SignalYield.SetBinContent(i+1,diff);
    h_SigPhoton_EBrest_SignalYield_StatOnly.SetBinContent(i+1,diff);

    if(_IsMC > 1) {
      diff =  h_SigPhoton_EBrest_Fit_BrecoTM.GetBinContent(i+1);
      diff -= h_EgBRF_AllBkgPrediction_tot_corr->GetBinContent(i+1);
      h_SigPhoton_EBrest_SignalYield_TrueSubTract.SetBinContent(i+1,diff);
    }

    double e = 0.0;
    e += h_SignalYield_vs_EgBRF_Cov_Stat->GetBinContent(i+1,i+1);
    e += h_SignalYield_vs_EgBRF_Cov_BBkgStat->GetBinContent(i+1,i+1);
    e += h_SignalYield_vs_EgBRF_Cov_AllHEEfficCorr->GetBinContent(i+1,i+1);
    e += h_SignalYield_vs_EgBRF_Cov_BrecoCorr->GetBinContent(i+1,i+1);
    e += h_SignalYield_vs_EgBRF_Cov_Pi0EtaBkgMCYield->GetBinContent(i+1,i+1);
    e += h_SignalYield_vs_EgBRF_Cov_Pi0EtaBkgRpi0Reta->GetBinContent(i+1,i+1);
    e += h_SignalYield_vs_EgBRF_Cov_Pi0EtaBkgRpi0CRetaC->GetBinContent(i+1,i+1);
    e += h_SignalYield_vs_EgBRF_Cov_Pi0EtaBkgmESFitYields->GetBinContent(i+1,i+1);
    e += h_SignalYield_vs_EgBRF_Cov_mESFixX1->GetBinContent(i+1,i+1);
    e += h_SignalYield_vs_EgBRF_Cov_mESFixX2->GetBinContent(i+1,i+1);
    e += h_SignalYield_vs_EgBRF_Cov_mESFixX3->GetBinContent(i+1,i+1);
    e += h_SignalYield_vs_EgBRF_Cov_mESFix_alphaP->GetBinContent(i+1,i+1);
    e += h_SignalYield_vs_EgBRF_Cov_mESFix_BBComb->GetBinContent(i+1,i+1);
    e += h_SignalYield_vs_EgBRF_Cov_mESFix_slope->GetBinContent(i+1,i+1);
    e += h_SignalYield_vs_EgBRF_Cov_mESFitBias->GetBinContent(i+1,i+1);
    e += h_SignalYield_vs_EgBRF_Cov_Pi0EtaBkgMggFitValSyst->GetBinContent(i+1,i+1);
    e += h_SignalYield_vs_EgBRF_Cov_Pi0EtaBkgLE_Base->GetBinContent(i+1,i+1);
    e += h_SignalYield_vs_EgBRF_Cov_Pi0EtaBkgULE->GetBinContent(i+1,i+1);
    e += h_SignalYield_vs_EgBRF_Cov_OmegaBkgAlpha->GetBinContent(i+1,i+1);
    e += h_SignalYield_vs_EgBRF_Cov_etaprimeBkgAlpha->GetBinContent(i+1,i+1);
    e += h_SignalYield_vs_EgBRF_Cov_BkgSLAlpha->GetBinContent(i+1,i+1);
    e += h_SignalYield_vs_EgBRF_Cov_BremsBkgMaterial->GetBinContent(i+1,i+1);
    e += h_SignalYield_vs_EgBRF_Cov_ElecBkgTrackIneffic->GetBinContent(i+1,i+1);
    e += h_SignalYield_vs_EgBRF_Cov_AntiNeutronsBkgAlpha->GetBinContent(i+1,i+1);
    e += h_SignalYield_vs_EgBRF_Cov_BDT->GetBinContent(i+1,i+1);
    e += h_SignalYield_vs_EgBRF_Cov_Pi0EtaVeto_MassLineShape->GetBinContent(i+1,i+1);
    e += h_SignalYield_vs_EgBRF_Cov_Pi0EtaVeto_2ndGammaMult->GetBinContent(i+1,i+1);
    e  = sqrt(e);
    h_SigPhoton_EBrest_SignalYield.SetBinError(i+1,e);

    e  = 0.0;
    e += h_SignalYield_vs_EgBRF_Cov_Stat->GetBinContent(i+1,i+1);
    e += h_SignalYield_vs_EgBRF_Cov_BBkgStat->GetBinContent(i+1,i+1);
    e  = sqrt(e);
    h_SigPhoton_EBrest_SignalYield_StatOnly.SetBinError(i+1,e);
  }


  cout << endl;
  cout << endl;
  cout << "========================================" << endl;
  cout << "Signal Yields results with syst. errors:" << endl;
  cout << "===========================================================================================================================" 
       << endl;
  cout << "Type                ";
  for(int i=0;i<Medium;i++) {
    double c = h_SigPhoton_EBrest_mESFit.GetBinCenter(i+1);
    double w = h_SigPhoton_EBrest_mESFit.GetBinWidth(i+1)*0.5;
    double R1,R2;
    R1 = c - w;
    R2 = c + w;
    sprintf(ytitle,"%.1f",R1);
    cout << "(" << ytitle << ",";
    sprintf(ytitle,"%.1f",R2);
    cout << ytitle << ")   ";
  }
  cout << endl;
  cout << "---------------------------------------------------------------------------------------------------------------------------"
       << endl;
  cout << "Yield               ";
  for(int i=0;i<Medium;i++) {
    sprintf(ytitle,"%.2f",h_SigPhoton_EBrest_SignalYield.GetBinContent(i+1));
    cout << ytitle
	 << "        ";
  }
  cout << endl;
  cout << "Err-Stat            ";
  for(int i=0;i<Medium;i++) {
    sprintf(ytitle,"%.2f",sqrt(h_SignalYield_vs_EgBRF_Cov_Stat->GetBinContent(i+1,i+1)));
    cout << ytitle
	 << "        ";
  }
  cout << endl;
  cout << "Err-BBkgStat        ";
  for(int i=0;i<Medium;i++) {
    sprintf(ytitle,"%.2f",sqrt(h_SignalYield_vs_EgBRF_Cov_BBkgStat->GetBinContent(i+1,i+1)));
    cout << ytitle
	 << "        ";
  }
  cout << endl;
  cout << "Err-HEgEffic        ";
  for(int i=0;i<Medium;i++) {
    sprintf(ytitle,"%.2f",sqrt(h_SignalYield_vs_EgBRF_Cov_AllHEEfficCorr->GetBinContent(i+1,i+1)));
    cout << ytitle
	 << "        ";
  }
  cout << endl;
  cout << "Err-BrecoCorr       ";
  for(int i=0;i<Medium;i++) {
    sprintf(ytitle,"%.2f",sqrt(h_SignalYield_vs_EgBRF_Cov_BrecoCorr->GetBinContent(i+1,i+1)));
    cout << ytitle
	 << "        ";
  }
  cout << endl;
  cout << "Err-Pi0EtaMCYield   ";
  for(int i=0;i<Medium;i++) {
    sprintf(ytitle,"%.2f",sqrt(h_SignalYield_vs_EgBRF_Cov_Pi0EtaBkgMCYield->GetBinContent(i+1,i+1)));
    cout << ytitle
	 << "        ";
  }
  cout << endl;
  cout << "Err-Rpi0Reta        ";
  for(int i=0;i<Medium;i++) {
    sprintf(ytitle,"%.2f",sqrt(h_SignalYield_vs_EgBRF_Cov_Pi0EtaBkgRpi0Reta->GetBinContent(i+1,i+1)));
    cout << ytitle
	 << "        ";
  }
  cout << endl;
  cout << "Err-Rpi0CRetaC      ";
  for(int i=0;i<Medium;i++) {
    sprintf(ytitle,"%.2f",sqrt(h_SignalYield_vs_EgBRF_Cov_Pi0EtaBkgRpi0CRetaC->GetBinContent(i+1,i+1)));
    cout << ytitle
	 << "        ";
  }
  cout << endl;
  cout << "Err-Pi0EtaFitYield  ";
  for(int i=0;i<Medium;i++) {
    sprintf(ytitle,"%.2f",sqrt(h_SignalYield_vs_EgBRF_Cov_Pi0EtaBkgmESFitYields->GetBinContent(i+1,i+1)));
    cout << ytitle
	 << "        ";
  }
  cout << endl;
  cout << "Err-mESFixX1        ";
  for(int i=0;i<Medium;i++) {
    sprintf(ytitle,"%.2f",sqrt(h_SignalYield_vs_EgBRF_Cov_mESFixX1->GetBinContent(i+1,i+1)));
    cout << ytitle
	 << "        ";
  }
  cout << endl;
  cout << "Err-mESFixX2        ";
  for(int i=0;i<Medium;i++) {
    sprintf(ytitle,"%.2f",sqrt(h_SignalYield_vs_EgBRF_Cov_mESFixX2->GetBinContent(i+1,i+1)));
    cout << ytitle
	 << "        ";
  }
  cout << endl;
  cout << "Err-mESFixX3        ";
  for(int i=0;i<Medium;i++) {
    sprintf(ytitle,"%.2f",sqrt(h_SignalYield_vs_EgBRF_Cov_mESFixX3->GetBinContent(i+1,i+1)));
    cout << ytitle
	 << "        ";
  }
  cout << endl;
  cout << "Err-mESFix_alphaP   ";
  for(int i=0;i<Medium;i++) {
    sprintf(ytitle,"%.2f",sqrt(h_SignalYield_vs_EgBRF_Cov_mESFix_alphaP->GetBinContent(i+1,i+1)));
    cout << ytitle
	 << "        ";
  }
  cout << endl;
  cout << "Err-mESFix_BBComb   ";
  for(int i=0;i<Medium;i++) {
    sprintf(ytitle,"%.2f",sqrt(h_SignalYield_vs_EgBRF_Cov_mESFix_BBComb->GetBinContent(i+1,i+1)));
    cout << ytitle
	 << "        ";
  }
  cout << endl;
  cout << "Err-mESFix_slope    ";
  for(int i=0;i<Medium;i++) {
    sprintf(ytitle,"%.2f",sqrt(h_SignalYield_vs_EgBRF_Cov_mESFix_slope->GetBinContent(i+1,i+1)));
    cout << ytitle
	 << "        ";
  }
  cout << endl;
  cout << "Err-mESFitBias      ";
  for(int i=0;i<Medium;i++) {
    sprintf(ytitle,"%.2f",sqrt(h_SignalYield_vs_EgBRF_Cov_mESFitBias->GetBinContent(i+1,i+1)));
    cout << ytitle
	 << "        ";
  }
  cout << endl;
  cout << "Err-MggFitSyst      ";
  for(int i=0;i<Medium;i++) {
    sprintf(ytitle,"%.2f",sqrt(h_SignalYield_vs_EgBRF_Cov_Pi0EtaBkgMggFitValSyst->GetBinContent(i+1,i+1)));
    cout << ytitle
	 << "        ";
  }
  cout << endl;
  cout << "Err-LE_Base         ";
  for(int i=0;i<Medium;i++) {
    sprintf(ytitle,"%.2f",sqrt(h_SignalYield_vs_EgBRF_Cov_Pi0EtaBkgLE_Base->GetBinContent(i+1,i+1)));
    cout << ytitle
	 << "        ";
  }
  cout << endl;
  cout << "Err-ULE             ";
  for(int i=0;i<Medium;i++) {
    sprintf(ytitle,"%.2f",sqrt(h_SignalYield_vs_EgBRF_Cov_Pi0EtaBkgULE->GetBinContent(i+1,i+1)));
    cout << ytitle
	 << "        ";
  }
  cout << endl;
  cout << "Err-Omega-Corr      ";
  for(int i=0;i<Medium;i++) {
    sprintf(ytitle,"%.2f",sqrt(h_SignalYield_vs_EgBRF_Cov_OmegaBkgAlpha->GetBinContent(i+1,i+1)));
    cout << ytitle
	 << "        ";
  }
  cout << endl;
  cout << "Err-eta'-Corr       ";
  for(int i=0;i<Medium;i++) {
    sprintf(ytitle,"%.2f",sqrt(h_SignalYield_vs_EgBRF_Cov_etaprimeBkgAlpha->GetBinContent(i+1,i+1)));
    cout << ytitle
	 << "        ";
  }
  cout << endl;
  cout << "Err-SL-Corr         ";
  for(int i=0;i<Medium;i++) {
    sprintf(ytitle,"%.2f",sqrt(h_SignalYield_vs_EgBRF_Cov_BkgSLAlpha->GetBinContent(i+1,i+1)));
    cout << ytitle
	 << "        ";
  }
  cout << endl;
  cout << "Err-Bremss-Syst     ";
  for(int i=0;i<Medium;i++) {
    sprintf(ytitle,"%.2f",sqrt(h_SignalYield_vs_EgBRF_Cov_BremsBkgMaterial->GetBinContent(i+1,i+1)));
    cout << ytitle
	 << "        ";
  }
  cout << endl;
  cout << "Err-Elec-Corr       ";
  for(int i=0;i<Medium;i++) {
    sprintf(ytitle,"%.2f",sqrt(h_SignalYield_vs_EgBRF_Cov_ElecBkgTrackIneffic->GetBinContent(i+1,i+1)));
    cout << ytitle
	 << "        ";
  }
  cout << endl;
  cout << "Err-antin0-Corr     ";
  for(int i=0;i<Medium;i++) {
    sprintf(ytitle,"%.2f",sqrt(h_SignalYield_vs_EgBRF_Cov_AntiNeutronsBkgAlpha->GetBinContent(i+1,i+1)));
    cout << ytitle
	 << "        ";
  }
  cout << endl;
  cout << "Err-BDT-syst        ";
  for(int i=0;i<Medium;i++) {
    sprintf(ytitle,"%.2f",sqrt(h_SignalYield_vs_EgBRF_Cov_BDT->GetBinContent(i+1,i+1)));
    cout << ytitle
	 << "        ";
  }
  cout << endl;
  cout << "Err-Pi0EtaVeto-LS   ";
  for(int i=0;i<Medium;i++) {
    sprintf(ytitle,"%.2f",sqrt(h_SignalYield_vs_EgBRF_Cov_Pi0EtaVeto_MassLineShape->GetBinContent(i+1,i+1)));
    cout << ytitle
	 << "        ";
  }
  cout << endl;
  cout << "Err-Pi0EtaVeto-Mult ";
  for(int i=0;i<Medium;i++) {
    sprintf(ytitle,"%.2f",sqrt(h_SignalYield_vs_EgBRF_Cov_Pi0EtaVeto_2ndGammaMult->GetBinContent(i+1,i+1)));
    cout << ytitle
	 << "        ";
  }
  cout << endl;
  cout << "---------------------------------------------------------------------------------------------------------------------------"
       << endl;
  cout << "Err-tot             ";
  for(int i=0;i<Medium;i++) {
    sprintf(ytitle,"%.2f",h_SigPhoton_EBrest_SignalYield.GetBinError(i+1));
    cout << ytitle
	 << "        ";
  }
  cout << endl;
  cout << "===========================================================================================================================" 
       << endl;
  cout << "Type                ";
  for(int i=Medium;i<Nbins_SigPhoton_EBrest_Fit;i++) {
    double c = h_SigPhoton_EBrest_mESFit.GetBinCenter(i+1);
    double w = h_SigPhoton_EBrest_mESFit.GetBinWidth(i+1)*0.5;
    double R1,R2;
    R1 = c - w;
    R2 = c + w;
    sprintf(ytitle,"%.1f",R1);
    cout << "(" << ytitle << ",";
    sprintf(ytitle,"%.1f",R2);
    cout << ytitle << ")   ";
  }
  cout << endl;
  cout << "---------------------------------------------------------------------------------------------------------------------------"
       << endl;
  cout << "Yield               ";
  for(int i=Medium;i<Nbins_SigPhoton_EBrest_Fit;i++) {
    sprintf(ytitle,"%.2f",h_SigPhoton_EBrest_SignalYield.GetBinContent(i+1));
    cout << ytitle
	 << "        ";
  }
  cout << endl;
  cout << "Err-Stat            ";
  for(int i=Medium;i<Nbins_SigPhoton_EBrest_Fit;i++) {
    sprintf(ytitle,"%.2f",sqrt(h_SignalYield_vs_EgBRF_Cov_Stat->GetBinContent(i+1,i+1)));
    cout << ytitle
	 << "        ";
  }
  cout << endl;
  cout << "Err-BBkgStat        ";
  for(int i=Medium;i<Nbins_SigPhoton_EBrest_Fit;i++) {
    sprintf(ytitle,"%.2f",sqrt(h_SignalYield_vs_EgBRF_Cov_BBkgStat->GetBinContent(i+1,i+1)));
    cout << ytitle
	 << "        ";
  }
  cout << endl;
  cout << "Err-HEgEffic        ";
  for(int i=Medium;i<Nbins_SigPhoton_EBrest_Fit;i++) {
    sprintf(ytitle,"%.2f",sqrt(h_SignalYield_vs_EgBRF_Cov_AllHEEfficCorr->GetBinContent(i+1,i+1)));
    cout << ytitle
	 << "        ";
  }
  cout << endl;
  cout << "Err-BrecoCorr       ";
  for(int i=Medium;i<Nbins_SigPhoton_EBrest_Fit;i++) {
    sprintf(ytitle,"%.2f",sqrt(h_SignalYield_vs_EgBRF_Cov_BrecoCorr->GetBinContent(i+1,i+1)));
    cout << ytitle
	 << "        ";
  }
  cout << endl;
  cout << "Err-Pi0EtaMCYield   ";
  for(int i=Medium;i<Nbins_SigPhoton_EBrest_Fit;i++) {
    sprintf(ytitle,"%.2f",sqrt(h_SignalYield_vs_EgBRF_Cov_Pi0EtaBkgMCYield->GetBinContent(i+1,i+1)));
    cout << ytitle
	 << "        ";
  }
  cout << endl;
  cout << "Err-Rpi0Reta        ";
  for(int i=Medium;i<Nbins_SigPhoton_EBrest_Fit;i++) {
    sprintf(ytitle,"%.2f",sqrt(h_SignalYield_vs_EgBRF_Cov_Pi0EtaBkgRpi0Reta->GetBinContent(i+1,i+1)));
    cout << ytitle
	 << "        ";
  }
  cout << endl;
  cout << "Err-Rpi0CRetaC      ";
  for(int i=Medium;i<Nbins_SigPhoton_EBrest_Fit;i++) {
    sprintf(ytitle,"%.2f",sqrt(h_SignalYield_vs_EgBRF_Cov_Pi0EtaBkgRpi0CRetaC->GetBinContent(i+1,i+1)));
    cout << ytitle
	 << "        ";
  }
  cout << endl;
  cout << "Err-Pi0EtaFitYield  ";
  for(int i=Medium;i<Nbins_SigPhoton_EBrest_Fit;i++) {
    sprintf(ytitle,"%.2f",sqrt(h_SignalYield_vs_EgBRF_Cov_Pi0EtaBkgmESFitYields->GetBinContent(i+1,i+1)));
    cout << ytitle
	 << "        ";
  }
  cout << endl;
  cout << "Err-mESFixX1        ";
  for(int i=Medium;i<Nbins_SigPhoton_EBrest_Fit;i++) {
    sprintf(ytitle,"%.2f",sqrt(h_SignalYield_vs_EgBRF_Cov_mESFixX1->GetBinContent(i+1,i+1)));
    cout << ytitle
	 << "        ";
  }
  cout << endl;
  cout << "Err-mESFixX2        ";
  for(int i=Medium;i<Nbins_SigPhoton_EBrest_Fit;i++) {
    sprintf(ytitle,"%.2f",sqrt(h_SignalYield_vs_EgBRF_Cov_mESFixX2->GetBinContent(i+1,i+1)));
    cout << ytitle
	 << "        ";
  }
  cout << endl;
  cout << "Err-mESFixX3        ";
  for(int i=Medium;i<Nbins_SigPhoton_EBrest_Fit;i++) {
    sprintf(ytitle,"%.2f",sqrt(h_SignalYield_vs_EgBRF_Cov_mESFixX3->GetBinContent(i+1,i+1)));
    cout << ytitle
	 << "        ";
  }
  cout << endl;
  cout << "Err-mESFix_alphaP   ";
  for(int i=Medium;i<Nbins_SigPhoton_EBrest_Fit;i++) {
    sprintf(ytitle,"%.2f",sqrt(h_SignalYield_vs_EgBRF_Cov_mESFix_alphaP->GetBinContent(i+1,i+1)));
    cout << ytitle
	 << "        ";
  }
  cout << endl;
  cout << "Err-mESFix_BBComb   ";
  for(int i=Medium;i<Nbins_SigPhoton_EBrest_Fit;i++) {
    sprintf(ytitle,"%.2f",sqrt(h_SignalYield_vs_EgBRF_Cov_mESFix_BBComb->GetBinContent(i+1,i+1)));
    cout << ytitle
	 << "        ";
  }
  cout << endl;
  cout << "Err-mESFix_slope    ";
  for(int i=Medium;i<Nbins_SigPhoton_EBrest_Fit;i++) {
    sprintf(ytitle,"%.2f",sqrt(h_SignalYield_vs_EgBRF_Cov_mESFix_slope->GetBinContent(i+1,i+1)));
    cout << ytitle
	 << "        ";
  }
  cout << endl;
  cout << "Err-mESFitBias      ";
  for(int i=Medium;i<Nbins_SigPhoton_EBrest_Fit;i++) {
    sprintf(ytitle,"%.2f",sqrt(h_SignalYield_vs_EgBRF_Cov_mESFitBias->GetBinContent(i+1,i+1)));
    cout << ytitle
	 << "        ";
  }
  cout << endl;
  cout << "Err-MggFitSyst      ";
  for(int i=Medium;i<Nbins_SigPhoton_EBrest_Fit;i++) {
    sprintf(ytitle,"%.2f",sqrt(h_SignalYield_vs_EgBRF_Cov_Pi0EtaBkgMggFitValSyst->GetBinContent(i+1,i+1)));
    cout << ytitle
	 << "        ";
  }
  cout << endl;
  cout << "Err-LE_Base         ";
  for(int i=Medium;i<Nbins_SigPhoton_EBrest_Fit;i++) {
    sprintf(ytitle,"%.2f",sqrt(h_SignalYield_vs_EgBRF_Cov_Pi0EtaBkgLE_Base->GetBinContent(i+1,i+1)));
    cout << ytitle
	 << "        ";
  }
  cout << endl;
  cout << "Err-ULE             ";
  for(int i=Medium;i<Nbins_SigPhoton_EBrest_Fit;i++) {
    sprintf(ytitle,"%.2f",sqrt(h_SignalYield_vs_EgBRF_Cov_Pi0EtaBkgULE->GetBinContent(i+1,i+1)));
    cout << ytitle
	 << "        ";
  }
  cout << endl;
  cout << "Err-Omega-Corr      ";
  for(int i=Medium;i<Nbins_SigPhoton_EBrest_Fit;i++) {
    sprintf(ytitle,"%.2f",sqrt(h_SignalYield_vs_EgBRF_Cov_OmegaBkgAlpha->GetBinContent(i+1,i+1)));
    cout << ytitle
	 << "        ";
  }
  cout << endl;
  cout << "Err-eta'-Corr       ";
  for(int i=Medium;i<Nbins_SigPhoton_EBrest_Fit;i++) {
    sprintf(ytitle,"%.2f",sqrt(h_SignalYield_vs_EgBRF_Cov_etaprimeBkgAlpha->GetBinContent(i+1,i+1)));
    cout << ytitle
	 << "        ";
  }
  cout << endl;
  cout << "Err-SL-Corr         ";
  for(int i=Medium;i<Nbins_SigPhoton_EBrest_Fit;i++) {
    sprintf(ytitle,"%.2f",sqrt(h_SignalYield_vs_EgBRF_Cov_BkgSLAlpha->GetBinContent(i+1,i+1)));
    cout << ytitle
	 << "        ";
  }
  cout << endl;
  cout << "Err-Bremss-Syst     ";
  for(int i=Medium;i<Nbins_SigPhoton_EBrest_Fit;i++) {
    sprintf(ytitle,"%.2f",sqrt(h_SignalYield_vs_EgBRF_Cov_BremsBkgMaterial->GetBinContent(i+1,i+1)));
    cout << ytitle
	 << "        ";
  }
  cout << endl;
  cout << "Err-Elec-Corr       ";
  for(int i=Medium;i<Nbins_SigPhoton_EBrest_Fit;i++) {
    sprintf(ytitle,"%.2f",sqrt(h_SignalYield_vs_EgBRF_Cov_ElecBkgTrackIneffic->GetBinContent(i+1,i+1)));
    cout << ytitle
	 << "        ";
  }
  cout << endl;
  cout << "Err-antin0-Corr     ";
  for(int i=Medium;i<Nbins_SigPhoton_EBrest_Fit;i++) {
    sprintf(ytitle,"%.2f",sqrt(h_SignalYield_vs_EgBRF_Cov_AntiNeutronsBkgAlpha->GetBinContent(i+1,i+1)));
    cout << ytitle
	 << "        ";
  }
  cout << endl;
  cout << "Err-BDT-syst        ";
  for(int i=Medium;i<Nbins_SigPhoton_EBrest_Fit;i++) {
    sprintf(ytitle,"%.2f",sqrt(h_SignalYield_vs_EgBRF_Cov_BDT->GetBinContent(i+1,i+1)));
    cout << ytitle
	 << "        ";
  }
  cout << endl;
  cout << "Err-Pi0EtaVeto-LS   ";
  for(int i=Medium;i<Nbins_SigPhoton_EBrest_Fit;i++) {
    sprintf(ytitle,"%.2f",sqrt(h_SignalYield_vs_EgBRF_Cov_Pi0EtaVeto_MassLineShape->GetBinContent(i+1,i+1)));
    cout << ytitle
	 << "        ";
  }
  cout << endl;
  cout << "Err-Pi0EtaVeto-Mult ";
  for(int i=Medium;i<Nbins_SigPhoton_EBrest_Fit;i++) {
    sprintf(ytitle,"%.2f",sqrt(h_SignalYield_vs_EgBRF_Cov_Pi0EtaVeto_2ndGammaMult->GetBinContent(i+1,i+1)));
    cout << ytitle
	 << "        ";
  }
  cout << endl;
  cout << "---------------------------------------------------------------------------------------------------------------------------"
       << endl;
  cout << "Err-tot             ";
  for(int i=Medium;i<Nbins_SigPhoton_EBrest_Fit;i++) {
    sprintf(ytitle,"%.2f",h_SigPhoton_EBrest_SignalYield.GetBinError(i+1));
    cout << ytitle
	 << "        ";
  }
  cout << endl;
  cout << "===========================================================================================================================" 
       << endl;
  cout << endl;

  TGraphErrors* gr_SigPhoton_EBrest_SignalYield = new TGraphErrors(h_SigPhoton_EBrest_SignalYield.GetXaxis()->GetNbins());
  gr_SigPhoton_EBrest_SignalYield->SetName("gr_SigPhoton_EBrest_SignalYield");
  gr_SigPhoton_EBrest_SignalYield->SetTitle(h_SigPhoton_EBrest_SignalYield.GetTitle());
  gr_SigPhoton_EBrest_SignalYield->GetXaxis()->SetTitle(h_SigPhoton_EBrest_SignalYield.GetXaxis()->GetTitle());
  gr_SigPhoton_EBrest_SignalYield->GetXaxis()->CenterTitle(true);
  gr_SigPhoton_EBrest_SignalYield->GetYaxis()->SetTitle(h_SigPhoton_EBrest_SignalYield.GetYaxis()->GetTitle());
  gr_SigPhoton_EBrest_SignalYield->GetYaxis()->CenterTitle(true);
  gr_SigPhoton_EBrest_SignalYield->SetLineColor(h_SigPhoton_EBrest_SignalYield.GetLineColor());
  gr_SigPhoton_EBrest_SignalYield->SetMarkerColor(h_SigPhoton_EBrest_SignalYield.GetMarkerColor());
  gr_SigPhoton_EBrest_SignalYield->SetMarkerStyle(h_SigPhoton_EBrest_SignalYield.GetMarkerStyle());
  gr_SigPhoton_EBrest_SignalYield->SetLineWidth(h_SigPhoton_EBrest_SignalYield.GetLineWidth());
  gr_SigPhoton_EBrest_SignalYield->GetXaxis()->SetTitleSize(TheSize);
  gr_SigPhoton_EBrest_SignalYield->GetXaxis()->SetLabelSize(TheSize);
  gr_SigPhoton_EBrest_SignalYield->GetYaxis()->SetTitleSize(TheSize);
  gr_SigPhoton_EBrest_SignalYield->GetYaxis()->SetLabelSize(TheSize);

  TGraphErrors* gr_SigPhoton_EBrest_SignalYield_StatOnly = new TGraphErrors(h_SigPhoton_EBrest_SignalYield_StatOnly.GetXaxis()->GetNbins());
  gr_SigPhoton_EBrest_SignalYield_StatOnly->SetName("gr_SigPhoton_EBrest_SignalYield_StatOnly");
  gr_SigPhoton_EBrest_SignalYield_StatOnly->SetTitle(h_SigPhoton_EBrest_SignalYield_StatOnly.GetTitle());
  gr_SigPhoton_EBrest_SignalYield_StatOnly->GetXaxis()->SetTitle(h_SigPhoton_EBrest_SignalYield_StatOnly.GetXaxis()->GetTitle());
  gr_SigPhoton_EBrest_SignalYield_StatOnly->GetXaxis()->CenterTitle(true);
  gr_SigPhoton_EBrest_SignalYield_StatOnly->GetYaxis()->SetTitle(h_SigPhoton_EBrest_SignalYield_StatOnly.GetYaxis()->GetTitle());
  gr_SigPhoton_EBrest_SignalYield_StatOnly->GetYaxis()->CenterTitle(true);
  gr_SigPhoton_EBrest_SignalYield_StatOnly->SetLineColor(h_SigPhoton_EBrest_SignalYield_StatOnly.GetLineColor());
  gr_SigPhoton_EBrest_SignalYield_StatOnly->SetMarkerColor(h_SigPhoton_EBrest_SignalYield_StatOnly.GetMarkerColor());
  gr_SigPhoton_EBrest_SignalYield_StatOnly->SetMarkerStyle(h_SigPhoton_EBrest_SignalYield_StatOnly.GetMarkerStyle());
  gr_SigPhoton_EBrest_SignalYield_StatOnly->SetLineWidth(h_SigPhoton_EBrest_SignalYield_StatOnly.GetLineWidth());
  gr_SigPhoton_EBrest_SignalYield_StatOnly->GetXaxis()->SetTitleSize(TheSize);
  gr_SigPhoton_EBrest_SignalYield_StatOnly->GetXaxis()->SetLabelSize(TheSize);
  gr_SigPhoton_EBrest_SignalYield_StatOnly->GetYaxis()->SetTitleSize(TheSize);
  gr_SigPhoton_EBrest_SignalYield_StatOnly->GetYaxis()->SetLabelSize(TheSize);

  for(int i=0;i<h_SigPhoton_EBrest_SignalYield.GetXaxis()->GetNbins();i++) {
    double v,e;
    double c,w;
    c = h_SigPhoton_EBrest_SignalYield.GetBinCenter(i+1);
    w = h_SigPhoton_EBrest_SignalYield.GetBinWidth(i+1)*0.5;

    v = h_SigPhoton_EBrest_SignalYield.GetBinContent(i+1);
    e = h_SigPhoton_EBrest_SignalYield.GetBinError(i+1);

    gr_SigPhoton_EBrest_SignalYield->SetPoint(i,c,v);
    gr_SigPhoton_EBrest_SignalYield->SetPointError(i,w,e);

    v = h_SigPhoton_EBrest_SignalYield_StatOnly.GetBinContent(i+1);
    e = h_SigPhoton_EBrest_SignalYield_StatOnly.GetBinError(i+1);

    gr_SigPhoton_EBrest_SignalYield_StatOnly->SetPoint(i,c,v);
    gr_SigPhoton_EBrest_SignalYield_StatOnly->SetPointError(i,w,e);
  }


  if(_IsMC > 1) {
    c1->Clear();
    Maximum = TMath::Max(h_SigPhoton_EBrest_SignalYield_TrueSubTract.GetMaximum(),
			 h_SigPhoton_EBrest_Fit_BrecoTM_SigGammaTM.GetMaximum());
    Minimum = TMath::Min(h_SigPhoton_EBrest_SignalYield_TrueSubTract.GetMinimum(),
			 h_SigPhoton_EBrest_Fit_BrecoTM_SigGammaTM.GetMinimum());
    porcent = 0.10;
    h_SigPhoton_EBrest_SignalYield_TrueSubTract.SetMaximum(Maximum + porcent*(Maximum-Minimum));
    h_SigPhoton_EBrest_SignalYield_TrueSubTract.SetMinimum(Minimum - porcent*(Maximum-Minimum));
    h_SigPhoton_EBrest_Fit_BrecoTM_SigGammaTM.SetMaximum(Maximum + porcent*(Maximum-Minimum));
    h_SigPhoton_EBrest_Fit_BrecoTM_SigGammaTM.SetMinimum(Minimum - porcent*(Maximum-Minimum));
    h_SigPhoton_EBrest_SignalYield_TrueSubTract.Draw();
    h_SigPhoton_EBrest_Fit_BrecoTM_SigGammaTM.Draw("same");
    c1->Print(EPSName.Data());
  }


  h_SigPhoton_EBrest_Fit_BrecoTM_SigGammaTM.SetStats(false);
  //Signal yield:
  c1->Clear();
  Maximum = -1.0e+20;
  Minimum =  1.0e+20;
  for(int i=0;i<h_SigPhoton_EBrest_SignalYield.GetXaxis()->GetNbins();i++) {
    double v,e;
    v = h_SigPhoton_EBrest_SignalYield.GetBinContent(i+1);
    e = h_SigPhoton_EBrest_SignalYield.GetBinError(i+1);
    if(Maximum < v+e) Maximum = v+e;
    if(Minimum > v-e) Minimum = v-e;
  }
  porcent = 0.20;
  h_SigPhoton_EBrest_SignalYield.SetMaximum(Maximum + porcent*(Maximum-Minimum));
  h_SigPhoton_EBrest_SignalYield.SetMinimum(Minimum - porcent*(Maximum-Minimum));
  gr_SigPhoton_EBrest_SignalYield->SetMaximum(Maximum + porcent*(Maximum-Minimum));
  gr_SigPhoton_EBrest_SignalYield->SetMinimum(Minimum - porcent*(Maximum-Minimum));
  h_SigPhoton_EBrest_SignalYield.Draw();
  //h_SigPhoton_EBrest_SignalYield_StatOnly.Draw("same");
  gr_SigPhoton_EBrest_SignalYield->Draw("P");
  gr_SigPhoton_EBrest_SignalYield_StatOnly->Draw("P");
  if(_IsMC > 1) h_SigPhoton_EBrest_Fit_BrecoTM_SigGammaTM.Draw("same");
  TLine* l_3 = new TLine(h_SigPhoton_EBrest_SignalYield.GetXaxis()->GetXmin(),0.0,
			 h_SigPhoton_EBrest_SignalYield.GetXaxis()->GetXmax(),0.0);
  l_3->SetLineColor(1);
  l_3->SetLineWidth(2);
  l_3->SetLineStyle(2);
  if(h_SigPhoton_EBrest_SignalYield.GetMinimum() < 0.0) l_3->Draw();
  TLine* l_sig1_del = new TLine(R_EgBRF_Sig[0],h_SigPhoton_EBrest_SignalYield.GetMaximum(),
				R_EgBRF_Sig[0],h_SigPhoton_EBrest_SignalYield.GetMinimum());
  l_sig1_del->SetLineColor(2);
  l_sig1_del->SetLineWidth(2);
  l_sig1_del->SetLineStyle(2);
  l_sig1_del->Draw();
  TLine* l_sig2_del = new TLine(R_EgBRF_Sig[1],h_SigPhoton_EBrest_SignalYield.GetMaximum(),
				R_EgBRF_Sig[1],h_SigPhoton_EBrest_SignalYield.GetMinimum());
  l_sig2_del->SetLineColor(2);
  l_sig2_del->SetLineWidth(2);
  l_sig2_del->SetLineStyle(2);
  l_sig2_del->Draw();
  TLegend* leg_sig = new TLegend(0.75,0.75,0.9,0.9);
  leg_sig->SetFillColor(10);
  if(_IsMC > 1) leg_sig->AddEntry(&h_SigPhoton_EBrest_Fit_BrecoTM_SigGammaTM,"MC","l");
  leg_sig->AddEntry(&h_SigPhoton_EBrest_SignalYield_StatOnly,  "Fit stat",     "l");
  leg_sig->AddEntry(&h_SigPhoton_EBrest_SignalYield,           "Fit stat+syst","l");
  leg_sig->Draw("same");
  c1->Print(EPSName.Data());


  if(_IsMC > 1) {
    //Signal yield:
    c1->Clear();
    h_SigPhoton_EBrest_SignalYield.Draw();
    gr_SigPhoton_EBrest_SignalYield->Draw("P");
    gr_SigPhoton_EBrest_SignalYield_StatOnly->Draw("P");
    h_SigPhoton_EBrest_Fit_BrecoTM_SigGammaTM.Draw("same");
    h_SigPhoton_EBrest_SignalYield_TrueSubTract.Draw("same");
    if(h_SigPhoton_EBrest_SignalYield.GetMinimum() < 0.0) l_3->Draw();
    l_sig1_del->Draw();
    l_sig2_del->Draw();
    leg_sig->Draw("same");
    c1->Print(EPSName.Data());
  }


  h_SigPhoton_EBrest_Fit_BrecoTM_SigGammaTM.SetStats(false);
  h_SigPhoton_EBrest_True_BrecoTM_SigGammaTM.SetStats(false);
  //Signal yield:
  c1->Clear();
  h_SigPhoton_EBrest_SignalYield.Draw();
  gr_SigPhoton_EBrest_SignalYield->Draw("P");
  gr_SigPhoton_EBrest_SignalYield_StatOnly->Draw("P");
  if(_IsMC > 1) {
    h_SigPhoton_EBrest_Fit_BrecoTM_SigGammaTM.Draw("same");
    h_SigPhoton_EBrest_True_BrecoTM_SigGammaTM.Draw("same");
  }
  if(h_SigPhoton_EBrest_SignalYield.GetMinimum() < 0.0) l_3->Draw();
  l_sig1_del->Draw();
  l_sig2_del->Draw();
  TLegend* leg_sig_v2 = new TLegend(0.75,0.75,0.9,0.9);
  leg_sig_v2->SetFillColor(10);
  if(_IsMC > 1) {
    leg_sig_v2->AddEntry(&h_SigPhoton_EBrest_Fit_BrecoTM_SigGammaTM,"MC","l");
    leg_sig_v2->AddEntry(&h_SigPhoton_EBrest_True_BrecoTM_SigGammaTM,"MC+true","l");
  }
  leg_sig_v2->AddEntry(&h_SigPhoton_EBrest_SignalYield_StatOnly,  "Fit stat",     "l");
  leg_sig_v2->AddEntry(&h_SigPhoton_EBrest_SignalYield,           "Fit stat+syst","l");
  leg_sig_v2->Draw("same");
  c1->Print(EPSName.Data());

  //Signal region Zoom:
  c1->Clear();
  Maximum = -1.0e+20;
  Minimum =  1.0e+20;
  for(int i=0;i<h_SigPhoton_EBrest_SignalYield.GetXaxis()->GetNbins();i++) {
    double v,e,c,w;

    c = h_SigPhoton_EBrest_SignalYield.GetBinCenter(i+1);
    w = h_SigPhoton_EBrest_SignalYield.GetBinWidth(i+1)*0.5;

    if(!(c-w + 1.0e-6 > R_EgBRF_Sig[0] && 
	 c+w - 1.0e-6 < R_EgBRF_Sig[1])) continue;

    v = h_SigPhoton_EBrest_SignalYield.GetBinContent(i+1);
    e = h_SigPhoton_EBrest_SignalYield.GetBinError(i+1);
    if(Maximum < v+e) Maximum = v+e;
    if(Minimum > v-e) Minimum = v-e;
  }
  porcent = 0.10;
  h_SigPhoton_EBrest_mESFit_REFSig.SetMaximum(Maximum + porcent*(Maximum-Minimum));
  h_SigPhoton_EBrest_mESFit_REFSig.SetMinimum(Minimum - porcent*(Maximum-Minimum));
  h_SigPhoton_EBrest_mESFit_REFSig.SetTitle(h_SigPhoton_EBrest_SignalYield.GetTitle());
  h_SigPhoton_EBrest_mESFit_REFSig.Draw();
  //h_SigPhoton_EBrest_SignalYield.Draw("same");
  //h_SigPhoton_EBrest_SignalYield_StatOnly.Draw("same");
  gr_SigPhoton_EBrest_SignalYield->Draw("P");
  gr_SigPhoton_EBrest_SignalYield_StatOnly->Draw("P");
  if(_IsMC > 1) h_SigPhoton_EBrest_Fit_BrecoTM_SigGammaTM.Draw("same");
  TLine* l_4 = new TLine(h_SigPhoton_EBrest_mESFit_REFSig.GetXaxis()->GetXmin(),0.0,
			 h_SigPhoton_EBrest_mESFit_REFSig.GetXaxis()->GetXmax(),0.0);
  l_4->SetLineColor(1);
  l_4->SetLineWidth(2);
  l_4->SetLineStyle(2);
  if(h_SigPhoton_EBrest_mESFit_REFSig.GetMinimum() < 0.0) l_4->Draw();
  if(_IsMC > 1) leg_sig->Draw("same");  
  c1->Print(EPSName.Data());


  //Signal region Zoom:
  c1->Clear();
  h_SigPhoton_EBrest_mESFit_REFSig.Draw();
  gr_SigPhoton_EBrest_SignalYield->Draw("P");
  gr_SigPhoton_EBrest_SignalYield_StatOnly->Draw("P");
  if(_IsMC > 1) {
    h_SigPhoton_EBrest_Fit_BrecoTM_SigGammaTM.Draw("same");
    h_SigPhoton_EBrest_True_BrecoTM_SigGammaTM.Draw("same");
  }
  if(h_SigPhoton_EBrest_mESFit_REFSig.GetMinimum() < 0.0) l_4->Draw();
  if(_IsMC > 1) leg_sig_v2->Draw("same");  
  c1->Print(EPSName.Data());


  //Low-E Bkg region Zoom:
  c1->Clear();
  Maximum = -1.0e+20;
  Minimum =  1.0e+20;
  for(int i=0;i<h_SigPhoton_EBrest_SignalYield.GetXaxis()->GetNbins();i++) {
    double v,e,c,w;

    c = h_SigPhoton_EBrest_SignalYield.GetBinCenter(i+1);
    w = h_SigPhoton_EBrest_SignalYield.GetBinWidth(i+1)*0.5;

    if(!(c-w + 1.0e-6 > R_EgBRF_Bkg[0] && 
	 c+w - 1.0e-6 < R_EgBRF_Bkg[1])) continue;

    v = h_SigPhoton_EBrest_SignalYield.GetBinContent(i+1);
    e = h_SigPhoton_EBrest_SignalYield.GetBinError(i+1);
    if(Maximum < v+e) Maximum = v+e;
    if(Minimum > v-e) Minimum = v-e;
  }
  porcent = 0.30;
  h_SigPhoton_EBrest_mESFit_REFBkg.SetMaximum(Maximum + porcent*(Maximum-Minimum));
  h_SigPhoton_EBrest_mESFit_REFBkg.SetMinimum(Minimum - porcent*(Maximum-Minimum));
  h_SigPhoton_EBrest_mESFit_REFBkg.SetTitle(h_SigPhoton_EBrest_SignalYield.GetTitle());
  h_SigPhoton_EBrest_mESFit_REFBkg.Draw();
  gr_SigPhoton_EBrest_SignalYield->Draw("P");
  gr_SigPhoton_EBrest_SignalYield_StatOnly->Draw("P");
  //h_SigPhoton_EBrest_SignalYield.Draw("same");
  //h_SigPhoton_EBrest_SignalYield_StatOnly.Draw("same");
  TLine* l_5 = new TLine(h_SigPhoton_EBrest_mESFit_REFBkg.GetXaxis()->GetXmin(),0.0,
			 h_SigPhoton_EBrest_mESFit_REFBkg.GetXaxis()->GetXmax(),0.0);
  l_5->SetLineColor(1);
  l_5->SetLineWidth(2);
  l_5->SetLineStyle(2);
  if(h_SigPhoton_EBrest_mESFit_REFBkg.GetMinimum() < 0.0) l_5->Draw();
  c1->Print(EPSName.Data());


  //High-E Bkg region Zoom:
  c1->Clear();
  Maximum = -1.0e+20;
  Minimum =  1.0e+20;
  for(int i=0;i<h_SigPhoton_EBrest_SignalYield.GetXaxis()->GetNbins();i++) {
    double v,e,c,w;

    c = h_SigPhoton_EBrest_SignalYield.GetBinCenter(i+1);
    w = h_SigPhoton_EBrest_SignalYield.GetBinWidth(i+1)*0.5;

    if(!(c-w + 1.0e-6 > R_EgBRF_Bkg2[0] && 
	 c+w - 1.0e-6 < R_EgBRF_Bkg2[1])) continue;

    v = h_SigPhoton_EBrest_SignalYield.GetBinContent(i+1);
    e = h_SigPhoton_EBrest_SignalYield.GetBinError(i+1);
    if(Maximum < v+e) Maximum = v+e;
    if(Minimum > v-e) Minimum = v-e;
  }
  porcent = 0.30;
  h_SigPhoton_EBrest_mESFit_REFBkg2.SetMaximum(Maximum + porcent*(Maximum-Minimum));
  h_SigPhoton_EBrest_mESFit_REFBkg2.SetMinimum(Minimum - porcent*(Maximum-Minimum));
  h_SigPhoton_EBrest_mESFit_REFBkg2.SetTitle(h_SigPhoton_EBrest_SignalYield.GetTitle());
  h_SigPhoton_EBrest_mESFit_REFBkg2.Draw();
  gr_SigPhoton_EBrest_SignalYield->Draw("P");
  gr_SigPhoton_EBrest_SignalYield_StatOnly->Draw("P");
  //h_SigPhoton_EBrest_SignalYield.Draw("same");
  //h_SigPhoton_EBrest_SignalYield_StatOnly.Draw("same");
  TLine* l_6 = new TLine(h_SigPhoton_EBrest_mESFit_REFBkg2.GetXaxis()->GetXmin(),0.0,
			 h_SigPhoton_EBrest_mESFit_REFBkg2.GetXaxis()->GetXmax(),0.0);
  l_6->SetLineColor(1);
  l_6->SetLineWidth(2);
  l_6->SetLineStyle(2);
  if(h_SigPhoton_EBrest_mESFit_REFBkg2.GetMinimum() < 0.0) l_6->Draw();
  c1->Print(EPSName.Data());

  c1->Clear();
  h_SignalYield_vs_EgBRF_Cov_Stat->Draw("colz");
  c1->Print(EPSName.Data());

  c1->Clear();
  h_SignalYield_vs_EgBRF_Cov_BBkgStat->Draw("colz");
  c1->Print(EPSName.Data());

  c1->Clear();
  h_SignalYield_vs_EgBRF_Cov_AllHEEfficCorr->Draw("colz");
  c1->Print(EPSName.Data());

  c1->Clear();
  h_SignalYield_vs_EgBRF_Cov_BrecoCorr->Draw("colz");
  c1->Print(EPSName.Data());

  c1->Clear();
  h_SignalYield_vs_EgBRF_Cov_Pi0EtaBkgMCYield->Draw("colz");
  c1->Print(EPSName.Data());

  c1->Clear();
  h_SignalYield_vs_EgBRF_Cov_Pi0EtaBkgRpi0Reta->Draw("colz");
  c1->Print(EPSName.Data());

  c1->Clear();
  h_SignalYield_vs_EgBRF_Cov_Pi0EtaBkgRpi0CRetaC->Draw("colz");
  c1->Print(EPSName.Data());

  c1->Clear();
  h_SignalYield_vs_EgBRF_Cov_Pi0EtaBkgmESFitYields->Draw("colz");
  c1->Print(EPSName.Data());

  c1->Clear();
  h_SignalYield_vs_EgBRF_Cov_mESFixX1->Draw("colz");
  c1->Print(EPSName.Data());

  c1->Clear();
  h_SignalYield_vs_EgBRF_Cov_mESFixX2->Draw("colz");
  c1->Print(EPSName.Data());

  c1->Clear();
  h_SignalYield_vs_EgBRF_Cov_mESFixX3->Draw("colz");
  c1->Print(EPSName.Data());

  c1->Clear();
  h_SignalYield_vs_EgBRF_Cov_mESFix_alphaP->Draw("colz");
  c1->Print(EPSName.Data());

  c1->Clear();
  h_SignalYield_vs_EgBRF_Cov_mESFix_BBComb->Draw("colz");
  c1->Print(EPSName.Data());

  c1->Clear();
  h_SignalYield_vs_EgBRF_Cov_mESFix_slope->Draw("colz");
  c1->Print(EPSName.Data());

  c1->Clear();
  h_SignalYield_vs_EgBRF_Cov_mESFitBias->Draw("colz");
  c1->Print(EPSName.Data());

  c1->Clear();
  h_SignalYield_vs_EgBRF_Cov_Pi0EtaBkgMggFitValSyst->Draw("colz");
  c1->Print(EPSName.Data());

  c1->Clear();
  h_SignalYield_vs_EgBRF_Cov_Pi0EtaBkgLE_Base->Draw("colz");
  c1->Print(EPSName.Data());

  c1->Clear();
  h_SignalYield_vs_EgBRF_Cov_Pi0EtaBkgULE->Draw("colz");
  c1->Print(EPSName.Data());

  c1->Clear();
  h_SignalYield_vs_EgBRF_Cov_OmegaBkgAlpha->Draw("colz");
  c1->Print(EPSName.Data());

  c1->Clear();
  h_SignalYield_vs_EgBRF_Cov_etaprimeBkgAlpha->Draw("colz");
  c1->Print(EPSName.Data());

  c1->Clear();
  h_SignalYield_vs_EgBRF_Cov_BkgSLAlpha->Draw("colz");
  c1->Print(EPSName.Data());

  c1->Clear();
  h_SignalYield_vs_EgBRF_Cov_BremsBkgMaterial->Draw("colz");
  c1->Print(EPSName.Data());

  c1->Clear();
  h_SignalYield_vs_EgBRF_Cov_ElecBkgTrackIneffic->Draw("colz");
  c1->Print(EPSName.Data());

  c1->Clear();
  h_SignalYield_vs_EgBRF_Cov_AntiNeutronsBkgAlpha->Draw("colz");
  c1->Print(EPSName.Data());

  c1->Clear();
  h_SignalYield_vs_EgBRF_Cov_BDT->Draw("colz");
  c1->Print(EPSName.Data());

  c1->Clear();
  h_SignalYield_vs_EgBRF_Cov_Pi0EtaVeto_MassLineShape->Draw("colz");
  c1->Print(EPSName.Data());

  c1->Clear();
  h_SignalYield_vs_EgBRF_Cov_Pi0EtaVeto_2ndGammaMult->Draw("colz");
  c1->Print(EPSName.Data());


  double NBkg_LowE_region[2];
  double NSig_LowE_region;
  double NBkg_HighE_region[2];
  double NSig_HighE_region;

  NBkg_LowE_region[0]  = 0.0;
  NBkg_LowE_region[1]  = 0.0;
  NBkg_HighE_region[0] = 0.0;
  NBkg_HighE_region[1] = 0.0;
  NSig_LowE_region     = 0.0;
  NSig_HighE_region    = 0.0;
  for(int i=0;i<h_SigPhoton_EBrest_SignalYield.GetXaxis()->GetNbins();i++) {
    double c = h_SigPhoton_EBrest_SignalYield.GetBinCenter(i+1);
    double w = h_SigPhoton_EBrest_SignalYield.GetBinWidth(i+1)*0.5;
    double Ri[2];
    Ri[0] = c-w + 1.0e-6;
    Ri[1] = c+w - 1.0e-6;
    double Ni = h_SigPhoton_EBrest_SignalYield.GetBinContent(i+1);
    double Ni_sig = h_SigPhoton_EBrest_Fit_BrecoTM_SigGammaTM.GetBinContent(i+1);

    //LowE region:
    if(Ri[0] > R_EgBRF_Bkg[0] && 
       Ri[1] < R_EgBRF_Bkg[1]) {
      NBkg_LowE_region[0] += Ni;
      NSig_LowE_region    += Ni_sig;
      for(int j=0;j<h_SigPhoton_EBrest_SignalYield.GetXaxis()->GetNbins();j++) {
	double cj = h_SigPhoton_EBrest_SignalYield.GetBinCenter(j+1);
	double wj = h_SigPhoton_EBrest_SignalYield.GetBinWidth(j+1)*0.5;
	double Rj[2];
	Rj[0] = cj-wj + 1.0e-6;
	Rj[1] = cj+wj - 1.0e-6;
	if(Rj[0] > R_EgBRF_Bkg[0] && 
	   Rj[1] < R_EgBRF_Bkg[1]) {
	  double cov = 0.0;
	  cov += h_SignalYield_vs_EgBRF_Cov_Stat->GetBinContent(i+1,j+1);
	  cov += h_SignalYield_vs_EgBRF_Cov_BBkgStat->GetBinContent(i+1,j+1);
	  cov += h_SignalYield_vs_EgBRF_Cov_AllHEEfficCorr->GetBinContent(i+1,j+1);
	  cov += h_SignalYield_vs_EgBRF_Cov_BrecoCorr->GetBinContent(i+1,j+1);
	  cov += h_SignalYield_vs_EgBRF_Cov_Pi0EtaBkgMCYield->GetBinContent(i+1,j+1);
	  cov += h_SignalYield_vs_EgBRF_Cov_Pi0EtaBkgRpi0Reta->GetBinContent(i+1,j+1);
	  cov += h_SignalYield_vs_EgBRF_Cov_Pi0EtaBkgRpi0CRetaC->GetBinContent(i+1,j+1);
	  cov += h_SignalYield_vs_EgBRF_Cov_Pi0EtaBkgmESFitYields->GetBinContent(i+1,j+1);
	  cov += h_SignalYield_vs_EgBRF_Cov_mESFixX1->GetBinContent(i+1,j+1);
	  cov += h_SignalYield_vs_EgBRF_Cov_mESFixX2->GetBinContent(i+1,j+1);
	  cov += h_SignalYield_vs_EgBRF_Cov_mESFixX3->GetBinContent(i+1,j+1);
	  cov += h_SignalYield_vs_EgBRF_Cov_mESFix_alphaP->GetBinContent(i+1,j+1);
	  cov += h_SignalYield_vs_EgBRF_Cov_mESFix_BBComb->GetBinContent(i+1,j+1);
	  cov += h_SignalYield_vs_EgBRF_Cov_mESFix_slope->GetBinContent(i+1,j+1);
	  cov += h_SignalYield_vs_EgBRF_Cov_mESFitBias->GetBinContent(i+1,j+1);
	  cov += h_SignalYield_vs_EgBRF_Cov_Pi0EtaBkgMggFitValSyst->GetBinContent(i+1,j+1);
	  cov += h_SignalYield_vs_EgBRF_Cov_Pi0EtaBkgLE_Base->GetBinContent(i+1,j+1);
	  cov += h_SignalYield_vs_EgBRF_Cov_Pi0EtaBkgULE->GetBinContent(i+1,j+1);
	  cov += h_SignalYield_vs_EgBRF_Cov_OmegaBkgAlpha->GetBinContent(i+1,j+1);
	  cov += h_SignalYield_vs_EgBRF_Cov_etaprimeBkgAlpha->GetBinContent(i+1,j+1);
	  cov += h_SignalYield_vs_EgBRF_Cov_BkgSLAlpha->GetBinContent(i+1,j+1);
	  cov += h_SignalYield_vs_EgBRF_Cov_BremsBkgMaterial->GetBinContent(i+1,j+1);
	  cov += h_SignalYield_vs_EgBRF_Cov_ElecBkgTrackIneffic->GetBinContent(i+1,j+1);
	  cov += h_SignalYield_vs_EgBRF_Cov_AntiNeutronsBkgAlpha->GetBinContent(i+1,j+1);
	  cov += h_SignalYield_vs_EgBRF_Cov_BDT->GetBinContent(i+1,j+1);
	  cov += h_SignalYield_vs_EgBRF_Cov_Pi0EtaVeto_MassLineShape->GetBinContent(i+1,j+1);
	  cov += h_SignalYield_vs_EgBRF_Cov_Pi0EtaVeto_2ndGammaMult->GetBinContent(i+1,j+1);
	  NBkg_LowE_region[1] += cov;
	}
      }
    }

    //HighE region:
    if(Ri[0] > R_EgBRF_Bkg2[0] && 
       Ri[1] < R_EgBRF_Bkg2[1]) {
      NBkg_HighE_region[0] += Ni;
      NSig_HighE_region    += Ni_sig;
      for(int j=0;j<h_SigPhoton_EBrest_SignalYield.GetXaxis()->GetNbins();j++) {
	double cj = h_SigPhoton_EBrest_SignalYield.GetBinCenter(j+1);
	double wj = h_SigPhoton_EBrest_SignalYield.GetBinWidth(j+1)*0.5;
	double Rj[2];
	Rj[0] = cj-wj + 1.0e-6;
	Rj[1] = cj+wj - 1.0e-6;
	if(Rj[0] > R_EgBRF_Bkg2[0] && 
	   Rj[1] < R_EgBRF_Bkg2[1]) {
	  double cov = 0.0;
	  cov += h_SignalYield_vs_EgBRF_Cov_Stat->GetBinContent(i+1,j+1);
	  cov += h_SignalYield_vs_EgBRF_Cov_BBkgStat->GetBinContent(i+1,j+1);
	  cov += h_SignalYield_vs_EgBRF_Cov_AllHEEfficCorr->GetBinContent(i+1,j+1);
	  cov += h_SignalYield_vs_EgBRF_Cov_BrecoCorr->GetBinContent(i+1,j+1);
	  cov += h_SignalYield_vs_EgBRF_Cov_Pi0EtaBkgMCYield->GetBinContent(i+1,j+1);
	  cov += h_SignalYield_vs_EgBRF_Cov_Pi0EtaBkgRpi0Reta->GetBinContent(i+1,j+1);
	  cov += h_SignalYield_vs_EgBRF_Cov_Pi0EtaBkgRpi0CRetaC->GetBinContent(i+1,j+1);
	  cov += h_SignalYield_vs_EgBRF_Cov_Pi0EtaBkgmESFitYields->GetBinContent(i+1,j+1);
	  cov += h_SignalYield_vs_EgBRF_Cov_mESFixX1->GetBinContent(i+1,j+1);
	  cov += h_SignalYield_vs_EgBRF_Cov_mESFixX2->GetBinContent(i+1,j+1);
	  cov += h_SignalYield_vs_EgBRF_Cov_mESFixX3->GetBinContent(i+1,j+1);
	  cov += h_SignalYield_vs_EgBRF_Cov_mESFix_alphaP->GetBinContent(i+1,j+1);
	  cov += h_SignalYield_vs_EgBRF_Cov_mESFix_BBComb->GetBinContent(i+1,j+1);
	  cov += h_SignalYield_vs_EgBRF_Cov_mESFix_slope->GetBinContent(i+1,j+1);
	  cov += h_SignalYield_vs_EgBRF_Cov_mESFitBias->GetBinContent(i+1,j+1);
	  cov += h_SignalYield_vs_EgBRF_Cov_Pi0EtaBkgMggFitValSyst->GetBinContent(i+1,j+1);
	  cov += h_SignalYield_vs_EgBRF_Cov_Pi0EtaBkgLE_Base->GetBinContent(i+1,j+1);
	  cov += h_SignalYield_vs_EgBRF_Cov_Pi0EtaBkgULE->GetBinContent(i+1,j+1);
	  cov += h_SignalYield_vs_EgBRF_Cov_OmegaBkgAlpha->GetBinContent(i+1,j+1);
	  cov += h_SignalYield_vs_EgBRF_Cov_etaprimeBkgAlpha->GetBinContent(i+1,j+1);
	  cov += h_SignalYield_vs_EgBRF_Cov_BkgSLAlpha->GetBinContent(i+1,j+1);
	  cov += h_SignalYield_vs_EgBRF_Cov_BremsBkgMaterial->GetBinContent(i+1,j+1);
	  cov += h_SignalYield_vs_EgBRF_Cov_ElecBkgTrackIneffic->GetBinContent(i+1,j+1);
	  cov += h_SignalYield_vs_EgBRF_Cov_AntiNeutronsBkgAlpha->GetBinContent(i+1,j+1);
	  cov += h_SignalYield_vs_EgBRF_Cov_BDT->GetBinContent(i+1,j+1);
	  cov += h_SignalYield_vs_EgBRF_Cov_Pi0EtaVeto_MassLineShape->GetBinContent(i+1,j+1);
	  cov += h_SignalYield_vs_EgBRF_Cov_Pi0EtaVeto_2ndGammaMult->GetBinContent(i+1,j+1);
	  NBkg_HighE_region[1] += cov;
	}
      }
    }    

  }
  NBkg_LowE_region[1]  = sqrt(NBkg_LowE_region[1]);
  NBkg_HighE_region[1] = sqrt(NBkg_HighE_region[1]);

  if(_IsMC > 1) {
    cout << endl;
    cout << endl;
    cout << "N events in (" << R_EgBRF_Bkg[0] << "," << R_EgBRF_Bkg[1] << ") GeV region = " 
	 << NBkg_LowE_region[0] << " +/- " << NBkg_LowE_region[1] 
	 << ", Nsignal = " << NSig_LowE_region
	 << "(Nmeas - Nexp)/sigma = " << (NBkg_LowE_region[0]-NSig_LowE_region)/NBkg_LowE_region[1] 
	 << ", prob = " << TMath::Prob(pow((NBkg_LowE_region[0]-NSig_LowE_region)/NBkg_LowE_region[1],2),2)
	 << endl;
    cout << "N events in (" << R_EgBRF_Bkg2[0] << "," << R_EgBRF_Bkg2[1] << ") GeV region = " 
	 << NBkg_HighE_region[0] << " +/- " << NBkg_HighE_region[1] 
	 << ", Nsignal = " << NSig_HighE_region
	 << "(Nmeas - Nexp)/sigma = " << (NBkg_HighE_region[0]-NSig_HighE_region)/NBkg_HighE_region[1] 
	 << ", prob = " << TMath::Prob(pow((NBkg_HighE_region[0]-NSig_HighE_region)/NBkg_HighE_region[1],2),2)
	 << endl;
    cout << endl;
    cout << endl;
  }
  else {
    cout << endl;
    cout << endl;
    cout << "N events in (" << R_EgBRF_Bkg[0] << "," << R_EgBRF_Bkg[1] << ") GeV region = " 
	 << NBkg_LowE_region[0] << " +/- " << NBkg_LowE_region[1] 
	 << endl;
    cout << "N events in (" << R_EgBRF_Bkg2[0] << "," << R_EgBRF_Bkg2[1] << ") GeV region = " 
	 << NBkg_HighE_region[0] << " +/- " << NBkg_HighE_region[1] 
	 << endl;
    cout << endl;
    cout << endl;
  }


  double X,Y;
  //Low-E Bkg region Zoom:
  c1->Clear();
  h_SigPhoton_EBrest_mESFit_REFBkg.Draw();
  gr_SigPhoton_EBrest_SignalYield->Draw("P");
  gr_SigPhoton_EBrest_SignalYield_StatOnly->Draw("P");
  //h_SigPhoton_EBrest_SignalYield.Draw("same");
  //h_SigPhoton_EBrest_SignalYield_StatOnly.Draw("same");
  if(h_SigPhoton_EBrest_mESFit_REFBkg.GetMinimum() < 0.0) l_5->Draw();
  porcent = 0.05;
  X  = h_SigPhoton_EBrest_mESFit_REFBkg.GetXaxis()->GetXmin();
  X += porcent*h_SigPhoton_EBrest_mESFit_REFBkg.GetXaxis()->GetXmax();
  X -= porcent*h_SigPhoton_EBrest_mESFit_REFBkg.GetXaxis()->GetXmin();
  porcent = 0.08;
  Y  = h_SigPhoton_EBrest_mESFit_REFBkg.GetMaximum();
  Y -= porcent*h_SigPhoton_EBrest_mESFit_REFBkg.GetMaximum();
  Y += porcent*h_SigPhoton_EBrest_mESFit_REFBkg.GetMinimum();
  sprintf(ytitle,"%.2f",NBkg_LowE_region[0]);
  HistName = TString("N^{tot} = ") + ytitle + TString(" #pm ");
  sprintf(ytitle,"%.2f",NBkg_LowE_region[1]);
  if(_IsMC > 1) {
    HistName += ytitle + TString(", N^{exp} = ");
    sprintf(ytitle,"%.2f",NSig_LowE_region);
    HistName += ytitle;
  }
  else HistName += ytitle;
  latex->DrawLatex(X,Y,HistName.Data());

  if(_IsMC > 1) {
    porcent = 0.15;
    Y  = h_SigPhoton_EBrest_mESFit_REFBkg.GetMaximum();
    Y -= porcent*h_SigPhoton_EBrest_mESFit_REFBkg.GetMaximum();
    Y += porcent*h_SigPhoton_EBrest_mESFit_REFBkg.GetMinimum();
    HistName  = TString("(N^{meas} - N^{exp})/#sigma = ");
    sprintf(ytitle,"%.3f",(NBkg_LowE_region[0]-NSig_LowE_region)/NBkg_LowE_region[1]);
    HistName += ytitle + TString(", prob = ");
    sprintf(ytitle,"%.3f",TMath::Prob(pow((NBkg_LowE_region[0]-NSig_LowE_region)/NBkg_LowE_region[1],2),2));
    HistName += ytitle;
    latex->DrawLatex(X,Y,HistName.Data());
  }
  c1->Print(EPSName.Data());

  //High-E Bkg region Zoom:
  c1->Clear();
  h_SigPhoton_EBrest_mESFit_REFBkg2.Draw();
  gr_SigPhoton_EBrest_SignalYield->Draw("P");
  gr_SigPhoton_EBrest_SignalYield_StatOnly->Draw("P");
  //h_SigPhoton_EBrest_SignalYield.Draw("same");
  //h_SigPhoton_EBrest_SignalYield_StatOnly.Draw("same");
  if(h_SigPhoton_EBrest_mESFit_REFBkg2.GetMinimum() < 0.0) l_6->Draw();
  porcent = 0.05;
  X  = h_SigPhoton_EBrest_mESFit_REFBkg2.GetXaxis()->GetXmin();
  X += porcent*h_SigPhoton_EBrest_mESFit_REFBkg2.GetXaxis()->GetXmax();
  X -= porcent*h_SigPhoton_EBrest_mESFit_REFBkg2.GetXaxis()->GetXmin();
  porcent = 0.08;
  Y  = h_SigPhoton_EBrest_mESFit_REFBkg2.GetMaximum();
  Y -= porcent*h_SigPhoton_EBrest_mESFit_REFBkg2.GetMaximum();
  Y += porcent*h_SigPhoton_EBrest_mESFit_REFBkg2.GetMinimum();
  sprintf(ytitle,"%.2f",NBkg_HighE_region[0]);
  HistName = TString("N^{tot} = ") + ytitle + TString(" #pm ");
  sprintf(ytitle,"%.2f",NBkg_HighE_region[1]);
  if(_IsMC > 1) {
    HistName += ytitle + TString(", N^{exp} = ");
    sprintf(ytitle,"%.2f",NSig_HighE_region);
    HistName += ytitle;
  }
  else HistName += ytitle;
  latex->DrawLatex(X,Y,HistName.Data());

  if(_IsMC > 1) {
    porcent = 0.15;
    Y  = h_SigPhoton_EBrest_mESFit_REFBkg2.GetMaximum();
    Y -= porcent*h_SigPhoton_EBrest_mESFit_REFBkg2.GetMaximum();
    Y += porcent*h_SigPhoton_EBrest_mESFit_REFBkg2.GetMinimum();
    HistName  = TString("(N^{meas} - N^{exp})/#sigma = ");
    sprintf(ytitle,"%.3f",(NBkg_HighE_region[0]-NSig_HighE_region)/NBkg_HighE_region[1]);
    HistName += ytitle + TString(", prob = ");
    sprintf(ytitle,"%.3f",TMath::Prob(pow((NBkg_HighE_region[0]-NSig_HighE_region)/NBkg_HighE_region[1],2),2));
    HistName += ytitle;
    latex->DrawLatex(X,Y,HistName.Data());
  }
  c1->Print(EPSName.Data());

  c1->Print(EPSNameC.Data());

  command = TString("rm ") + Name_mES;
  gSystem->Exec(command.Data());

  TString ROOTFile = TString(outputFile) + TString(".root");
  TFile f_out(ROOTFile.Data(),"RECREATE");
  h_SigPhoton_EBrest_SignalYield.Write();
  h_SigPhoton_EBrest_SignalYield_StatOnly.Write();

  gr_SigPhoton_EBrest_SignalYield->Write();
  gr_SigPhoton_EBrest_SignalYield_StatOnly->Write();

  if(_IsMC > 1) {
    h_SigPhoton_EBrest_Fit_BrecoTM_SigGammaTM.Write();
    h_SigPhoton_EBrest_True_BrecoTM_SigGammaTM.Write();
  }

  h_SignalYield_vs_EgBRF_Cov_Stat->Write();
  h_SignalYield_vs_EgBRF_Cov_BBkgStat->Write();
  h_SignalYield_vs_EgBRF_Cov_AllHEEfficCorr->Write();
  h_SignalYield_vs_EgBRF_Cov_BrecoCorr->Write();
  h_SignalYield_vs_EgBRF_Cov_Pi0EtaBkgMCYield->Write();
  h_SignalYield_vs_EgBRF_Cov_Pi0EtaBkgRpi0Reta->Write();
  h_SignalYield_vs_EgBRF_Cov_Pi0EtaBkgRpi0CRetaC->Write();
  h_SignalYield_vs_EgBRF_Cov_Pi0EtaBkgmESFitYields->Write();
  h_SignalYield_vs_EgBRF_Cov_mESFixX1->Write();
  h_SignalYield_vs_EgBRF_Cov_mESFixX2->Write();
  h_SignalYield_vs_EgBRF_Cov_mESFixX3->Write();
  h_SignalYield_vs_EgBRF_Cov_mESFix_alphaP->Write();
  h_SignalYield_vs_EgBRF_Cov_mESFix_BBComb->Write();
  h_SignalYield_vs_EgBRF_Cov_mESFix_slope->Write();
  h_SignalYield_vs_EgBRF_Cov_mESFitBias->Write();
  h_SignalYield_vs_EgBRF_Cov_Pi0EtaBkgMggFitValSyst->Write();
  h_SignalYield_vs_EgBRF_Cov_Pi0EtaBkgLE_Base->Write();
  h_SignalYield_vs_EgBRF_Cov_Pi0EtaBkgULE->Write();
  h_SignalYield_vs_EgBRF_Cov_OmegaBkgAlpha->Write();
  h_SignalYield_vs_EgBRF_Cov_etaprimeBkgAlpha->Write();
  h_SignalYield_vs_EgBRF_Cov_BkgSLAlpha->Write();
  h_SignalYield_vs_EgBRF_Cov_BremsBkgMaterial->Write();
  h_SignalYield_vs_EgBRF_Cov_ElecBkgTrackIneffic->Write();
  h_SignalYield_vs_EgBRF_Cov_AntiNeutronsBkgAlpha->Write();
  h_SignalYield_vs_EgBRF_Cov_BDT->Write();
  h_SignalYield_vs_EgBRF_Cov_Pi0EtaVeto_MassLineShape->Write();
  h_SignalYield_vs_EgBRF_Cov_Pi0EtaVeto_2ndGammaMult->Write();
  f_out.Close();

  return;

}
//============================================================================

