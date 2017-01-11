#include "Riostream.h"
#include "TF1.h"
#include "TFile.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TGraphErrors.h"
#include "TGraphAsymmErrors.h"
#include "TGraphBentErrors.h"
#include "TCanvas.h"
#include "TLegend.h"

void MomFromDedx( Char_t *inputFileName){

  // For each bin of the X axis of a 2D histogram,
  //  build the distribution, for the Y axis variable,
  //  of most probable value.
  // This distribution is the native resolution 
  //  on variable X as measured by variable Y.
  //  Parametrise the curve: dEdx per layer MPV = f(momentum)

  //
  //  .L MomFromDedx.C++
  //  MomFromDedx("histos.root")
  //
  // IRB 2014/03/10

  // Control variables
  Int_t minimalEntriesForFit = 20;
  Int_t MPVmethod = 1; // 0=highest bin, 1=Landau fit, 2=Gaussian fit
  
  // Note: function ranges will be defined later, from histogram axis
  TF1 *fXvsYModel = new TF1( "fXvsYModel", "[1]+[0]/sqrt(x)",.1, 1.);

  
  Int_t canvasWidth=800;

  // Open file and load histogram
  TFile *inputFile = new TFile( inputFileName, "UPDATE");
  if( inputFile->IsOpen()==0 ) {
    cout << "File name is wrong !" << " --> EXIT" << endl;
    return;
  }
  TH2F *h2 = (TH2F*)inputFile->Get("dEdxVsMom");

  if( h2==NULL ) {
    cout << "Histo name is wrong or does not exist !" << " --> EXIT" << endl;
    return;
  }
  cout << "Got histo with name " << h2->GetName() << endl;
  
// fXvsYModel->SetRange( h2->GetYaxis()->GetXmin(), h2->GetYaxis()->GetXmax());
//  fXvsYModel->SetRange( 70000., h2->GetYaxis()->GetXmax());
  fXvsYModel->SetRange( 600000., h2->GetYaxis()->GetXmax());
 
  
  TH1F *hchi2 = new TH1F("hchi2","chi2 distribution of Edeposit distribution fit",100,0.,100.);
  
  // --------------------------------------------------------------
  // Display original histo
  TCanvas *co = new TCanvas("co","Original 2D histo", 0, 0, canvasWidth, canvasWidth*3/4);
  co->SetGrid(1,1); 
 
  cout<< "Number of bins along Y (Edeposit) = "<<h2->GetNbinsY()<<endl;

  h2->DrawClone("colz");
  co->Update();


  // Prepare variables
  TCanvas *cp = new TCanvas("cp","Projection", 5, 5, canvasWidth, canvasWidth*3/4);  
  cp->cd();
  cp->SetGrid(1,1); 
  TH1F *h1Yprojection;
  TH1F *h1YprojectionToKeep;  
  Int_t nBinsX = h2->GetNbinsX();
  Int_t nBinsY = h2->GetNbinsY();
  
  cout << "Original histogram " << h2->GetTitle() << " has " << nBinsX << " bins in X.(momentum)" << endl << endl;

  
  Int_t nLimits = 3;
  Double_t probabilities[6] = { 0.00135, 0.02275, 0.15865, 0.84135, 0.97725, 0.99865 };
  if( nLimits>3 ) nLimits=3; // maxim nb of intervals is 3
  Double_t *quantiles = new Double_t[2*nLimits];
  Double_t *confidences = new Double_t[nLimits];
  Char_t **labels = new Char_t*[nLimits];
  Char_t name[20], title[150];

  Double_t *xValues = new Double_t[nBinsX];
  Double_t *yValues = new Double_t[nBinsX];
  Double_t **xLimits = new Double_t*[2*nLimits];
  Double_t **yLimits = new Double_t*[2*nLimits];
  for( Int_t iLimit = 0; iLimit<2*nLimits; iLimit++) {
    xLimits[iLimit] = new Double_t[nBinsX];
    yLimits[iLimit] = new Double_t[nBinsX];
  }
  for( Int_t iLimit = 0; iLimit<nLimits; iLimit++) {
    confidences[iLimit] = -probabilities[iLimit]+probabilities[2*nLimits-1-iLimit];
    labels[iLimit] = new Char_t[100];
    sprintf( labels[iLimit], "%.2f %%", confidences[iLimit]*100);
  }  
  
  TGraphAsymmErrors *gXvsY;
  
  // --------------------------------------------------------------
  // For each slice in momentum, fit the MPV of the Edeposit Landau distribution:

	Double_t chi2fit = -1.;
//  for( Int_t iBinX = 1; iBinX < nBinsX; iBinX++ ) { // Loop on X bins
  for( Int_t iBinX = 1; iBinX < nBinsX+1; iBinX++ ) { // Loop on X bins
    cout << " -- Considering momentum bin " << iBinX << endl;
    sprintf( name, "pX%d", iBinX);
	h1Yprojection = (TH1F*)h2->ProjectionY( name, iBinX, iBinX, "E");
//	h1Yprojection->Rebin(4);  (attention a sumw2)
	xValues[iBinX] = h2->GetXaxis()->GetBinCenter( iBinX);
    cout << "    x = " << xValues[iBinX] << endl;
    if( h1Yprojection->GetEntries() >= minimalEntriesForFit ) {
// Display Y projection distribution for Momentum (X) bin ~ in the middle: 		
      if( 47 < iBinX && iBinX < 49 ) {
        cout << "Keeping projection on Edeposit " << iBinX << " for display." << endl;
		  h1YprojectionToKeep = (TH1F*)h2->ProjectionY( name, iBinX, iBinX, "E");        
		  sprintf( title, "Projection on Y for X=%.2f (bin %d) - for display", xValues[iBinX], iBinX);
        h1YprojectionToKeep->SetTitle(title);
        h1YprojectionToKeep->SetYTitle(h2->GetYaxis()->GetTitle());  
      }
      h1Yprojection->GetQuantiles( 2*nLimits, quantiles, probabilities);
      // Define the most probable value according to given method
		chi2fit = 0.;
      switch (MPVmethod) {
        case 1:
          h1Yprojection->Fit("landau","Q");
          yValues[iBinX] = h1Yprojection->GetFunction( "landau" )->GetParameter(1);
		  chi2fit = h1Yprojection->GetFunction( "landau" )->GetChisquare();
          break;
        case 2:
          h1Yprojection->Fit("gaus","QL");
          yValues[iBinX] = h1Yprojection->GetFunction( "gaus" )->GetParameter(1);
		  chi2fit = h1Yprojection->GetFunction( "gaus" )->GetChisquare();			  
          break;
        case 0:
          yValues[iBinX] = h1Yprojection->GetBinCenter( h1Yprojection->GetMaximumBin() );
        default:
          break;
      }
      cout << " number of entries = " << h1Yprojection->GetEntries() << endl;
      cout << "    y = " << yValues[iBinX] << endl;
      cout << "   chi2 of the fit = " << chi2fit << endl;
      cout << endl;
      hchi2->Fill(chi2fit);
      
      for( Int_t iLimit = 0; iLimit<nLimits; iLimit++) {
        xLimits[2*iLimit][iBinX] = h2->GetXaxis()->GetBinWidth(iBinX)/2.;
        xLimits[2*iLimit+1][iBinX] = h2->GetXaxis()->GetBinWidth(iBinX)/2.;
        yLimits[2*iLimit][iBinX] = yValues[iBinX] - quantiles[iLimit];
        yLimits[2*iLimit+1][iBinX] = quantiles[2*nLimits-1-iLimit] - yValues[iBinX];
      }
    } //end if projection integral not zero  (or > given value)
    else {
      yValues[iBinX] = 0.;
      for( Int_t iLimit = 0; iLimit<nLimits; iLimit++) {
        xLimits[2*iLimit][iBinX] = 0.;
        xLimits[2*iLimit+1][iBinX] = 0.;
        yLimits[2*iLimit][iBinX] = 0.;
        yLimits[2*iLimit+1][iBinX] = 0.;
      }      
      cout << "    y = " << yValues[iBinX] << ", null limits (only " << h1Yprojection->GetEntries() << " entries)." << endl;

    }
    
  } // end loop on X bins
  cp->cd();
  h1YprojectionToKeep->Draw("");
  cp->Update();
  //return;

  TCanvas *cc = new TCanvas("cc","Chi2 from Landau fit on Edeposit", 0, 0, canvasWidth, canvasWidth*3/4);
  cc->cd();
  hchi2->Draw("");
  cc->Update();

  // --------------------------------------------------------------

  
  cout << endl << "Drawing Momentum (X) versus dE/dx (Y)" << endl;
  
  TCanvas *cv = new TCanvas("cv","Most probable Momentum (X) for a given dE/dx (Y)", 10, 10, canvasWidth, canvasWidth*3/4);  
  cv->cd();
  cv->SetGrid(1,1); 
  
  TH1F *hframeXvsY = new TH1F("hframeXvsY","Most probable Momentum (X) for a given dE/dx (Y)" , nBinsY, h2->GetYaxis()->GetXmin(), h2->GetYaxis()->GetXmax());
  hframeXvsY->SetXTitle(h2->GetYaxis()->GetTitle());  
  hframeXvsY->SetYTitle(h2->GetXaxis()->GetTitle());
  hframeXvsY->SetMinimum(h2->GetXaxis()->GetXmin());
  hframeXvsY->SetMaximum(h2->GetXaxis()->GetXmax());
  hframeXvsY->SetStats(0);
  
  hframeXvsY->Draw();  
// yLimits[4] et [5] -> quantiles of 1sigma on each side
  gXvsY = new TGraphAsymmErrors( nBinsX, yValues, xValues, yLimits[4], yLimits[5], xLimits[0], xLimits[1] );
  gXvsY->SetName( "gXvsY");
  gXvsY->SetTitle( "Most probale X versus Y");
  gXvsY->SetMarkerStyle(20);
  gXvsY->SetMarkerSize(1.);
  
  gXvsY->Draw("P");
  fXvsYModel->FixParameter(1,0.);
  gXvsY->Fit("fXvsYModel","R");
  Float_t param0=0.,param1=0.;
  param0 = gXvsY->GetFunction("fXvsYModel")->GetParameter(0);
  param1 = gXvsY->GetFunction("fXvsYModel")->GetParameter(1);
  cout<<endl<<endl<<endl;
  cout<<"--------------- parametrisation = "<< param1 <<" + "<< param0 << "/ sqr(Edeposit)" <<endl;
  
  cv->Update();
    
        

}
