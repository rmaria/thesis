void Pdftest(void)
{

  double TheSize  = 0.06;
  char ytitle[100];

  double m;
  double m0 = 134.0;
  double sL = 6.0;
  double sR = 8.0;
  double aL = 1.5;
  double aR = 1.8;
  double oL = 5;
  double oR = 5;

  double Rm[2];
  Rm[0] =  80.0;
  Rm[1] = 250.0;

  int Nbins = 1000;
  TH1D* h_doubleCBShape = new TH1D("h_doubleCBShape","The doubleCBShape",
			 Nbins,
			 Rm[0],Rm[1]);
  h_doubleCBShape->SetXTitle("m (MeV/c^{2})");
  h_doubleCBShape->GetXaxis()->CenterTitle(true);
  sprintf(ytitle,"Events / (%.3f MeV/c^{2})",h_doubleCBShape->GetBinWidth(0));
  h_doubleCBShape->SetYTitle(ytitle);
  h_doubleCBShape->GetYaxis()->CenterTitle(true);
  h_doubleCBShape->SetLineColor(2);
  h_doubleCBShape->SetLineWidth(2);
  h_doubleCBShape->SetMinimum(0.0);
  h_doubleCBShape->GetXaxis()->SetTitleSize(TheSize);
  h_doubleCBShape->GetXaxis()->SetLabelSize(TheSize);
  h_doubleCBShape->GetYaxis()->SetTitleSize(TheSize);
  h_doubleCBShape->GetYaxis()->SetLabelSize(TheSize);
  h_doubleCBShape->SetStats(false);

  TH1D* h_bifurcatedCBShape = new TH1D("h_bifurcatedCBShape","The bifurcatedCBShape",
				       Nbins,
				       Rm[0],Rm[1]);
  h_bifurcatedCBShape->SetXTitle("m (MeV/c^{2})");
  h_bifurcatedCBShape->GetXaxis()->CenterTitle(true);
  sprintf(ytitle,"Events / (%.3f MeV/c^{2})",h_bifurcatedCBShape->GetBinWidth(0));
  h_bifurcatedCBShape->SetYTitle(ytitle);
  h_bifurcatedCBShape->GetYaxis()->CenterTitle(true);
  h_bifurcatedCBShape->SetLineColor(kGreen+2);
  h_bifurcatedCBShape->SetLineWidth(2);
  h_bifurcatedCBShape->SetMinimum(0.0);
  h_bifurcatedCBShape->GetXaxis()->SetTitleSize(TheSize);
  h_bifurcatedCBShape->GetXaxis()->SetLabelSize(TheSize);
  h_bifurcatedCBShape->GetYaxis()->SetTitleSize(TheSize);
  h_bifurcatedCBShape->GetYaxis()->SetLabelSize(TheSize);
  h_bifurcatedCBShape->SetStats(false);

  TH1D* h_gauss = new TH1D("h_gauss","The gauss",
			 Nbins,
			 Rm[0],Rm[1]);
  h_gauss->SetXTitle("m (MeV/c^{2})");
  h_gauss->GetXaxis()->CenterTitle(true);
  sprintf(ytitle,"Events / (%.3f MeV/c^{2})",h_gauss->GetBinWidth(0));
  h_gauss->SetYTitle(ytitle);
  h_gauss->GetYaxis()->CenterTitle(true);
  h_gauss->SetLineColor(4);
  h_gauss->SetLineWidth(2);
  h_gauss->SetMinimum(0.0);
  h_gauss->GetXaxis()->SetTitleSize(TheSize);
  h_gauss->GetXaxis()->SetLabelSize(TheSize);
  h_gauss->GetYaxis()->SetTitleSize(TheSize);
  h_gauss->GetYaxis()->SetLabelSize(TheSize);
  h_gauss->SetStats(false);

  TH1D* h_diff = new TH1D("h_diff","The diff",
			 Nbins,
			 Rm[0],Rm[1]);
  h_diff->SetXTitle("m (MeV/c^{2})");
  h_diff->GetXaxis()->CenterTitle(true);
  sprintf(ytitle,"Events / (%.3f MeV/c^{2})",h_diff->GetBinWidth(0));
  h_diff->SetYTitle(ytitle);
  h_diff->GetYaxis()->CenterTitle(true);
  h_diff->SetLineColor(2);
  h_diff->SetLineWidth(2);
  h_diff->GetXaxis()->SetTitleSize(TheSize);
  h_diff->GetXaxis()->SetLabelSize(TheSize);
  h_diff->GetYaxis()->SetTitleSize(TheSize);
  h_diff->GetYaxis()->SetLabelSize(TheSize);
  h_diff->SetStats(false);

  TH1D* h_diff2 = new TH1D("h_diff2","The diff2",
			   Nbins,
			   Rm[0],Rm[1]);
  h_diff2->SetXTitle("m (MeV/c^{2})");
  h_diff2->GetXaxis()->CenterTitle(true);
  sprintf(ytitle,"Events / (%.3f MeV/c^{2})",h_diff2->GetBinWidth(0));
  h_diff2->SetYTitle(ytitle);
  h_diff2->GetYaxis()->CenterTitle(true);
  h_diff2->SetLineColor(4);
  h_diff2->SetLineWidth(2);
  h_diff2->GetXaxis()->SetTitleSize(TheSize);
  h_diff2->GetXaxis()->SetLabelSize(TheSize);
  h_diff2->GetYaxis()->SetTitleSize(TheSize);
  h_diff2->GetYaxis()->SetLabelSize(TheSize);
  h_diff2->SetStats(false);

  double min = 1.0e+20;
  for(int i=0;i<Nbins;i++) {
    m = Rm[0] + (i+0.5)*(Rm[1]-Rm[0])/Nbins;
    double val = doubleCBShape(m,m0,sL,sR,aL,aR,oL,oR);
    double diff = val;
    h_doubleCBShape->SetBinContent(i+1,val);
    val = gaussian(m,m0,sL);
    diff -= val;
    h_gauss->SetBinContent(i+1,val);
    h_diff->SetBinContent(i+1,diff);

    val = BifurcatedCBShape(m,m0,sL,sR,aL,oL);
    double diff2 = val;
    h_bifurcatedCBShape->SetBinContent(i+1,val);
    diff2 -= gaussian(m,m0,sL);
    h_diff2->SetBinContent(i+1,diff2);
  }

  cout << endl;
  cout << "RooDoubleCBShape: Numeric integral from hist  = " << h_doubleCBShape->Integral("width") << endl;
  cout << "RooDoubleCBShape: Analytic integral from func = " 
       << doubleCBShape_integral(m0,sL,sR,aL,aR,oL,oR,Rm[0],Rm[1])
       << endl;
  cout << endl;
  cout << "RooBifurcatedCBShape: Numeric integral from hist  = " << h_bifurcatedCBShape->Integral("width") << endl;
  cout << "RooBifurcatedCBShape: Analytic integral from func = " 
       << bifurcatedCBShape_integral(m0,sL,sR,aL,oL,Rm[0],Rm[1])
       << endl;
  cout << endl;

  TString EPSName;
  EPSName = TString("Plots/Pdf_test") + TString(".eps");
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
  c1->SetBottomMargin(0.12);
  c1->SetLeftMargin(0.12);

  c1->Print(EPSNameO.Data());

  c1->Clear();
  c1->Divide(1,2);
  c1->cd(1);
  gPad->SetFillColor(10);
  gPad->SetFrameFillColor(10);
  gPad->SetTickx(1);
  gPad->SetTicky(1);
  gPad->SetBottomMargin(0.12);
  gPad->SetLeftMargin(0.12);
  h_doubleCBShape->Draw("l");
  h_gauss->Draw("lsame");
  h_bifurcatedCBShape->Draw("lsame");
  double Max = TMath::Max(h_doubleCBShape->GetMaximum(),
			  h_gauss->GetMaximum());
  
  TLine* l1 = new TLine(m0,0.0,m0,Max);
  l1->SetLineColor(1);
  l1->SetLineWidth(2);
  l1->SetLineStyle(2);

  TLine* l2 = new TLine(m0-sL*aL,0.0,m0-sL*aL,
			h_doubleCBShape->GetBinContent(h_doubleCBShape->FindBin(m0-sL*aL)));
  l2->SetLineColor(1);
  l2->SetLineWidth(2);
  l2->SetLineStyle(2);

  TLine* l3 = new TLine(m0+sR*aR,0.0,m0+sR*aR,
			h_doubleCBShape->GetBinContent(h_doubleCBShape->FindBin(m0+sR*aR)));
  l3->SetLineColor(1);
  l3->SetLineWidth(2);
  l3->SetLineStyle(2);

  l1->Draw();
  l2->Draw();
  l3->Draw();
  c1->cd(2);
  gPad->SetFillColor(10);
  gPad->SetFrameFillColor(10);
  gPad->SetTickx(1);
  gPad->SetTicky(1);
  gPad->SetBottomMargin(0.12);
  gPad->SetLeftMargin(0.12);
  h_diff->Draw("l");
  h_diff2->Draw("lsame");
  c1->Print(EPSName.Data());

  c1->Print(EPSNameC.Data());

}
//==================================================================================
double doubleCBShape(double m,
		     double m0,
		     double sL,
		     double sR,
		     double aL,
		     double aR,
		     double oL,
		     double oR)
{

  if(aL < 0 || aR < 0) {
    cout << endl;
    cout << "RooDoubleCBShape::alphaL and alphaR parameters need to be positive." << endl;
    cout << "RooDoubleCBShape::Exiting now!!!" << endl;
    cout << endl;
    assert(false);
  }

  double Dm = m-m0;
  if(Dm < 0) {
    double t = fabs(Dm/sL);
    if(t <= aL) {
      return exp(-0.5*t*t);
    }
    else {
      double a = TMath::Power(oL/aL,oL)*exp(-0.5*aL*aL);
      double b = oL/aL - aL;
      return a/TMath::Power(b+t,oL);
    }
  }
  else {
    double t = fabs(Dm/sR);
    if(t <= aR) {
      return exp(-0.5*t*t);
    }
    double a = TMath::Power(oR/aR,oR)*exp(-0.5*aR*aR);
    double b = oR/aR - aR;
    return a/TMath::Power(b+t,oR);
  }

}
//==================================================================================
double gaussian(double m,
		double m0,
		double s)
{

  double Dm = m-m0;
  double t = fabs(Dm/s);
  return exp(-0.5*t*t);

  return 0.0;

}
//==================================================================================
double doubleCBShape_integral(double m0,
			      double sL,
			      double sR,
			      double aL,
			      double aR,
			      double oL,
			      double oR,
			      double m_min,
			      double m_max)
{

  if(aL < 0 || aR < 0) {
    cout << endl;
    cout << "PeakPi0EtaPdf::alphaL and alphaR parameters need to be positive." << endl;
    cout << "PeakPi0EtaPdf::Exiting now!!!" << endl;
    cout << endl;
    assert(false);
  }

  if(m_min > m_max) {
    cout << endl;
    cout << "PeakPi0EtaPdf::m_min is higher then m_max." << endl;
    cout << "PeakPi0EtaPdf::Exiting now!!!" << endl;
    cout << endl;
    assert(false);
  }

  if(m_min <= m0 - aL*sL) {
    if(m_max <= m0 - aL*sL) {
      double integral = 0.0;
      integral += TMath::Power((oL/aL) - aL + fabs((m_max-m0)/sL),-(oL-1));
      integral -= TMath::Power((oL/aL) - aL + fabs((m_min-m0)/sL),-(oL-1));
      integral *= TMath::Power(oL/aL,oL)*exp(-0.5*aL*aL)*sL/((oL-1));

      return integral;
    }
    else if(m_max > m0 - aL*sL && 
	    m_max <= m0) {
      double integral  = 0.0;
      double integral1 = 0.0;
      double integral2 = 0.0;

      integral1 += TMath::Power((oL/aL) - aL + fabs(aL),-(oL-1));
      integral1 -= TMath::Power((oL/aL) - aL + fabs((m_min-m0)/sL),-(oL-1));
      integral1 *= TMath::Power(oL/aL,oL)*exp(-0.5*aL*aL)*sL/(oL-1);

      integral2 += TMath::Erf(aL/sqrt(2.0));
      integral2 -= TMath::Erf(fabs((m_max-m0)/sL)/sqrt(2.0));
      integral2 *= sL*sqrt(TMath::Pi()/2.0);

      integral += integral1;
      integral += integral2;

      return integral;
    }
    else if(m_max > m0 && 
	    m_max <= m0 + aR*sR) {
      double integral  = 0.0;
      double integral1 = 0.0;
      double integral2 = 0.0;
      double integral3 = 0.0;

      integral1 += TMath::Power((oL/aL) - aL + fabs(aL),-(oL-1));
      integral1 -= TMath::Power((oL/aL) - aL + fabs((m_min-m0)/sL),-(oL-1));
      integral1 *= TMath::Power(oL/aL,oL)*exp(-0.5*aL*aL)*sL/(oL-1);

      integral2 += TMath::Erf(aL/sqrt(2.0));
      integral2 *= sL*sqrt(TMath::Pi()/2.0);

      integral3 += TMath::Erf(((m_max-m0)/sR)/sqrt(2.0));
      integral3 *= sR*sqrt(TMath::Pi()/2.0);

      integral += integral1;
      integral += integral2;
      integral += integral3;

      return integral;
    }
    else if(m_max > m0 + aR*sR) {
      double integral  = 0.0;
      double integral1 = 0.0;
      double integral2 = 0.0;
      double integral3 = 0.0;
      double integral4 = 0.0;

      integral1 += TMath::Power((oL/aL) - aL + fabs(aL),-(oL-1));
      integral1 -= TMath::Power((oL/aL) - aL + fabs((m_min-m0)/sL),-(oL-1));
      integral1 *= TMath::Power(oL/aL,oL)*exp(-0.5*aL*aL)*sL/(oL-1);

      integral2 += TMath::Erf(aL/sqrt(2.0));
      integral2 *= sL*sqrt(TMath::Pi()/2.0);

      integral3 += TMath::Erf(aR/sqrt(2.0));
      integral3 *= sR*sqrt(TMath::Pi()/2.0);

      integral4 += TMath::Power((oR/aR) - aR + fabs((m_max-m0)/sR),-(oR-1));
      integral4 -= TMath::Power((oR/aR) - aR + fabs(aR),-(oR-1));
      integral4 *= TMath::Power(oR/aR,oR)*exp(-0.5*aR*aR)*sR/(-(oR-1));

      integral += integral1;
      integral += integral2;
      integral += integral3;
      integral += integral4;

      return integral;
    }
  }
  else if(m_min > m0 - aL*sL && 
	  m_min <= m0) {
    if(m_max > m0 - aL*sL && 
       m_max <= m0) {
      double integral  = 0.0;
      double integral1 = 0.0;

      integral1 += TMath::Erf(fabs((m_min-m0)/sL)/sqrt(2.0));
      integral1 -= TMath::Erf(fabs((m_max-m0)/sL)/sqrt(2.0));
      integral1 *= sL*sqrt(TMath::Pi()/2.0);

      integral += integral1;

      return integral;
    }
    if(m_max > m0 && 
       m_max <= m0 + aR*sR) {
      double integral  = 0.0;
      double integral1 = 0.0;
      double integral2 = 0.0;

      integral1 += TMath::Erf(fabs((m_min-m0)/sL)/sqrt(2.0));
      integral1 *= sL*sqrt(TMath::Pi()/2.0);      

      integral2 += TMath::Erf(fabs((m_max-m0)/sR)/sqrt(2.0));
      integral2 *= sR*sqrt(TMath::Pi()/2.0);

      integral += integral1;
      integral += integral2;

      return integral;
    }
    if(m_max > m0 + aR*sR) {
      double integral  = 0.0;
      double integral1 = 0.0;
      double integral2 = 0.0;
      double integral3 = 0.0;

      integral1 += TMath::Erf(fabs((m_min-m0)/sL)/sqrt(2.0));
      integral1 *= sL*sqrt(TMath::Pi()/2.0);

      integral2 += TMath::Erf(fabs(aR/sqrt(2.0)));
      integral2 *= sR*sqrt(TMath::Pi()/2.0);

      integral3 += TMath::Power((oR/aR) - aR + fabs((m_max-m0)/sR),-(oR-1));
      integral3 -= TMath::Power((oR/aR) - aR + fabs(aR),-(oR-1));
      integral3 *= TMath::Power(oR/aR,oR)*exp(-0.5*aR*aR)*sR/(-(oR-1));

      integral += integral1;
      integral += integral2;
      integral += integral3;

      return integral;
    }
  }
  else if(m_min > m0 && 
	  m_min <= m0 + aR*sR) {
    if(m_max > m0 && 
       m_max <= m0 + aR*sR) {
      double integral  = 0.0;
      double integral1 = 0.0;

      integral1 += TMath::Erf(fabs((m_max-m0)/sR)/sqrt(2.0));
      integral1 -= TMath::Erf(fabs((m_min-m0)/sR)/sqrt(2.0));
      integral1 *= sR*sqrt(TMath::Pi()/2.0);

      integral += integral1;

      return integral;
    }
    else if(m_max > m0 + aR*sR) {
      double integral  = 0.0;
      double integral1 = 0.0;
      double integral2 = 0.0;

      integral1 += TMath::Erf(fabs(aR/sqrt(2.0)));
      integral1 -= TMath::Erf(fabs((m_min-m0)/sR)/sqrt(2.0));
      integral1 *= sR*sqrt(TMath::Pi()/2.0);

      integral2 += TMath::Power((oR/aR) - aR + fabs((m_max-m0)/sR),-(oR-1));
      integral2 -= TMath::Power((oR/aR) - aR + fabs(aR),-(oR-1));
      integral2 *= TMath::Power(oR/aR,oR)*exp(-0.5*aR*aR)*sR/(-(oR-1));

      integral += integral1;
      integral += integral2;

      return integral;
    }
  }
  else if(m_min > m0 + aR*sR) {
    if(m_max > m0 + aR*sR) {
      double integral  = 0.0;
      double integral1 = 0.0;

      integral1 += TMath::Power((oR/aR) - aR + fabs((m_max-m0)/sR),-(oR-1));
      integral1 -= TMath::Power((oR/aR) - aR + fabs((m_min-m0)/sR),-(oR-1));
      integral1 *= TMath::Power(oR/aR,oR)*exp(-0.5*aR*aR)*sR/(-(oR-1));

      integral += integral1;

      return integral;
    }
  }

  return -1.0;

}
//==================================================================================
double BifurcatedCBShape(double m,
			 double m0,
			 double sL,
			 double sR,
			 double a,
			 double o)
{

  if(a < 0 || sL < 0 || sR < 0) {
    cout << endl;
    cout << "RooBifurcatedCBShape::alpha or sigmaR or sigmaL parameters need to be positive." << endl;
    cout << "RooBifurcatedCBShape::Exiting now!!!" << endl;
    cout << endl;
    assert(false);
  }


  double Dm = m-m0;
  if(Dm < 0) {
    double t = fabs(Dm/sL);
    if(t <= a) {
      return exp(-0.5*t*t);
    }
    else {
      double ap = TMath::Power(o/a,o)*exp(-0.5*a*a);
      double bp = o/a - a;
      return ap/TMath::Power(bp + t,o);
    }
  }
  else {
    double t = fabs(Dm/sR);
    return exp(-0.5*t*t);
  }

  return -1.0;

}
//==================================================================================
double bifurcatedCBShape_integral(double m0,
				  double sL,
				  double sR,
				  double a,
				  double o,
				  double m_min,
				  double m_max)
{

  if(a < 0 || sR < 0 || sL < 0) {
    cout << endl;
    cout << "RooBifurcatedCBShape::alpha  or sigmaR or sigmaL parameters need to be positive." << endl;
    cout << "RooBifurcatedCBShape::Exiting now!!!" << endl;
    cout << endl;
    assert(false);
  }

  if(m_min > m_max) {
    cout << endl;
    cout << "RooBifurcatedCBShape::m_min is higher then m_max." << endl;
    cout << "RooBifurcatedCBShape::Exiting now!!!" << endl;
    cout << endl;
    assert(false);
  }

  if(m_min <= m0 - a*sL) {
    if(m_max <= m0 - a*sL) {

      double integral = 0.0;
      integral += TMath::Power((o/a) - a + fabs((m_max-m0)/sL),-(o-1));
      integral -= TMath::Power((o/a) - a + fabs((m_min-m0)/sL),-(o-1));
      integral *= TMath::Power(o/a,o)*exp(-0.5*a*a)*sL/((o-1));

      return integral;
    }
    else if(m_max > m0 - a*sL && 
	    m_max <= m0) {

      double integral  = 0.0;
      double integral1 = 0.0;
      double integral2 = 0.0;

      integral1 += TMath::Power((o/a) - a + fabs(a),-(o-1));
      integral1 -= TMath::Power((o/a) - a + fabs((m_min-m0)/sL),-(o-1));
      integral1 *= TMath::Power(o/a,o)*exp(-0.5*a*a)*sL/(o-1);

      integral2 += TMath::Erf(a/sqrt(2.0));
      integral2 -= TMath::Erf(fabs((m_max-m0)/sL)/sqrt(2.0));
      integral2 *= sL*sqrt(TMath::Pi()/2.0);

      integral += integral1;
      integral += integral2;

      return integral;
    }
    else if(m_max > m0) {
      double integral  = 0.0;
      double integral1 = 0.0;
      double integral2 = 0.0;
      double integral3 = 0.0;

      integral1 += TMath::Power((o/a) - a + fabs(a),-(o-1));
      integral1 -= TMath::Power((o/a) - a + fabs((m_min-m0)/sL),-(o-1));
      integral1 *= TMath::Power(o/a,o)*exp(-0.5*a*a)*sL/(o-1);

      integral2 += TMath::Erf(a/sqrt(2.0));
      integral2 *= sL*sqrt(TMath::Pi()/2.0);

      integral3 += TMath::Erf(((m_max-m0)/sR)/sqrt(2.0));
      integral3 *= sR*sqrt(TMath::Pi()/2.0);

      integral += integral1;
      integral += integral2;
      integral += integral3;

      return integral;
    }
  }
  else if(m_min > m0 - a*sL && 
	  m_min <= m0) {
    if(m_max > m0 - a*sL && 
       m_max <= m0) {
      double integral  = 0.0;
      double integral1 = 0.0;

      integral1 += TMath::Erf(fabs((m_min-m0)/sL)/sqrt(2.0));
      integral1 -= TMath::Erf(fabs((m_max-m0)/sL)/sqrt(2.0));
      integral1 *= sL*sqrt(TMath::Pi()/2.0);

      integral += integral1;

      return integral;
    }
    if(m_max > m0) {
      double integral  = 0.0;
      double integral1 = 0.0;
      double integral2 = 0.0;

      integral1 += TMath::Erf(fabs((m_min-m0)/sL)/sqrt(2.0));
      integral1 *= sL*sqrt(TMath::Pi()/2.0);      

      integral2 += TMath::Erf(fabs((m_max-m0)/sR)/sqrt(2.0));
      integral2 *= sR*sqrt(TMath::Pi()/2.0);

      integral += integral1;
      integral += integral2;

      return integral;
    }
  }
  else if(m_min > m0) {
    if(m_max > m0) {
      double integral  = 0.0;
      double integral1 = 0.0;

      integral1 += TMath::Erf(fabs((m_max-m0)/sR)/sqrt(2.0));
      integral1 -= TMath::Erf(fabs((m_min-m0)/sR)/sqrt(2.0));
      integral1 *= sR*sqrt(TMath::Pi()/2.0);

      integral += integral1;

      return integral;
    }
  }

  return -1.0;

}
//==================================================================================
