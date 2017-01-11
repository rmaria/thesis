//==================================================================================
void Pdftest_PolyAndExpDecayV3(void)
{

  TString HistTitle,HistName;

  double TheSize  = 0.06;
  char ytitle[100];

  double m;
  double rho   = 0.884286;
  double alpha = 0.282634;
  double m0    = 5.27;
  double mt    = 5.28970;

  double Rm[2];
  Rm[0] = 5.22;
  Rm[1] = 5.29;

 double Rm_g[2];
  Rm_g[0] = 5.22;
  Rm_g[1] = 5.29;

  int Nbins = 10000;
  TH1D* h_polyanddecay = new TH1D("h_polyanddecay","The PolyAndDecay pdf",
				  Nbins,
				  Rm[0],Rm[1]);
  h_polyanddecay->SetXTitle("m (GeV/c^{2})");
  h_polyanddecay->GetXaxis()->CenterTitle(true);
  sprintf(ytitle,"Events / (%.3f MeV/c^{2})",h_polyanddecay->GetBinWidth(0)*1.0e+3);
  h_polyanddecay->SetYTitle(ytitle);
  h_polyanddecay->GetYaxis()->CenterTitle(true);
  h_polyanddecay->SetLineColor(2);
  h_polyanddecay->SetLineWidth(2);
  //h_polyanddecay->SetMinimum(0.0);
  h_polyanddecay->GetXaxis()->SetTitleSize(TheSize);
  h_polyanddecay->GetXaxis()->SetLabelSize(TheSize);
  h_polyanddecay->GetYaxis()->SetTitleSize(TheSize);
  h_polyanddecay->GetYaxis()->SetLabelSize(TheSize);
  h_polyanddecay->SetStats(false);

  for(int i=0;i<Nbins;i++) {
    m = Rm[0] + (i+0.5)*(Rm[1]-Rm[0])/Nbins;
    double val = polyanddecay(m,rho,alpha,m0,mt,Rm_g[0],Rm_g[1]);
    h_polyanddecay->SetBinContent(i+1,val);
  }

  cout << endl;
  cout << "RooDoubleCBShape: Numeric integral from hist  = " << h_polyanddecay->Integral("width") << endl;
  cout << "RooDoubleCBShape: Analytic integral from func = " 
       << polyanddecay_integral(rho,alpha,m0,mt,Rm_g[0],Rm_g[1],Rm[0],Rm[1])
       << endl;
  cout << endl;

  TString EPSName;
  EPSName = TString("Plots/Pdf_test_PolyAndExpDecayV3") + TString(".eps");
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
  h_polyanddecay->Draw("l");
  c1->Print(EPSName.Data());

  c1->Print(EPSNameC.Data());

  return;

}
//==================================================================================
double polyanddecay(double m,
		    double rho,
		    double alpha,
		    double m0,
		    double mt,
		    double R_min,
		    double R_max)
{

  if(m >= mt) return 0.0;

  double m_min = m0;
  //double m_max = R_max;
  double m_max = mt;
  if(m_min < R_min) m_min = R_min;
  double w = m_max - m_min;

  double a2 = rho*TMath::Cos(0.5*TMath::Pi()*(alpha+1));
  double a3 = 0.5*(rho*TMath::Sin(0.5*TMath::Pi()*(alpha+1)) - 1.0);

  double a1 = -1.0 - a2 - a3;

  if(m >= m0) {
    double x  = -1 + 2*(m - m_min)/(m_max - m_min);
    double x2 = x*x;

    double sum = 0.0;
    sum += 1;
    sum += a1*x;
    sum += a2*(2.0*x2 - 1.0);
    sum += a3*x*(4.0*x2 - 3.0);

    return sum;
  }
  else {
    double p2_m0   = 1.0 - a1 + a2 - a3;
    double p2_m0_p = (2.0/w)*(a1 - 4.0*a2 + 9.0*a3);
    

    return p2_m0*TMath::Exp(-(p2_m0_p/p2_m0)*(m0-m));
  }

  return -1.0;

}
//==================================================================================
double polyanddecay_integral(double rho,
			     double alpha,
			     double m0,
			     double mt,
			     double R_min_g,
			     double R_max_g,
			     double R_min_l,
			     double R_max_l)
{

  double m_min = m0;
  //double m_max = R_max_g;
  double m_max = mt;
  if(m_min < R_min_g) m_min = R_min_g;
  double w = m_max - m_min;

  double m_min_int = R_min_l;
  double m_max_int = R_max_l;

  double a2 = rho*TMath::Cos(0.5*TMath::Pi()*(alpha+1));
  double a3 = 0.5*(rho*TMath::Sin(0.5*TMath::Pi()*(alpha+1)) - 1.0);
  double a1 = -1.0 - a2 - a3;

  if(m_min_int > m_max_int) {
    cout << endl;
    cout << "RooPolyAndDecay::m_min is higher then m_max." << endl;
    cout << "RooPolyAndDecay::Exiting now!!!" << endl;
    cout << endl;
    assert(false);
  }

  if(m_max_int < m0) {
    //integrate polynomial-decay:
    double tot_integral = 0.0;

    double p2_m0   = 1.0 - a1 + a2 - a3;
    double p2_m0_p = (2.0/w)*(a1 - 4.0*a2 + 9.0*a3);

    double limit_a = m_min_int;
    double limit_b = m_max_int;

    double int_a = (pow(p2_m0,2)/p2_m0_p)*TMath::Exp(-(p2_m0_p/p2_m0)*(m0 - limit_a));
    double int_b = (pow(p2_m0,2)/p2_m0_p)*TMath::Exp(-(p2_m0_p/p2_m0)*(m0 - limit_b));
    tot_integral += int_b - int_a;

    return tot_integral;

  }
  else if(m_min_int < m0 && m_max_int >= m0) {
    //partially integrate polynominal and polynominal-decay:
    double tot_integral = 0.0;

    double p2_m0   = 1.0 - a1 + a2 - a3;
    double p2_m0_p = (2.0/w)*(a1 - 4.0*a2 + 9.0*a3);

    double limit_a = m_min_int;
    double limit_b = m0;

    double int_a = (pow(p2_m0,2)/p2_m0_p)*TMath::Exp(-(p2_m0_p/p2_m0)*(m0 - limit_a));
    double int_b = (pow(p2_m0,2)/p2_m0_p)*TMath::Exp(-(p2_m0_p/p2_m0)*(m0 - limit_b));
    tot_integral += int_b - int_a;

    limit_a = m0;
    limit_b = m_max_int;
    if(m_max_int >= mt) limit_b = mt;

    double x_a  = -1 + 2*(limit_a - m0)/(mt - m0);
    double x_b  = -1 + 2*(limit_b - m0)/(mt - m0);

    int_a = ((mt - m0)/2.0)*(x_a + 0.5*a1*x_a*x_a + a2*x_a*((2./3.)*x_a*x_a - 1.0) + a3*x_a*x_a*(x_a*x_a - (3./2.)));
    int_b = ((mt - m0)/2.0)*(x_b + 0.5*a1*x_b*x_b + a2*x_b*((2./3.)*x_b*x_b - 1.0) + a3*x_b*x_b*(x_b*x_b - (3./2.)));

    tot_integral += int_b - int_a;

    return tot_integral;
  }
  else if(m_min_int >= m0) {
    //integrate polynominal:
    double tot_integral = 0.0;

    double limit_a = m_min_int;
    double limit_b = m_max_int;
    if(m_max_int >= mt) limit_b = mt;

    double x_a  = -1 + 2*(limit_a - m0)/(mt - m0);
    double x_b  = -1 + 2*(limit_b - m0)/(mt - m0);

    double int_a, int_b;

    int_a = ((mt - m0)/2.0)*(x_a + 0.5*a1*x_a*x_a + a2*x_a*((2./3.)*x_a*x_a - 1.0) + a3*x_a*x_a*(x_a*x_a - (3./2.)));
    int_b = ((mt - m0)/2.0)*(x_b + 0.5*a1*x_b*x_b + a2*x_b*((2./3.)*x_b*x_b - 1.0) + a3*x_b*x_b*(x_b*x_b - (3./2.)));

    tot_integral += int_b - int_a;

    return tot_integral;
  }

  return -1.0;

}
//==================================================================================
double f_ellipse(double a2,
		 int up)
{

  if(TMath::Abs(a2) > 1.0) return -1.0;

  double b   =  1.0;
  double a   =  0.5;
  double a20 =  0.0;
  double a30 = -0.5;

  if(up != 1 && up != 0) {
    cout << endl;
    cout << "The up parameter can only have the values 0 and 1." << endl;
    cout << "Exiting now!!!" << endl;
    cout << endl;
    assert(false);
  }

  if(up == 1) return a30 + a*sqrt(1.0 - pow((a2-a20)/b,2));
  else        return a30 - a*sqrt(1.0 - pow((a2-a20)/b,2));

}
//==================================================================================
double theline(double a2)
{

  if(TMath::Abs(a2) > 1.0) return -1.0;

  double alpha = -3.0/8.0;

  return alpha*(a2+1.0) + 0.5;

}
//==================================================================================



