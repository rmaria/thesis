//==================================================================================
void Pdftest_PolyAndDecay(void)
{

  double TheSize  = 0.06;
  char ytitle[100];

  double m;
  double a1 = -2.29045e-01;
  double a2 = -3.21106e-01;
  double a3 = -1.15583e-01;
  double order = 2.0;
  double m0 = 5.27;

  double gamma = 10;

  double order1 = 2.0;
  double order2 = 5.0;

  double Rm[2];
  Rm[0] = 5.22;
  Rm[1] = 5.29;

 double Rm_g[2];
  Rm_g[0] = 5.22;
  Rm_g[1] = 5.29;

  int Nbins = 100000;
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
  h_polyanddecay->SetMinimum(0.0);
  h_polyanddecay->GetXaxis()->SetTitleSize(TheSize);
  h_polyanddecay->GetXaxis()->SetLabelSize(TheSize);
  h_polyanddecay->GetYaxis()->SetTitleSize(TheSize);
  h_polyanddecay->GetYaxis()->SetLabelSize(TheSize);
  h_polyanddecay->SetStats(false);

  TH1D* h_Ratio_polyanddecay = new TH1D("h_Ratio_polyanddecay","The PolyAndDecay pdf Ratio",
					Nbins,
					Rm[0],m0);
  h_Ratio_polyanddecay->SetXTitle("m (GeV/c^{2})");
  h_Ratio_polyanddecay->GetXaxis()->CenterTitle(true);
  sprintf(ytitle,"Events / (%.3f MeV/c^{2})",h_Ratio_polyanddecay->GetBinWidth(0)*1.0e+3);
  h_Ratio_polyanddecay->SetYTitle(ytitle);
  h_Ratio_polyanddecay->GetYaxis()->CenterTitle(true);
  h_Ratio_polyanddecay->SetLineColor(2);
  h_Ratio_polyanddecay->SetLineWidth(2);
  //h_Ratio_polyanddecay->SetMinimum(0.0);
  h_Ratio_polyanddecay->GetXaxis()->SetTitleSize(TheSize);
  h_Ratio_polyanddecay->GetXaxis()->SetLabelSize(TheSize);
  h_Ratio_polyanddecay->GetYaxis()->SetTitleSize(TheSize);
  h_Ratio_polyanddecay->GetYaxis()->SetLabelSize(TheSize);
  h_Ratio_polyanddecay->SetStats(false);

  for(int i=0;i<Nbins;i++) {
    m = Rm[0] + (i+0.5)*(Rm[1]-Rm[0])/Nbins;
    double val = polyanddecay(m,a1,a2,a3,order,gamma,m0,Rm_g[0],Rm_g[1]);
    h_polyanddecay->SetBinContent(i+1,val);
  }

  for(int i=0;i<Nbins;i++) {
    m = Rm[0] + (i+0.5)*(m0-Rm[0])/Nbins;
    double val = Ratio_polyanddecay(m,a1,a2,a3,order1,order2,gamma,m0,Rm_g[0],Rm_g[1]);
    h_Ratio_polyanddecay->SetBinContent(i+1,val);
  }

  cout << endl;
  cout << "RooDoubleCBShape: Numeric integral from hist  = " << h_polyanddecay->Integral("width") << endl;
  cout << "RooDoubleCBShape: Analytic integral from func = " 
       << polyanddecay_integral(a1,a2,a3,order,gamma,m0,Rm_g[0],Rm_g[1],Rm[0],Rm[1])
       << endl;
  cout << endl;

  cout << "Value at m0      is " << polyanddecay(m0,     a1,a2,a3,order,gamma,m0,Rm_g[0],Rm_g[1]) << endl;
  cout << "Value at Rm_g[0] is " << polyanddecay(Rm_g[0],a1,a2,a3,order,gamma,m0,Rm_g[0],Rm_g[1]) << endl;
  cout << "The Delta Value at m0 and Rm_g[0] is " 
       << polyanddecay(m0,a1,a2,a3,order,gamma,m0,Rm_g[0],Rm_g[1]) - 
    polyanddecay(Rm_g[0],a1,a2,a3,order,gamma,m0,Rm_g[0],Rm_g[1])
       << endl;

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
  h_polyanddecay->Draw("l");
  c1->Print(EPSName.Data());

  c1->Clear();
  h_Ratio_polyanddecay->Draw("l");
  c1->Print(EPSName.Data());

  c1->Print(EPSNameC.Data());

  return;

}
//==================================================================================
double polyanddecay(double m,
		    double a1,
		    double a2,
		    double a3,
		    double order,
		    double gamma,
		    double m0,
		    double R_min,
		    double R_max)
{

  if(order < 0) {
    cout << endl;
    cout << "RooPolyAndDecay::order parameter need to be positive." << endl;
    cout << "RooPolyAndDecay::order = " << order << endl;
    cout << "RooPolyAndDecay::Exiting now!!!" << endl;
    cout << endl;
    assert(false);
  }

  double m_min = m0;
  double m_max = R_max;
  if(m_min < R_min) m_min = R_min;
  double w = m_max - m_min;

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

    return p2_m0/TMath::Power(1 + (p2_m0_p/(order*p2_m0))*(m0 - m)*gamma,order);
  }

  return -1.0;

}
//==================================================================================
double Ratio_polyanddecay(double m,
			  double a1,
			  double a2,
			  double a3,
			  double order1,
			  double order2,
			  double gamma,
			  double m0,
			  double R_min,
			  double R_max)
{

  if(order1 < 0 || order2 < 0) {
    cout << endl;
    cout << "RooPolyAndDecay::order1 or order2 parameters need to be positive." << endl;
    cout << "RooPolyAndDecay::order1 = " << order1 << endl;
    cout << "RooPolyAndDecay::order2 = " << order2 << endl;
    cout << "RooPolyAndDecay::Exiting now!!!" << endl;
    cout << endl;
    assert(false);
  }

  if(m > m0) {
    cout << endl;
    cout << "RooPolyAndDecay::m variable is higher than m0." << endl;
    cout << "RooPolyAndDecay::m = " << m << endl;
    cout << "RooPolyAndDecay::Exiting now!!!" << endl;
    cout << endl;
    assert(false);
  }

  double m_min = m0;
  double m_max = R_max;
  if(m_min < R_min) m_min = R_min;
  double w = m_max - m_min;

  double p2_m0   = 1.0 - a1 + a2 - a3;
  double p2_m0_p = (2.0/w)*(a1 - 4.0*a2 + 9.0*a3);

  double f1 = 1.0 + (p2_m0_p/(order1*p2_m0))*(m0 - m)*gamma;
  f1        = TMath::Power(f1,order1);
  double f2 = 1.0 + (p2_m0_p/(order2*p2_m0))*(m0 - m)*gamma;
  f2        = TMath::Power(f2,order2);
  
  return f2/f1;

}
//==================================================================================
double polyanddecay_integral(double a1,
			     double a2,
			     double a3,
			     double order,
			     double gamma,
			     double m0,
			     double R_min_g,
			     double R_max_g,
			     double R_min_l,
			     double R_max_l)
{

  if(order < 0) {
    cout << endl;
    cout << "RooPolyAndDecay::order parameter need to be positive." << endl;
    cout << "RooPolyAndDecay::Exiting now!!!" << endl;
    cout << endl;
    assert(false);
  }

  double m_min = m0;
  double m_max = R_max_g;
  if(m_min < R_min_g) m_min = R_min_g;
  double w = m_max - m_min;

  double m_min_int = R_min_l;
  double m_max_int = R_max_l;

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

    double factor = (order*pow(p2_m0,2))/(gamma*p2_m0_p*(order-1));

    double int_a = factor/TMath::Power(1.0 + ((p2_m0_p*gamma)/(order*p2_m0))*(m0 - limit_a),order-1);
    double int_b = factor/TMath::Power(1.0 + ((p2_m0_p*gamma)/(order*p2_m0))*(m0 - limit_b),order-1);
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

    double factor = (order*pow(p2_m0,2))/(gamma*p2_m0_p*(order-1));

    double int_a = factor/TMath::Power(1.0 + ((p2_m0_p*gamma)/(order*p2_m0))*(m0 - limit_a),order-1);
    double int_b = factor/TMath::Power(1.0 + ((p2_m0_p*gamma)/(order*p2_m0))*(m0 - limit_b),order-1);
    tot_integral += int_b - int_a;

    limit_a = m0;
    limit_b = m_max_int;

    double x_a  = -1 + 2*(limit_a - m0)/(R_max_g - m0);
    double x_b  = -1 + 2*(limit_b - m0)/(R_max_g - m0);

    int_a = ((R_max_g - m0)/2.0)*(x_a + 0.5*a1*x_a*x_a + a2*x_a*((2./3.)*x_a*x_a - 1.0) + a3*x_a*x_a*(x_a*x_a - (3./2.)));
    int_b = ((R_max_g - m0)/2.0)*(x_b + 0.5*a1*x_b*x_b + a2*x_b*((2./3.)*x_b*x_b - 1.0) + a3*x_b*x_b*(x_b*x_b - (3./2.)));

    tot_integral += int_b - int_a;

    return tot_integral;
  }
  else if(m_min_int >= m0) {
    //integrate polynominal:
    double tot_integral = 0.0;

    double limit_a = m_min_int;
    double limit_b = m_max_int;

    double x_a  = -1 + 2*(limit_a - m0)/(R_max_g - m0);
    double x_b  = -1 + 2*(limit_b - m0)/(R_max_g - m0);

    double int_a, int_b;

    int_a = ((R_max_g - m0)/2.0)*(x_a + 0.5*a1*x_a*x_a + a2*x_a*((2./3.)*x_a*x_a - 1.0) + a3*x_a*x_a*(x_a*x_a - (3./2.)));
    int_b = ((R_max_g - m0)/2.0)*(x_b + 0.5*a1*x_b*x_b + a2*x_b*((2./3.)*x_b*x_b - 1.0) + a3*x_b*x_b*(x_b*x_b - (3./2.)));

    tot_integral += int_b - int_a;

    return tot_integral;
  }


  return -1.0;

}
//==================================================================================


