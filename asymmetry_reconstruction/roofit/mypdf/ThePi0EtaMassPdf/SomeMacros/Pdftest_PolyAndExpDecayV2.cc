//==================================================================================
void Pdftest_PolyAndExpDecayV2(void)
{

  TString HistTitle,HistName;

  double TheSize  = 0.06;
  char ytitle[100];

  double m;
  double a2 =  0.0;
  double a3 =  0.4;
  double m0 = 5.27;
  double mt = 5.28970;

  double R_a2[2];
  R_a2[0] = -1.1;
  R_a2[1] =  1.1;
  double R_a3[2];
  R_a3[0] = -1.1;
  R_a3[1] =  0.6;

  double Rm[2];
  Rm[0] = 5.22;
  Rm[1] = 5.29;

 double Rm_g[2];
  Rm_g[0] = 5.22;
  Rm_g[1] = 5.29;

  int Nbins = 1000;
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
    double val = polyanddecay(m,a2,a3,m0,mt,Rm_g[0],Rm_g[1]);
    h_polyanddecay->SetBinContent(i+1,val);
  }

  cout << endl;
  cout << "RooDoubleCBShape: Numeric integral from hist  = " << h_polyanddecay->Integral("width") << endl;
  cout << "RooDoubleCBShape: Analytic integral from func = " 
       << polyanddecay_integral(a2,a3,m0,mt,Rm_g[0],Rm_g[1],Rm[0],Rm[1])
       << endl;
  cout << endl;

  int Nbins_a = 100;
  TH2D* h_polyanddecay_pos = new TH2D("h_polyanddecay_pos","The a_{2} and a_{3} where the funtion is >= 0 for x in (-1,1)",
				      Nbins_a,R_a2[0],R_a2[1],
				      Nbins_a,R_a3[0],R_a3[1]);
  h_polyanddecay_pos->SetXTitle("a_{2}");
  h_polyanddecay_pos->GetXaxis()->CenterTitle(true);
  h_polyanddecay_pos->SetYTitle("a_{3}");
  h_polyanddecay_pos->GetYaxis()->CenterTitle(true);
  h_polyanddecay_pos->SetLineColor(4);
  h_polyanddecay_pos->SetLineWidth(2);
  h_polyanddecay_pos->SetMinimum(0.0);
  h_polyanddecay_pos->GetXaxis()->SetTitleSize(TheSize);
  h_polyanddecay_pos->GetXaxis()->SetLabelSize(TheSize);
  h_polyanddecay_pos->GetYaxis()->SetTitleSize(TheSize);
  h_polyanddecay_pos->GetYaxis()->SetLabelSize(TheSize);
  h_polyanddecay_pos->SetStats(false);

  TH2D* h_polyanddecay_2der_pos = new TH2D("h_polyanddecay_2der_pos","The a_{2} and a_{3} where the 2nd funtion derivative is >= 0 for x in (-1,1)",
					   Nbins_a,R_a2[0],R_a2[1],
					   Nbins_a,R_a3[0],R_a3[1]);
  h_polyanddecay_2der_pos->SetXTitle("a_{2}");
  h_polyanddecay_2der_pos->GetXaxis()->CenterTitle(true);
  h_polyanddecay_2der_pos->SetYTitle("a_{3}");
  h_polyanddecay_2der_pos->GetYaxis()->CenterTitle(true);
  h_polyanddecay_2der_pos->SetLineColor(4);
  h_polyanddecay_2der_pos->SetLineWidth(2);
  h_polyanddecay_2der_pos->SetMinimum(0.0);
  h_polyanddecay_2der_pos->GetXaxis()->SetTitleSize(TheSize);
  h_polyanddecay_2der_pos->GetXaxis()->SetLabelSize(TheSize);
  h_polyanddecay_2der_pos->GetYaxis()->SetTitleSize(TheSize);
  h_polyanddecay_2der_pos->GetYaxis()->SetLabelSize(TheSize);
  h_polyanddecay_2der_pos->SetStats(false);

  TH2D* h_polyanddecay_f_and_2der_pos = new TH2D("h_polyanddecay_f_and_2der_pos","The a_{2} and a_{3} where the f and the 2nd funtion derivative is >= 0 for x in (-1,1)",
						 Nbins_a,R_a2[0],R_a2[1],
						 Nbins_a,R_a3[0],R_a3[1]);
  h_polyanddecay_f_and_2der_pos->SetXTitle("a_{2}");
  h_polyanddecay_f_and_2der_pos->GetXaxis()->CenterTitle(true);
  h_polyanddecay_f_and_2der_pos->SetYTitle("a_{3}");
  h_polyanddecay_f_and_2der_pos->GetYaxis()->CenterTitle(true);
  h_polyanddecay_f_and_2der_pos->SetLineColor(4);
  h_polyanddecay_f_and_2der_pos->SetLineWidth(2);
  h_polyanddecay_f_and_2der_pos->SetMinimum(0.0);
  h_polyanddecay_f_and_2der_pos->GetXaxis()->SetTitleSize(TheSize);
  h_polyanddecay_f_and_2der_pos->GetXaxis()->SetLabelSize(TheSize);
  h_polyanddecay_f_and_2der_pos->GetYaxis()->SetTitleSize(TheSize);
  h_polyanddecay_f_and_2der_pos->GetYaxis()->SetLabelSize(TheSize);
  h_polyanddecay_f_and_2der_pos->SetStats(false);

  for(int ia2=0;ia2<Nbins_a;ia2++) {
    double a2_val = R_a2[0] + (ia2+0.5)*(R_a2[1]-R_a2[0])/Nbins_a;
    for(int ia3=0;ia3<Nbins_a;ia3++) {
      double a3_val = R_a3[0] + (ia3+0.5)*(R_a3[1]-R_a3[0])/Nbins_a;

      int N = 100;
      bool fIsPos = true;
      for(int i=0;i<N;i++) {
	double x = -1.0 + (i+0.5)*(2.0)/N;
	double a1_val = -1 - a2_val - a3_val;
	double f = 1 + a1_val*x + a2_val*(2.0*x*x - 1.0) + a3_val*x*(4*x*x - 3.0);
	if(f < 0) {
	  fIsPos = false;
	  break;
	}
      }

      bool f2derIsPos = true;
      for(int i=0;i<N;i++) {
	double x = -1.0 + (i+0.5)*(2.0)/N;
	double f = 4.0*a2_val + 24.0*a3_val*x;

	if(f > 0) {
	  f2derIsPos = false;
	  break;
	}
      }

      if(fIsPos) h_polyanddecay_pos->SetBinContent(ia2+1,ia3+1, 1.0);
      else       h_polyanddecay_pos->SetBinContent(ia2+1,ia3+1,-1.0);


      if(f2derIsPos) h_polyanddecay_2der_pos->SetBinContent(ia2+1,ia3+1,-1.0);
      else           h_polyanddecay_2der_pos->SetBinContent(ia2+1,ia3+1, 1.0);

      if(fIsPos && !f2derIsPos) h_polyanddecay_f_and_2der_pos->SetBinContent(ia2+1,ia3+1, 1.0);
      else                      h_polyanddecay_f_and_2der_pos->SetBinContent(ia2+1,ia3+1,-1.0);

    }
  }

  double R_a2_plots[2];
  R_a2_plots[0] = -1.0;
  R_a2_plots[1] =  1.0;
  double R_a2_plots_2[2];
  R_a2_plots_2[0] = -1.0;
  R_a2_plots_2[1] =  3./5.;
  int Npoints_h = 15;
  int Npoints_v = 10;

  //Intra-ellipse:
  TH1D* h_polyanddecay_intraellipse[100][100];
  TH1D* h_polyanddecay_outellipse_good[100][100];
  TH1D* h_polyanddecay_outellipse_bad_up[100][100];
  TH1D* h_polyanddecay_outellipse_bad_dn[100][100];

  double a2_points[100];
  double a2_points_v2[100];
  double a3_points_int_ellipse[100][100];
  double a3_points_out_ellipse_good[100][100];
  double a3_points_out_ellipse_bad_up[100][100];
  double a3_points_out_ellipse_bad_dn[100][100];

  for(int ih=0;ih<Npoints_h;ih++) {
    double a2_val = R_a2_plots[0] + (ih+0.5)*(R_a2_plots[1]-R_a2_plots[0])/Npoints_h;
    a2_points[ih] = a2_val;
    double R_a3_intra_ellipse[2];
    R_a3_intra_ellipse[0] = f_ellipse(a2_val,0);
    R_a3_intra_ellipse[1] = f_ellipse(a2_val,1);
    
    double R_a3_out_ellipse_bad_up[2];
    R_a3_out_ellipse_bad_up[1] = 0.5;
    if(a2_val <= 3./5.) R_a3_out_ellipse_bad_up[0] = theline(a2_val);
    else                R_a3_out_ellipse_bad_up[0] = f_ellipse(a2_val,1);
    double R_a3_out_ellipse_bad_dn[2];
    R_a3_out_ellipse_bad_dn[0] = -1.0;
    R_a3_out_ellipse_bad_dn[1] = f_ellipse(a2_val,0);

    for(int iv=0;iv<Npoints_v;iv++) {
      double a3_val;
      TString ThePoint;

      a3_val = R_a3_intra_ellipse[0] + (iv+0.5)*(R_a3_intra_ellipse[1]-R_a3_intra_ellipse[0])/Npoints_v;
      a3_points_int_ellipse[ih][iv] = a3_val;
      sprintf(ytitle,"%.2f",a2_val);
      ThePoint  = TString("(") + TString(ytitle) + TString(",");
      sprintf(ytitle,"%.2f",a3_val);
      ThePoint += TString(ytitle) + TString(")");

      HistName  = TString("h_polyanddecay_intraellipse_a2_") + long(ih+1) + TString("_a3_") + long(iv+1);
      HistTitle = TString("pdf for (a_{2},a_{3}) = ") + ThePoint;
      h_polyanddecay_intraellipse[ih][iv] = new TH1D(HistName.Data(),
						     HistTitle.Data(),
						     Nbins,
						     Rm[0],Rm[1]);
      h_polyanddecay_intraellipse[ih][iv]->SetXTitle("m (GeV/c^{2})");
      h_polyanddecay_intraellipse[ih][iv]->GetXaxis()->CenterTitle(true);
      sprintf(ytitle,"Events / (%.3f MeV/c^{2})",h_polyanddecay_intraellipse[ih][iv]->GetBinWidth(0)*1.0e+3);
      h_polyanddecay_intraellipse[ih][iv]->SetYTitle(ytitle);
      h_polyanddecay_intraellipse[ih][iv]->GetYaxis()->CenterTitle(true);
      h_polyanddecay_intraellipse[ih][iv]->SetLineColor(2);
      h_polyanddecay_intraellipse[ih][iv]->SetLineWidth(2);
      //h_polyanddecay_intraellipse[ih][iv]->SetMinimum(0.0);
      h_polyanddecay_intraellipse[ih][iv]->GetXaxis()->SetTitleSize(TheSize);
      h_polyanddecay_intraellipse[ih][iv]->GetXaxis()->SetLabelSize(TheSize);
      h_polyanddecay_intraellipse[ih][iv]->GetYaxis()->SetTitleSize(TheSize);
      h_polyanddecay_intraellipse[ih][iv]->GetYaxis()->SetLabelSize(TheSize);
      h_polyanddecay_intraellipse[ih][iv]->SetStats(false);

      for(int i=0;i<Nbins;i++) {
	m = Rm[0] + (i+0.5)*(Rm[1]-Rm[0])/Nbins;
	double val = polyanddecay(m,a2_val,a3_val,m0,mt,Rm_g[0],Rm_g[1]);
	h_polyanddecay_intraellipse[ih][iv]->SetBinContent(i+1,val);
      }
      
      a3_val = R_a3_out_ellipse_bad_up[0] + (iv+0.5)*(R_a3_out_ellipse_bad_up[1]-R_a3_out_ellipse_bad_up[0])/Npoints_v;
      a3_points_out_ellipse_bad_up[ih][iv] = a3_val;
      sprintf(ytitle,"%.2f",a2_val);
      ThePoint  = TString("(") + TString(ytitle) + TString(",");
      sprintf(ytitle,"%.2f",a3_val);
      ThePoint += TString(ytitle) + TString(")");

      HistName  = TString("h_polyanddecay_outellipse_bad_up_a2_") + long(ih+1) + TString("_a3_") + long(iv+1);
      HistTitle = TString("pdf for (a_{2},a_{3}) = ") + ThePoint;
      h_polyanddecay_outellipse_bad_up[ih][iv] = new TH1D(HistName.Data(),
							  HistTitle.Data(),
							  Nbins,
							  Rm[0],Rm[1]);
      h_polyanddecay_outellipse_bad_up[ih][iv]->SetXTitle("m (GeV/c^{2})");
      h_polyanddecay_outellipse_bad_up[ih][iv]->GetXaxis()->CenterTitle(true);
      sprintf(ytitle,"Events / (%.3f MeV/c^{2})",h_polyanddecay_outellipse_bad_up[ih][iv]->GetBinWidth(0)*1.0e+3);
      h_polyanddecay_outellipse_bad_up[ih][iv]->SetYTitle(ytitle);
      h_polyanddecay_outellipse_bad_up[ih][iv]->GetYaxis()->CenterTitle(true);
      h_polyanddecay_outellipse_bad_up[ih][iv]->SetLineColor(2);
      h_polyanddecay_outellipse_bad_up[ih][iv]->SetLineWidth(2);
      //h_polyanddecay_outellipse_bad_up[ih][iv]->SetMinimum(0.0);
      h_polyanddecay_outellipse_bad_up[ih][iv]->GetXaxis()->SetTitleSize(TheSize);
      h_polyanddecay_outellipse_bad_up[ih][iv]->GetXaxis()->SetLabelSize(TheSize);
      h_polyanddecay_outellipse_bad_up[ih][iv]->GetYaxis()->SetTitleSize(TheSize);
      h_polyanddecay_outellipse_bad_up[ih][iv]->GetYaxis()->SetLabelSize(TheSize);
      h_polyanddecay_outellipse_bad_up[ih][iv]->SetStats(false);

      for(int i=0;i<Nbins;i++) {
	m = Rm[0] + (i+0.5)*(Rm[1]-Rm[0])/Nbins;
	double val = polyanddecay(m,a2_val,a3_val,m0,mt,Rm_g[0],Rm_g[1]);
	h_polyanddecay_outellipse_bad_up[ih][iv]->SetBinContent(i+1,val);
      }

      a3_val = R_a3_out_ellipse_bad_dn[0] + (iv+0.5)*(R_a3_out_ellipse_bad_dn[1]-R_a3_out_ellipse_bad_dn[0])/Npoints_v;
      a3_points_out_ellipse_bad_dn[ih][iv] = a3_val;
      sprintf(ytitle,"%.2f",a2_val);
      ThePoint  = TString("(") + TString(ytitle) + TString(",");
      sprintf(ytitle,"%.2f",a3_val);
      ThePoint += TString(ytitle) + TString(")");

      HistName  = TString("h_polyanddecay_outellipse_bad_dn_a2_") + long(ih+1) + TString("_a3_") + long(iv+1);
      HistTitle = TString("pdf for (a_{2},a_{3}) = ") + ThePoint;
      h_polyanddecay_outellipse_bad_dn[ih][iv] = new TH1D(HistName.Data(),
							  HistTitle.Data(),
							  Nbins,
							  Rm[0],Rm[1]);
      h_polyanddecay_outellipse_bad_dn[ih][iv]->SetXTitle("m (GeV/c^{2})");
      h_polyanddecay_outellipse_bad_dn[ih][iv]->GetXaxis()->CenterTitle(true);
      sprintf(ytitle,"Events / (%.3f MeV/c^{2})",h_polyanddecay_outellipse_bad_dn[ih][iv]->GetBinWidth(0)*1.0e+3);
      h_polyanddecay_outellipse_bad_dn[ih][iv]->SetYTitle(ytitle);
      h_polyanddecay_outellipse_bad_dn[ih][iv]->GetYaxis()->CenterTitle(true);
      h_polyanddecay_outellipse_bad_dn[ih][iv]->SetLineColor(2);
      h_polyanddecay_outellipse_bad_dn[ih][iv]->SetLineWidth(2);
      //h_polyanddecay_outellipse_bad_dn[ih][iv]->SetMinimum(0.0);
      h_polyanddecay_outellipse_bad_dn[ih][iv]->GetXaxis()->SetTitleSize(TheSize);
      h_polyanddecay_outellipse_bad_dn[ih][iv]->GetXaxis()->SetLabelSize(TheSize);
      h_polyanddecay_outellipse_bad_dn[ih][iv]->GetYaxis()->SetTitleSize(TheSize);
      h_polyanddecay_outellipse_bad_dn[ih][iv]->GetYaxis()->SetLabelSize(TheSize);
      h_polyanddecay_outellipse_bad_dn[ih][iv]->SetStats(false);

      for(int i=0;i<Nbins;i++) {
	m = Rm[0] + (i+0.5)*(Rm[1]-Rm[0])/Nbins;
	double val = polyanddecay(m,a2_val,a3_val,m0,mt,Rm_g[0],Rm_g[1]);
	h_polyanddecay_outellipse_bad_dn[ih][iv]->SetBinContent(i+1,val);
      }
    }
  }


  for(int ih=0;ih<Npoints_h;ih++) {
    double a2_val = R_a2_plots_2[0] + (ih+0.5)*(R_a2_plots_2[1]-R_a2_plots_2[0])/Npoints_h;
    a2_points_v2[ih] = a2_val;
    double R_a3_out_ellipse_good[2];
    R_a3_out_ellipse_good[0] = f_ellipse(a2_val,1);
    R_a3_out_ellipse_good[1] = theline(a2_val);

    for(int iv=0;iv<Npoints_v;iv++) {
      double a3_val;
      TString ThePoint;

      a3_val = R_a3_out_ellipse_good[0] + (iv+0.5)*(R_a3_out_ellipse_good[1]-R_a3_out_ellipse_good[0])/Npoints_v;
      a3_points_out_ellipse_good[ih][iv] = a3_val;
      sprintf(ytitle,"%.2f",a2_val);
      ThePoint  = TString("(") + TString(ytitle) + TString(",");
      sprintf(ytitle,"%.2f",a3_val);
      ThePoint += TString(ytitle) + TString(")");

      HistName  = TString("h_polyanddecay_outellipse_good_a2_") + long(ih+1) + TString("_a3_") + long(iv+1);
      HistTitle = TString("pdf for (a_{2},a_{3}) = ") + ThePoint;
      h_polyanddecay_outellipse_good[ih][iv] = new TH1D(HistName.Data(),
							HistTitle.Data(),
							Nbins,
							Rm[0],Rm[1]);
      h_polyanddecay_outellipse_good[ih][iv]->SetXTitle("m (GeV/c^{2})");
      h_polyanddecay_outellipse_good[ih][iv]->GetXaxis()->CenterTitle(true);
      sprintf(ytitle,"Events / (%.3f MeV/c^{2})",h_polyanddecay_outellipse_good[ih][iv]->GetBinWidth(0)*1.0e+3);
      h_polyanddecay_outellipse_good[ih][iv]->SetYTitle(ytitle);
      h_polyanddecay_outellipse_good[ih][iv]->GetYaxis()->CenterTitle(true);
      h_polyanddecay_outellipse_good[ih][iv]->SetLineColor(2);
      h_polyanddecay_outellipse_good[ih][iv]->SetLineWidth(2);
      //h_polyanddecay_outellipse_good[ih][iv]->SetMinimum(0.0);
      h_polyanddecay_outellipse_good[ih][iv]->GetXaxis()->SetTitleSize(TheSize);
      h_polyanddecay_outellipse_good[ih][iv]->GetXaxis()->SetLabelSize(TheSize);
      h_polyanddecay_outellipse_good[ih][iv]->GetYaxis()->SetTitleSize(TheSize);
      h_polyanddecay_outellipse_good[ih][iv]->GetYaxis()->SetLabelSize(TheSize);
      h_polyanddecay_outellipse_good[ih][iv]->SetStats(false);

      for(int i=0;i<Nbins;i++) {
	m = Rm[0] + (i+0.5)*(Rm[1]-Rm[0])/Nbins;
	double val = polyanddecay(m,a2_val,a3_val,m0,mt,Rm_g[0],Rm_g[1]);
	h_polyanddecay_outellipse_good[ih][iv]->SetBinContent(i+1,val);
      }

    }
  }

  TString EPSName;
  EPSName = TString("Plots/Pdf_test_PolyAndExpDecayV2") + TString(".eps");
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
  c1->Divide(2,2);
  c1->cd(1);
  gPad->SetFillColor(10);
  gPad->SetFrameFillColor(10);
  gPad->SetTickx(1);
  gPad->SetTicky(1);
  gPad->SetBottomMargin(0.12);
  gPad->SetLeftMargin(0.12);
  h_polyanddecay_pos->Draw("col");
  c1->cd(2);
  gPad->SetFillColor(10);
  gPad->SetFrameFillColor(10);
  gPad->SetTickx(1);
  gPad->SetTicky(1);
  gPad->SetBottomMargin(0.12);
  gPad->SetLeftMargin(0.12);
  h_polyanddecay_2der_pos->Draw("col");
  c1->cd(3);
  gPad->SetFillColor(10);
  gPad->SetFrameFillColor(10);
  gPad->SetTickx(1);
  gPad->SetTicky(1);
  gPad->SetBottomMargin(0.12);
  gPad->SetLeftMargin(0.12);
  h_polyanddecay_f_and_2der_pos->Draw("col");
  TLine* l1 = new TLine(h_polyanddecay_f_and_2der_pos->GetXaxis()->GetXmin(),-1.0,
			h_polyanddecay_f_and_2der_pos->GetXaxis()->GetXmax(),-1.0);
  l1->SetLineColor(4);
  l1->SetLineWidth(2);
  l1->SetLineStyle(2);
  l1->Draw();
  TLine* l2 = new TLine(h_polyanddecay_f_and_2der_pos->GetXaxis()->GetXmin(),0.5,
			h_polyanddecay_f_and_2der_pos->GetXaxis()->GetXmax(),0.5);
  l2->SetLineColor(4);
  l2->SetLineWidth(2);
  l2->SetLineStyle(2);
  l2->Draw();
  TLine* l3 = new TLine(-1.0,h_polyanddecay_f_and_2der_pos->GetYaxis()->GetXmin(),
			-1.0,h_polyanddecay_f_and_2der_pos->GetYaxis()->GetXmax());
  l3->SetLineColor(4);
  l3->SetLineWidth(2);
  l3->SetLineStyle(2);
  l3->Draw();
  TLine* l4 = new TLine(1.0,h_polyanddecay_f_and_2der_pos->GetYaxis()->GetXmin(),
			1.0,h_polyanddecay_f_and_2der_pos->GetYaxis()->GetXmax());
  l4->SetLineColor(4);
  l4->SetLineWidth(2);
  l4->SetLineStyle(2);
  l4->Draw();

  double t = 1.0/8.0;
  TLine* l5 = new TLine(0.0,h_polyanddecay_f_and_2der_pos->GetYaxis()->GetXmin(),
			0.0,t);
  l5->SetLineColor(4);
  l5->SetLineWidth(2);
  l5->SetLineStyle(2);
  l5->Draw();

  double a = (0.5-t)/(-1.0-0.0);
  double b = 0.5 - a*(-1.0);

  TLine* l6 = new TLine(-1.0,0.5,
			1.0,a*1.0+b);
  l6->SetLineColor(4);
  l6->SetLineWidth(2);
  l6->SetLineStyle(2);
  l6->Draw();

  TLine* l7 = new TLine(-1.0,1.0/6.0,
			0.0,0.0);
  l7->SetLineColor(4);
  l7->SetLineWidth(2);
  l7->SetLineStyle(2);
  l7->Draw();

  TLine* l8 = new TLine(-1.0,-1.0/6.0,
			0.0,0.0);
  l8->SetLineColor(4);
  l8->SetLineWidth(2);
  l8->SetLineStyle(2);
  l8->Draw();

  TEllipse* e = new TEllipse(0.0,-0.5,1.0,0.5);
  e->SetLineColor(1);
  e->SetLineWidth(2);
  e->SetLineStyle(2);
  e->SetFillStyle(0);
  e->Draw("same");

  TEllipse* p = new TEllipse(a2,a3,4.0e-2,4.0e-2);
  p->SetLineColor(1);
  p->SetLineWidth(2);
  p->SetFillColor(1);
  p->Draw("same");

  c1->cd(4);
  gPad->SetFillColor(10);
  gPad->SetFrameFillColor(10);
  gPad->SetTickx(1);
  gPad->SetTicky(1);
  gPad->SetBottomMargin(0.12);
  gPad->SetLeftMargin(0.12);
  h_polyanddecay->Draw("l");
  c1->Print(EPSName.Data());

  c1->Clear();
  h_polyanddecay_f_and_2der_pos->Draw("col");

  TF1 f_up_ellip("f_up_ellip","f_ellipse(x,[0])",-1.0,1.0);
  f_up_ellip.SetParameter(0,1);
  f_up_ellip.SetLineColor(4);
  f_up_ellip.SetLineWidth(2);

  TF1 f_dn_ellip("f_up_ellip","f_ellipse(x,[0])",-1.0,1.0);
  f_dn_ellip.SetParameter(0,0);
  f_dn_ellip.SetLineColor(4);
  f_dn_ellip.SetLineWidth(2);

  TF1 f_line("f_line","theline(x)",-1.0,3.0/5.0);
  f_line.SetParameter(0,0);
  f_line.SetLineColor(4);
  f_line.SetLineWidth(2);

  f_up_ellip.Draw("same");
  f_dn_ellip.Draw("same");
  f_line.Draw("same");

  TLine* l_ref1 = new TLine(h_polyanddecay_f_and_2der_pos->GetXaxis()->GetXmin(),-0.5,
			    h_polyanddecay_f_and_2der_pos->GetXaxis()->GetXmax(),-0.5);
  l_ref1->SetLineColor(1);
  l_ref1->SetLineWidth(2);
  l_ref1->SetLineStyle(2);
  l_ref1->Draw();

  TLine* l_ref2 = new TLine(0.0,h_polyanddecay_f_and_2der_pos->GetYaxis()->GetXmin(),
			   0.0,h_polyanddecay_f_and_2der_pos->GetYaxis()->GetXmax());
  l_ref2->SetLineColor(1);
  l_ref2->SetLineWidth(2);
  l_ref2->SetLineStyle(2);
  l_ref2->Draw();
  c1->Print(EPSName.Data());

  double wh = (R_a2[1]-R_a2[0])/80.0;
  double wv = (R_a3[1]-R_a3[0])/80.0;

  //The intra-ellipse points:
  for(int ih=0;ih<Npoints_h;ih++) {
    for(int iv=0;iv<Npoints_v;iv++) {
      c1->Clear();
      c1->Divide(2,2);
      c1->cd(1);
      gPad->SetFillColor(10);
      gPad->SetFrameFillColor(10);
      gPad->SetTickx(1);
      gPad->SetTicky(1);
      gPad->SetBottomMargin(0.12);
      gPad->SetLeftMargin(0.12); 
      h_polyanddecay_f_and_2der_pos->Draw("col");
      f_up_ellip.Draw("same");
      f_dn_ellip.Draw("same");
      f_line.Draw("same");
      l_ref1->Draw();
      l_ref2->Draw();

      TEllipse p_tmp(a2_points[ih],a3_points_int_ellipse[ih][iv],wh,wv);
      p_tmp.SetLineColor(1);
      p_tmp.SetLineWidth(2);
      p_tmp.SetFillColor(1);
      p_tmp.Draw("same");
      c1->cd(2);
      gPad->SetFillColor(10);
      gPad->SetFrameFillColor(10);
      gPad->SetTickx(1);
      gPad->SetTicky(1);
      gPad->SetBottomMargin(0.12);
      gPad->SetLeftMargin(0.12);
      h_polyanddecay_intraellipse[ih][iv]->Draw("l");
      c1->Print(EPSName.Data());
    }
  }

  //The out-ellipse good points:
  for(int ih=0;ih<Npoints_h;ih++) {
    for(int iv=0;iv<Npoints_v;iv++) {
      c1->Clear();
      c1->Divide(2,2);
      c1->cd(1);
      gPad->SetFillColor(10);
      gPad->SetFrameFillColor(10);
      gPad->SetTickx(1);
      gPad->SetTicky(1);
      gPad->SetBottomMargin(0.12);
      gPad->SetLeftMargin(0.12); 
      h_polyanddecay_f_and_2der_pos->Draw("col");
      f_up_ellip.Draw("same");
      f_dn_ellip.Draw("same");
      f_line.Draw("same");
      l_ref1->Draw();
      l_ref2->Draw();
      TEllipse p_tmp(a2_points_v2[ih],a3_points_out_ellipse_good[ih][iv],wh,wv);
      p_tmp.SetLineColor(1);
      p_tmp.SetLineWidth(2);
      p_tmp.SetFillColor(1);
      p_tmp.Draw("same");
      c1->cd(2);
      gPad->SetFillColor(10);
      gPad->SetFrameFillColor(10);
      gPad->SetTickx(1);
      gPad->SetTicky(1);
      gPad->SetBottomMargin(0.12);
      gPad->SetLeftMargin(0.12);
      h_polyanddecay_outellipse_good[ih][iv]->Draw("l");
      c1->Print(EPSName.Data());
    }
  }

  //The out-ellipse bad up points:
  for(int ih=0;ih<Npoints_h;ih++) {
    for(int iv=0;iv<Npoints_v;iv++) {
      c1->Clear();
      c1->Divide(2,2);
      c1->cd(1);
      gPad->SetFillColor(10);
      gPad->SetFrameFillColor(10);
      gPad->SetTickx(1);
      gPad->SetTicky(1);
      gPad->SetBottomMargin(0.12);
      gPad->SetLeftMargin(0.12); 
      h_polyanddecay_f_and_2der_pos->Draw("col");
      f_up_ellip.Draw("same");
      f_dn_ellip.Draw("same");
      f_line.Draw("same");
      l_ref1->Draw();
      l_ref2->Draw();
      TEllipse p_tmp(a2_points[ih],a3_points_out_ellipse_bad_up[ih][iv],wh,wv);
      p_tmp.SetLineColor(1);
      p_tmp.SetLineWidth(2);
      p_tmp.SetFillColor(1);
      p_tmp.Draw("same");
      c1->cd(2);
      gPad->SetFillColor(10);
      gPad->SetFrameFillColor(10);
      gPad->SetTickx(1);
      gPad->SetTicky(1);
      gPad->SetBottomMargin(0.12);
      gPad->SetLeftMargin(0.12);
      h_polyanddecay_outellipse_bad_up[ih][iv]->Draw("l");
      c1->Print(EPSName.Data());
    }
  }

  //The out-ellipse bad dn points:
  for(int ih=0;ih<Npoints_h;ih++) {
    for(int iv=0;iv<Npoints_v;iv++) {
      c1->Clear();
      c1->Divide(2,2);
      c1->cd(1);
      gPad->SetFillColor(10);
      gPad->SetFrameFillColor(10);
      gPad->SetTickx(1);
      gPad->SetTicky(1);
      gPad->SetBottomMargin(0.12);
      gPad->SetLeftMargin(0.12); 
      h_polyanddecay_f_and_2der_pos->Draw("col");
      f_up_ellip.Draw("same");
      f_dn_ellip.Draw("same");
      f_line.Draw("same");
      l_ref1->Draw();
      l_ref2->Draw();
      TEllipse p_tmp(a2_points[ih],a3_points_out_ellipse_bad_dn[ih][iv],wh,wv);
      p_tmp.SetLineColor(1);
      p_tmp.SetLineWidth(2);
      p_tmp.SetFillColor(1);
      p_tmp.Draw("same");
      c1->cd(2);
      gPad->SetFillColor(10);
      gPad->SetFrameFillColor(10);
      gPad->SetTickx(1);
      gPad->SetTicky(1);
      gPad->SetBottomMargin(0.12);
      gPad->SetLeftMargin(0.12);
      h_polyanddecay_outellipse_bad_dn[ih][iv]->Draw("l");
      c1->Print(EPSName.Data());
    }
  }

  c1->Print(EPSNameC.Data());

  return;

}
//==================================================================================
double polyanddecay(double m,
		    double a2,
		    double a3,
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
double polyanddecay_integral(double a2,
			     double a3,
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



