/*****************************************************************************
 * Project: BaBar detector at the SLAC PEP-II B-factory
 * Package: RooDKDalitz
 *    File: $Id: DmixJohnsonSU.cc,v 1.1 2010/12/19 13:58:19 gcasa Exp $
 * Authors:
 *   Fernando Martinez-Vidal martinef@slac.stanford.edu, IFIC-Valencia
 * History:
 *     - Apr 20: creation
 * Copyright (C) 2007 IFIC-Valencia
 *****************************************************************************/
// rmaria -> modified for the D0 physics

// -- CLASS DESCRIPTION [PDF] --

//#include "RooFit.h"

#include "Riostream.h"
#include "RooAbsReal.h" 
#include "RooAbsCategory.h" 
#include "RooMath.h"
#include "TMath.h"
#include "RooDoubleCBShape.hh"


//ClassImp(RooDoubleCBShape)

RooDoubleCBShape::RooDoubleCBShape(){} //the default constructor

RooDoubleCBShape::RooDoubleCBShape(const char *name, const char *title,
				   RooAbsReal& _m,
				   RooAbsReal& _m0,
				   RooAbsReal& _sL,
				   RooAbsReal& _sR,
				   RooAbsReal& _aL,
				   RooAbsReal& _aR,
				   RooAbsReal& _oL,
				   RooAbsReal& _oR):
  RooAbsPdf(name, title), 
  m("m",  "Dependent",this, _m),
  m0("m0","m0",       this, _m0),
  sL("sL","sL",       this, _sL),
  sR("sR","sR",       this, _sR),
  aL("aL","aL",       this, _aL),
  aR("aR","aR",       this, _aR),
  oL("oL","oL",       this, _oL),
  oR("oR","oR",       this, _oR)
{
}

RooDoubleCBShape::RooDoubleCBShape(const RooDoubleCBShape& other, const char* name) :
  RooAbsPdf(other, name),
  m ("m", this,other.m), 
  m0("m0",this,other.m0),
  sL("sL",this,other.sL),
  sR("sR",this,other.sR),
  aL("aL",this,other.aL),
  aR("aR",this,other.aR),
  oL("oL",this,other.oL),
  oR("oR",this,other.oR)
{
}

Double_t RooDoubleCBShape::evaluate() const
{

  if(!(aL >= 0 && aR >= 0 && sL >= 0 && sR >= 0 && oL >=0 && oR >= 0)) {
    cout << endl;
    cout << "RooDoubleCBShape::aL, aR, sL, sR, oL and oR parameters need to be positive." << endl;
    cout << "RooDoubleCBShape::aL = " << aL << endl;
    cout << "RooDoubleCBShape::aR = " << aR << endl;
    cout << "RooDoubleCBShape::sL = " << sL << endl;
    cout << "RooDoubleCBShape::sR = " << sR << endl;
    cout << "RooDoubleCBShape::oL = " << oL << endl;
    cout << "RooDoubleCBShape::oR = " << oR << endl;
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


Int_t RooDoubleCBShape::getAnalyticalIntegral(RooArgSet& allVars, RooArgSet& analVars, const char* /*rangeName*/) const
{

  if(matchArgs(allVars,analVars,m)) return 1;

  return 0;
}

Double_t RooDoubleCBShape::analyticalIntegral(Int_t code, const char* rangeName) const
{

  if(!(aL >= 0 && aR >= 0 && sL >= 0 && sR >= 0)) {
    cout << endl;
    cout << "RooDoubleCBShape::aL,aR, oL and oR parameters need to be positive." << endl;
    cout << "RooDoubleCBShape::Exiting now!!!" << endl;
    cout << endl;
    assert(false);
  }


  double m_min = m.min(rangeName);
  double m_max = m.max(rangeName);

  if(m_min > m_max) {
    cout << endl;
    cout << "RooDoubleCBShape::m_min is higher then m_max." << endl;
    cout << "RooDoubleCBShape::Exiting now!!!" << endl;
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

