/*****************************************************************************
 * Project: BaBar detector at the SLAC PEP-II B-factory
 * Package: RooDKDalitz
 *    File: $Id: DmixJohnsonSU.rdl,v 1.1 2010/12/19 13:58:19 gcasa Exp $
 *  See .cc file
 *****************************************************************************/

#ifndef ROODOUBLECBSHAPE
#define ROODOUBLECBSHAPE


#include "RooAbsPdf.h"
#include "RooRealProxy.h"
#include "RooCategoryProxy.h"
#include "RooAbsReal.h"
#include "RooAbsCategory.h"
#include "RooAbsGenContext.h"

class RooDoubleCBShape : public RooAbsPdf {

protected:

  RooRealProxy m;
  RooRealProxy m0;
  RooRealProxy sL;
  RooRealProxy sR;
  RooRealProxy aL;
  RooRealProxy aR;
  RooRealProxy oL;
  RooRealProxy oR;
  
  Double_t evaluate() const;

public:
  RooDoubleCBShape() ;
  RooDoubleCBShape(const char *name, const char *title,
		   RooAbsReal& _m,
		   RooAbsReal& _m0,
		   RooAbsReal& _sL,
		   RooAbsReal& _sR,
		   RooAbsReal& _aL,
		   RooAbsReal& _aR,
		   RooAbsReal& _oL,
		   RooAbsReal& _oR);
  RooDoubleCBShape(const RooDoubleCBShape& other, const char* name=0);
  virtual TObject* clone(const char* newname) const { return new RooDoubleCBShape(*this,newname); }
  inline virtual ~RooDoubleCBShape() { }

  virtual Int_t    getAnalyticalIntegral(RooArgSet& allVars, RooArgSet& analVars, const char* rangeName=0) const; 
  virtual Double_t analyticalIntegral(Int_t code, const char* rangeName=0) const; 
  

  //ClassDef(RooDoubleCBShape,0)

    };

#endif

