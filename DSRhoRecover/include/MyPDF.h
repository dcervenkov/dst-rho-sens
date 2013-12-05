#ifndef MYPDF
#define MYPDF

#include "RooAbsPdf.h"
#include "RooRealProxy.h"
#include "RooCategoryProxy.h"
#include "RooAbsReal.h"
#include "RooAbsCategory.h"

class MyPDF : public RooAbsPdf {
public:
    MyPDF() {} ;
    MyPDF(const char* name, const char* title,
          RooAbsReal& _xp,
          RooAbsReal& _x0,
          RooAbsReal& _xt,
          RooAbsReal& _yp,
          RooAbsReal& _y0,
          RooAbsReal& _yt,
          RooAbsReal& _xbp,
          RooAbsReal& _xb0,
          RooAbsReal& _xbt,
          RooAbsReal& _ybp,
          RooAbsReal& _yb0,
          RooAbsReal& _ybt,
          RooAbsReal& _xpe,
          RooAbsReal& _x0e,
          RooAbsReal& _xte,
          RooAbsReal& _ype,
          RooAbsReal& _y0e,
          RooAbsReal& _yte,
          RooAbsReal& _xbpe,
          RooAbsReal& _xb0e,
          RooAbsReal& _xbte,
          RooAbsReal& _ybpe,
          RooAbsReal& _yb0e,
          RooAbsReal& _ybte,
          RooAbsReal& _rp,
          RooAbsReal& _r0,
          RooAbsReal& _rt,
          RooAbsReal& _sp,
          RooAbsReal& _s0,
          RooAbsReal& _st,
          RooAbsReal& _phiw
         );
    MyPDF(const MyPDF& other, const char* name=0) ;
    virtual TObject* clone(const char* newname) const {
        return new MyPDF(*this,newname);
    }
    inline virtual ~MyPDF() { }

    Int_t getAnalyticalIntegral(RooArgSet& allVars, RooArgSet& analVars, const char* rangeName=0) const ;
    Double_t analyticalIntegral(Int_t code, const char* rangeName=0) const ;

protected:

    RooRealProxy xp ;
    RooRealProxy x0 ;
    RooRealProxy xt ;
    RooRealProxy yp ;
    RooRealProxy y0 ;
    RooRealProxy yt ;
    RooRealProxy xbp ;
    RooRealProxy xb0 ;
    RooRealProxy xbt ;
    RooRealProxy ybp ;
    RooRealProxy yb0 ;
    RooRealProxy ybt ;

    RooRealProxy xpe;
    RooRealProxy x0e;
    RooRealProxy xte;
    RooRealProxy ype;
    RooRealProxy y0e;
    RooRealProxy yte;
    RooRealProxy xbpe;
    RooRealProxy xb0e;
    RooRealProxy xbte;
    RooRealProxy ybpe;
    RooRealProxy yb0e;
    RooRealProxy ybte;

    RooRealProxy rp;
    RooRealProxy r0;
    RooRealProxy rt;

    RooRealProxy sp;
    RooRealProxy s0;
    RooRealProxy st;

    RooRealProxy phiw;

    Double_t evaluate() const ;

};

#endif
