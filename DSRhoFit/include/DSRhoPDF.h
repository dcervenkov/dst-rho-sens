/*****************************************************************************
 * Project: RooFit                                                           *
 *                                                                           *
  * This code was autogenerated by RooClassFactory                            *
 *****************************************************************************/

#ifndef DSRHOPDF
#define DSRHOPDF

#include "RooAbsPdf.h"
#include "RooRealProxy.h"
#include "RooCategoryProxy.h"
#include "RooAbsReal.h"
#include "RooAbsCategory.h"

class DSRhoPDF : public RooAbsPdf
{
public:
    DSRhoPDF() {} ;
    DSRhoPDF(const char *name, const char *title, const char* _type,
             RooAbsReal& _tht,
             RooAbsReal& _thb,
             RooAbsReal& _phit,
             RooAbsReal& _dt,
             RooAbsReal& _ap,
             RooAbsReal& _apa,
             RooAbsReal& _a0,
             RooAbsReal& _ata,
             RooAbsReal& _phiw,
             RooAbsReal& _rp,
             RooAbsReal& _r0,
             RooAbsReal& _rt,
             RooAbsReal& _sp,
             RooAbsReal& _s0,
             RooAbsReal& _st);
    DSRhoPDF(const DSRhoPDF& other, const char* name=0) ;
    virtual TObject* clone(const char* newname) const
    {
        return new DSRhoPDF(*this,newname);
    }
    inline virtual ~DSRhoPDF() { }

    Int_t getAnalyticalIntegral(RooArgSet& allVars, RooArgSet& analVars, const char* rangeName=0) const ;
    Double_t analyticalIntegral(Int_t code, const char* rangeName=0) const ;
    void setType(Int_t new_type) {type = new_type;} ;
    Int_t getType() const {return type;};

protected:

    RooRealProxy tht ;
    RooRealProxy thb ;
    RooRealProxy phit ;
    RooRealProxy dt ;
    RooRealProxy ap ;
    RooRealProxy apa ;
    RooRealProxy a0 ;
    RooRealProxy ata ;
    RooRealProxy phiw ;
    RooRealProxy rp ;
    RooRealProxy r0 ;
    RooRealProxy rt ;
    RooRealProxy sp ;
    RooRealProxy s0 ;
    RooRealProxy st ;

    Int_t type;

    Double_t evaluate() const ;

private:
    //ClassDef(DSRhoPDF,1) // Your description goes here...
};

#endif
