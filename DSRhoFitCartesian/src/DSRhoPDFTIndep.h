/*****************************************************************************
 * Project: RooFit                                                           *
 *                                                                           *
  * This code was autogenerated by RooClassFactory                            *
 *****************************************************************************/

#ifndef DSRHOPDFTINDEP
#define DSRHOPDFTINDEP

#include "RooAbsPdf.h"
#include "RooRealProxy.h"
#include "RooCategoryProxy.h"
#include "RooAbsReal.h"
#include "RooAbsCategory.h"

class DSRhoPDFTIndep : public RooAbsPdf
{
public:
    DSRhoPDFTIndep() {} ;
    DSRhoPDFTIndep(const char *name, const char *title,
             RooAbsReal& _tht,
             RooAbsReal& _thb,
             RooAbsReal& _phit,
             RooAbsReal& _ap,
             RooAbsReal& _apa,
             RooAbsReal& _a0,
             RooAbsReal& _ata);
    DSRhoPDFTIndep(const DSRhoPDFTIndep& other, const char* name=0) ;
    virtual TObject* clone(const char* newname) const
    {
        return new DSRhoPDFTIndep(*this,newname);
    }
    inline virtual ~DSRhoPDFTIndep() { }

    Int_t getAnalyticalIntegral(RooArgSet& allVars, RooArgSet& analVars, const char* rangeName=0) const ;
    Double_t analyticalIntegral(Int_t code, const char* rangeName=0) const ;

protected:

    RooRealProxy tht ;
    RooRealProxy thb ;
    RooRealProxy phit ;
    RooRealProxy ap ;
    RooRealProxy apa ;
    RooRealProxy a0 ;
    RooRealProxy ata ;

    Double_t evaluate() const ;

private:
    //ClassDef(DSRhoPDF,1) // Your description goes here...
};

#endif
