/*****************************************************************************
 * Project: RooFit                                                           *
 *                                                                           *
 * This code was autogenerated by RooClassFactory                            *
 *****************************************************************************/

// Your description goes here...

#include "Riostream.h"

#include "RooAbsReal.h"
#include "RooAbsCategory.h"
#include <math.h>
#include "TMath.h"

#include "DSRhoPDF.h"
#include "Constants.h"

//ClassImp(DSRhoPDF)

DSRhoPDF::DSRhoPDF(const char *name, const char *title, const char *_type,
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
                   RooAbsReal& _st) :
    RooAbsPdf(name,title),
    tht("tht","tht",this,_tht),
    thb("thb","thb",this,_thb),
    phit("phit","phit",this,_phit),
    dt("dt","dt",this,_dt),
    ap("ap","ap",this,_ap),
    apa("apa","apa",this,_apa),
    a0("a0","a0",this,_a0),
    ata("ata","ata",this,_ata),
    phiw("phiw","phiw",this,_phiw),
    rp("rp","rp",this,_rp),
    r0("r0","r0",this,_r0),
    rt("rt","rt",this,_rt),
    sp("sp","sp",this,_sp),
    s0("s0","s0",this,_s0),
    st("st","st",this,_st)
{
    if     (strcmp(_type,"a") == 0)
        type = 1;
    else if(strcmp(_type,"ab") == 0)
        type = 2;
    else if(strcmp(_type,"b") == 0)
        type = 3;
    else if(strcmp(_type,"bb") == 0)
        type = 4;
    else
        printf("ERROR: unknown _type: %s\n",_type);

}


DSRhoPDF::DSRhoPDF(const DSRhoPDF& other, const char* name) :
    RooAbsPdf(other,name),
    tht("tht",this,other.tht),
    thb("thb",this,other.thb),
    phit("phit",this,other.phit),
    dt("dt",this,other.dt),
    ap("ap",this,other.ap),
    apa("apa",this,other.apa),
    a0("a0",this,other.a0),
    ata("ata",this,other.ata),
    phiw("phiw",this,other.phiw),
    rp("rp",this,other.rp),
    r0("r0",this,other.r0),
    rt("rt",this,other.rt),
    sp("sp",this,other.sp),
    s0("s0",this,other.s0),
    st("st",this,other.st)
{
    type = other.type;
}



Double_t DSRhoPDF::evaluate() const
{
    Int_t phiw_sign = 1;

    Double_t a0a = 0;
    Double_t at = sqrt(1-ap*ap-a0*a0);

    Double_t ap0r = ap*a0*cos(-apa+a0a);
    Double_t ap0i = ap*a0*sin(-apa+a0a);
    Double_t a0tr = a0*at*cos(-a0a+ata);
    Double_t a0ti = a0*at*sin(-a0a+ata);
    Double_t aptr = ap*at*cos(-apa+ata);
    Double_t apti = ap*at*sin(-apa+ata);

    Double_t At2 = 0;
    Double_t Ap2 = 0;
    Double_t A02 = 0;
    Double_t Ap0r = 0;
    Double_t A0ti = 0;
    Double_t Apti = 0;

    if(type < 1 || type > 4)
    {
        printf("ERROR: undefined type: %i",type);
        return 0;
    }

    if(type == 2 || type == 3)
        phiw_sign = -1;


    if(type == 1 || type == 2)
    {
        At2 = at*at*((1+rt*rt)+(1-rt*rt)*cos(Bfreq*dt)+2*rt*sin(phiw_sign*phiw-st)*sin(Bfreq*dt));
        Ap2 = ap*ap*((1+rp*rp)+(1-rp*rp)*cos(Bfreq*dt)-2*rp*sin(phiw_sign*phiw-sp)*sin(Bfreq*dt));
        A02 = a0*a0*((1+r0*r0)+(1-r0*r0)*cos(Bfreq*dt)-2*r0*sin(phiw_sign*phiw-s0)*sin(Bfreq*dt));

        Ap0r =    ap0r*(1+rp*r0*cos(sp-s0))+ap0i*rp*r0*sin(sp-s0)+\
                          (ap0r*(1-rp*r0*cos(sp-s0))-ap0i*rp*r0*sin(sp-s0))*cos(Bfreq*dt)-\
                          (ap0r*(rp*sin(phiw_sign*phiw-sp)+r0*sin(phiw_sign*phiw-s0))+\
                           ap0i*(rp*cos(phiw_sign*phiw-sp)-r0*cos(phiw_sign*phiw-s0)))*sin(Bfreq*dt);

        A0ti =    a0ti*(1-r0*rt*cos(s0-st))+a0tr*r0*rt*sin(s0-st)+\
                      (a0ti*(1+r0*rt*cos(s0-st))-a0tr*r0*rt*sin(s0-st))*cos(Bfreq*dt)-\
                      (a0ti*(r0*sin(phiw_sign*phiw-s0)-rt*sin(phiw_sign*phiw-st))-\
                       a0tr*(r0*cos(phiw_sign*phiw-s0)+rt*cos(phiw_sign*phiw-st)))*sin(Bfreq*dt);

        Apti =    apti*(1-rp*rt*cos(sp-st))+aptr*rp*rt*sin(sp-st)+\
                          (apti*(1+rp*rt*cos(sp-st))-aptr*rp*rt*sin(sp-st))*cos(Bfreq*dt)-\
                          (apti*(rp*sin(phiw_sign*phiw-sp)-rt*sin(phiw_sign*phiw-st))-\
                           aptr*(rp*cos(phiw_sign*phiw-sp)+rt*cos(phiw_sign*phiw-st)))*sin(Bfreq*dt);
    }

    /// Writing this again explicitly with the changed sign of sin(Bfreq*dt) and cos(Bfreq*dt) is safer than changing
    /// Bfreq or dt to exploit sin(Bfreq*dt + PI) = -sin(Bfreq*dt) because of numerical problems that can cause.
    if(type == 3 || type == 4)
    {
        At2 = at*at*((1+rt*rt)+(1-rt*rt)*(-1)*cos(Bfreq*dt)+2*rt*sin(phiw_sign*phiw-st)*(-1)*sin(Bfreq*dt));
        Ap2 = ap*ap*((1+rp*rp)+(1-rp*rp)*(-1)*cos(Bfreq*dt)-2*rp*sin(phiw_sign*phiw-sp)*(-1)*sin(Bfreq*dt));
        A02 = a0*a0*((1+r0*r0)+(1-r0*r0)*(-1)*cos(Bfreq*dt)-2*r0*sin(phiw_sign*phiw-s0)*(-1)*sin(Bfreq*dt));

        Ap0r =    ap0r*(1+rp*r0*cos(sp-s0))+ap0i*rp*r0*sin(sp-s0)+\
                          (ap0r*(1-rp*r0*cos(sp-s0))-ap0i*rp*r0*sin(sp-s0))*(-1)*cos(Bfreq*dt)-\
                          (ap0r*(rp*sin(phiw_sign*phiw-sp)+r0*sin(phiw_sign*phiw-s0))+\
                           ap0i*(rp*cos(phiw_sign*phiw-sp)-r0*cos(phiw_sign*phiw-s0)))*(-1)*sin(Bfreq*dt);

        A0ti =    a0ti*(1-r0*rt*cos(s0-st))+a0tr*r0*rt*sin(s0-st)+\
                      (a0ti*(1+r0*rt*cos(s0-st))-a0tr*r0*rt*sin(s0-st))*(-1)*cos(Bfreq*dt)-\
                      (a0ti*(r0*sin(phiw_sign*phiw-s0)-rt*sin(phiw_sign*phiw-st))-\
                       a0tr*(r0*cos(phiw_sign*phiw-s0)+rt*cos(phiw_sign*phiw-st)))*(-1)*sin(Bfreq*dt);

        Apti =    apti*(1-rp*rt*cos(sp-st))+aptr*rp*rt*sin(sp-st)+\
                          (apti*(1+rp*rt*cos(sp-st))-aptr*rp*rt*sin(sp-st))*(-1)*cos(Bfreq*dt)-\
                          (apti*(rp*sin(phiw_sign*phiw-sp)-rt*sin(phiw_sign*phiw-st))-\
                           aptr*(rp*cos(phiw_sign*phiw-sp)+rt*cos(phiw_sign*phiw-st)))*(-1)*sin(Bfreq*dt);
    }



    Double_t value =    exp(-Bgamma*fabs(dt))*(Ap2*2*sin(tht)*sin(tht)*sin(tht)*sin(thb)*sin(thb)*sin(thb)*sin(phit)*sin(phit)+\
                        At2*2*cos(tht)*cos(tht)*sin(tht)*sin(thb)*sin(thb)*sin(thb)+\
                        A02*4*sin(tht)*sin(tht)*sin(tht)*cos(thb)*cos(thb)*sin(thb)*cos(phit)*cos(phit)+\
                        sqrt(2)*Ap0r*sin(tht)*sin(tht)*sin(tht)*sin(2*thb)*sin(thb)*sin(2*phit)-\
                        sqrt(2)*A0ti*sin(2*tht)*sin(tht)*sin(2*thb)*sin(thb)*cos(phit)-\
                        2*Apti*sin(2*tht)*sin(tht)*sin(thb)*sin(thb)*sin(thb)*sin(phit));

    return value ;
}



Int_t DSRhoPDF::getAnalyticalIntegral(RooArgSet& allVars, RooArgSet& analVars, const char* /*rangeName*/) const
{
    // LIST HERE OVER WHICH VARIABLES ANALYTICAL INTEGRATION IS SUPPORTED,
    // ASSIGN A NUMERIC CODE FOR EACH SUPPORTED (SET OF) PARAMETERS
    // THE EXAMPLE BELOW ASSIGNS CODE 1 TO INTEGRATION OVER VARIABLE X
    // YOU CAN ALSO IMPLEMENT MORE THAN ONE ANALYTICAL INTEGRAL BY REPEATING THE matchArgs
    // EXPRESSION MULTIPLE TIMES

    // if (matchArgs(allVars,analVars,x)) return 1 ;

    if(matchArgs(allVars,analVars,tht,thb,phit,dt)) return 1;
    if(matchArgs(allVars,analVars,tht,thb,phit)) return 2;
    if(matchArgs(allVars,analVars,tht,thb)) return 3;
    if(matchArgs(allVars,analVars,tht,phit)) return 4;
    if(matchArgs(allVars,analVars,thb,phit)) return 5;
    if(matchArgs(allVars,analVars,tht)) return 6;
    if(matchArgs(allVars,analVars,thb)) return 7;
    if(matchArgs(allVars,analVars,phit)) return 8;

    return 0 ;
}



Double_t DSRhoPDF::analyticalIntegral(Int_t code, const char* rangeName) const
{
    // RETURN ANALYTICAL INTEGRAL DEFINED BY RETURN CODE ASSIGNED BY getAnalyticalIntegral
    // THE MEMBER FUNCTION x.min(rangeName) AND x.max(rangeName) WILL RETURN THE INTEGRATION
    // BOUNDARIES FOR EACH OBSERVABLE x

    // assert(code==1) ;
    // return (x.max(rangeName)-x.min(rangeName)) ;

    Int_t phiw_sign = 1;

    Double_t a0a = 0;
    Double_t at = sqrt(1-ap*ap-a0*a0);

    Double_t ap0r = ap*a0*cos(-apa+a0a);
    Double_t ap0i = ap*a0*sin(-apa+a0a);
//    Double_t a0tr = a0*at*cos(-a0a+ata);
//    Double_t a0ti = a0*at*sin(-a0a+ata);
    Double_t aptr = ap*at*cos(-apa+ata);
    Double_t apti = ap*at*sin(-apa+ata);

    Double_t At2 = 0;
    Double_t Ap2 = 0;
    Double_t A02 = 0;
    Double_t Ap0r = 0;
    Double_t Apti = 0;

    if(type < 1 || type > 4)
    {
        printf("ERROR: undefined type: %i",type);
        return 0;
    }

    if(type == 2 || type == 3)
        phiw_sign = -1;

    if(type == 1 || type == 2)
    {
        At2 = at*at*((1+rt*rt)+(1-rt*rt)*cos(Bfreq*dt)+2*rt*sin(phiw_sign*phiw-st)*sin(Bfreq*dt));
        Ap2 = ap*ap*((1+rp*rp)+(1-rp*rp)*cos(Bfreq*dt)-2*rp*sin(phiw_sign*phiw-sp)*sin(Bfreq*dt));
        A02 = a0*a0*((1+r0*r0)+(1-r0*r0)*cos(Bfreq*dt)-2*r0*sin(phiw_sign*phiw-s0)*sin(Bfreq*dt));

        Ap0r =    ap0r*(1+rp*r0*cos(sp-s0))+ap0i*rp*r0*sin(sp-s0)+\
                          (ap0r*(1-rp*r0*cos(sp-s0))-ap0i*rp*r0*sin(sp-s0))*cos(Bfreq*dt)-\
                          (ap0r*(rp*sin(phiw_sign*phiw-sp)+r0*sin(phiw_sign*phiw-s0))+\
                           ap0i*(rp*cos(phiw_sign*phiw-sp)-r0*cos(phiw_sign*phiw-s0)))*sin(Bfreq*dt);

        Apti =    apti*(1-rp*rt*cos(sp-st))+aptr*rp*rt*sin(sp-st)+\
                          (apti*(1+rp*rt*cos(sp-st))-aptr*rp*rt*sin(sp-st))*cos(Bfreq*dt)-\
                          (apti*(rp*sin(phiw_sign*phiw-sp)-rt*sin(phiw_sign*phiw-st))-\
                           aptr*(rp*cos(phiw_sign*phiw-sp)+rt*cos(phiw_sign*phiw-st)))*sin(Bfreq*dt);
    }

    /// Writing this again explicitly with the changed sign of sin(Bfreq*dt) and cos(Bfreq*dt) is safer than changing
    /// Bfreq or dt to exploit sin(Bfreq*dt + PI) = -sin(Bfreq*dt) because of numerical problems that can cause.
    if(type == 3 || type == 4)
    {
        At2 = at*at*((1+rt*rt)+(1-rt*rt)*(-1)*cos(Bfreq*dt)+2*rt*sin(phiw_sign*phiw-st)*(-1)*sin(Bfreq*dt));
        Ap2 = ap*ap*((1+rp*rp)+(1-rp*rp)*(-1)*cos(Bfreq*dt)-2*rp*sin(phiw_sign*phiw-sp)*(-1)*sin(Bfreq*dt));
        A02 = a0*a0*((1+r0*r0)+(1-r0*r0)*(-1)*cos(Bfreq*dt)-2*r0*sin(phiw_sign*phiw-s0)*(-1)*sin(Bfreq*dt));

        Ap0r =    ap0r*(1+rp*r0*cos(sp-s0))+ap0i*rp*r0*sin(sp-s0)+\
                          (ap0r*(1-rp*r0*cos(sp-s0))-ap0i*rp*r0*sin(sp-s0))*(-1)*cos(Bfreq*dt)-\
                          (ap0r*(rp*sin(phiw_sign*phiw-sp)+r0*sin(phiw_sign*phiw-s0))+\
                           ap0i*(rp*cos(phiw_sign*phiw-sp)-r0*cos(phiw_sign*phiw-s0)))*(-1)*sin(Bfreq*dt);

        Apti =    apti*(1-rp*rt*cos(sp-st))+aptr*rp*rt*sin(sp-st)+\
                          (apti*(1+rp*rt*cos(sp-st))-aptr*rp*rt*sin(sp-st))*(-1)*cos(Bfreq*dt)-\
                          (apti*(rp*sin(phiw_sign*phiw-sp)-rt*sin(phiw_sign*phiw-st))-\
                           aptr*(rp*cos(phiw_sign*phiw-sp)+rt*cos(phiw_sign*phiw-st)))*(-1)*sin(Bfreq*dt);
    }


    switch(code)
    {
    case 1: // Int[g,{tht,thb,phit,dt}]
        if(type == 1 || type == 2)
            return 128.*PI/9.*(1/Bgamma+Bgamma/(Bfreq*Bfreq+Bgamma*Bgamma)+(ap*ap*rp*rp+a0*a0*r0*r0+at*at*rt*rt)*(1/Bgamma-Bgamma/(Bfreq*Bfreq+Bgamma*Bgamma)));
        else
            return 128.*PI/9.*(1/Bgamma-Bgamma/(Bfreq*Bfreq+Bgamma*Bgamma)+(ap*ap*rp*rp+a0*a0*r0*r0+at*at*rt*rt)*(1/Bgamma+Bgamma/(Bfreq*Bfreq+Bgamma*Bgamma)));

    case 2: // Int[g,{tht,thb,phit}]
        return 64.*PI/9.*exp(-Bgamma*TMath::Abs(dt))*(Ap2+At2+A02);

    case 3: // Int[g,{tht,thb}]
        return 32./9.*exp(-Bgamma*TMath::Abs(dt))*(At2 + 2*A02*cos(phit)*cos(phit) + 2*Ap2*sin(phit)*sin(phit));

    case 4: // Int[g,{tht,phit}]
        return 16.*PI/3.*exp(-Bgamma*TMath::Abs(dt))*(2*A02*cos(thb)*cos(thb) + (Ap2 + At2)*sin(thb)*sin(thb))*sin(thb);

    case 5: // Int[g,{thb,phit}]
        return 16.*PI/3.*exp(-Bgamma*TMath::Abs(dt))*(2*At2*cos(tht)*cos(tht) + (A02 + Ap2)*sin(tht)*sin(tht))*sin(tht);

    case 6: // Int[g,{tht}]
        return 8./3.*exp(-Bgamma*TMath::Abs(dt))*(4*A02*cos(phit)*cos(phit)*cos(thb)*cos(thb) + \
       (At2 + 2*Ap2*sin(phit)*sin(phit))*sin(thb)*sin(thb) + sqrt(2)*Ap0r*sin(2*phit)*sin(2*thb))*sin(thb);

    case 7: // Int[g,{thb}]
        return 16./3.*exp(-Bgamma*TMath::Abs(dt))*(At2*cos(tht)*cos(tht) + A02*cos(phit)*cos(phit)*sin(tht)*sin(tht) + \
       sin(phit)*(Ap2*sin(phit)*sin(tht)*sin(tht) - Apti*sin(2*tht)))*sin(tht);

    case 8: // Int[g,{phit}]
        return 4*PI*exp(-Bgamma*TMath::Abs(dt))*(2*At2*cos(tht)*cos(tht)*sin(thb)*sin(thb) +
       (2*A02*cos(thb)*cos(thb) + Ap2*sin(thb)*sin(thb))*sin(tht)*sin(tht))*(sin(tht)*sin(thb));

    default:
        return 0;
    }

}

