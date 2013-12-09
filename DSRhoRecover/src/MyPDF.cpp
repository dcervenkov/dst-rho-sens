#include "Riostream.h"

#include "RooAbsReal.h"
#include "RooAbsCategory.h"
#include <math.h>
#include "TMath.h"

#include "MyPDF.h"
#include "Constants.h"

MyPDF::MyPDF(const char* name, const char* title,
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
             RooAbsReal& _phiw) :
    RooAbsPdf(name,title),
    xp("xp","xp",this,_xp),
    x0("x0","x0",this,_x0),
    xt("xt","xt",this,_xt),
    yp("yp","yp",this,_yp),
    y0("y0","y0",this,_y0),
    yt("yt","yt",this,_yt),
    xbp("xbp","xbp",this,_xbp),
    xb0("xb0","xb0",this,_xb0),
    xbt("xbt","xbt",this,_xbt),
    ybp("ybp","ybp",this,_ybp),
    yb0("yb0","yb0",this,_yb0),
    ybt("ybt","ybt",this,_ybt),
    xpe("xpe","xpe",this,_xpe),
    x0e("x0e","x0e",this,_x0e),
    xte("xte","xte",this,_xte),
    ype("ype","ype",this,_ype),
    y0e("y0e","y0e",this,_y0e),
    yte("yte","yte",this,_yte),
    xbpe("xbpe","xbpe",this,_xbpe),
    xb0e("xb0e","xb0e",this,_xb0e),
    xbte("xbte","xbte",this,_xbte),
    ybpe("ybpe","ybpe",this,_ybpe),
    yb0e("yb0e","yb0e",this,_yb0e),
    ybte("ybte","ybte",this,_ybte),
    rp("rp","rp",this,_rp),
    r0("r0","r0",this,_r0),
    rt("rt","rt",this,_rt),
    sp("sp","sp",this,_sp),
    s0("s0","s0",this,_s0),
    st("st","st",this,_st),
    phiw("phiw","phiw",this,_phiw) {
}


MyPDF::MyPDF(const MyPDF& other, const char* name) :
    RooAbsPdf(other,name),
    xp("xp",this,other.xp),
    x0("x0",this,other.x0),
    xt("xt",this,other.xt),
    yp("yp",this,other.yp),
    y0("y0",this,other.y0),
    yt("yt",this,other.yt),
    xbp("xbp",this,other.xbp),
    xb0("xb0",this,other.xb0),
    xbt("xbt",this,other.xbt),
    ybp("ybp",this,other.ybp),
    yb0("yb0",this,other.yb0),
    ybt("ybt",this,other.ybt),
    xpe("xpe",this,other.xpe),
    x0e("x0e",this,other.x0e),
    xte("xte",this,other.xte),
    ype("ype",this,other.ype),
    y0e("y0e",this,other.y0e),
    yte("yte",this,other.yte),
    xbpe("xbpe",this,other.xbpe),
    xb0e("xb0e",this,other.xb0e),
    xbte("xbte",this,other.xbte),
    ybpe("ybpe",this,other.ybpe),
    yb0e("yb0e",this,other.yb0e),
    ybte("ybte",this,other.ybte),
    rp("rp",this,other.rp),
    r0("r0",this,other.r0),
    rt("rt",this,other.rt),
    sp("sp",this,other.sp),
    s0("s0",this,other.s0),
    st("st",this,other.st),
    phiw("phiw",this,other.phiw) {
}



Double_t MyPDF::evaluate() const {
    Double_t mxp = rp*cos(-phiw+sp);
    Double_t mx0 = r0*cos(-phiw+s0);
    Double_t mxt = rt*cos(-phiw+st);
    Double_t myp = rp*sin(-phiw+sp);
    Double_t my0 = r0*sin(-phiw+s0);
    Double_t myt = rt*sin(-phiw+st);

    Double_t mxbp = rp*cos(+phiw+sp);
    Double_t mxb0 = r0*cos(+phiw+s0);
    Double_t mxbt = rt*cos(+phiw+st);
    Double_t mybp = rp*sin(+phiw+sp);
    Double_t myb0 = r0*sin(+phiw+s0);
    Double_t mybt = rt*sin(+phiw+st);

    Double_t gp = 1./(2*PI*xpe*ype)*exp(-0.5*((xp-mxp)*(xp-mxp)/(xpe*xpe) + (yp-myp)*(yp-myp)/(ype*ype)));
    Double_t g0 = 1./(2*PI*x0e*y0e)*exp(-0.5*((x0-mx0)*(x0-mx0)/(x0e*x0e) + (y0-my0)*(y0-my0)/(y0e*y0e)));
    Double_t gt = 1./(2*PI*xte*yte)*exp(-0.5*((xt-mxt)*(xt-mxt)/(xte*xte) + (yt-myt)*(yt-myt)/(yte*yte)));

    Double_t gbp = 1./(2*PI*xbpe*ybpe)*exp(-0.5*((xbp-mxbp)*(xbp-mxbp)/(xbpe*xbpe) + (ybp-mybp)*(ybp-mybp)/(ybpe*ybpe)));
    Double_t gb0 = 1./(2*PI*xb0e*yb0e)*exp(-0.5*((xb0-mxb0)*(xb0-mxb0)/(xb0e*xb0e) + (yb0-myb0)*(yb0-myb0)/(yb0e*yb0e)));
    Double_t gbt = 1./(2*PI*xbte*ybte)*exp(-0.5*((xbt-mxbt)*(xbt-mxbt)/(xbte*xbte) + (ybt-mybt)*(ybt-mybt)/(ybte*ybte)));

    Double_t value = gp*g0*gt*gbp*gb0*gbt;

//    printf("returning val = %f\n",value);
//    printf("gs: %f %f %f %f %f %f\n",gp,g0,gt,gbp,gb0,gbt);
//    printf("xs: %f %f %f %f %f %f\n",(double)xp,(double)x0,(double)xt,(double)xbp,(double)xb0,(double)xbt);
//    printf("mxs: %f %f %f %f %f %f\n",(double)mxp,(double)mx0,(double)mxt,(double)mxbp,(double)mxb0,(double)mxbt);
//    printf("xes: %f %f %f %f %f %f\n",(double)xpe,(double)x0e,(double)xte,(double)xbpe,(double)xb0e,(double)xbte);
//    printf("term: %f\n",-0.5*((xb0-mxb0)*(xb0-mxb0)/(xb0e*xb0e) + (yb0-myb0)*(yb0-myb0)/(yb0e*yb0e)));

    return value;
}



Int_t MyPDF::getAnalyticalIntegral(RooArgSet& allVars, RooArgSet& analVars, const char* /*rangeName*/) const {
    // LIST HERE OVER WHICH VARIABLES ANALYTICAL INTEGRATION IS SUPPORTED,
    // ASSIGN A NUMERIC CODE FOR EACH SUPPORTED (SET OF) PARAMETERS
    // THE EXAMPLE BELOW ASSIGNS CODE 1 TO INTEGRATION OVER VARIABLE X
    // YOU CAN ALSO IMPLEMENT MORE THAN ONE ANALYTICAL INTEGRAL BY REPEATING THE matchArgs
    // EXPRESSION MULTIPLE TIMES

    // if (matchArgs(allVars,analVars,x)) return 1 ;
//    RooArgSet intSet(static_cast<RooAbsArg>(xp),static_cast<RooAbsArg>(x0),static_cast<RooAbsArg>(xt),\
//                     static_cast<RooAbsArg>(yp),static_cast<RooAbsArg>(y0),static_cast<RooAbsArg>(yt));
    RooArgSet intSet(xp.arg(),x0.arg(),xt.arg(),yp.arg(),y0.arg(),yt.arg());
    intSet.add(xbp.arg());
    intSet.add(xb0.arg());
    intSet.add(xbt.arg());
    intSet.add(ybp.arg());
    intSet.add(yb0.arg());
    intSet.add(ybt.arg());

    intSet.add(xpe.arg());
    intSet.add(x0e.arg());
    intSet.add(xte.arg());
    intSet.add(ype.arg());
    intSet.add(y0e.arg());
    intSet.add(yte.arg());

    intSet.add(xbpe.arg());
    intSet.add(xb0e.arg());
    intSet.add(xbte.arg());
    intSet.add(ybpe.arg());
    intSet.add(yb0e.arg());
    intSet.add(ybte.arg());

    if (matchArgs(allVars,analVars,intSet)) return 1 ;

    return 0;
}



Double_t MyPDF::analyticalIntegral(Int_t code, const char* rangeName) const {
    // RETURN ANALYTICAL INTEGRAL DEFINED BY RETURN CODE ASSIGNED BY getAnalyticalIntegral
    // THE MEMBER FUNCTION x.min(rangeName) AND x.max(rangeName) WILL RETURN THE INTEGRATION
    // BOUNDARIES FOR EACH OBSERVABLE x
    // assert(code==1) ;
    // return (x.max(rangeName)-x.min(rangeName)) ;
    if(code == 1) return 1;
}

