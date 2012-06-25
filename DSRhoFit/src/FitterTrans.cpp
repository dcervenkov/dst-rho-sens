#include "RooRealVar.h"
#include "RooFormulaVar.h"
#include "RooCategory.h"
#include "RooGenericPdf.h"
#include "Constants.h"
#include "RooDataSet.h"
#include "RooRandom.h"
#include "RooChi2Var.h"
#include "RooSimultaneous.h"
#include "Minuit2/Minuit2Minimizer.h"
#include "TPluginManager.h"
#include "TMath.h"
#include "TIterator.h"

#include "DSRhoPDF.h"
#include "FitterTrans.h"

#include "TH1D.h"
#include "TH2F.h"
#include "TCanvas.h"
#include "TFile.h"

FitterTrans::FitterTrans(RooDataSet* outer_dataSet, Double_t* outer_par_input)
{
    gPluginMgr = new TPluginManager;
    gPluginMgr->AddHandler("ROOT::Math::Minimizer", "Minuit2", "Minuit2Minimizer", "Minuit2", "Minuit2Minimizer(const char *)");

    dataSet = outer_dataSet;
    for(int i = 0; i < 11; i++)
        par_input[i] = outer_par_input[i];

    chi2Var = 0;
    result = 0;

    tht_bins = 30;
    thb_bins = 30;
    phit_bins = 30;
    dt_bins = 40;

    thb = new RooRealVar("thb","thb",0,PI);
    tht = new RooRealVar("tht","tht",0,PI);
    phit = new RooRealVar("phit","phit",-PI,PI);
    dt = new RooRealVar("dt","dt",-6,6);
    decType = new RooCategory("decType","decType");
    decType->defineType("a",1);
    decType->defineType("ab",2);
    decType->defineType("b",3);
    decType->defineType("bb",4);
    gamma = new RooRealVar("gamma","gamma",2.83);

    ap = new RooRealVar("ap","ap",par_input[0],0.1,0.4);
    apa = new RooRealVar("apa","apa",par_input[1],0.3,0.7);
    apr = new RooFormulaVar("apr","ap*cos(apa)",RooArgSet(*ap,*apa));
    api = new RooFormulaVar("api","ap*sin(apa)",RooArgSet(*ap,*apa));
    a0 = new RooRealVar("a0","a0",par_input[2],0.8,1);
    a0a = new RooRealVar("a0a","a0a",0);
    a0r = new RooFormulaVar("a0r","a0*cos(a0a)",RooArgSet(*a0,*a0a));
    a0i = new RooFormulaVar("a0i","a0*sin(a0a)",RooArgSet(*a0,*a0a));
    at = new RooFormulaVar("at","sqrt(1-ap*ap-a0*a0)",RooArgSet(*ap,*a0));
    ata = new RooRealVar("ata","ata",par_input[3],2.7,3);
    atr = new RooFormulaVar("atr","at*cos(ata)",RooArgSet(*at,*ata));
    ati = new RooFormulaVar("ati","at*sin(ata)",RooArgSet(*at,*ata));

    ap0r = new RooFormulaVar("ap0r","ap*a0*cos(-apa+a0a)",RooArgSet(*ap,*apa,*a0,*a0a));
    a0ti = new RooFormulaVar("a0ti","a0*at*sin(-a0a+ata)",RooArgSet(*a0,*a0a,*at,*ata));
    apti = new RooFormulaVar("apti","ap*at*sin(-apa+ata)",RooArgSet(*ap,*apa,*at,*ata));

    /// Time-dep additional vars

    dm = new RooRealVar("dm","dm",0.507e12);
    phiw = new RooRealVar("phiw","phiw",par_input[4],1.7,1.9);

    rp = new RooRealVar("rp","rp",par_input[5],0,0.1);
    r0 = new RooRealVar("r0","r0",par_input[6],0,0.1);
    rt = new RooRealVar("rt","rt",par_input[7],0,0.1); /// eq. (100) in BN419 approximates this

    /// s is strong phase; delta_polarization in BN419
    sp = new RooRealVar("sp","sp",par_input[8],0,2*PI);
    s0 = new RooRealVar("s0","s0",par_input[9],0,2*PI);
    st = new RooRealVar("st","st",par_input[10],0,2*PI);

    ap0i = new RooFormulaVar("ap0i","ap*a0*sin(-apa+a0a)",RooArgSet(*ap,*apa,*a0,*a0a));
    a0tr = new RooFormulaVar("a0tr","a0*at*cos(-a0a+ata)",RooArgSet(*a0,*a0a,*at,*ata));
    aptr = new RooFormulaVar("aptr","ap*at*cos(-apa+ata)",RooArgSet(*ap,*apa,*at,*ata));

    At2_a = new RooFormulaVar("At2_a","at*at*((1+rt*rt)+(1-rt*rt)*cos(dm*dt)+2*rt*sin(phiw-st)*sin(dm*dt))",RooArgSet(*at,*rt,*dm,*dt,*phiw,*st));
    Ap2_a = new RooFormulaVar("Ap2_a","ap*ap*((1+rp*rp)+(1-rp*rp)*cos(dm*dt)-2*rp*sin(phiw-sp)*sin(dm*dt))",RooArgSet(*ap,*rp,*dm,*dt,*phiw,*sp));
    A02_a = new RooFormulaVar("A02_a","a0*a0*((1+r0*r0)+(1-r0*r0)*cos(dm*dt)-2*r0*sin(phiw-s0)*sin(dm*dt))",RooArgSet(*a0,*r0,*dm,*dt,*phiw,*s0));

    Ap0r_a = new RooFormulaVar("Ap0r_a","ap0r*(1+rp*r0*cos(sp-s0))+ap0i*rp*r0*sin(sp-s0)+\
                              (ap0r*(1-rp*r0*cos(sp-s0))-ap0i*rp*r0*sin(sp-s0))*cos(dm*dt)-\
                              (ap0r*(rp*sin(phiw-sp)+r0*sin(phiw-s0))+\
                               ap0i*(rp*cos(phiw-sp)-r0*cos(phiw-s0)))*sin(dm*dt)",RooArgSet(*ap0r,*ap0i,*rp,*r0,*sp,*s0,*dm,*dt,*phiw));

    A0ti_a = new RooFormulaVar("A0ti_a","a0ti*(1-r0*rt*cos(s0-st))+a0tr*r0*rt*sin(s0-st)+\
                              (a0ti*(1+r0*rt*cos(s0-st))-a0tr*r0*rt*sin(s0-st))*cos(dm*dt)-\
                              (a0ti*(r0*sin(phiw-s0)-rt*sin(phiw-st))-\
                               a0tr*(r0*cos(phiw-s0)+rt*cos(phiw-st)))*sin(dm*dt)",RooArgSet(*a0ti,*a0tr,*r0,*rt,*s0,*st,*dm,*dt,*phiw));

    Apti_a = new RooFormulaVar("Apti_a","apti*(1-rp*rt*cos(sp-st))+aptr*rp*rt*sin(sp-st)+\
                              (apti*(1+rp*rt*cos(sp-st))-aptr*rp*rt*sin(sp-st))*cos(dm*dt)-\
                              (apti*(rp*sin(phiw-sp)-rt*sin(phiw-st))-\
                               aptr*(rp*cos(phiw-sp)+rt*cos(phiw-st)))*sin(dm*dt)",RooArgSet(*apti,*aptr,*rp,*rt,*sp,*st,*dm,*dt,*phiw));


    At2_ab = new RooFormulaVar("At2_ab","at*at*((1+rt*rt)+(1-rt*rt)*cos(dm*dt)+2*rt*sin(-phiw-st)*sin(dm*dt))",RooArgSet(*at,*rt,*dm,*dt,*phiw,*st));
    Ap2_ab = new RooFormulaVar("Ap2_ab","ap*ap*((1+rp*rp)+(1-rp*rp)*cos(dm*dt)-2*rp*sin(-phiw-sp)*sin(dm*dt))",RooArgSet(*ap,*rp,*dm,*dt,*phiw,*sp));
    A02_ab = new RooFormulaVar("A02_ab","a0*a0*((1+r0*r0)+(1-r0*r0)*cos(dm*dt)-2*r0*sin(-phiw-s0)*sin(dm*dt))",RooArgSet(*a0,*r0,*dm,*dt,*phiw,*s0));

    Ap0r_ab = new RooFormulaVar("Ap0r_ab","ap0r*(1+rp*r0*cos(sp-s0))+ap0i*rp*r0*sin(sp-s0)+\
                              (ap0r*(1-rp*r0*cos(sp-s0))-ap0i*rp*r0*sin(sp-s0))*cos(dm*dt)-\
                              (ap0r*(rp*sin(-phiw-sp)+r0*sin(-phiw-s0))+\
                               ap0i*(rp*cos(-phiw-sp)-r0*cos(-phiw-s0)))*sin(dm*dt)",RooArgSet(*ap0r,*ap0i,*rp,*r0,*sp,*s0,*dm,*dt,*phiw));

    A0ti_ab = new RooFormulaVar("A0ti_ab","a0ti*(1-r0*rt*cos(s0-st))+a0tr*r0*rt*sin(s0-st)+\
                              (a0ti*(1+r0*rt*cos(s0-st))-a0tr*r0*rt*sin(s0-st))*cos(dm*dt)-\
                              (a0ti*(r0*sin(-phiw-s0)-rt*sin(-phiw-st))-\
                               a0tr*(r0*cos(-phiw-s0)+rt*cos(-phiw-st)))*sin(dm*dt)",RooArgSet(*a0ti,*a0tr,*r0,*rt,*s0,*st,*dm,*dt,*phiw));

    Apti_ab = new RooFormulaVar("Apti_ab","apti*(1-rp*rt*cos(sp-st))+aptr*rp*rt*sin(sp-st)+\
                              (apti*(1+rp*rt*cos(sp-st))-aptr*rp*rt*sin(sp-st))*cos(dm*dt)-\
                              (apti*(rp*sin(-phiw-sp)-rt*sin(-phiw-st))-\
                               aptr*(rp*cos(-phiw-sp)+rt*cos(-phiw-st)))*sin(dm*dt)",RooArgSet(*apti,*aptr,*rp,*rt,*sp,*st,*dm,*dt,*phiw));


    At2_b = new RooFormulaVar("At2_b","at*at*((1+rt*rt)+(1-rt*rt)*(-1)*cos(dm*dt)+2*rt*sin(phiw-st)*(-1)*sin(dm*dt))",RooArgSet(*at,*rt,*dm,*dt,*phiw,*st));
    Ap2_b = new RooFormulaVar("Ap2_b","ap*ap*((1+rp*rp)+(1-rp*rp)*(-1)*cos(dm*dt)-2*rp*sin(phiw-sp)*(-1)*sin(dm*dt))",RooArgSet(*ap,*rp,*dm,*dt,*phiw,*sp));
    A02_b = new RooFormulaVar("A02_b","a0*a0*((1+r0*r0)+(1-r0*r0)*(-1)*cos(dm*dt)-2*r0*sin(phiw-s0)*(-1)*sin(dm*dt))",RooArgSet(*a0,*r0,*dm,*dt,*phiw,*s0));

    Ap0r_b = new RooFormulaVar("Ap0r_b","ap0r*(1+rp*r0*cos(sp-s0))+ap0i*rp*r0*sin(sp-s0)+\
                              (ap0r*(1-rp*r0*cos(sp-s0))-ap0i*rp*r0*sin(sp-s0))*(-1)*cos(dm*dt)-\
                              (ap0r*(rp*sin(phiw-sp)+r0*sin(phiw-s0))+\
                               ap0i*(rp*cos(phiw-sp)-r0*cos(phiw-s0)))*(-1)*sin(dm*dt)",RooArgSet(*ap0r,*ap0i,*rp,*r0,*sp,*s0,*dm,*dt,*phiw));

    A0ti_b = new RooFormulaVar("A0ti_b","a0ti*(1-r0*rt*cos(s0-st))+a0tr*r0*rt*sin(s0-st)+\
                              (a0ti*(1+r0*rt*cos(s0-st))-a0tr*r0*rt*sin(s0-st))*(-1)*cos(dm*dt)-\
                              (a0ti*(r0*sin(phiw-s0)-rt*sin(phiw-st))-\
                               a0tr*(r0*cos(phiw-s0)+rt*cos(phiw-st)))*(-1)*sin(dm*dt)",RooArgSet(*a0ti,*a0tr,*r0,*rt,*s0,*st,*dm,*dt,*phiw));

    Apti_b = new RooFormulaVar("Apti_b","apti*(1-rp*rt*cos(sp-st))+aptr*rp*rt*sin(sp-st)+\
                              (apti*(1+rp*rt*cos(sp-st))-aptr*rp*rt*sin(sp-st))*(-1)*cos(dm*dt)-\
                              (apti*(rp*sin(phiw-sp)-rt*sin(phiw-st))-\
                               aptr*(rp*cos(phiw-sp)+rt*cos(phiw-st)))*(-1)*sin(dm*dt)",RooArgSet(*apti,*aptr,*rp,*rt,*sp,*st,*dm,*dt,*phiw));


    At2_bb = new RooFormulaVar("At2_bb","at*at*((1+rt*rt)+(1-rt*rt)*(-1)*cos(dm*dt)+2*rt*sin(-phiw-st)*(-1)*sin(dm*dt))",RooArgSet(*at,*rt,*dm,*dt,*phiw,*st));
    Ap2_bb = new RooFormulaVar("Ap2_bb","ap*ap*((1+rp*rp)+(1-rp*rp)*(-1)*cos(dm*dt)-2*rp*sin(-phiw-sp)*(-1)*sin(dm*dt))",RooArgSet(*ap,*rp,*dm,*dt,*phiw,*sp));
    A02_bb = new RooFormulaVar("A02_bb","a0*a0*((1+r0*r0)+(1-r0*r0)*(-1)*cos(dm*dt)-2*r0*sin(-phiw-s0)*(-1)*sin(dm*dt))",RooArgSet(*a0,*r0,*dm,*dt,*phiw,*s0));

    Ap0r_bb = new RooFormulaVar("Ap0r_bb","ap0r*(1+rp*r0*cos(sp-s0))+ap0i*rp*r0*sin(sp-s0)+\
                              (ap0r*(1-rp*r0*cos(sp-s0))-ap0i*rp*r0*sin(sp-s0))*(-1)*cos(dm*dt)-\
                              (ap0r*(rp*sin(-phiw-sp)+r0*sin(-phiw-s0))+\
                               ap0i*(rp*cos(-phiw-sp)-r0*cos(-phiw-s0)))*(-1)*sin(dm*dt)",RooArgSet(*ap0r,*ap0i,*rp,*r0,*sp,*s0,*dm,*dt,*phiw));

    A0ti_bb = new RooFormulaVar("A0ti_bb","a0ti*(1-r0*rt*cos(s0-st))+a0tr*r0*rt*sin(s0-st)+\
                              (a0ti*(1+r0*rt*cos(s0-st))-a0tr*r0*rt*sin(s0-st))*(-1)*cos(dm*dt)-\
                              (a0ti*(r0*sin(-phiw-s0)-rt*sin(-phiw-st))-\
                               a0tr*(r0*cos(-phiw-s0)+rt*cos(-phiw-st)))*(-1)*sin(dm*dt)",RooArgSet(*a0ti,*a0tr,*r0,*rt,*s0,*st,*dm,*dt,*phiw));

    Apti_bb = new RooFormulaVar("Apti_bb","apti*(1-rp*rt*cos(sp-st))+aptr*rp*rt*sin(sp-st)+\
                              (apti*(1+rp*rt*cos(sp-st))-aptr*rp*rt*sin(sp-st))*(-1)*cos(dm*dt)-\
                              (apti*(rp*sin(-phiw-sp)-rt*sin(-phiw-st))-\
                               aptr*(rp*cos(-phiw-sp)+rt*cos(-phiw-st)))*(-1)*sin(dm*dt)",RooArgSet(*apti,*aptr,*rp,*rt,*sp,*st,*dm,*dt,*phiw));

    /// a decays are favored, b corresponding suppressed; a is B0 -> D*- + rho+
    formula_a ="exp(-gamma*abs(dt))*(Ap2_a*2*sin(tht)*sin(tht)*sin(tht)*sin(thb)*sin(thb)*sin(thb)*sin(phit)*sin(phit)+\
                        At2_a*2*cos(tht)*cos(tht)*sin(tht)*sin(thb)*sin(thb)*sin(thb)+\
                        A02_a*4*sin(tht)*sin(tht)*sin(tht)*cos(thb)*cos(thb)*sin(thb)*cos(phit)*cos(phit)+\
                        sqrt(2)*Ap0r_a*sin(tht)*sin(tht)*sin(tht)*sin(2*thb)*sin(thb)*sin(2*phit)-\
                        sqrt(2)*A0ti_a*sin(2*tht)*sin(tht)*sin(2*thb)*sin(thb)*cos(phit)-\
                        2*Apti_a*sin(2*tht)*sin(tht)*sin(thb)*sin(thb)*sin(thb)*sin(phit))";


    formula_ab="exp(-gamma*abs(dt))*(Ap2_ab*2*sin(tht)*sin(tht)*sin(tht)*sin(thb)*sin(thb)*sin(thb)*sin(phit)*sin(phit)+\
                        At2_ab*2*cos(tht)*cos(tht)*sin(tht)*sin(thb)*sin(thb)*sin(thb)+\
                        A02_ab*4*sin(tht)*sin(tht)*sin(tht)*cos(thb)*cos(thb)*sin(thb)*cos(phit)*cos(phit)+\
                        sqrt(2)*Ap0r_ab*sin(tht)*sin(tht)*sin(tht)*sin(2*thb)*sin(thb)*sin(2*phit)-\
                        sqrt(2)*A0ti_ab*sin(2*tht)*sin(tht)*sin(2*thb)*sin(thb)*cos(phit)-\
                        2*Apti_ab*sin(2*tht)*sin(tht)*sin(thb)*sin(thb)*sin(thb)*sin(phit))";


    formula_b="exp(-gamma*abs(dt))*(Ap2_b*2*sin(tht)*sin(tht)*sin(tht)*sin(thb)*sin(thb)*sin(thb)*sin(phit)*sin(phit)+\
                        At2_b*2*cos(tht)*cos(tht)*sin(tht)*sin(thb)*sin(thb)*sin(thb)+\
                        A02_b*4*sin(tht)*sin(tht)*sin(tht)*cos(thb)*cos(thb)*sin(thb)*cos(phit)*cos(phit)+\
                        sqrt(2)*Ap0r_b*sin(tht)*sin(tht)*sin(tht)*sin(2*thb)*sin(thb)*sin(2*phit)-\
                        sqrt(2)*A0ti_b*sin(2*tht)*sin(tht)*sin(2*thb)*sin(thb)*cos(phit)-\
                        2*Apti_b*sin(2*tht)*sin(tht)*sin(thb)*sin(thb)*sin(thb)*sin(phit))";


    formula_bb="exp(-gamma*abs(dt))*(Ap2_bb*2*sin(tht)*sin(tht)*sin(tht)*sin(thb)*sin(thb)*sin(thb)*sin(phit)*sin(phit)+\
                        At2_bb*2*cos(tht)*cos(tht)*sin(tht)*sin(thb)*sin(thb)*sin(thb)+\
                        A02_bb*4*sin(tht)*sin(tht)*sin(tht)*cos(thb)*cos(thb)*sin(thb)*cos(phit)*cos(phit)+\
                        sqrt(2)*Ap0r_bb*sin(tht)*sin(tht)*sin(tht)*sin(2*thb)*sin(thb)*sin(2*phit)-\
                        sqrt(2)*A0ti_bb*sin(2*tht)*sin(tht)*sin(2*thb)*sin(thb)*cos(phit)-\
                        2*Apti_bb*sin(2*tht)*sin(tht)*sin(thb)*sin(thb)*sin(thb)*sin(phit))";

    //formula_bb="3*dt";

    varSet_a = new RooArgSet(*Ap2_a,*At2_a,*A02_a,*Ap0r_a,*A0ti_a,*Apti_a,*tht,*thb,*phit);
    varSet_a->add(*dt);
    varSet_a->add(*gamma);
    varSet_b = new RooArgSet(*Ap2_b,*At2_b,*A02_b,*Ap0r_b,*A0ti_b,*Apti_b,*tht,*thb,*phit);
    varSet_b->add(*dt);
    varSet_b->add(*gamma);
    varSet_ab = new RooArgSet(*Ap2_ab,*At2_ab,*A02_ab,*Ap0r_ab,*A0ti_ab,*Apti_ab,*tht,*thb,*phit);
    varSet_ab->add(*dt);
    varSet_ab->add(*gamma);
    varSet_bb = new RooArgSet(*Ap2_bb,*At2_bb,*A02_bb,*Ap0r_bb,*A0ti_bb,*Apti_bb,*tht,*thb,*phit);
    varSet_bb->add(*dt);
    varSet_bb->add(*gamma);

    pdf_a = new DSRhoPDF("pdf_a","pdf_a","a",*tht,*thb,*phit,*dt,*ap,*apa,*a0,*ata,*phiw,*rp,*r0,*rt,*sp,*s0,*st);
    pdf_b = new DSRhoPDF("pdf_b","pdf_b","b",*tht,*thb,*phit,*dt,*ap,*apa,*a0,*ata,*phiw,*rp,*r0,*rt,*sp,*s0,*st);
    pdf_ab = new DSRhoPDF("pdf_ab","pdf_ab","ab",*tht,*thb,*phit,*dt,*ap,*apa,*a0,*ata,*phiw,*rp,*r0,*rt,*sp,*s0,*st);
    pdf_bb = new DSRhoPDF("pdf_bb","pdf_bb","bb",*tht,*thb,*phit,*dt,*ap,*apa,*a0,*ata,*phiw,*rp,*r0,*rt,*sp,*s0,*st);

    simPdf = new RooSimultaneous("simPdf","simPdf",*decType);
    simPdf->addPdf(*pdf_a,"a");
    simPdf->addPdf(*pdf_ab,"ab");
    simPdf->addPdf(*pdf_b,"b");
    simPdf->addPdf(*pdf_bb,"bb");

    parameters = new RooArgSet(*ap,*apa,*a0,*ata,*phiw,*rp,*r0,*rt);
    parameters->add(*sp);
    parameters->add(*s0);
    parameters->add(*st);

    varSet = new RooArgSet(*tht,*thb,*phit,*ap,*a0,*at,*ap0r,*a0ti,*apti);

    /// numFitParameters holds # of NON-constant fit parameters
    fitParameters = new RooArgSet(*ap,*apa,*a0,*ata);
    numFitParameters = (parameters->selectByAttrib("Constant",kFALSE))->getSize();

//    pdfFormula =   "ap*ap*2*sin(tht)*sin(tht)*sin(tht)*sin(thb)*sin(thb)*sin(thb)*sin(phit)*sin(phit)+\
//                                at*at*2*cos(tht)*cos(tht)*sin(tht)*sin(thb)*sin(thb)*sin(thb)+\
//                                a0*a0*4*sin(tht)*sin(tht)*sin(tht)*cos(thb)*cos(thb)*sin(thb)*cos(phit)*cos(phit)+\
//                                sqrt(2)*ap0r*sin(tht)*sin(tht)*sin(tht)*sin(2*thb)*sin(thb)*sin(2*phit)-\
//                                sqrt(2)*a0ti*sin(2*tht)*sin(tht)*sin(2*thb)*sin(thb)*cos(phit)-\
//                                2*apti*sin(2*tht)*sin(tht)*sin(thb)*sin(thb)*sin(thb)*sin(phit)";
//
//    pdf = new RooGenericPdf("pdf","Generic PDF",pdfFormula,*varSet);
//
//    myPdf_a = new DSRhoPDF("myPdf","myPdf","b",*tht,*thb,*phit,*dt,*ap,*apa,*a0,*ata,*phiw,*rp,*r0,*rt,*sp,*s0,*st);
//    myPdf_b = new DSRhoPDF("myPdf","myPdf","b",*tht,*thb,*phit,*dt,*ap,*apa,*a0,*ata,*phiw,*rp,*r0,*rt,*sp,*s0,*st);
}

FitterTrans::~FitterTrans()
{
    delete gPluginMgr;
    delete thb;
    delete tht;
    delete phit;
    delete dt;
    delete decType;
    delete gamma;
    delete ap;
    delete apa;
    delete apr;
    delete api;
    delete a0;
    delete a0a;
    delete a0r;
    delete a0i;
    delete at;
    delete ata;
    delete atr;
    delete ati;
    delete ap0r;
    delete a0ti;
    delete apti;
    delete dm;
    delete phiw;
    delete rp;
    delete r0;
    delete rt;
    delete sp;
    delete s0;
    delete st;
    delete ap0i;
    delete a0tr;
    delete aptr;
    delete At2_a;
    delete Ap2_a;
    delete A02_a;
    delete Ap0r_a;
    delete A0ti_a;
    delete Apti_a;
    delete At2_ab;
    delete Ap2_ab;
    delete A02_ab;
    delete Ap0r_ab;
    delete A0ti_ab;
    delete Apti_ab;
    delete At2_b;
    delete Ap2_b;
    delete A02_b;
    delete Ap0r_b;
    delete A0ti_b;
    delete Apti_b;
    delete At2_bb;
    delete Ap2_bb;
    delete A02_bb;
    delete Ap0r_bb;
    delete A0ti_bb;
    delete Apti_bb;
    delete varSet_a;
    delete varSet_b;
    delete varSet_ab;
    delete varSet_bb;
    delete pdf_a;
    delete pdf_b;
    delete pdf_ab;
    delete pdf_bb;
    delete simPdf;
    delete parameters;
    delete fitParameters;
}

Int_t FitterTrans::Fit()
{
    //TPluginManager* gPluginMgr = new TPluginManager;
    //gPluginMgr->AddHandler("ROOT::Math::Minimizer", "Minuit2", "Minuit2Minimizer", "Minuit2", "Minuit2Minimizer(const char *)");
    //result = pdf->fitTo(*dataSet,RooFit::Save(),RooFit::Timer(true),RooFit::Minimizer("Minuit2"));//,RooFit::NumCPU(2));
    result = simPdf->fitTo(*dataSet,RooFit::Save(),RooFit::Timer(true),RooFit::Minimizer("Minuit2"));//,RooFit::NumCPU(2));
}

void FitterTrans::CreateBinnedDataSet(const char* type)
{
    tht->setBins(tht_bins);
    thb->setBins(thb_bins);
    phit->setBins(phit_bins);
    dt->setBins(dt_bins);

    TString cut = "decType==decType::";
    cut += type;

    RooRandom::randomGenerator()->SetSeed(0);
    RooDataSet* dataSet_reduced = (RooDataSet*)dataSet->reduce(cut);
    binnedNumEntries = dataSet_reduced->numEntries();

    /// Create a binned dataSet which is needed for chi2 calculation
    dataSet_binned = new RooDataHist("dataSet_binned","dataSet_binned",RooArgSet(*tht,*thb,*phit,*dt),*dataSet_reduced);

    //dataSet_binned = new RooDataHist("dataSet_binned","dataSet_binned",RooArgSet(*tht,*thb,*phit,*dt),*pdf_a);
    //RooDataHist* dataSet_binned = pdf->generateBinned(RooArgSet(var1,var2,var3),dataSet->numEntries(),kFALSE);
}

Int_t FitterTrans::ComputeChi2(const char* type)
{
    //if(dataSet_binned == NULL)
        CreateBinnedDataSet(type);

    numFitParameters = (parameters->selectByAttrib("Constant",kFALSE))->getSize();

    if      (strcmp(type,"a") == 0)   chi2Var = new RooChi2Var("chi2Var","chi2Var",*pdf_a,*dataSet_binned);
    else if (strcmp(type,"b") == 0)   chi2Var = new RooChi2Var("chi2Var","chi2Var",*pdf_b,*dataSet_binned);
    else if (strcmp(type,"ab") == 0)  chi2Var = new RooChi2Var("chi2Var","chi2Var",*pdf_ab,*dataSet_binned);
    else if (strcmp(type,"bb") == 0)  chi2Var = new RooChi2Var("chi2Var","chi2Var",*pdf_bb,*dataSet_binned);

	//chi2Var = new RooChi2Var("chi2Var","chi2Var",*simPdf,*dataSet_binned);

	RooRealVar* ndof     = new RooRealVar("ndof","number of degrees of freedom",0);
	RooRealVar* chi2red  = new RooRealVar("chi2red","reduced chi^2",0);
	RooRealVar* prob     = new RooRealVar("prob","prob(chi2,ndof)",0);

	ndof->setVal(dataSet_binned->numEntries()-numFitParameters);
	chi2red->setVal(chi2Var->getVal()/ndof->getVal()) ;
	prob->setVal(TMath::Prob(chi2Var->getVal(),static_cast<int>(ndof->getVal())));

	delete dataSet_binned;

	printf("%s:\tchi2 = %f\n%s:\tndof = %f\n%s:\tchi2red = %f\n%s:\tprob = %f\n\n",type,chi2Var->getVal(),type,ndof->getVal(),type,chi2red->getVal(),type,prob->getVal());

}


Double_t FitterTrans::GetChi2(const char* type)
{
    Double_t mychi2 = 0;
    Double_t n = 0;
    Double_t v = 0;

//    TH1D* h_dchi2 = new TH1D("dchi2","dchi2",100,0,100000);


    //if(dataSet_binned == NULL)
        CreateBinnedDataSet(type);

    DSRhoPDF* pdf = 0;

    if(strcmp(type,"a") == 0)       pdf = pdf_a;
    else if(strcmp(type,"b") == 0)  pdf = pdf_b;
    else if(strcmp(type,"ab") == 0) pdf = pdf_ab;
    else if(strcmp(type,"bb") == 0) pdf = pdf_bb;

    Double_t binVolume = tht->getBinWidth(0)*thb->getBinWidth(0)*phit->getBinWidth(0)*dt->getBinWidth(0);
    Int_t numBins = dataSet_binned->numEntries();

    Int_t numVPrecise = 0;

    /// Cycle through the centers of all bins
    /// I'm getting width of the first bin, because all bins are of equal width
    for(*tht = tht->getMin()+tht->getBinWidth(0)/2; tht->getVal() < tht->getMax(); tht->setVal(tht->getVal()+tht->getBinWidth(0)))
    {
        for(*thb = thb->getMin()+thb->getBinWidth(0)/2; thb->getVal() < thb->getMax(); thb->setVal(thb->getVal()+thb->getBinWidth(0)))
        {
            for(*phit = phit->getMin()+phit->getBinWidth(0)/2; phit->getVal() < phit->getMax(); phit->setVal(phit->getVal()+phit->getBinWidth(0)))
            {
                for(*dt = dt->getMin()+dt->getBinWidth(0)/2; dt->getVal() < dt->getMax(); dt->setVal(dt->getVal()+dt->getBinWidth(0)))
                {
                    /// Weight is actually the bin content
                    n = dataSet_binned->weight(RooArgSet(*tht,*thb,*phit,*dt),0);
                    if(n < 5)
                    {
                        numBins--;
                        continue;
                    }

                    v = pdf->getVal(RooArgSet(*tht,*thb,*phit,*dt))*binVolume*binnedNumEntries;

                    if(((n-v)*(n-v)/v) > 10)
                    {
                        v = GetVPrecise(pdf);
                        numVPrecise++;
                    }

//                    printf("%.10f\t%.10f\n",v,GetVPrecise(pdf));

                    mychi2 += (n-v)*(n-v)/v;
//                    h_dchi2->Fill((n-v)*(n-v)/v);
                }
            }
        }
    }

    printf("# VPrecise called: %i\n",numVPrecise);

//    TCanvas* c1 = new TCanvas("c1","c1",800,600);
//    c1->SetLogy();
//    h_dchi2->Draw();

    delete dataSet_binned;

    //printf("binVolume = %f\n",binVolume);
    printf("%s: numEntries = %i\n",type,binnedNumEntries);
    printf("%s: numBins = %i\n",type,numBins);
    printf("%s: mychi2 = %f\n",type,mychi2);
    printf("%s: mychi2red = %f\n",type,mychi2/numBins);
    printf("%s: prob = %.10f\n\n",type,TMath::Prob(mychi2,numBins));
    return mychi2;
}

Double_t FitterTrans::GetVPrecise(DSRhoPDF* pdf)
{
    /// The pdf seems to be varying quite rapidly at some places, so the approximation of constant pdf in a voxel is sometimes bad.
    /// This function gives a better approximation of a pdf value in a voxel by averaging through multiple points inside it.

    Double_t v = 0;

    Int_t tht_subbins = 3;
    Int_t thb_subbins = 3;
    Int_t phit_subbins = 3;
    Int_t dt_subbins = 3;

    Double_t binVolume = tht->getBinWidth(0)*thb->getBinWidth(0)*phit->getBinWidth(0)*dt->getBinWidth(0);

    Double_t tht_binmin = tht->getVal() - tht->getBinWidth(0) /2;
    Double_t tht_binmax = tht->getVal() + tht->getBinWidth(0) /2;
    Double_t thb_binmin = thb->getVal() - thb->getBinWidth(0) /2;
    Double_t thb_binmax = thb->getVal() + thb->getBinWidth(0) /2;
    Double_t phit_binmin = phit->getVal() - phit->getBinWidth(0) /2;
    Double_t phit_binmax = phit->getVal() + phit->getBinWidth(0) /2;
    Double_t dt_binmin = dt->getVal() - dt->getBinWidth(0) /2;
    Double_t dt_binmax = dt->getVal() + dt->getBinWidth(0) /2;

    Int_t numPasses = 0;
    Int_t pass = 0;

    /// Using *_binmax - 0.001 because when one is at a boundary of e.g. thb, thb->getVal() < thb_binmax is never violated, even though it should be equal.
    for(*tht = tht_binmin+tht->getBinWidth(0)/(2*tht_subbins); tht->getVal() < tht_binmax-0.001; tht->setVal(tht->getVal()+tht->getBinWidth(0)/tht_subbins))
    {
        for(*thb = thb_binmin+thb->getBinWidth(0)/(2*thb_subbins); thb->getVal() < thb_binmax-0.001; thb->setVal(thb->getVal()+thb->getBinWidth(0)/thb_subbins))
        {
            for(*phit = phit_binmin+phit->getBinWidth(0)/(2*phit_subbins); phit->getVal() < phit_binmax-0.001; phit->setVal(phit->getVal()+phit->getBinWidth(0)/phit_subbins))
            {
                for(*dt = dt_binmin+dt->getBinWidth(0)/(2*dt_subbins); dt->getVal() < dt_binmax-0.001; dt->setVal(dt->getVal()+dt->getBinWidth(0)/dt_subbins))
                {
                    v += pdf->getVal(RooArgSet(*tht,*thb,*phit,*dt))*binVolume*binnedNumEntries;
                }
            }
        }
    }

    /// These lines return the variables to their initial states, before calling GetVPrecise()
    tht->setVal(tht_binmin + tht->getBinWidth(0)/2);
    thb->setVal(thb_binmin + thb->getBinWidth(0)/2);
    phit->setVal(phit_binmin + phit->getBinWidth(0)/2);
    dt->setVal(dt_binmin + dt->getBinWidth(0)/2);

    v = v/(tht_subbins*thb_subbins*phit_subbins*dt_subbins);

    return v;
}

Double_t FitterTrans::SaveChi2Maps(const char* type)
{
    /// Create a histogram from the pdf with the expected number of events with no statistical fluctuation
    //RooDataHist* pdf_binned = pdf->generateBinned(RooArgSet(var1,var2,var3),numEvents,kTRUE);

    Double_t mychi2 = 0;
    Double_t dchi2 = 0;
    Double_t n = 0;
    Double_t v = 0;

    RooRealVar* var1 = tht;
    RooRealVar* var2 = thb;
    RooRealVar* var3 = phit;
    RooRealVar* var4 = dt;
    Int_t var1_bins = tht_bins;
    Int_t var2_bins = thb_bins;
    Int_t var3_bins = phit_bins;
    Int_t var4_bins = dt_bins;

    DSRhoPDF* pdf = 0;

    if(strcmp(type,"a") == 0)       pdf = pdf_a;
    else if(strcmp(type,"b") == 0)  pdf = pdf_b;
    else if(strcmp(type,"ab") == 0) pdf = pdf_ab;
    else if(strcmp(type,"bb") == 0) pdf = pdf_bb;

    //if(dataSet_binned == NULL)
        CreateBinnedDataSet(type);

    Double_t binVolume = tht->getBinWidth(0)*thb->getBinWidth(0)*phit->getBinWidth(0)*dt->getBinWidth(0);

    TH1F* h1_chi2 = new TH1F("h1_chi2","h1_chi2",100,var1->getMin(),var1->getMax());
    TH2F* h2_chi2_1 = new TH2F("h2_chi2_1","h2_chi2_1",var1_bins,var1->getMin(),var1->getMax(),var2_bins,var2->getMin(),var2->getMax());
    TH2F* h2_chi2_2 = new TH2F("h2_chi2_2","h2_chi2_2",var1_bins,var1->getMin(),var1->getMax(),var3_bins,var3->getMin(),var3->getMax());
    TH2F* h2_chi2_3 = new TH2F("h2_chi2_3","h2_chi2_3",var1_bins,var1->getMin(),var1->getMax(),var4_bins,var4->getMin(),var4->getMax());
    TH2F* h2_chi2_4 = new TH2F("h2_chi2_4","h2_chi2_4",var2_bins,var2->getMin(),var2->getMax(),var3_bins,var3->getMin(),var3->getMax());
    TH2F* h2_chi2_5 = new TH2F("h2_chi2_5","h2_chi2_5",var2_bins,var2->getMin(),var2->getMax(),var4_bins,var4->getMin(),var4->getMax());
    TH2F* h2_chi2_6 = new TH2F("h2_chi2_6","h2_chi2_6",var3_bins,var3->getMin(),var3->getMax(),var4_bins,var4->getMin(),var4->getMax());

    Int_t vars_bins[4] = {tht_bins,thb_bins,phit_bins,dt_bins};
    RooRealVar* vars[4] = {tht,thb,phit,dt};

    TString name;
    TH1F* h1_resid[4];

    for(int i = 0; i < 4; i++)
    {
        name = "h1_resid_";
        name += i+1;
        h1_resid[i] = new TH1F(name,name,vars_bins[i],vars[i]->getMin(),vars[i]->getMax());
    }


    TH1F* h1_resid_temp = new TH1F("resid_test","resid_test",vars_bins[0],vars[0]->getMin(),vars[0]->getMax());

    TString cut = "decType==decType::";
    cut += type;

    RooRandom::randomGenerator()->SetSeed(0);
    RooDataSet* dataSet_reduced = (RooDataSet*)dataSet->reduce(cut);
    Int_t temp_binnedNumEntries = dataSet_reduced->numEntries();

    /// Create a binned dataSet which is needed for chi2 calculation
    RooDataHist* temp_dataSet_binned = new RooDataHist("temp_dataSet_binned","temp_dataSet_binned",RooArgSet(*tht),*dataSet_reduced);

    for(*var1 = var1->getMin()+var1->getBinWidth(0)/2; var1->getVal() < var1->getMax(); var1->setVal(var1->getVal()+var1->getBinWidth(0)))
    {
        n = temp_dataSet_binned->weight(RooArgSet(*var1),0);
        v = pdf->getVal(RooArgSet(*var1,*var2,*var3,*var4))*binVolume*binnedNumEntries;
    }



    /// Cycle through the centers of all bins
    /// I'm getting width of the first bin, because all bins are of equal width
    for(*var1 = var1->getMin()+var1->getBinWidth(0)/2; var1->getVal() < var1->getMax(); var1->setVal(var1->getVal()+var1->getBinWidth(0)))
    {
        for(*var2 = var2->getMin()+var2->getBinWidth(0)/2; var2->getVal() < var2->getMax(); var2->setVal(var2->getVal()+var2->getBinWidth(0)))
        {
            for(*var3 = var3->getMin()+var3->getBinWidth(0)/2; var3->getVal() < var3->getMax(); var3->setVal(var3->getVal()+var3->getBinWidth(0)))
            {
                for(*var4 = var4->getMin()+var4->getBinWidth(0)/2; var4->getVal() < var4->getMax(); var4->setVal(var4->getVal()+var4->getBinWidth(0)))
                {

                    /// Weight is actually the bin content
                    n = dataSet_binned->weight(RooArgSet(*var1,*var2,*var3,*var4),0);
                    //if(n == 0) continue;
                    v = pdf->getVal(RooArgSet(*var1,*var2,*var3,*var4))*binVolume*binnedNumEntries;

                    if(((n-v)*(n-v)/v) > 10)
                    {
                        v = GetVPrecise(pdf);
                    }

                    dchi2 = (n-v)*(n-v)/v;

                    h1_chi2->Fill(dchi2);
                    h2_chi2_1->Fill(var1->getVal(),var2->getVal(),dchi2);
                    h2_chi2_2->Fill(var1->getVal(),var3->getVal(),dchi2);
                    h2_chi2_3->Fill(var1->getVal(),var4->getVal(),dchi2);
                    h2_chi2_4->Fill(var2->getVal(),var3->getVal(),dchi2);
                    h2_chi2_5->Fill(var2->getVal(),var4->getVal(),dchi2);
                    h2_chi2_6->Fill(var3->getVal(),var4->getVal(),dchi2);
                    mychi2 += dchi2;

                    for(int i = 0; i < 4; i++)
                        h1_resid[i]->Fill(vars[i]->getVal(),(n-v)/v);
                }
            }
        }
    }

    delete dataSet_binned;

    TFile* file = new TFile("plots/chi2maps.root","RECREATE");
    TCanvas* c2 = new TCanvas("c2","c2",600,600);
    TString path;

    for(int i = 0; i < 4; i++)
    {
        h1_resid[i]->GetXaxis()->SetTitle(vars[i]->GetName());
        h1_resid[i]->Draw();
        h1_resid[i]->Write();
        path = "plots/residual_";
        path += vars[i]->GetName();
        path += ".png";
        c2->SaveAs(path);
    }


    //c2->SetLogy(kTRUE);
    h1_chi2->GetXaxis()->SetTitle("dchi");
    h1_chi2->GetYaxis()->SetTitle("num bins");
    h1_chi2->Draw();
    h1_chi2->Write();
    c2->SaveAs("plots/chi2_delta.png");

    c2->SetLogy(kFALSE);

    h2_chi2_1->SetOption("colz");
    //hdchi2->SetMinimum(0);
    //hdchi2->SetMaximum(100);
    h2_chi2_1->SetStats(kFALSE);
    h2_chi2_1->GetXaxis()->SetTitle(var1->GetName());
    h2_chi2_1->GetYaxis()->SetTitle(var2->GetName());
    h2_chi2_1->Draw();
    h2_chi2_1->Write();
    path = "plots/chi2map_";
    path += var1->GetName();
    path += "_";
    path += var2->GetName();
    path += ".png";
    c2->SaveAs(path);

    h2_chi2_2->SetOption("colz");
    //hdchi2->SetMinimum(0);
    //hdchi2->SetMaximum(100);
    h2_chi2_2->SetStats(kFALSE);
    h2_chi2_2->GetXaxis()->SetTitle(var1->GetName());
    h2_chi2_2->GetYaxis()->SetTitle(var3->GetName());
    h2_chi2_2->Draw();
    h2_chi2_2->Write();
    path = "plots/chi2map_";
    path += var1->GetName();
    path += "_";
    path += var3->GetName();
    path += ".png";
    c2->SaveAs(path);

    h2_chi2_3->SetOption("colz");
    //hdchi2->SetMinimum(0);
    //hdchi2->SetMaximum(100);
    h2_chi2_3->SetStats(kFALSE);
    h2_chi2_3->GetXaxis()->SetTitle(var1->GetName());
    h2_chi2_3->GetYaxis()->SetTitle(var4->GetName());
    h2_chi2_3->Draw();
    h2_chi2_3->Write();
    path = "plots/chi2map_";
    path += var1->GetName();
    path += "_";
    path += var4->GetName();
    path += ".png";
    c2->SaveAs(path);

    h2_chi2_4->SetOption("colz");
    //hdchi2->SetMinimum(0);
    //hdchi2->SetMaximum(100);
    h2_chi2_4->SetStats(kFALSE);
    h2_chi2_4->GetXaxis()->SetTitle(var2->GetName());
    h2_chi2_4->GetYaxis()->SetTitle(var3->GetName());
    h2_chi2_4->Draw();
    h2_chi2_4->Write();
    path = "plots/chi2map_";
    path += var2->GetName();
    path += "_";
    path += var3->GetName();
    path += ".png";
    c2->SaveAs(path);

    h2_chi2_5->SetOption("colz");
    //hdchi2->SetMinimum(0);
    //hdchi2->SetMaximum(100);
    h2_chi2_5->SetStats(kFALSE);
    h2_chi2_5->GetXaxis()->SetTitle(var2->GetName());
    h2_chi2_5->GetYaxis()->SetTitle(var4->GetName());
    h2_chi2_5->Draw();
    h2_chi2_5->Write();
    path = "plots/chi2map_";
    path += var2->GetName();
    path += "_";
    path += var4->GetName();
    path += ".png";
    c2->SaveAs(path);

    h2_chi2_6->SetOption("colz");
    //h2_chi2_6->SetMinimum(0);
    //h2_chi2_6->SetMaximum(100000);
    h2_chi2_6->SetStats(kFALSE);
    h2_chi2_6->GetXaxis()->SetTitle(var3->GetName());
    h2_chi2_6->GetYaxis()->SetTitle(var4->GetName());
    h2_chi2_6->Draw();
    h2_chi2_6->Write();
    path = "plots/chi2map_";
    path += var3->GetName();
    path += "_";
    path += var4->GetName();
    path += ".png";
    c2->SaveAs(path);

    file->Close();

    delete h1_chi2;
    delete h2_chi2_1;
    delete h2_chi2_2;
    delete h2_chi2_3;
    delete h2_chi2_4;
    delete h2_chi2_5;
    delete h2_chi2_6;
    delete c2;

    return mychi2;
}

void FitterTrans::GetRecoveredParameters(Int_t& numParameters, Double_t** recoveredParameters)
{
    RooRealVar* chi2red  = new RooRealVar("chi2red","reduced chi^2",0);
    chi2red->setVal(chi2Var->getVal()/(dataSet_binned->numEntries()-numFitParameters));

    numParameters = 17;
    Double_t* parameters = new Double_t[numParameters];

    parameters[0] = chi2red->getVal();
    parameters[1] = ap->getVal();
    parameters[2] = ap->getError();
    parameters[3] = apa->getVal();
    parameters[4] = apa->getError();
    parameters[5] = a0->getVal();
    parameters[6] = a0->getError();
    parameters[7] = a0a->getVal();
    parameters[8] = a0a->getError();
    parameters[9] = at->getVal();
    parameters[10] = at->getPropagatedError(*result);
    parameters[11] = ata->getVal();
    parameters[12] = ata->getError();
    parameters[13] = par_input[0];
    parameters[14] = par_input[1];
    parameters[15] = par_input[2];
    parameters[16] = par_input[3];

    *recoveredParameters = parameters;
}

RooDataHist* FitterTrans::GetBinnedDataSet()
{
    if(dataSet_binned == NULL)
        CreateBinnedDataSet("a");

    return dataSet_binned;
}

void FitterTrans::FixAllParameters()
{
    RooRealVar* rooPar = 0;
    TIterator* parIter = parameters->createIterator();
    while(rooPar = (RooRealVar*)parIter->Next())
        rooPar->setConstant();

}

void FitterTrans::FixParameter(const char* par)
{
    RooRealVar* rooPar = 0;
    rooPar = (RooRealVar*)parameters->find(par);
    if(rooPar != 0)
        rooPar->setConstant();
}

void FitterTrans::FreeParameter(const char* par)
{
    RooRealVar* rooPar = 0;
    rooPar = (RooRealVar*)parameters->find(par);
    if(rooPar != 0)
        rooPar->setConstant(kFALSE);
}

