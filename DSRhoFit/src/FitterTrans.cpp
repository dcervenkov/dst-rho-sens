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

FitterTrans::FitterTrans(RooDataSet* outer_dataSet, Double_t* outer_par_input)
{
    gPluginMgr = new TPluginManager;
    gPluginMgr->AddHandler("ROOT::Math::Minimizer", "Minuit2", "Minuit2Minimizer", "Minuit2", "Minuit2Minimizer(const char *)");

    dataSet = outer_dataSet;
    for(int i = 0; i < 11; i++)
        par_input[i] = outer_par_input[i];

    chi2Var = 0;
    result = 0;

    tht_bins = 20;
    thb_bins = 20;
    phit_bins = 40;
    dt_bins = 80;

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
    apa = new RooRealVar("apa","apa",par_input[1],0,2*PI);
    apr = new RooFormulaVar("apr","ap*cos(apa)",RooArgSet(*ap,*apa));
    api = new RooFormulaVar("api","ap*sin(apa)",RooArgSet(*ap,*apa));
    a0 = new RooRealVar("a0","a0",par_input[2],0.8,1);
    a0a = new RooRealVar("a0a","a0a",0);
    a0r = new RooFormulaVar("a0r","a0*cos(a0a)",RooArgSet(*a0,*a0a));
    a0i = new RooFormulaVar("a0i","a0*sin(a0a)",RooArgSet(*a0,*a0a));
    at = new RooFormulaVar("at","sqrt(1-ap*ap-a0*a0)",RooArgSet(*ap,*a0));
    ata = new RooRealVar("ata","ata",par_input[3],0,2*PI);
    atr = new RooFormulaVar("atr","at*cos(ata)",RooArgSet(*at,*ata));
    ati = new RooFormulaVar("ati","at*sin(ata)",RooArgSet(*at,*ata));

    ap0r = new RooFormulaVar("ap0r","ap*a0*cos(-apa+a0a)",RooArgSet(*ap,*apa,*a0,*a0a));
    a0ti = new RooFormulaVar("a0ti","a0*at*sin(-a0a+ata)",RooArgSet(*a0,*a0a,*at,*ata));
    apti = new RooFormulaVar("apti","ap*at*sin(-apa+ata)",RooArgSet(*ap,*apa,*at,*ata));

    /// Time-dep additional vars

    dm = new RooRealVar("dm","dm",0.507e12);
    phiw = new RooRealVar("phiw","phiw",par_input[4],0,2*PI);

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

    pdf_a = new RooGenericPdf("pdf_a","pdf_a",formula_a,*varSet_a);
    pdf_b = new RooGenericPdf("pdf_b","pdf_b",formula_b,*varSet_b);
    pdf_ab = new RooGenericPdf("pdf_ab","pdf_ab",formula_ab,*varSet_ab);
    pdf_bb = new RooGenericPdf("pdf_bb","pdf_bb",formula_bb,*varSet_bb);

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
    numFitParameters = (fitParameters->selectByAttrib("Constant",kFALSE))->getSize();

    pdfFormula =   "ap*ap*2*sin(tht)*sin(tht)*sin(tht)*sin(thb)*sin(thb)*sin(thb)*sin(phit)*sin(phit)+\
                                at*at*2*cos(tht)*cos(tht)*sin(tht)*sin(thb)*sin(thb)*sin(thb)+\
                                a0*a0*4*sin(tht)*sin(tht)*sin(tht)*cos(thb)*cos(thb)*sin(thb)*cos(phit)*cos(phit)+\
                                sqrt(2)*ap0r*sin(tht)*sin(tht)*sin(tht)*sin(2*thb)*sin(thb)*sin(2*phit)-\
                                sqrt(2)*a0ti*sin(2*tht)*sin(tht)*sin(2*thb)*sin(thb)*cos(phit)-\
                                2*apti*sin(2*tht)*sin(tht)*sin(thb)*sin(thb)*sin(thb)*sin(phit)";

    pdf = new RooGenericPdf("pdf","Generic PDF",pdfFormula,*varSet);

    myPdf_a = new DSRhoPDF("myPdf","myPdf","b",*tht,*thb,*phit,*dt,*ap,*apa,*a0,*ata,*phiw,*rp,*r0,*rt,*sp,*s0,*st);
    myPdf_b = new DSRhoPDF("myPdf","myPdf","b",*tht,*thb,*phit,*dt,*ap,*apa,*a0,*ata,*phiw,*rp,*r0,*rt,*sp,*s0,*st);
}

FitterTrans::~FitterTrans()
{
    //dtor
}

Int_t FitterTrans::Fit()
{
    TPluginManager* gPluginMgr = new TPluginManager;
    gPluginMgr->AddHandler("ROOT::Math::Minimizer", "Minuit2", "Minuit2Minimizer", "Minuit2", "Minuit2Minimizer(const char *)");
    //result = pdf->fitTo(*dataSet,RooFit::Save(),RooFit::Timer(true),RooFit::Minimizer("Minuit2"));//,RooFit::NumCPU(2));
    result = simPdf->fitTo(*dataSet,RooFit::Save(),RooFit::Timer(true),RooFit::Minimizer("Minuit2"));//,RooFit::NumCPU(2));
}

void FitterTrans::CreateBinnedDataSet()
{
    tht->setBins(tht_bins);
    thb->setBins(thb_bins);
    phit->setBins(phit_bins);
    dt->setBins(dt_bins);

    RooRandom::randomGenerator()->SetSeed(0);

    /// Create a binned dataSet which is needed for chi2 calculation
    dataSet_binned = new RooDataHist("dataSet_binned","dataSet_binned",RooArgSet(*tht,*thb,*phit),*dataSet);
    //RooDataHist* dataSet_binned = pdf->generateBinned(RooArgSet(var1,var2,var3),dataSet->numEntries(),kFALSE);
}

Int_t FitterTrans::ComputeChi2()
{
    if(dataSet_binned == NULL)
        CreateBinnedDataSet();

	chi2Var = new RooChi2Var("chi2Var","chi2Var",*pdf,*dataSet_binned);

	RooRealVar* ndof     = new RooRealVar("ndof","number of degrees of freedom",0);
	RooRealVar* chi2red  = new RooRealVar("chi2red","reduced chi^2",0);
	RooRealVar* prob     = new RooRealVar("prob","prob(chi2,ndof)",0);

	ndof->setVal(dataSet_binned->numEntries()-numFitParameters);
	chi2red->setVal(chi2Var->getVal()/ndof->getVal()) ;
	prob->setVal(TMath::Prob(chi2Var->getVal(),static_cast<int>(ndof->getVal())));

	printf("chi2 = %f\nndof = %f\nchi2red = %f\nprob = %f\n",chi2Var->getVal(),ndof->getVal(),chi2red->getVal(),prob->getVal());
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
        CreateBinnedDataSet();

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

