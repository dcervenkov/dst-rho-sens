#include <stdio.h>
//#include <stdlib.h>

#include "TROOT.h"
#include "TStyle.h"
#include "TApplication.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TNtuple.h"
#include "TTree.h"
#include "TRandom1.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TStopwatch.h"
#include "TMath.h"

#include "RooFit.h"
#include "RooRealVar.h"
#include "RooPlot.h"
#include "RooDataSet.h"
#include "RooDataHist.h"
#include "RooArgSet.h"
#include "RooGenericPdf.h"
#include "RooFitResult.h"
#include "RooChi2Var.h"
#include "RooGaussian.h"
#include "RooProdPdf.h"
#include "RooNLLVar.h"
#include "RooMinuit.h"
#include "RooPolynomial.h"
#include "RooAddPdf.h"


#include "DSRhoFit.h"
#include "Constants.h"
#include "ASSERT.h"

#define DEBUG
//#define VERBOSE
//#define GRAPHIC
#define HELICITY
//#define TRANSVERSITY

#ifdef HELICITY
TH1D *hChi = new TH1D("hChi", "Chi distribution", 50, 0, 2*PI);
TH1D *hThA = new TH1D("hThA", "Theta A distribution", 50, 0, PI);
TH1D *hThB = new TH1D("hThB", "Theta B distribution", 50, 0, PI);
TH2D *hG = new TH2D("hG","Gamma distribution", 30, 0, PI, 30, 0, PI);
#endif

#ifdef TRANSVERSITY
TH1D *hPhiT = new TH1D("hPhiT", "Phi T distribution", 50, -PI, PI);
TH1D *hThT = new TH1D("hThT", "Theta T distribution", 50, 0, PI);
TH1D *hThB = new TH1D("hThB", "Theta B distribution", 50, 0, PI);
#endif


int main(int argc, char* argv[])
{
    #ifdef GRAPHIC
    TApplication* rootapp = new TApplication("DSRhoSens",&argc,argv); //For graphic-output apps only
    /**
     * When using TApplication it removes arguments it "handles" from
     * the argument array. E.g. -b, -x, -q, --help, <dir>, <file>, etc.
     * For more info read TApplication's GetOptions function help.
     * The solution is to use rootapp->Argc() and rootapp->Argv(i).
     * The next few lines are for compatibility of GRAPHIC vs. non-GRAPHIC -
     * they recreate the original argc and argv when using GRAPHIC.
     **/
    argc = rootapp->Argc();
    for(int i = 0; i < argc; i++)
    {
        argv[i] = rootapp->Argv(i);
    }
    #endif

    inputFile = argv[1];
    outputFile = argv[2];

    TStopwatch timer;
    timer.Start();

    #ifdef HELICITY
    RooRealVar tha("tha","tha",0,PI);
    RooRealVar thb("thb","thb",0,PI);
    RooRealVar chi("chi","chi",0,2*PI);

    RooDataSet* dataSet = new RooDataSet("data","data",RooArgList(tha,thb,chi));
    #endif

    #ifdef TRANSVERSITY
    RooRealVar phit("phit","phit",-PI,PI);
    RooRealVar tht("tht","tht",0,PI);
    RooRealVar thb("thb","thb",0,PI);
    RooRealVar ap2("ap2","ap2",0,1);
    RooDataSet* dataSet = new RooDataSet("data","data",RooArgSet(phit,tht,thb));
    RooGenericPdf* pdf = new RooGenericPdf("pdf","Generic PDF","ap2*sin(phit)*sin(phit)*sin(tht)*sin(tht)*sin(tht)*sin(thb)*sin(thb)*sin(thb)",RooArgSet(phit,tht,thb,ap2));
    #endif

    dataSet = RooDataSet::read(inputFile,RooArgList(tha,thb,chi));
    Fit(dataSet,tha,thb,chi);

    timer.Stop();
    timer.Print();

    #ifdef GRAPHIC
//    TCanvas* cc = new TCanvas("cc","Canvas",800,600);

        #ifdef HELICITY
//    	RooPlot* frame1 = chi.frame();
//        RooPlot* frame2 = tha.frame();
//        RooPlot* frame3 = thb.frame();
//        binnedDataSet->plotOn(frame1);
//        pdf->plotOn(frame1);
//        binnedDataSet->plotOn(frame2);
//        pdf->plotOn(frame2);
//        binnedDataSet->plotOn(frame3);
//        pdf->plotOn(frame3);
//        cc->Divide(2,2);
//        cc->cd(1);
//        frame1->Draw();
//        cc->cd(2);
//    	frame2->Draw();
//    	cc->cd(3);
//    	frame3->Draw();
//    	printf("chi1red = %f\nchi2red = %f\nchi3red = %f\n",frame1->chiSquare(5),frame2->chiSquare(5),frame3->chiSquare(5));
        #endif

        #ifdef TRANSVERSITY
        cc->Divide(2,2);
        cc->cd(1);
        hPhiT->Draw();
        cc->cd(2);
        hThT->Draw();
        cc->cd(3);
        hThB->Draw();
        cc->cd(4);

        TCanvas* c1 = new TCanvas("c1","c1");
        pdf->fitTo(*dataSet);
        RooPlot* xframe = phit.frame();
        dataSet->plotOn(xframe);
        pdf->plotOn(xframe);
        xframe->Draw();

        TCanvas* c2 = new TCanvas("c2","c2");
        RooPlot* yframe = tht.frame();
        dataSet->plotOn(yframe);
        pdf->plotOn(yframe);
        yframe->Draw();

        TCanvas* c3 = new TCanvas("c3","c3");
        RooPlot* zframe = thb.frame();
        dataSet->plotOn(zframe);
        pdf->plotOn(zframe);
        zframe->Draw();
        #endif

    printf("\nProgram execution has finished.\n");
    rootapp->Run(); //For graphic-output apps only
    #endif
    return 0;
}

int Fit(RooDataSet* dataSet, RooRealVar& tha, RooRealVar& thb, RooRealVar& chi)
{

    /**
     * The f one gets from fitting with the following PDF is actually (|H+|^2 + |H-|^2)
     * viz formula (35) in BN419. Equvivalentely (1-f) = |H0|^2.
     **/
//    RooRealVar f("f","fraction",0.5,0,1);
//    RooRealVar c1("c1","c1 var",0.05,0,0.5);
//    RooRealVar c2("c2","c2 var",0.05,0,0.5);

    //RooRealVar hp("hp","hp",0.152,0,0.3);
    RooRealVar hp("hp","hp",0.152);
    //RooRealVar hpa("hpa","hpa",1.47,0,2*PI);
    RooRealVar hpa("hpa","hpa",1.47);
    RooFormulaVar hpr("hpr","hp*cos(hpa)",RooArgSet(hp,hpa));
    RooFormulaVar hpi("hpi","hp*sin(hpa)",RooArgSet(hp,hpa));
    //RooRealVar h0("h0","h0",0.936,0.8,1);
    RooRealVar h0("h0","h0",0.936);
    //RooRealVar h0a("h0a","h0a",0,0,2*PI);
    RooRealVar h0a("h0a","h0a",0);
    RooFormulaVar h0r("h0r","h0*cos(h0a)",RooArgSet(h0,h0a));
    RooFormulaVar h0i("h0i","h0*sin(h0a)",RooArgSet(h0,h0a));
    RooFormulaVar hm("hm","sqrt(1-hp*hp-h0*h0)",RooArgSet(hp,h0));
    //RooRealVar hma("hma","hma",0.19,0,2*PI);
    RooRealVar hma("hma","hma",0.19);
    RooFormulaVar hmr("hmr","hm*cos(hma)",RooArgSet(hm,hma));
    RooFormulaVar hmi("hmi","hm*sin(hma)",RooArgSet(hm,hma));

    RooFormulaVar hptr("hptr","hp*hm*cos(hpa-hma)",RooArgSet(hp,hpa,hm,hma));
    RooFormulaVar hpti("hpti","hp*hm*sin(hpa-hma)",RooArgSet(hp,hpa,hm,hma));

    /**
     * hat: h - helicity, a - addition, t - transverse
     * hst: s - subtraction
     */
    RooFormulaVar hatr("hatr","hp*cos(hpa) + hm*cos(hma)",RooArgSet(hp,hpa,hm,hma));
    RooFormulaVar hati("hati","hp*sin(hpa) + hm*sin(hma)",RooArgSet(hp,hpa,hm,hma));
    RooFormulaVar hat("hat","sqrt(hatr*hatr + hati*hati)",RooArgSet(hatr,hati));
    RooFormulaVar hata("hata","atan2(hati,hatr)",RooArgSet(hati,hatr));

    RooFormulaVar hstr("hstr","hp*cos(hpa) - hm*cos(hma)",RooArgSet(hp,hpa,hm,hma));
    RooFormulaVar hsti("hsti","hp*sin(hpa) - hm*sin(hma)",RooArgSet(hp,hpa,hm,hma));
    RooFormulaVar hst("hst","sqrt(hstr*hstr + hsti*hsti)",RooArgSet(hstr,hsti));
    RooFormulaVar hsta("hsta","atan2(hsti,hstr)",RooArgSet(hsti,hstr));

    RooArgSet varSet(tha,thb,hp,hpa,hpr,hpi,h0,h0a,h0r);
    varSet.add(h0i);
    varSet.add(hm);
    varSet.add(hma);
    varSet.add(hmr);
    varSet.add(hmi);
    varSet.add(hat);
    varSet.add(hata);
    varSet.add(hst);
    varSet.add(hsta);
    varSet.add(chi);

    const char* pdfFormula =   "(hp*hp+hm*hm)*sin(tha)*sin(tha)*sin(tha)*sin(thb)*sin(thb)*sin(thb)+\
                                4*h0*h0*cos(tha)*cos(tha)*sin(tha)*cos(thb)*cos(thb)*sin(thb)+\
                                2*(hp*hm*cos(hpa-hma)*cos(2*chi)-hp*hm*sin(hpa-hma)*sin(2*chi))*sin(tha)*sin(tha)*sin(tha)*sin(thb)*sin(thb)*sin(thb)+\
                                (hat*h0*cos(hata-h0a)*cos(chi)-hst*h0*sin(hsta-h0a)*sin(chi))*sin(2*tha)*sin(tha)*sin(2*thb)*sin(thb)";

    RooGenericPdf* pdf = new RooGenericPdf("pdf","Generic PDF",pdfFormula,varSet);

    //RooGenericPdf* pdf_int_chi = new RooGenericPdf("pdf_int_chi","pdf_int_chi","f*sin(tha)*sin(tha)*sin(tha)*sin(thb)*sin(thb)*sin(thb)+(1-f)*4*cos(tha)*cos(tha)*sin(tha)*cos(thb)*cos(thb)*sin(thb)",RooArgSet(tha,thb,f));
    //RooGenericPdf* pdf_int_chi_thb = new RooGenericPdf("pdf_int_chi","pdf_int_chi","f*(4/3)*sin(tha)*sin(tha)*sin(tha)+(1-f)*(2/3)*4*cos(tha)*cos(tha)*sin(tha)",RooArgSet(tha,f));
    //RooGenericPdf* pdf_int_tha_thb = new RooGenericPdf("pdf_int_tha_thb","pdf_int_tha_thb","1+2*(c1*cos(2*chi)-c2*sin(2*chi))",RooArgSet(c1,c2,chi));


    //RooDataSet* data = pdf->generate(RooArgSet(tha,thb),10000);

//    RooFitResult* result = pdf_int_chi->fitTo(*dataSet,RooFit::Save(),RooFit::Timer(true));//,RooFit::NumCPU(2));
//    result->Print();
//    RooPlot* frame1 = tha.frame(RooFit::Bins(100));
//    dataSet->plotOn(frame1);
//    pdf_int_chi->plotOn(frame1);
//    pdf_int_chi->paramOn(frame1);
//    frame1->Draw();
//    printf("(int_chi) chi^2 = %f\n",frame1->chiSquare(1));

//    RooFormulaVar h0fit("h0fit","sqrt(1-f)",f);
//    h0 = h0fit;
//    h0fit.Print();
//    h0.Print();
//    h0.setConstant(true);
    //c1 = 0.0138;
    //c2 = 0.0462;
    //c1.setConstant(true);
    //c2.setConstant(true);

//    RooFitResult* result2 = pdf_int_tha_thb->fitTo(*dataSet,RooFit::Save(),RooFit::Timer(true));
//    result2->Print();
//    c1.Print();
//    c2.Print();
//    h0.Print();
    //printf("h0 = %f\n",h0.getVal());
//    c1.setConstant(true);
//    c2.setConstant(true);
//
//	RooPlot* frame2 = chi.frame(RooFit::Bins(100));
//	dataSet->plotOn(frame2);
//	pdf_int_tha_thb->plotOn(frame2);
//	pdf_int_tha_thb->paramOn(frame2);
//	frame2->Draw();
//	printf("(int_tha_thb) chi^2 = %f\n",frame2->chiSquare(2));

	//RooRealVar c1err("c1err","c1err",c1.getError());
	//RooRealVar c2err("c2err","c2err",c2.getError());
	//RooGaussian* c1_constraint = new RooGaussian("c1_constraint","c1_constraint",hptr,c1,c1err);
	//RooGaussian* c2_constraint = new RooGaussian("c2_constraint","c2_constraint",hpti,c2,c2err);

    //RooFitResult* result3 = pdf->fitTo(*dataSet,RooFit::ExternalConstraints(RooArgSet(*c1_constraint,*c2_constraint)),RooFit::Save(),RooFit::Timer(true));
    //RooFitResult* result3 = pdf->fitTo(*dataSet,RooFit::Save(),RooFit::Timer(true));
    //result3->Print();

    hp.Print();
    hpa.Print();
    h0.Print();
    h0a.Print();
    hm.Print();
    hma.Print();
//    c1.Print();
//    c2.Print();
//    hptr.Print();
//    hpti.Print();

    chi.setBins(40);
	tha.setBins(20);
	thb.setBins(20);

	RooDataHist* binnedDataSet = new RooDataHist("binnedDataSet","binnedDataSet",RooArgSet(chi,tha,thb),*dataSet);
	RooChi2Var chi2Var("chi2Var","chi2Var",*pdf,*binnedDataSet);

	RooRealVar* _chi2     = new RooRealVar("chi2","chi^2",0) ;
	RooRealVar* _ndof     = new RooRealVar("ndof","number of degrees of freedom",0) ;
	RooRealVar* _chi2red  = new RooRealVar("chi2red","reduced chi^2",0) ;
	RooRealVar* _prob     = new RooRealVar("prob","prob(chi2,ndof)",0) ;

	_chi2->setVal(chi2Var.getVal()) ;
	//RooArgSet* floatPars = (RooArgSet*) fitParams()->selectByAttrib("Constant",kFALSE);
	// Should use the above line instead of 5
	_ndof->setVal(binnedDataSet->numEntries()-5) ;
	_chi2red->setVal(_chi2->getVal()/_ndof->getVal()) ;
	_prob->setVal(TMath::Prob(_chi2->getVal(),static_cast<int>(_ndof->getVal())));

	printf("chi2 = %f\nndof = %f\nchi2red = %f\nprob = %f\n",_chi2->getVal(),_ndof->getVal(),_chi2red->getVal(),_prob->getVal());

    FILE* pFile;
    pFile = fopen (outputFile,"w");
    if (pFile!=NULL)
    {
        fprintf(pFile,"%f\n",_chi2->getVal());
        fclose (pFile);
    }
    else
        printf("ERROR: couldn't open file %s for writing!\n",outputFile);

    return 0;
}