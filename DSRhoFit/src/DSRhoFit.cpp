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
#include "RooRandom.h"

#include "DSRhoFit.h"
#include "Constants.h"
#include "ASSERT.h"

#define DEBUG
//#define VERBOSE
//#define GRAPHIC
//#define HELICITY
#define TRANSVERSITY

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
    RooDataSet* dataSet = new RooDataSet("data","data",RooArgList(tha,thb,chi));
    dataSet = RooDataSet::read(inputFile,RooArgList(tha,thb,chi));
    FitHel(dataSet,tha,thb,chi);
    #endif

    #ifdef TRANSVERSITY
    RooDataSet* dataSet = new RooDataSet("data","data",RooArgSet(tht,thb,phit));
    dataSet = RooDataSet::read(inputFile,RooArgList(tht,thb,phit));
    FitTrans(dataSet,tht,thb,phit);
    #endif



    timer.Stop();
    timer.Print();

    #ifdef GRAPHIC
//    TCanvas* cc = new TCanvas("cc","Canvas",800,600);

        #ifdef HELICITY
//    	RooPlot* frame1 = chi.frame();
//        RooPlot* frame2 = tha.frame();
//        RooPlot* frame3 = thb.frame();
//        dataSet_binned->plotOn(frame1);
//        pdf->plotOn(frame1);
//        dataSet_binned->plotOn(frame2);
//        pdf->plotOn(frame2);
//        dataSet_binned->plotOn(frame3);
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

int FitHel(RooDataSet* dataSet, RooRealVar& tha, RooRealVar& thb, RooRealVar& chi)
{

    /**
     * The f one gets from fitting with the following PDF is actually (|H+|^2 + |H-|^2)
     * viz formula (35) in BN419. Equvivalentely (1-f) = |H0|^2.
     **/
    RooRealVar f("f","fraction",0.5,0,1);
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

    RooGenericPdf* pdf_int_chi = new RooGenericPdf("pdf_int_chi","pdf_int_chi","f*sin(tha)*sin(tha)*sin(tha)*sin(thb)*sin(thb)*sin(thb)+(1-f)*4*cos(tha)*cos(tha)*sin(tha)*cos(thb)*cos(thb)*sin(thb)",RooArgSet(tha,thb,f));
    //RooGenericPdf* pdf_int_chi_thb = new RooGenericPdf("pdf_int_chi","pdf_int_chi","f*(4/3)*sin(tha)*sin(tha)*sin(tha)+(1-f)*(2/3)*4*cos(tha)*cos(tha)*sin(tha)",RooArgSet(tha,f));
    //RooGenericPdf* pdf_int_tha_thb = new RooGenericPdf("pdf_int_tha_thb","pdf_int_tha_thb","1+2*(c1*cos(2*chi)-c2*sin(2*chi))",RooArgSet(c1,c2,chi));

    chi.setBins(chi_bins);
	tha.setBins(tha_bins);
	thb.setBins(thb_bins);

    //printf("(int_thb_chi) chi^2 = %f\n",frame1->chiSquare(5));



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


//	printf("(int_tha_chi) chi^2 = %f\n",frame2->chiSquare(5));


//	printf("(int_tha_thb) chi^2 = %f\n",frame3->chiSquare(5));

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

//    tha = PI/2;
//    thb = PI/2;
//    chi = PI/4;
//    printf("*** pdf value = %f\n",pdf->getVal());


    RooRandom::randomGenerator()->SetSeed(0);

	RooDataHist* dataSet_binned = new RooDataHist("dataSet_binned","dataSet_binned",RooArgSet(chi,tha,thb),*dataSet);
    //RooDataHist* dataSet_binned = pdf->generateBinned(RooArgSet(tha,thb,chi),dataSet->numEntries(),kFALSE);
    RooDataHist* pdf_binned = pdf->generateBinned(RooArgSet(tha,thb,chi),dataSet->numEntries(),kTRUE);

	RooChi2Var chi2Var("chi2Var","chi2Var",*pdf,*dataSet_binned);

	RooRealVar* _chi2     = new RooRealVar("chi2","chi^2",0) ;
	RooRealVar* _ndof     = new RooRealVar("ndof","number of degrees of freedom",0) ;
	RooRealVar* _chi2red  = new RooRealVar("chi2red","reduced chi^2",0) ;
	RooRealVar* _prob     = new RooRealVar("prob","prob(chi2,ndof)",0) ;

	_chi2->setVal(chi2Var.getVal());
	//RooArgSet* floatPars = (RooArgSet*) fitParams()->selectByAttrib("Constant",kFALSE);
	// Should use the above line instead of 5
	_ndof->setVal(dataSet_binned->numEntries()-5) ;
	_chi2red->setVal(_chi2->getVal()/_ndof->getVal()) ;
	_prob->setVal(TMath::Prob(_chi2->getVal(),static_cast<int>(_ndof->getVal())));

	printf("chi2 = %f\nndof = %f\nchi2red = %f\nprob = %f\n",_chi2->getVal(),_ndof->getVal(),_chi2red->getVal(),_prob->getVal());
    //WriteToFile(_chi2->getVal());


//    RooFitResult* result = pdf_int_chi->fitTo(*dataSet,RooFit::Save(),RooFit::Timer(true));//,RooFit::NumCPU(2));
//    result->Print();

    TFile* file = new TFile("plots/frames.root","RECREATE");
    TCanvas* c1 = new TCanvas("c1","c1",1000,800);

    RooPlot* frame1 = tha.frame();
    pdf_binned->plotOn(frame1,RooFit::MarkerColor(2),RooFit::Name("pdf"));
    dataSet_binned->plotOn(frame1,RooFit::MarkerColor(3),RooFit::Name("gen data"));
    pdf->plotOn(frame1);
    //dataSet->plotOn(frame1,RooFit::Name("data"));
    frame1->Draw();
    frame1->Write();
    c1->SaveAs("plots/frame1.png");

    RooPlot* frame2 = thb.frame();
    pdf_binned->plotOn(frame2,RooFit::MarkerColor(2),RooFit::Name("pdf"));
    dataSet_binned->plotOn(frame2,RooFit::MarkerColor(3),RooFit::Name("gen data"));
    //dataSet->plotOn(frame2,RooFit::Name("data"));
	pdf->plotOn(frame2);
	frame2->Draw();
    c1->SaveAs("plots/frame2.png");

    RooPlot* frame3 = chi.frame();
    pdf_binned->plotOn(frame3,RooFit::MarkerColor(2),RooFit::Name("pdf"));
    dataSet_binned->plotOn(frame3,RooFit::MarkerColor(3),RooFit::Name("gen data"));
    //dataSet->plotOn(frame3,RooFit::Name("data"));
	pdf->plotOn(frame3);
	frame3->Draw();
    c1->SaveAs("plots/frame3.png");

    file->Close();

    TH2* h2 = dynamic_cast<TH2*>(dataSet_binned->createHistogram("h2",tha,RooFit::YVar(thb)));
    h2->SetOption("colz");
    h2->Draw();
    c1->SaveAs("plots/tha_thb_data.png");

    h2 = dynamic_cast<TH2*>(pdf_binned->createHistogram("h2",tha,RooFit::YVar(thb)));
    h2->SetOption("colz");
    h2->Draw();
    c1->SaveAs("plots/tha_thb_pdf.png");

    h2 = dynamic_cast<TH2*>(pdf_binned->createHistogram("h2",tha,RooFit::YVar(chi)));
    h2->SetOption("colz");
    h2->Draw();
    c1->SaveAs("plots/tha_chi_pdf.png");

    h2 = dynamic_cast<TH2*>(dataSet_binned->createHistogram("h2",tha,RooFit::YVar(chi)));
    h2->SetOption("colz");
    h2->Draw();
    c1->SaveAs("plots/tha_chi_data.png");

    h2 = dynamic_cast<TH2*>(pdf_binned->createHistogram("h2",thb,RooFit::YVar(chi)));
    h2->SetOption("colz");
    h2->Draw();
    c1->SaveAs("plots/thb_chi_pdf.png");

    h2 = dynamic_cast<TH2*>(dataSet_binned->createHistogram("h2",thb,RooFit::YVar(chi)));
    h2->SetOption("colz");
    h2->Draw();
    c1->SaveAs("plots/thb_chi_data.png");


    GetChi2(dataSet_binned,pdf_binned);

    return 0;
}

int FitTrans(RooDataSet* dataSet, RooRealVar& tht, RooRealVar& thb, RooRealVar& phit)
{

    RooRealVar f("f","fraction",0.5,0,1);
//    RooRealVar c1("c1","c1 var",0.05,0,0.5);
//    RooRealVar c2("c2","c2 var",0.05,0,0.5);

    RooRealVar ap("ap","ap",0.275);
    RooRealVar apa("apa","apa",0.57);
    RooFormulaVar apr("apr","ap*cos(apa)",RooArgSet(ap,apa));
    RooFormulaVar api("api","ap*sin(apa)",RooArgSet(ap,apa));
    RooRealVar a0("a0","a0",0.936);
    RooRealVar a0a("a0a","a0a",0);
    RooFormulaVar a0r("a0r","a0*cos(a0a)",RooArgSet(a0,a0a));
    RooFormulaVar a0i("a0i","a0*sin(a0a)",RooArgSet(a0,a0a));
    RooFormulaVar at("at","sqrt(1-ap*ap-a0*a0)",RooArgSet(ap,a0));
    RooRealVar ata("ata","ata",2.84);
    RooFormulaVar atr("atr","at*cos(ata)",RooArgSet(at,ata));
    RooFormulaVar ati("ati","at*sin(ata)",RooArgSet(at,ata));

    RooFormulaVar ap0r("ap0r","ap*a0*cos(-apa+a0a)",RooArgSet(ap,apa,a0,a0a));
    RooFormulaVar a0ti("a0ti","a0*at*sin(-a0a+ata)",RooArgSet(a0,a0a,at,ata));
    RooFormulaVar apti("apti","ap*at*sin(-apa+ata)",RooArgSet(ap,apa,at,ata));

    RooArgSet varSet(tht,thb,phit,ap,apa,a0,a0a,at,ata);
    varSet.add(ap0r);
    varSet.add(a0ti);
    varSet.add(apti);

    const char* pdfFormula =   "ap*ap*2*sin(tht)*sin(tht)*sin(tht)*sin(thb)*sin(thb)*sin(thb)*sin(phit)*sin(phit)+\
                                at*at*2*cos(tht)*cos(tht)*sin(tht)*sin(thb)*sin(thb)*sin(thb)+\
                                a0*a0*4*sin(tht)*sin(tht)*sin(tht)*cos(thb)*cos(thb)*sin(thb)*cos(phit)*cos(phit)+\
                                sqrt(2)*ap0r*sin(tht)*sin(tht)*sin(tht)*sin(2*thb)*sin(thb)*sin(2*phit)-\
                                sqrt(2)*a0ti*sin(2*tht)*sin(tht)*sin(2*thb)*sin(thb)*cos(phit)-\
                                2*apti*sin(2*tht)*sin(tht)*sin(thb)*sin(thb)*sin(thb)*sin(phit)";

    RooGenericPdf* pdf = new RooGenericPdf("pdf","Generic PDF",pdfFormula,varSet);

    //RooGenericPdf* pdf_int_phit = new RooGenericPdf("pdf_int_phit","pdf_int_phit","f*sin(tha)*sin(tha)*sin(tha)*sin(thb)*sin(thb)*sin(thb)+(1-f)*4*cos(tha)*cos(tha)*sin(tha)*cos(thb)*cos(thb)*sin(thb)",RooArgSet(tht,thb,phit));
    //RooGenericPdf* pdf_int_chi_thb = new RooGenericPdf("pdf_int_chi","pdf_int_chi","f*(4/3)*sin(tha)*sin(tha)*sin(tha)+(1-f)*(2/3)*4*cos(tha)*cos(tha)*sin(tha)",RooArgSet(tha,f));
    RooGenericPdf* pdf_int_thb_phit = new RooGenericPdf("pdf_int_tht_phit","pdf_int_tht_phit","(1-f)*sin(tht)*sin(tht)*sin(tht)+f*2*cos(tht)*cos(tht)*sin(tht)",RooArgSet(f,tht,thb,phit));

    tht.setBins(tht_bins);
    thb.setBins(thb_bins);
    phit.setBins(phit_bins);

    //printf("(int_thb_chi) chi^2 = %f\n",frame1->chiSquare(5));
    //printf("(int_tha_chi) chi^2 = %f\n",frame2->chiSquare(5));
    //printf("(int_tha_thb) chi^2 = %f\n",frame3->chiSquare(5));


    RooRandom::randomGenerator()->SetSeed(0);

	//RooDataHist* dataSet_binned = new RooDataHist("dataSet_binned","dataSet_binned",RooArgSet(tht,thb,phit),*dataSet);
    RooDataHist* dataSet_binned = pdf->generateBinned(RooArgSet(tht,thb,phit),dataSet->numEntries(),kFALSE);
    RooDataHist* pdf_binned = pdf->generateBinned(RooArgSet(tht,thb,phit),dataSet->numEntries(),kTRUE);

	RooChi2Var chi2Var("chi2Var","chi2Var",*pdf,*dataSet_binned);

	RooRealVar* _chi2     = new RooRealVar("chi2","chi^2",0) ;
	RooRealVar* _ndof     = new RooRealVar("ndof","number of degrees of freedom",0) ;
	RooRealVar* _chi2red  = new RooRealVar("chi2red","reduced chi^2",0) ;
	RooRealVar* _prob     = new RooRealVar("prob","prob(chi2,ndof)",0) ;

	_chi2->setVal(chi2Var.getVal());
	//RooArgSet* floatPars = (RooArgSet*) fitParams()->selectByAttrib("Constant",kFALSE);
	// Should use the above line instead of 5
	_ndof->setVal(dataSet_binned->numEntries()-5) ;
	_chi2red->setVal(_chi2->getVal()/_ndof->getVal()) ;
	_prob->setVal(TMath::Prob(_chi2->getVal(),static_cast<int>(_ndof->getVal())));

	printf("chi2 = %f\nndof = %f\nchi2red = %f\nprob = %f\n",_chi2->getVal(),_ndof->getVal(),_chi2red->getVal(),_prob->getVal());
    //WriteToFile(_chi2->getVal());


//    RooFitResult* result = pdf_int_chi->fitTo(*dataSet,RooFit::Save(),RooFit::Timer(true));//,RooFit::NumCPU(2));
//    result->Print();

    TFile* file = new TFile("plots/frames.root","RECREATE");
    TCanvas* c1 = new TCanvas("c1","c1",800,600);

    RooPlot* frame1 = tht.frame();
    pdf_binned->plotOn(frame1,RooFit::MarkerColor(2),RooFit::Name("pdf"));
    dataSet_binned->plotOn(frame1,RooFit::MarkerColor(3),RooFit::Name("gen data"));
    pdf->plotOn(frame1);
    //dataSet->plotOn(frame1,RooFit::Name("data"));
    frame1->Draw();
    frame1->Write();
    c1->SaveAs("plots/frame1.png");

    RooPlot* frame2 = thb.frame();
    pdf_binned->plotOn(frame2,RooFit::MarkerColor(2),RooFit::Name("pdf"));
    dataSet_binned->plotOn(frame2,RooFit::MarkerColor(3),RooFit::Name("gen data"));
    //dataSet->plotOn(frame2,RooFit::Name("data"));
	pdf->plotOn(frame2);
	frame2->Draw();
    c1->SaveAs("plots/frame2.png");

    RooPlot* frame3 = phit.frame();
    pdf_binned->plotOn(frame3,RooFit::MarkerColor(2),RooFit::Name("pdf"));
    dataSet_binned->plotOn(frame3,RooFit::MarkerColor(3),RooFit::Name("gen data"));
    //dataSet->plotOn(frame3,RooFit::Name("data"));
	pdf->plotOn(frame3);
	frame3->Draw();
    c1->SaveAs("plots/frame3.png");

    file->Close();

    TH2* h2 = dynamic_cast<TH2*>(dataSet_binned->createHistogram("h2",tht,RooFit::YVar(thb)));
    h2->SetOption("colz");
    h2->Draw();
    c1->SaveAs("plots/tht_thb_data.png");

    h2 = dynamic_cast<TH2*>(pdf_binned->createHistogram("h2",tht,RooFit::YVar(thb)));
    h2->SetOption("colz");
    h2->Draw();
    c1->SaveAs("plots/tht_thb_pdf.png");

    h2 = dynamic_cast<TH2*>(pdf_binned->createHistogram("h2",tht,RooFit::YVar(phit)));
    h2->SetOption("colz");
    h2->Draw();
    c1->SaveAs("plots/tht_phit_pdf.png");

    h2 = dynamic_cast<TH2*>(dataSet_binned->createHistogram("h2",tht,RooFit::YVar(phit)));
    h2->SetOption("colz");
    h2->Draw();
    c1->SaveAs("plots/tht_phit_data.png");

    h2 = dynamic_cast<TH2*>(pdf_binned->createHistogram("h2",thb,RooFit::YVar(phit)));
    h2->SetOption("colz");
    h2->Draw();
    c1->SaveAs("plots/thb_phit_pdf.png");

    h2 = dynamic_cast<TH2*>(dataSet_binned->createHistogram("h2",thb,RooFit::YVar(phit)));
    h2->SetOption("colz");
    h2->Draw();
    c1->SaveAs("plots/thb_phit_data.png");


    GetChi2Trans(dataSet_binned,pdf_binned);

    return 0;
}

Double_t GetChi2(RooDataHist* data, RooDataHist* pdf)
{
    Double_t mychi2 = 0;
    Double_t n = 0;
    Double_t v = 0;

    TH1F* hdchi = new TH1F("hdchi","hdchi",100,0,40);
    TH2F* hdchi2 = new TH2F("hdchi2","hdchi2",tha_bins,tha.getMin(),tha.getMax(),chi_bins,chi.getMin(),chi.getMax());

    /// I'm getting width of the first bin, because all bins are of equal width
    for(tha = tha.getBinWidth(0)/2; tha.getVal() < tha.getMax(); tha.setVal(tha.getVal()+tha.getBinWidth(0)))
    {
        for(thb = thb.getBinWidth(0)/2; thb.getVal() < thb.getMax(); thb.setVal(thb.getVal()+thb.getBinWidth(0)))
        {
            for(chi = chi.getBinWidth(0)/2; chi.getVal() < chi.getMax(); chi.setVal(chi.getVal()+chi.getBinWidth(0)))
            {
                n = data->weight(RooArgSet(tha,thb,chi),0);
                v = pdf->weight(RooArgSet(tha,thb,chi),0);
                hdchi->Fill((n-v)*(n-v)/v);
                hdchi2->Fill(tha.getVal(),chi.getVal(),(n-v)*(n-v)/v);
                mychi2 += (n-v)*(n-v)/v;
                //printf("pdf = %0.3f\tdata = %0.3f\n",pdf->weight(),data->weight());
            }
        }
    }

//    TCanvas* c2 = new TCanvas("c2","c2");
//    c2->Divide(2);
//    c2->cd(1);
//    hdchi->Draw();
//    c2->cd(2);
    hdchi2->SetOption("colz");
    //hdchi2->SetMinimum(0);
    //hdchi2->SetMaximum(100);
    hdchi2->SetStats(kFALSE);
    hdchi2->SaveAs(outputFile);
//    hdchi2->Draw();
//    c2->SaveAs("plots/dchi.root");
//    c2->SaveAs("plots/dchi.png");
    printf("my chi2 = %f\n",mychi2);
    return 3; //mychi2;
}

Double_t GetChi2Trans(RooDataHist* data, RooDataHist* pdf)
{
    Double_t mychi2 = 0;
    Double_t n = 0;
    Double_t v = 0;

    TH1F* hdchi = new TH1F("hdchi","hdchi",100,0,40);
    TH2F* hdchi2 = new TH2F("hdchi2","hdchi2",tht_bins,tht.getMin(),tht.getMax(),phit_bins,phit.getMin(),phit.getMax());

    /// I'm getting width of the first bin, because all bins are of equal width
    for(tht = tht.getBinWidth(0)/2; tht.getVal() < tht.getMax(); tht.setVal(tht.getVal()+tht.getBinWidth(0)))
    {
        for(thb = thb.getBinWidth(0)/2; thb.getVal() < thb.getMax(); thb.setVal(thb.getVal()+thb.getBinWidth(0)))
        {
            for(phit = phit.getMin()+phit.getBinWidth(0)/2; phit.getVal() < phit.getMax(); phit.setVal(phit.getVal()+phit.getBinWidth(0)))
            {
                n = data->weight(RooArgSet(tht,thb,phit),0);
                v = pdf->weight(RooArgSet(tht,thb,phit),0);
                hdchi->Fill((n-v)*(n-v)/v);
                hdchi2->Fill(tht.getVal(),phit.getVal(),(n-v)*(n-v)/v);
                mychi2 += (n-v)*(n-v)/v;
                //printf("pdf = %0.3f\tdata = %0.3f\n",pdf->weight(),data->weight());
            }
        }
    }

    TCanvas* c2 = new TCanvas("c2","c2",1000,600);
    c2->Divide(2);
    c2->cd(1);
    hdchi->Draw();
    c2->cd(2);
    hdchi2->SetOption("colz");
    hdchi2->SetMinimum(0);
    hdchi2->SetMaximum(1000);
    hdchi2->SetStats(kFALSE);
    //hdchi2->SaveAs(outputFile);
    hdchi2->Draw();
    c2->SaveAs("plots/dchi.root");
    c2->SaveAs("plots/dchi.png");
    printf("my chi2 = %f\n",mychi2);
    return 3; //mychi2;
}

void WriteToFile(Double_t var1)//, Double_t var2, Double_t var3, Double_t var4)
{
    FILE* pFile;
    pFile = fopen (outputFile,"w");
    if (pFile!=NULL)
    {
        //fprintf(pFile,"%f %f %f %f\n",var1,var2,var3,var4);
        fprintf(pFile,"%f\n",var1);
        fclose (pFile);
    }
    else
        printf("ERROR: couldn't open file %s for writing!\n",outputFile);

    return;
}