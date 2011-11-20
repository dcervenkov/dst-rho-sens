#include <stdio.h>
#include <stdlib.h>

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
#include "RooGlobalFunc.h"
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

const Int_t tha_bins = 20;
const Int_t thb_bins = 20;
const Int_t chi_bins = 40;
RooRealVar tha("tha","tha",0,PI);
RooRealVar thb("thb","thb",0,PI);
RooRealVar chi("chi","chi",0,2*PI);

const Int_t tht_bins = 20;
const Int_t phit_bins = 40;
RooRealVar tht("tht","tht",0,PI);
RooRealVar phit("phit","phit",-PI,PI);

char* inputFile = 0;
char* outputFile = 0;

Double_t par1_input = 0.275;
Double_t par2_input = 0.57;
Double_t par3_input = 0.936;
Double_t par4_input = 2.84;

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

    if(argc < 3)
    {
        #ifdef HELICITY
        printf("Usage: DSRhoFit inputFile outputFile [hp] [hpa] [h0] [hma]\n");
        #endif
        #ifdef TRANSVERSITY
        printf("Usage: DSRhoFit inputFile outputFile [ap] [apa] [a0] [ata]\n");
        #endif
        return 1;
    }

    inputFile = argv[1];
    outputFile = argv[2];

    par1_input = atof(argv[3]);
    par2_input = atof(argv[4]);
    par3_input = atof(argv[5]);
    par4_input = atof(argv[6]);

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

    RooRealVar hp("hp","hp",par1_input,0,0.3);
    RooRealVar hpa("hpa","hpa",par2_input,0,2*PI);
    RooFormulaVar hpr("hpr","hp*cos(hpa)",RooArgSet(hp,hpa));
    RooFormulaVar hpi("hpi","hp*sin(hpa)",RooArgSet(hp,hpa));
    RooRealVar h0("h0","h0",par3_input,0.8,1);
    RooRealVar h0a("h0a","h0a",0);
    RooFormulaVar h0r("h0r","h0*cos(h0a)",RooArgSet(h0,h0a));
    RooFormulaVar h0i("h0i","h0*sin(h0a)",RooArgSet(h0,h0a));
    RooFormulaVar hm("hm","sqrt(1-hp*hp-h0*h0)",RooArgSet(hp,h0));
    RooRealVar hma("hma","hma",par4_input,0,2*PI);
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

    RooFitResult* result = pdf->fitTo(*dataSet,RooFit::Save(),RooFit::Timer(true));//,RooFit::NumCPU(2));
    result->Print();

    chi.setBins(chi_bins);
	tha.setBins(tha_bins);
	thb.setBins(thb_bins);

    RooRandom::randomGenerator()->SetSeed(0);

	RooDataHist* dataSet_binned = new RooDataHist("dataSet_binned","dataSet_binned",RooArgSet(chi,tha,thb),*dataSet);
    //RooDataHist* dataSet_binned = pdf->generateBinned(RooArgSet(tha,thb,chi),dataSet->numEntries(),kFALSE);
    //RooDataHist* pdf_binned = pdf->generateBinned(RooArgSet(tha,thb,chi),dataSet->numEntries(),kTRUE);

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

//    Int_t numParameters = 13;
//	Double_t recoveredParameters[13] = {_chi2red->getVal(),hp.getVal(),hp.getError(),hpa.getVal(),hpa.getError(),
//                                                   h0.getVal(),h0.getError(),h0a.getVal(),h0a.getError(),hm.getVal(),
//                                                   hm.getPropagatedError(*result),hma.getVal(),hma.getError()};
//
//    WriteToFile(numParameters, recoveredParameters);

    SavePlots(dataSet,pdf,tht,thb,phit);

//    GetChi2(dataSet_binned,pdf_binned);

    return 0;
}

int FitTrans(RooDataSet* dataSet, RooRealVar& tht, RooRealVar& thb, RooRealVar& phit)
{
    RooRealVar ap("ap","ap",par1_input,0.1,0.4);
    RooRealVar apa("apa","apa",par2_input,0,2*PI);
    RooFormulaVar apr("apr","ap*cos(apa)",RooArgSet(ap,apa));
    RooFormulaVar api("api","ap*sin(apa)",RooArgSet(ap,apa));
    RooRealVar a0("a0","a0",par3_input,0.8,1);
    RooRealVar a0a("a0a","a0a",0);
    RooFormulaVar a0r("a0r","a0*cos(a0a)",RooArgSet(a0,a0a));
    RooFormulaVar a0i("a0i","a0*sin(a0a)",RooArgSet(a0,a0a));
    RooFormulaVar at("at","sqrt(1-ap*ap-a0*a0)",RooArgSet(ap,a0));
    RooRealVar ata("ata","ata",par4_input,0,2*PI);
    RooFormulaVar atr("atr","at*cos(ata)",RooArgSet(at,ata));
    RooFormulaVar ati("ati","at*sin(ata)",RooArgSet(at,ata));

//    ap.setConstant();
//    apa.setConstant();
//    a0.setConstant();
//    ata.setConstant();

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

    RooFitResult* result = pdf->fitTo(*dataSet,RooFit::Save(),RooFit::Timer(true));//,RooFit::NumCPU(2));
    result->Print();

    tht.setBins(tht_bins);
    thb.setBins(thb_bins);
    phit.setBins(phit_bins);

    RooRandom::randomGenerator()->SetSeed(0);

	RooDataHist* dataSet_binned = new RooDataHist("dataSet_binned","dataSet_binned",RooArgSet(tht,thb,phit),*dataSet);
    //RooDataHist* dataSet_binned = pdf->generateBinned(RooArgSet(tht,thb,phit),dataSet->numEntries(),kFALSE);
    //RooDataHist* pdf_binned = pdf->generateBinned(RooArgSet(tht,thb,phit),dataSet->numEntries(),kTRUE);

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

//	Int_t numParameters = 13;
//	Double_t recoveredParameters[13] = {_chi2red->getVal(),ap.getVal(),ap.getError(),apa.getVal(),apa.getError(),
//                                                   a0.getVal(),a0.getError(),a0a.getVal(),a0a.getError(),at.getVal(),
//                                                   at.getPropagatedError(*result),ata.getVal(),ata.getError()};
//
//    WriteToFile(numParameters, recoveredParameters);

    SavePlots(dataSet,pdf,tht,thb,phit);

//    GetChi2Trans(dataSet_binned,pdf_binned);

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

    TCanvas* c2 = new TCanvas("c2","c2");
    c2->Divide(2);
    c2->cd(1);
    hdchi->Draw();
    c2->cd(2);
    hdchi2->SetOption("colz");
    //hdchi2->SetMinimum(0);
    //hdchi2->SetMaximum(100);
    hdchi2->SetStats(kFALSE);
//    hdchi2->SaveAs(outputFile);
    hdchi2->Draw();
    c2->SaveAs("plots/dchi.root");
    c2->SaveAs("plots/dchi.png");
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
    //hdchi2->SetMinimum(0);
    //hdchi2->SetMaximum(1000);
    hdchi2->SetStats(kFALSE);
    //hdchi2->SaveAs(outputFile);
    hdchi2->Draw();
    c2->SaveAs("plots/dchi.root");
    c2->SaveAs("plots/dchi.png");
    printf("my chi2 = %f\n",mychi2);
    return 3; //mychi2;
}

void WriteToFile(Int_t numEntries, Double_t* vars)
{
    FILE* pFile;
    pFile = fopen (outputFile,"w");
    if (pFile!=NULL)
    {
        for(Int_t i = 0; i < numEntries; i++, vars++)
        {
            fprintf(pFile,"%f ",*vars);
        }
        //fprintf(pFile,"\n");
        fclose (pFile);
    }
    else
        printf("ERROR: couldn't open file %s for writing!\n",outputFile);

    return;
}

void SavePlots(RooDataSet* dataSet, RooGenericPdf* pdf, const RooRealVar& var1, const RooRealVar& var2, const RooRealVar& var3)
{
    TFile* file = new TFile("plots/projections.root","RECREATE");
    TCanvas* c1 = new TCanvas("c1","c1",800,600);
    TString path;
    RooDataHist* pdf_binned = pdf->generateBinned(RooArgSet(tht,thb,phit),dataSet->numEntries(),kTRUE);

    RooPlot* frame1 = var1.frame();
    dataSet->plotOn(frame1,RooFit::Name("data"));
    pdf->plotOn(frame1);
    frame1->Draw();
    frame1->Write();
    path = "plots/proj_";
    path += var1.GetName();
    path += ".png";
    c1->SaveAs(path);

    RooPlot* frame2 = var2.frame();
    dataSet->plotOn(frame2,RooFit::Name("data"));
    pdf->plotOn(frame2);
    frame2->Draw();
    frame2->Write();
    path = "plots/proj_";
    path += var2.GetName();
    path += ".png";
    c1->SaveAs(path);

    RooPlot* frame3 = var3.frame();
    dataSet->plotOn(frame3,RooFit::Name("data"));
    pdf->plotOn(frame3);
    frame3->Draw();
    frame3->Write();
    path = "plots/proj_";
    path += var3.GetName();
    path += ".png";
    c1->SaveAs(path);


    TH2* h2_1_pdf = dynamic_cast<TH2*>(pdf_binned->createHistogram("h2_1_pdf",var1,RooFit::YVar(var2)));
    h2_1_pdf->SetOption("colz");
    h2_1_pdf->SetStats(kFALSE);
    h2_1_pdf->SetMinimum(0);
    h2_1_pdf->Draw();
    h2_1_pdf->Write();
    path = "plots/proj_";
    path += var1.GetName();
    path += "_";
    path += var2.GetName();
    path += "_pdf.png";
    c1->SaveAs(path);

    TH2* h2_1_data = dynamic_cast<TH2*>(dataSet->createHistogram("h2_1_data",var1,RooFit::YVar(var2)));
    h2_1_data->SetOption("colz");
    h2_1_data->SetStats(kFALSE);
    h2_1_data->SetMinimum(0);
    h2_1_data->SetMaximum(h2_1_pdf->GetMaximum());
    h2_1_data->Draw();
    h2_1_data->Write();
    path = "plots/proj_";
    path += var1.GetName();
    path += "_";
    path += var2.GetName();
    path += "_data.png";
    c1->SaveAs(path);

    TH2* h2_2_pdf = dynamic_cast<TH2*>(pdf_binned->createHistogram("h2_2_pdf",var1,RooFit::YVar(var3)));
    h2_2_pdf->SetOption("colz");
    h2_2_pdf->SetStats(kFALSE);
    h2_2_pdf->SetMinimum(0);
    h2_2_pdf->Draw();
    h2_2_pdf->Write();
    path = "plots/proj_";
    path += var1.GetName();
    path += "_";
    path += var3.GetName();
    path += "_pdf.png";
    c1->SaveAs(path);

    TH2* h2_2_data = dynamic_cast<TH2*>(dataSet->createHistogram("h2_2_data",var1,RooFit::YVar(var3)));
    h2_2_data->SetOption("colz");
    h2_2_data->SetStats(kFALSE);
    h2_2_data->SetMinimum(0);
    h2_2_data->SetMaximum(h2_2_pdf->GetMaximum());
    h2_2_data->Draw();
    h2_2_data->Write();
    path = "plots/proj_";
    path += var1.GetName();
    path += "_";
    path += var3.GetName();
    path += "_data.png";
    c1->SaveAs(path);

    TH2* h2_3_pdf = dynamic_cast<TH2*>(pdf_binned->createHistogram("h2_3_pdf",var2,RooFit::YVar(var3)));
    h2_3_pdf->SetOption("colz");
    h2_3_pdf->SetStats(kFALSE);
    h2_3_pdf->SetMinimum(0);
    h2_3_pdf->Draw();
    h2_3_pdf->Write();
    path = "plots/proj_";
    path += var2.GetName();
    path += "_";
    path += var3.GetName();
    path += "_pdf.png";
    c1->SaveAs(path);

    TH2* h2_3_data = dynamic_cast<TH2*>(dataSet->createHistogram("h2_3_data",var2,RooFit::YVar(var3)));
    h2_3_data->SetOption("colz");
    h2_3_data->SetStats(kFALSE);
    h2_3_data->SetMinimum(0);
    h2_3_data->SetMaximum(h2_3_pdf->GetMaximum());
    h2_3_data->Draw();
    h2_3_data->Write();
    path = "plots/proj_";
    path += var2.GetName();
    path += "_";
    path += var3.GetName();
    path += "_data.png";
    c1->SaveAs(path);

    file->Close();

}
