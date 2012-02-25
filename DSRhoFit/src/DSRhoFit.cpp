#include <stdio.h>
#include <stdlib.h>
#include <math.h>

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
#include "RooCategory.h"
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
#include "RooExponential.h"

#include "DSRhoFit.h"
#include "Constants.h"
#include "ASSERT.h"
#include "FitterTrans.h"

#define DEBUG
//#define VERBOSE
//#define GRAPHIC
//#define HELICITY
#define TRANSVERSITY

const Int_t var1_bins = 20;
const Int_t var2_bins = 20;
const Int_t var3_bins = 40;
const Int_t dt_bins = 40;

char* inputFile = 0;
char* outputFile = 0;

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

    if(argc != 16)
    {
        printf("ERROR: Wrong number of arguments.\n");
        #ifdef HELICITY
        printf("Usage: DSRhoFit inputFile outputFile hp hpa h0 hma [doFit] [doPlot]\n");
        #endif
        #ifdef TRANSVERSITY
        printf("Usage: DSRhoFit inputFile outputFile ap apa a0 ata phiw rp r0 rt sp s0 st doFit doPlot\n");
        #endif
        return 1;
    }

    TStopwatch timer;
    timer.Start();

    /// This is so I have to change only the next block if I change the
    /// ordering, etc. of arguments
    inputFile = argv[1];
    outputFile = argv[2];
    Int_t numPars = argc - 5;
    Double_t par_input[numPars];
    for(Int_t i = 0; i < numPars; i++)
        par_input[i] = atof(argv[i+3]);

    Bool_t doFit = atoi(argv[argc-2]);
    Bool_t doPlot = atoi(argv[argc-1]);

    RooRealVar tha("tha","tha",0,PI);
    RooRealVar thb("thb","thb",0,PI);
    RooRealVar chi("chi","chi",0,2*PI);
    RooRealVar tht("tht","tht",0,PI);
    RooRealVar phit("phit","phit",-PI,PI);
    RooRealVar dt("dt","dt",-10,10);
    RooCategory decType("decType","decType");
    decType.defineType("a",1);
    decType.defineType("ab",2);
    decType.defineType("b",3);
    decType.defineType("bb",4);

    #ifdef HELICITY
    RooDataSet* dataSet = new RooDataSet("data","data",RooArgList(tha,thb,chi,dt));
    dataSet = RooDataSet::read(inputFile,RooArgList(tha,thb,chi,dt));
    ConvertTransToHel(par_input);
    ProcessHel(dataSet,tha,thb,chi,dt,par_input,doFit,doPlot);
    #endif

    #ifdef TRANSVERSITY
    RooDataSet* dataSet = new RooDataSet("data","data",RooArgSet(tht,thb,phit,dt,decType));
    dataSet = RooDataSet::read(inputFile,RooArgList(tht,thb,phit,dt,decType));
    ProcessTrans(dataSet,par_input,doFit,doPlot);
    #endif

    timer.Stop();
    timer.Print();

    #ifdef GRAPHIC
    printf("\nProgram execution has finished.\n");
    rootapp->Run();
    #endif
    return 0;
}

int ProcessHel(RooDataSet* dataSet, RooRealVar& tha, RooRealVar& thb, RooRealVar& chi, RooRealVar& dt, Double_t* par_input, Bool_t doFit, Bool_t doPlot)
{
    RooRealVar hp("hp","hp",par_input[0],0,0.4);
    RooRealVar hpa("hpa","hpa",par_input[1],0,2*PI);
    RooFormulaVar hpr("hpr","hp*cos(hpa)",RooArgSet(hp,hpa));
    RooFormulaVar hpi("hpi","hp*sin(hpa)",RooArgSet(hp,hpa));
    RooRealVar h0("h0","h0",par_input[2],0.8,1);
    RooRealVar h0a("h0a","h0a",0);
    RooFormulaVar h0r("h0r","h0*cos(h0a)",RooArgSet(h0,h0a));
    RooFormulaVar h0i("h0i","h0*sin(h0a)",RooArgSet(h0,h0a));
    RooFormulaVar hm("hm","sqrt(1-hp*hp-h0*h0)",RooArgSet(hp,h0));
    RooRealVar hma("hma","hma",par_input[3],-PI,PI);
    RooFormulaVar hmr("hmr","hm*cos(hma)",RooArgSet(hm,hma));
    RooFormulaVar hmi("hmi","hm*sin(hma)",RooArgSet(hm,hma));

    if(doFit == kFALSE)
    {
        hp.setConstant();
        hpa.setConstant();
        h0.setConstant();
        hma.setConstant();
    }

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

    RooArgSet varSet(tha,thb,chi,hp,hpa,h0,h0a,hm,hma);
    varSet.add(hat);
    varSet.add(hata);
    varSet.add(hst);
    varSet.add(hsta);

    /// numFitParameters holds # of NON-constant fit parameters
    RooArgSet fitParameters(hp,hpa,h0,hma);
    Int_t numFitParameters = (fitParameters.selectByAttrib("Constant",kFALSE))->getSize();

    const char* pdfFormula =   "(hp*hp+hm*hm)*sin(tha)*sin(tha)*sin(tha)*sin(thb)*sin(thb)*sin(thb)+\
                                4*h0*h0*cos(tha)*cos(tha)*sin(tha)*cos(thb)*cos(thb)*sin(thb)+\
                                2*(hp*hm*cos(hpa-hma)*cos(2*chi)-hp*hm*sin(hpa-hma)*sin(2*chi))*sin(tha)*sin(tha)*sin(tha)*sin(thb)*sin(thb)*sin(thb)+\
                                (hat*h0*cos(hata-h0a)*cos(chi)-hst*h0*sin(hsta-h0a)*sin(chi))*sin(2*tha)*sin(tha)*sin(2*thb)*sin(thb)";

    RooGenericPdf* pdf = new RooGenericPdf("pdf","Generic PDF",pdfFormula,varSet);

    RooFitResult* result = 0;
    if(doFit == kTRUE)
    {
        result = pdf->fitTo(*dataSet,RooFit::Save(),RooFit::Timer(true));//,RooFit::NumCPU(2));
        result->Print();
    }

    /// This probably needs to be way over here for the fit few lines up to be unbinned
	tha.setBins(var1_bins);
	thb.setBins(var2_bins);
    chi.setBins(var3_bins);
    dt.setBins(dt_bins);

    RooRandom::randomGenerator()->SetSeed(0);

    /// Create a binned dataSet which is needed for chi2 calculation
	RooDataHist* dataSet_binned = new RooDataHist("dataSet_binned","dataSet_binned",RooArgSet(tha,thb,chi),*dataSet);
    //RooDataHist* dataSet_binned = pdf->generateBinned(RooArgSet(tha,thb,chi),dataSet->numEntries(),kFALSE);

	RooChi2Var chi2Var("chi2Var","chi2Var",*pdf,*dataSet_binned);

	RooRealVar* chi2     = new RooRealVar("chi2","chi^2",0) ;
	RooRealVar* ndof     = new RooRealVar("ndof","number of degrees of freedom",0) ;
	RooRealVar* chi2red  = new RooRealVar("chi2red","reduced chi^2",0) ;
	RooRealVar* prob     = new RooRealVar("prob","prob(chi2,ndof)",0) ;

	chi2->setVal(chi2Var.getVal());
	ndof->setVal(dataSet_binned->numEntries()-numFitParameters);
	chi2red->setVal(chi2->getVal()/ndof->getVal()) ;
	prob->setVal(TMath::Prob(chi2->getVal(),static_cast<int>(ndof->getVal())));

	printf("chi2 = %f\nndof = %f\nchi2red = %f\nprob = %f\n",chi2->getVal(),ndof->getVal(),chi2red->getVal(),prob->getVal());

    if(doFit == kTRUE)
    {
        Int_t numParameters = 17;
        Double_t recoveredParameters[17] = {chi2red->getVal(),hp.getVal(),hp.getError(),hpa.getVal(),hpa.getError(),
                                                       h0.getVal(),h0.getError(),h0a.getVal(),h0a.getError(),hm.getVal(),
                                                       hm.getPropagatedError(*result),hma.getVal(),hma.getError(),
                                                       par_input[0],par_input[1],par_input[2],par_input[3]};
        WriteToFile(numParameters, recoveredParameters, outputFile);
    }

    if(doPlot == kTRUE)
    {
        SavePlots(dataSet,pdf,tha,thb,chi,dt);
        SaveChi2Maps(dataSet_binned,dataSet->numEntries(),pdf,tha,thb,chi);
    }

    return 0;
}

int ProcessTrans(RooDataSet* dataSet, Double_t* par_input, Bool_t doFit, Bool_t doPlot)
{
    FitterTrans* fitter = new FitterTrans(dataSet,par_input);

    if(doFit)
    {
        fitter->FixAllParameters();
        fitter->FreeParameter("s0");
        fitter->Fit();
        fitter->ComputeChi2();

        Int_t numParameters = 0;
        Double_t* recoveredParameters = 0;

        //fitter->GetRecoveredParameters(numParameters,&recoveredParameters);

        //WriteToFile(numParameters,recoveredParameters,outputFile);
    }
    else
        fitter->ComputeChi2();

    if(doPlot == kTRUE)
    {
        SavePlots(dataSet,fitter->GetPdf(),*(fitter->GetTht()),*(fitter->GetThb()),*(fitter->GetPhit()),*(fitter->GetDt()));
        SaveChi2Maps(fitter->GetBinnedDataSet(),dataSet->numEntries(),fitter->GetPdf(),*(fitter->GetTht()),*(fitter->GetThb()),*(fitter->GetPhit()));
    }

    return 0;
}

Double_t SaveChi2Maps(RooDataHist* data_binned, Int_t numEvents, RooGenericPdf* pdf, RooRealVar var1, RooRealVar var2, RooRealVar var3)
{
    /// Create a histogram from the pdf with the expected number of events with no statistical fluctuation
    RooDataHist* pdf_binned = pdf->generateBinned(RooArgSet(var1,var2,var3),numEvents,kTRUE);

    Double_t mychi2 = 0;
    Double_t dchi2 = 0;
    Double_t n = 0;
    Double_t v = 0;

    TH1F* h1_chi2 = new TH1F("h1_chi2","h1_chi2",100,0,100);
    TH2F* h2_chi2_1 = new TH2F("h2_chi2_1","h2_chi2_1",var1_bins,var1.getMin(),var1.getMax(),var2_bins,var2.getMin(),var2.getMax());
    TH2F* h2_chi2_2 = new TH2F("h2_chi2_2","h2_chi2_2",var1_bins,var1.getMin(),var1.getMax(),var3_bins,var3.getMin(),var3.getMax());
    TH2F* h2_chi2_3 = new TH2F("h2_chi2_3","h2_chi2_3",var2_bins,var2.getMin(),var2.getMax(),var3_bins,var3.getMin(),var3.getMax());

    /// Cycle through the centers of all bins
    /// I'm getting width of the first bin, because all bins are of equal width
    for(var1 = var1.getMin()+var1.getBinWidth(0)/2; var1.getVal() < var1.getMax(); var1.setVal(var1.getVal()+var1.getBinWidth(0)))
    {
        for(var2 = var2.getMin()+var2.getBinWidth(0)/2; var2.getVal() < var2.getMax(); var2.setVal(var2.getVal()+var2.getBinWidth(0)))
        {
            for(var3 = var3.getMin()+var3.getBinWidth(0)/2; var3.getVal() < var3.getMax(); var3.setVal(var3.getVal()+var3.getBinWidth(0)))
            {
                /// Weight is actually the bin content
                n = data_binned->weight(RooArgSet(var1,var2,var3),0);
                v = pdf_binned->weight(RooArgSet(var1,var2,var3),0);
                dchi2 = (n-v)*(n-v)/v;
                if(dchi2 > h1_chi2->GetXaxis()->GetXmax()-1)
                    dchi2 = h1_chi2->GetXaxis()->GetXmax()-1;
                h1_chi2->Fill(dchi2);
                h2_chi2_1->Fill(var1.getVal(),var2.getVal(),dchi2);
                h2_chi2_2->Fill(var1.getVal(),var3.getVal(),dchi2);
                h2_chi2_3->Fill(var2.getVal(),var3.getVal(),dchi2);
                mychi2 += dchi2;
            }
        }
    }

    TFile* file = new TFile("plots/chi2maps.root","RECREATE");
    TCanvas* c2 = new TCanvas("c2","c2",600,600);
    TString path;

    c2->SetLogy(kTRUE);
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
    h2_chi2_1->GetXaxis()->SetTitle(var1.GetName());
    h2_chi2_1->GetYaxis()->SetTitle(var2.GetName());
    h2_chi2_1->Draw();
    h2_chi2_1->Write();
    path = "plots/chi2map_";
    path += var1.GetName();
    path += "_";
    path += var2.GetName();
    path += ".png";
    c2->SaveAs(path);

    h2_chi2_2->SetOption("colz");
    //hdchi2->SetMinimum(0);
    //hdchi2->SetMaximum(100);
    h2_chi2_2->SetStats(kFALSE);
    h2_chi2_2->GetXaxis()->SetTitle(var1.GetName());
    h2_chi2_2->GetYaxis()->SetTitle(var3.GetName());
    h2_chi2_2->Draw();
    h2_chi2_2->Write();
    path = "plots/chi2map_";
    path += var1.GetName();
    path += "_";
    path += var3.GetName();
    path += ".png";
    c2->SaveAs(path);

    h2_chi2_3->SetOption("colz");
    //hdchi2->SetMinimum(0);
    //hdchi2->SetMaximum(100);
    h2_chi2_3->SetStats(kFALSE);
    h2_chi2_3->GetXaxis()->SetTitle(var2.GetName());
    h2_chi2_3->GetYaxis()->SetTitle(var3.GetName());
    h2_chi2_3->Draw();
    h2_chi2_3->Write();
    path = "plots/chi2map_";
    path += var2.GetName();
    path += "_";
    path += var3.GetName();
    path += ".png";
    c2->SaveAs(path);

    file->Close();

    delete h1_chi2;
    delete h2_chi2_1;
    delete h2_chi2_2;
    delete h2_chi2_3;
    delete c2;

    return mychi2;
}

Double_t GetChi2(RooDataHist* data_binned, Int_t numEvents, RooGenericPdf* pdf, RooRealVar var1, RooRealVar var2, RooRealVar var3)
{
    /// Create a histogram from the pdf with the expected number of events with no statistical fluctuation
    RooDataHist* pdf_binned = pdf->generateBinned(RooArgSet(var1,var2,var3),numEvents,kTRUE);

    Double_t mychi2 = 0;
    Double_t n = 0;
    Double_t v = 0;

    /// Cycle through the centers of all bins
    /// I'm getting width of the first bin, because all bins are of equal width
    for(var1 = var1.getMin()+var1.getBinWidth(0)/2; var1.getVal() < var1.getMax(); var1.setVal(var1.getVal()+var1.getBinWidth(0)))
    {
        for(var2 = var2.getMin()+var2.getBinWidth(0)/2; var2.getVal() < var2.getMax(); var2.setVal(var2.getVal()+var2.getBinWidth(0)))
        {
            for(var3 = var3.getMin()+var3.getBinWidth(0)/2; var3.getVal() < var3.getMax(); var3.setVal(var3.getVal()+var3.getBinWidth(0)))
            {
                /// Weight is actually the bin content
                n = data_binned->weight(RooArgSet(var1,var2,var3),0);
                v = pdf_binned->weight(RooArgSet(var1,var2,var3),0);
                mychi2 += (n-v)*(n-v)/v;
            }
        }
    }

    return mychi2;
}


void WriteToFile(Int_t numEntries, Double_t* vars, char* file)
{
    FILE* pFile;
    pFile = fopen (file,"w");
    if (pFile == NULL)
    {
        printf("ERROR: couldn't open file %s for writing!\n",file);
        return;
    }

    for(Int_t i = 0; i < numEntries; i++, vars++)
    {
        fprintf(pFile,"%f ",*vars);
    }
    fprintf(pFile,"\n");
    fclose (pFile);

    return;
}

void SavePlots(RooDataSet* dataSet, RooGenericPdf* pdf, const RooRealVar& var1, const RooRealVar& var2, const RooRealVar& var3, RooRealVar& dt)
{
    /// Directory and format of the saved plots
    const TString dir = "plots/";
    const TString format = ".png";
    TString path;
    TString name;

    path = dir + "projections.root";
    TFile* file = new TFile(path,"RECREATE");

    TCanvas* c1 = new TCanvas("c1","c1",800,600);
    RooPlot* frame = 0;

    /// Create a binned pdf with the same number of events as the data, so that the 2d plots of pdf
    /// and data are the same scale
    RooDataHist* pdf_binned = pdf->generateBinned(RooArgSet(var1,var2,var3),dataSet->numEntries(),kTRUE);

    /// This is a quick and dirty solution to be able to loop through all variables (except dt, which is treated separetely)
    const int numVars = 3;
    RooRealVar vars[numVars] = {var1,var2,var3};

    /// Saving simple projections on each of the variables
    for(int i = 0; i < numVars; i++)
    {
        name = "proj_";
        name += vars[i].GetName();
        frame = vars[i].frame();
        dataSet->plotOn(frame,RooFit::Name("data"));
        //pdf->plotOn(frame);
        frame->SetName(name);
        frame->Draw();
        frame->Write();
        path = dir + name + format;
        c1->SaveAs(path);
        delete frame;
    }

    /// The next 2 lines enable getting category items' names and therefore reduced datasets in a loop
    const RooArgSet* args = dataSet->get();
    const RooCategory* cat = (RooCategory*)args->find("decType");
    RooDataSet* datacut;

    RooFormulaVar absdt("absdt","abs(dt)",RooArgSet(dt));
    RooRealVar gamma("gamma","gamma",-10,10);
    RooExponential myexp("myexp","myexp",absdt,gamma);

    /// Saving dt plots for all 4 decay types
    for(int i = 1; i <= 4; i++)
    {
        frame = dt.frame();
        TString type = (char*)cat->lookupType(i)->GetName();
        TString cut = "decType==decType::" + type;
        name = "proj_" + (dt.GetName() + ("_" + type));
        datacut = (RooDataSet*)dataSet->reduce(dt,cut);
        datacut->plotOn(frame,RooFit::Name("data"));
        if(i == 1)
        {
            myexp.fitTo(*datacut,RooFit::Range(-2,2,kTRUE),RooFit::Minimizer("Minuit2"));
            myexp.plotOn(frame);
        }
        frame->SetName(name);
        frame->Draw();
        frame->Write();
        path = dir + name + format;
        c1->SaveAs(path);
        delete frame;
    }

    TH2* h2_pdf = 0;
    TH2* h2_data = 0;

    /// Saving projections of both data and pdf on 2 dimensions
    for(int i = 0; i < numVars; i++)
    {
        for(int j = i+1; j < numVars; j++)
        {
            name = "proj_";
            name += vars[i].GetName();
            name += "_";
            name += vars[j].GetName();
            h2_pdf = dynamic_cast<TH2*>(pdf_binned->createHistogram(name + "_pdf",vars[i],RooFit::YVar(vars[j])));
            h2_pdf->SetOption("colz");
            h2_pdf->SetStats(kFALSE);
            h2_pdf->SetMinimum(0);
            h2_data = dynamic_cast<TH2*>(dataSet->createHistogram(name + "_data",vars[i],RooFit::YVar(vars[j])));
            h2_data->SetOption("colz");
            h2_data->SetStats(kFALSE);
            h2_data->SetMinimum(0);
            h2_data->SetMaximum(h2_pdf->GetMaximum());
            h2_pdf->Draw();
            h2_pdf->Write();
            path = dir + name + "_pdf" + format;
            c1->SaveAs(path);
            h2_data->Draw();
            h2_data->Write();
            path = dir + name + "_data" + format;
            c1->SaveAs(path);
            delete h2_pdf;
            delete h2_data;
        }
    }

    file->Close();
}

void ConvertTransToHel(Double_t* par_input)
{
    RooRealVar ap("ap","ap",par_input[0]);
    RooRealVar apa("apa","apa",par_input[1]);
    RooFormulaVar apr("apr","ap*cos(apa)",RooArgSet(ap,apa));
    RooFormulaVar api("api","ap*sin(apa)",RooArgSet(ap,apa));
    RooRealVar a0("a0","a0",par_input[2]);
    RooFormulaVar at("at","sqrt(1-ap*ap-a0*a0)",RooArgSet(ap,a0));
    RooRealVar ata("ata","ata",par_input[3]);
    RooFormulaVar atr("atr","at*cos(ata)",RooArgSet(at,ata));
    RooFormulaVar ati("ati","at*sin(ata)",RooArgSet(at,ata));

    RooFormulaVar hpr("hpr","(apr + atr)/sqrt(2)",RooArgSet(apr,atr));
    RooFormulaVar hpi("hpi","(api + ati)/sqrt(2)",RooArgSet(api,ati));
    RooFormulaVar hp("hp","sqrt(hpr*hpr+hpi*hpi)",RooArgSet(hpr,hpi));
    RooFormulaVar hpa("hpa","atan2(hpi,hpr)",RooArgSet(hpr,hpi));

    RooFormulaVar hmr("hmr","(apr - atr)/sqrt(2)",RooArgSet(apr,atr));
    RooFormulaVar hmi("hmi","(api - ati)/sqrt(2)",RooArgSet(api,ati));
    RooFormulaVar hm("hm","sqrt(hmr*hmr+hmi*hmi)",RooArgSet(hmr,hmi));
    RooFormulaVar hma("hma","atan2(hmi,hmr)",RooArgSet(hmr,hmi));

    printf("original trans:\t");
    printf("par[0] = %f\tpar[1] = %f\tpar[2] = %f\tpar[3] = %f\n",par_input[0],par_input[1],par_input[2],par_input[3]);

    par_input[0] = Round(hp.getVal(),2);

    if(hpa.getVal() < 0)
        par_input[1] = Round(hpa.getVal()+2*PI,2);
    else
        par_input[1] = Round(hpa.getVal(),2);

    /// par_input[2] = par_input[2]; because a0 = h0

    par_input[3] = Round(hma.getVal(),2);

    printf("converted hel:\t");
    printf("par[0] = %f\tpar[1] = %f\tpar[2] = %f\tpar[3] = %f\n",par_input[0],par_input[1],par_input[2],par_input[3]);
}

void ConvertHelToTrans(Double_t* par_input)
{
    RooRealVar hp("hp","hp",par_input[0]);
    RooRealVar hpa("hpa","hpa",par_input[1]);
    RooFormulaVar hpr("hpr","hp*cos(hpa)",RooArgSet(hp,hpa));
    RooFormulaVar hpi("hpi","hp*sin(hpa)",RooArgSet(hp,hpa));
    RooRealVar h0("h0","h0",par_input[2]);
    RooFormulaVar hm("hm","sqrt(1-hp*hp-h0*h0)",RooArgSet(hp,h0));
    RooRealVar hma("hma","hma",par_input[3]);
    RooFormulaVar hmr("hmr","hm*cos(hma)",RooArgSet(hm,hma));
    RooFormulaVar hmi("hmi","hm*sin(hma)",RooArgSet(hm,hma));

    RooFormulaVar apr("apr","(hpr + hmr)/sqrt(2)",RooArgSet(hpr,hmr));
    RooFormulaVar api("api","(hpi + hmi)/sqrt(2)",RooArgSet(hpi,hmi));
    RooFormulaVar ap("ap","sqrt(apr*apr+api*api)",RooArgSet(apr,api));
    RooFormulaVar apa("apa","atan2(api,apr)",RooArgSet(apr,api));

    RooFormulaVar atr("atr","(hpr - hmr)/sqrt(2)",RooArgSet(hpr,hmr));
    RooFormulaVar ati("ati","(hpi - hmi)/sqrt(2)",RooArgSet(hpi,hmi));
    RooFormulaVar at("at","sqrt(atr*atr+ati*ati)",RooArgSet(atr,ati));
    RooFormulaVar ata("ata","atan2(ati,atr)",RooArgSet(atr,ati));

    printf("original hel:\t");
    printf("par[0] = %f\tpar[1] = %f\tpar[2] = %f\tpar[3] = %f\n",par_input[0],par_input[1],par_input[2],par_input[3]);

    par_input[0] = Round(ap.getVal(),2);

    if(apa.getVal() < 0)
        par_input[1] = Round(apa.getVal()+2*PI,2);
    else
        par_input[1] = Round(apa.getVal(),2);

    /// par_input[2] = par_input[2]; because a0 = h0

    if(ata.getVal() < 0)
        par_input[3] = Round(ata.getVal()+2*PI,2);
    else
        par_input[3] = Round(ata.getVal(),2);

    printf("converted trans:\t");
    printf("par[0] = %f\tpar[1] = %f\tpar[2] = %f\tpar[3] = %f\n",par_input[0],par_input[1],par_input[2],par_input[3]);
}

Double_t Round(Double_t number, Int_t digits)
{
    number = number * pow(10,digits);
    if(fmod(number,1)>0.5)
        number = ceil(number);
    else
        number = floor(number);

    return number/pow(10,digits);
}
