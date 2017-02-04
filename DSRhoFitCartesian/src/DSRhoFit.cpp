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
#include "TComplex.h"

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
#include "FitterTransTIndep.h"

#include <unistd.h>

#define DEBUG
//#define VERBOSE
//#define GRAPHIC

char* inputFile = 0;
char* outputFile = 0;

int main(int argc, char* argv[]) {
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
    for(int i = 0; i < argc; i++) {
        argv[i] = rootapp->Argv(i);
    }
#endif

    if(argc != 21) {
        printf("ERROR: Wrong number of arguments.\n");
        printf("Usage: DSRhoFit inputFile outputFile ap apa a0 ata xp x0 xt yp y0 yt xpb x0b xtb ypb y0b ytb doFit doPlot\n");
        return 85;
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

    Int_t doFit = atoi(argv[argc-2]);
    Int_t doPlot = atoi(argv[argc-1]);

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

    if(doFit == 3) {
        ConvertBetweenHelAndTrans(par_input);
    } else if(doFit == 4) {
        FitterTrans* fitter = new FitterTrans(par_input);
        fitter->GenerateDataSet(100000);
        fitter->GetDataSet()->write(inputFile);
    } else {
        //ConvertBetweenHelAndTrans(par_input);
        FitterTrans* fitter = new FitterTrans(par_input);
        fitter->ReadDataSet(inputFile);
        ProcessTrans(fitter,doFit,doPlot);
    }

    timer.Stop();
    timer.Print();

#ifdef GRAPHIC
    printf("\nProgram execution has finished.\n");
    rootapp->Run();
#endif
    return 0;
}


int ProcessTrans(FitterTrans* fitter, Int_t doFit, Int_t doPlot) {
    if(doFit) {
        fitter->FixAllParameters();
        fitter->FreeParameter("ap");
        fitter->FreeParameter("apa");
        fitter->FreeParameter("a0");
        fitter->FreeParameter("ata");
        fitter->FreeParameter("xp");
        fitter->FreeParameter("x0");
        fitter->FreeParameter("xt");
        fitter->FreeParameter("yp");
        fitter->FreeParameter("y0");
        fitter->FreeParameter("yt");
        fitter->FreeParameter("xpb");
        fitter->FreeParameter("x0b");
        fitter->FreeParameter("xtb");
        fitter->FreeParameter("ypb");
        fitter->FreeParameter("y0b");
        fitter->FreeParameter("ytb");
        fitter->Fit();
        fitter->SaveParameters(outputFile);
    }

    if(doPlot == kTRUE) {
        //SaveChi2Maps(fitter->GetBinnedDataSet(),dataSet->numEntries(),fitter->GetPdf(),*(fitter->GetTht()),*(fitter->GetThb()),*(fitter->GetPhit()));
        //Double_t mychi2 = fitter->SaveChi2Maps("a");
        //fitter->SaveResiduals();
        //fitter->SaveNllPlot("yt");
        //printf("mychi2 from SaveChi2Maps = %f\n",mychi2);
        SavePlots(fitter->GetDataSet(),fitter->GetPdf(),fitter->GetPdfBar(),*(fitter->GetTht()),*(fitter->GetThb()),*(fitter->GetPhit()),*(fitter->GetDt()));
    }

    return 0;
}



Double_t SaveChi2Maps(RooDataHist* data_binned, Int_t numEvents, RooGenericPdf* pdf, RooRealVar var1, RooRealVar var2, RooRealVar var3) {
    /// Create a histogram from the pdf with the expected number of events with no statistical fluctuation
    RooDataHist* pdf_binned = pdf->generateBinned(RooArgSet(var1,var2,var3),numEvents,kTRUE);

    const Int_t var1_bins = 30;
    const Int_t var2_bins = 30;
    const Int_t var3_bins = 30;
    const Int_t dt_bins = 50;

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
    for(var1 = var1.getMin()+var1.getBinWidth(0)/2; var1.getVal() < var1.getMax(); var1.setVal(var1.getVal()+var1.getBinWidth(0))) {
        for(var2 = var2.getMin()+var2.getBinWidth(0)/2; var2.getVal() < var2.getMax(); var2.setVal(var2.getVal()+var2.getBinWidth(0))) {
            for(var3 = var3.getMin()+var3.getBinWidth(0)/2; var3.getVal() < var3.getMax(); var3.setVal(var3.getVal()+var3.getBinWidth(0))) {
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

void SavePlots(RooDataSet* dataSet, DSRhoPDF* pdf, DSRhoPDF* pdf_bar, const RooRealVar& var1, const RooRealVar& var2, const RooRealVar& var3, RooRealVar& dt) {
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
    //printf("Generation of binned dataset starting...\n");
    //RooDataHist* pdf_binned = pdf->generateBinned(RooArgSet(var1,var2,var3,dt),dataSet->numEntries(),kTRUE);
    //printf("Generation of binned dataset finished\n");

    /// This is a quick and dirty solution to be able to loop through all variables (except dt, which is treated separetely)
    const int numVars = 3;
    RooRealVar vars[numVars] = {var1,var2,var3};

    /// Saving simple projections on each of the variables
    for(int i = 0; i < numVars; i++) {
        name = "proj_";
        name += vars[i].GetName();
        frame = vars[i].frame();
        dataSet->plotOn(frame,RooFit::Name("data"));
        //if(i == 1)
        pdf->plotOn(frame,RooFit::Project(RooArgSet(var1,var2,var3,dt)));
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

    /// Saving dt plots for all 4 decay types
    for(int i = 1; i <= 4; i++) {
        frame = dt.frame();
        TString type = (char*)cat->lookupType(i)->GetName();
        TString cut = "decType==decType::" + type;
        name = "proj_" + (dt.GetName() + ("_" + type));
        datacut = (RooDataSet*)dataSet->reduce(dt,cut);
        datacut->plotOn(frame,RooFit::Name("data"));

        pdf->setType(i);
        pdf_bar->setType(i);
//        if(i == 3)
//            pdf->setType(4);
//        else if(i == 4)
//            pdf->setType(3);
//        else if(i == 1)
//            pdf->setType(2);
//        else if(i == 2)
//            pdf->setType(1);

        if (i == 1 || i == 4) {
        	pdf->plotOn(frame,RooFit::Project(RooArgSet(var1,var2,var3)));
        } else {
        	pdf_bar->plotOn(frame,RooFit::Project(RooArgSet(var1,var2,var3)));
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
//    for(int i = 0; i < numVars; i++)
//    {
//        for(int j = i+1; j < numVars; j++)
//        {
//            name = "proj_";
//            name += vars[i].GetName();
//            name += "_";
//            name += vars[j].GetName();
//            h2_pdf = dynamic_cast<TH2*>(pdf_binned->createHistogram(name + "_pdf",vars[i],RooFit::YVar(vars[j])));
//            h2_pdf->SetOption("colz");
//            h2_pdf->SetStats(kFALSE);
//            h2_pdf->SetMinimum(0);
//            h2_data = dynamic_cast<TH2*>(dataSet->createHistogram(name + "_data",vars[i],RooFit::YVar(vars[j])));
//            h2_data->SetOption("colz");
//            h2_data->SetStats(kFALSE);
//            h2_data->SetMinimum(0);
//            h2_data->SetMaximum(h2_pdf->GetMaximum());
//            h2_pdf->Draw();
//            h2_pdf->Write();
//            path = dir + name + "_pdf" + format;
//            c1->SaveAs(path);
//            h2_data->Draw();
//            h2_data->Write();
//            path = dir + name + "_data" + format;
//            c1->SaveAs(path);
//            delete h2_pdf;
//            delete h2_data;
//        }
//    }

    file->Close();
}

void ConvertBetweenHelAndTrans(Double_t* par_input) {
    /// The variables are named as if converting parameters from helicity to transversity
    /// but the transformation is symmetric and can be used to convert from transversity
    /// to helicity as well.

    TComplex hp(par_input[0],par_input[1],true);
    TComplex h0(par_input[2],0,true);
    TComplex hm(sqrt(1-hp.Rho2()-h0.Rho2()),par_input[3],true);

    TComplex hRhop(par_input[5],par_input[8],true);
    TComplex hRho0(par_input[6],par_input[9],true);
    TComplex hRhom(par_input[7],par_input[10],true);

    TComplex hps = hRhop * hp;
    TComplex h0s = hRho0 * h0;
    TComplex hms = hRhom * hm;

    TComplex ap = (hp + hm)/sqrt(2);
    TComplex a0 = h0;
    TComplex at = (hp - hm)/sqrt(2);

    TComplex aps = (hps + hms)/sqrt(2);
    TComplex a0s = h0s;
    TComplex ats = (hps - hms)/sqrt(2);

    TComplex tRhop = aps/ap;
    TComplex tRho0 = a0s/a0;
    TComplex tRhot = ats/at;

    printf("original hel:\t");
    printf("hp = %.3f\thpa = %.2f\th0 = %.3f\thma = %.2f\tphiw = %.4f\trp = %.3f\tr0 = %.3f\trm = %.3f\tsp = %.3f\ts0 = %.3f\tsm = %.3f\n",\
           par_input[0],par_input[1],par_input[2],par_input[3],par_input[4],par_input[5],par_input[6],par_input[7],par_input[8],par_input[9],par_input[10]);

    par_input[0] = ap.Rho();
    par_input[1] = ap.Theta();
    par_input[2] = a0.Rho();
    par_input[3] = at.Theta();
    par_input[5] = tRhop.Rho();
    par_input[6] = tRho0.Rho();
    par_input[7] = tRhot.Rho();
    par_input[8] = tRhop.Theta();
    par_input[9] = tRho0.Theta();
    par_input[10] = tRhot.Theta();

    printf("converted trans:");
    printf("ap = %.3f\tapa = %.2f\ta0 = %.3f\tata = %.2f\tphiw = %.4f\trp = %.3f\tr0 = %.3f\trt = %.3f\tsp = %.3f\ts0 = %.3f\tst = %.3f\n",\
           par_input[0],par_input[1],par_input[2],par_input[3],par_input[4],par_input[5],par_input[6],par_input[7],par_input[8],par_input[9],par_input[10]);

}

Double_t Round(Double_t number, Int_t digits) {
    number = number * pow(10,digits);
    if(fmod(number,1)>0.5)
        number = ceil(number);
    else
        number = floor(number);

    return number/pow(10,digits);
}
