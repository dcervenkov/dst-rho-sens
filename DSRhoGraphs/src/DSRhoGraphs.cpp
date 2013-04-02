#include <stdio.h>

#include "TROOT.h"
#include "TStyle.h"
#include "TApplication.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TTree.h"
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
#include "RooArgSet.h"

#include "ObservablesCollection.h"
#include "DSRhoGraphs.h"

#define GRAPHIC

int main(int argc, char* argv[]){

    #ifdef GRAPHIC
    TApplication* rootapp = new TApplication("DSRhoGraphs",&argc,argv); //For graphic-output apps only
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
    ObservablesCollection c;
    RooDataSet* dataset = RooDataSet::read("data/fit5.res",c.CreateArgList());
    c.BindToDataset(dataset);
    c.CreateResidualsAndPulls(dataset);
    Int_t plotMap[] = {1,4,2,5,3,6};
    RooPlot* frame;
    TCanvas* c_amp_pulls = new TCanvas("c_amp_pulls","Amplitude pulls",1280,800);
    c_amp_pulls->Divide(3,2);
    for (int i = 0; i < 6; i++){
        c_amp_pulls->cd(plotMap[i]);
        frame = c.pulls[i]->frame(RooFit::AutoRange(*dataset));
        dataset->plotOn(frame,RooFit::DrawOption("BX"));
        frame->Draw("BX");
    }
    TCanvas* c_amp_residuals = new TCanvas("c_amp_residuals","Amplitude residuals",1280,800);
    c_amp_residuals->Divide(3,2);
    for (int i = 0; i < 6; i++){
        c_amp_residuals->cd(plotMap[i]);
        frame = c.residuals[i]->frame(RooFit::AutoRange(*dataset));
        dataset->plotOn(frame,RooFit::DrawOption("BX"));
        frame->Draw("BX");
    }
    TCanvas* c_amp_errors = new TCanvas("c_amp_errors","Amplitude errors",1280,800);
    c_amp_errors->Divide(3,2);
    for (int i = 0; i < 6; i++){
        c_amp_errors->cd(plotMap[i]);
        frame = c.errors[i]->frame(RooFit::AutoRange(*dataset));
        dataset->plotOn(frame,RooFit::DrawOption("BX"));
        frame->Draw("BX");
    }
    #ifdef GRAPHIC
    printf("\nProgram execution has finished.\n");
    rootapp->Run();
    #endif
    return 0;
}

