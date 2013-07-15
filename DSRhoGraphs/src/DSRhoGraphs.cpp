#include <stdio.h>

#include "TROOT.h"
#include "TEnv.h"
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
#include "RooArgList.h"
#include "RooGaussian.h"

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
//    RooDataSet* dataset = RooDataSet::read("data/fit7a.res",c.CreateArgList());
//    RooDataSet* dataset = RooDataSet::read("data/f",c.CreateArgList());
    RooDataSet* dataset = RooDataSet::read("data/fit_gen7.res",c.CreateArgList());
    c.BindToDataset(dataset);
    c.AdjustInputSForPeriodicity(dataset);
    c.CreateResidualsAndPulls(dataset);

    Int_t plotMap[] = {1,4,2,5,3,6};
    RooGaussian* gaus[13];
    RooRealVar mean("mean","mean",0,-1,1);
    RooRealVar sigma("sigma","sigma",1,-10,10);;
    RooPlot* frame;
    TCanvas* c_amp_pulls = new TCanvas("c_amp_pulls","Amplitude pulls",1280,800);
    c_amp_pulls->Divide(3,2);
    for (int i = 0; i < 6; i++){
        c_amp_pulls->cd(plotMap[i]);
        frame = c.pulls[i]->frame(RooFit::AutoRange(*dataset));
        if(frame->GetXaxis()->GetXmin() == frame->GetXaxis()->GetXmax()){
            frame->GetXaxis()->SetLimits(-1,1);
            dataset->plotOn(frame,RooFit::DrawOption("BX"));
        } else {
            dataset->plotOn(frame,RooFit::DrawOption("BX"));
            gaus[i] = new RooGaussian((TString("gaus") + TString(i)).Data(),(TString("gaus") + TString(i)).Data(),*(c.pulls[i]),mean,sigma);
            gaus[i]->fitTo(*dataset);
            gaus[i]->paramOn(frame);
            gaus[i]->plotOn(frame);
        }
        frame->Draw();
    }
    TCanvas* c_amp_residuals = new TCanvas("c_amp_residuals","Amplitude residuals",1280,800);
    c_amp_residuals->Divide(3,2);
    for (int i = 0; i < 6; i++){
        c_amp_residuals->cd(plotMap[i]);
        frame = c.residuals[i]->frame(RooFit::AutoRange(*dataset));
        if(frame->GetXaxis()->GetXmin() == frame->GetXaxis()->GetXmax()){
            frame->GetXaxis()->SetLimits(-1,1);
        }
        dataset->plotOn(frame,RooFit::DrawOption("BX"));
        frame->Draw("BX");
    }
    TCanvas* c_amp_errors = new TCanvas("c_amp_errors","Amplitude errors",1280,800);
    c_amp_errors->Divide(3,2);
    for (int i = 0; i < 6; i++){
        c_amp_errors->cd(plotMap[i]);
        frame = c.errors[i]->frame(RooFit::AutoRange(*dataset));
        if(frame->GetXaxis()->GetXmin() == frame->GetXaxis()->GetXmax()){
            frame->GetXaxis()->SetLimits(-1,1);
        }
        dataset->plotOn(frame,RooFit::DrawOption("BX"));
        frame->Draw("BX");
    }

    TCanvas* c_phiw = new TCanvas("c_phiw","Weak phase",1280,400);
    c_phiw->Divide(3);
    c_phiw->cd(1);
    frame = c.pull_phiw->frame(RooFit::AutoRange(*dataset));
    if(frame->GetXaxis()->GetXmin() == frame->GetXaxis()->GetXmax()){
        frame->GetXaxis()->SetLimits(-1,1);
        dataset->plotOn(frame,RooFit::DrawOption("BX"));
    } else {
        dataset->plotOn(frame,RooFit::DrawOption("BX"));
        gaus[6] = new RooGaussian("gaus_phiw","gaus_phiw",*(c.pull_phiw),mean,sigma);
        gaus[6]->fitTo(*dataset);
        gaus[6]->paramOn(frame);
        gaus[6]->plotOn(frame);
    }
    frame->Draw();

    c_phiw->cd(2);
    frame = c.residual_phiw->frame(RooFit::AutoRange(*dataset));
    if(frame->GetXaxis()->GetXmin() == frame->GetXaxis()->GetXmax()){
            frame->GetXaxis()->SetLimits(-1,1);
        }
    dataset->plotOn(frame,RooFit::DrawOption("BX"));
    frame->Draw();

    c_phiw->cd(3);
    frame = c.phiwe->frame(RooFit::AutoRange(*dataset));
    if(frame->GetXaxis()->GetXmin() == frame->GetXaxis()->GetXmax()){
            frame->GetXaxis()->SetLimits(-1,1);
        }
    dataset->plotOn(frame,RooFit::DrawOption("BX"));
    frame->Draw();
//
    TCanvas* c_rs_pulls = new TCanvas("c_rs_pulls","r and strong phase pulls",1280,800);
    c_rs_pulls->Divide(3,2);
    for (int i = 7; i < 13; i++){
        c_rs_pulls->cd(i-6);
        frame = c.pulls[i]->frame(RooFit::AutoRange(*dataset));
        if(frame->GetXaxis()->GetXmin() == frame->GetXaxis()->GetXmax()){
            frame->GetXaxis()->SetLimits(-1,1);
            dataset->plotOn(frame,RooFit::DrawOption("BX"));
        } else {
            dataset->plotOn(frame,RooFit::DrawOption("BX"));
            gaus[i] = new RooGaussian((TString("gaus") + TString(i)).Data(),(TString("gaus") + TString(i)).Data(),*(c.pulls[i]),mean,sigma);
            gaus[i]->fitTo(*dataset);
            gaus[i]->paramOn(frame);
            gaus[i]->plotOn(frame);
        }
        frame->Draw();
    }
    TCanvas* c_rs_residuals = new TCanvas("c_rs_residuals","r and strong phase residuals",1280,800);
    c_rs_residuals->Divide(3,2);
    for (int i = 7; i < 13; i++){
        c_rs_residuals->cd(i-6);
        frame = c.residuals[i]->frame(RooFit::AutoRange(*dataset));
        if(frame->GetXaxis()->GetXmin() == frame->GetXaxis()->GetXmax()){
            frame->GetXaxis()->SetLimits(-1,1);
        }
        dataset->plotOn(frame,RooFit::DrawOption("BX"));
        frame->Draw("BX");
    }
    TCanvas* c_rs_errors = new TCanvas("c_rs_errors","r and strong phase errors",1280,800);
    c_rs_errors->Divide(3,2);
    for (int i = 7; i < 13; i++){
        c_rs_errors->cd(i-6);
        frame = c.errors[i]->frame(RooFit::AutoRange(*dataset));
        if(frame->GetXaxis()->GetXmin() == frame->GetXaxis()->GetXmax()){
            frame->GetXaxis()->SetLimits(-1,1);
        }
        dataset->plotOn(frame,RooFit::DrawOption("BX"));
        frame->Draw("BX");
    }

    gEnv->SetValue("Canvas.PrintDirectory","plots");
    c_amp_pulls->SaveAs(".gif");
    c_amp_residuals->SaveAs(".gif");
    c_amp_errors->SaveAs(".gif");
    c_phiw->SaveAs(".gif");
    c_rs_pulls->SaveAs(".gif");
    c_rs_residuals->SaveAs(".gif");
    c_rs_errors->SaveAs(".gif");

    #ifdef GRAPHIC
    printf("\nProgram execution has finished.\n");
    rootapp->Run();
    #endif
    return 0;
}

