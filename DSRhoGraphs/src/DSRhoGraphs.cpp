#include <stdio.h>
#include <vector>

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
#include "TString.h"

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

int main(int argc, char* argv[]) {

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
    for(int i = 0; i < argc; i++) {
        argv[i] = rootapp->Argv(i);
    }
#endif // GRAPHIC

    if(argc != 2 || argc != 3){
        printf("ERROR: Wrong number of arguments.\n");
        printf("Usage: DSRhoGraph successfulResultsFile [allResultsFile]\n");
    }

    gStyle->SetOptStat(0);
    gStyle->SetPaintTextFormat("3.2g");
    gStyle->SetMarkerSize(2);
    gEnv->SetValue("Canvas.PrintDirectory","plots");

    ObservablesCollection c;
    RooDataSet* dataset = RooDataSet::read(argv[1],c.CreateArgList());
    c.BindToDataset(dataset);
    c.AdjustInputSForPeriodicity(dataset);
    c.CreateResidualsAndPulls(dataset);

    CreateGeneralPlots(dataset,c);

    /// If 2 arguments were supplied, create maps as well
    if(argc == 3){
        ObservablesCollection c_all;
        RooDataSet* dataset_all = RooDataSet::read(argv[2],c_all.CreateArgList());
        c_all.BindToDataset(dataset_all);
        c_all.AdjustInputSForPeriodicity(dataset_all);
        c_all.CreateResidualsAndPulls(dataset_all);

        std::vector<TH2D*> maps_ok = CreateMaps(dataset);
        std::vector<TH2D*> maps_all = CreateMaps(dataset_all);
        std::vector<TH2D*> maps_div = CreateMaps(dataset);

        TCanvas* c_maps = new TCanvas("c_maps","c_maps",1000,1000);
        c_maps->Divide(3,3);
        for(int i = 0; i < maps_ok.size(); i++) {
            c_maps->cd(i+1);
            maps_ok[i]->Draw("TEXT");
            c_maps->cd(i+1+3);
            maps_all[i]->Draw("TEXT");
            c_maps->cd(i+1+6);
            maps_div[i]->Divide(maps_all[i]);
            maps_div[i]->Draw("TEXT");
        }
        c_maps->SaveAs(".gif");
    }

#ifdef GRAPHIC
    printf("\nProgram execution has finished.\n");
    rootapp->Run();
#endif
    return 0;
}

std::vector<TH2D*> CreateMaps(const RooDataSet* const dataset) {
    const RooArgSet* args = dataset->get();
    RooArgList vars(*args->find("spi"),*args->find("s0i"),*args->find("sti"));
    std::vector<TH2D*> maps((vars.getSize()*(vars.getSize()-1))/2);
    TString name;
    int counter = 0;
    for(int j = 0; j < vars.getSize(); j++) {
        for(int i = 0; i < j; i++) {
            name = "map_";
            name += vars[i].GetName();
            name += "_";
            name += vars[j].GetName();
            maps[counter] = new TH2D(name,name,6,(dynamic_cast<RooRealVar*>(&vars[i]))->getMin(),(dynamic_cast<RooRealVar*>(&vars[i]))->getMax() \
                                     ,6,(dynamic_cast<RooRealVar*>(&vars[j]))->getMin(),(dynamic_cast<RooRealVar*>(&vars[j]))->getMax());
            maps[counter]->GetXaxis()->SetTitle(vars[i].GetName());
            maps[counter]->GetYaxis()->SetTitle(vars[j].GetName());
            counter++;
        }
    }
    for(int i = 0; i < dataset->sumEntries(); i++) {
        args = dataset->get(i);
        counter = 0;
        for(int j = 0; j < vars.getSize(); j++) {
            for(int i = 0; i < j; i++) {
                maps[counter++]->Fill((Double_t)((dynamic_cast<RooRealVar*>(&vars[i]))->getValV()),(Double_t)((dynamic_cast<RooRealVar*>(&vars[j]))->getValV()),1);
//                printf("%s: %f, %s: %f, w: %i\n",(dynamic_cast<RooRealVar*>(&vars[i]))->GetName(),(Double_t)((dynamic_cast<RooRealVar*>(&vars[i]))->getValV()),\
//                       (dynamic_cast<RooRealVar*>(&vars[j]))->GetName(),(Double_t)((dynamic_cast<RooRealVar*>(&vars[j]))->getValV()),1);
            }
        }
    }
    return maps;

}

void CreateGeneralPlots(RooDataSet* const dataset, const ObservablesCollection c) {
    Int_t plotMap[] = {1,4,2,5,3,6};
    RooGaussian* gaus[13];
    RooRealVar mean("mean","mean",0,-1,1);
    RooRealVar sigma("sigma","sigma",1,-10,10);;
    RooPlot* frame;

    TCanvas* c_amp_pulls = new TCanvas("c_amp_pulls","Amplitude pulls",1280,800);
    c_amp_pulls->Divide(3,2);
    for(int i = 0; i < 6; i++) {
        c_amp_pulls->cd(plotMap[i]);
        frame = c.pulls[i]->frame(RooFit::AutoRange(*dataset));
        if(frame->GetXaxis()->GetXmin() == frame->GetXaxis()->GetXmax()) {
            frame->GetXaxis()->SetLimits(frame->GetXaxis()->GetXmax()-1,frame->GetXaxis()->GetXmax()+1);
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
    for(int i = 0; i < 6; i++) {
        c_amp_residuals->cd(plotMap[i]);
        frame = c.residuals[i]->frame(RooFit::AutoRange(*dataset));
        if(frame->GetXaxis()->GetXmin() == frame->GetXaxis()->GetXmax()) {
            frame->GetXaxis()->SetLimits(frame->GetXaxis()->GetXmax()-1,frame->GetXaxis()->GetXmax()+1);
        }
        dataset->plotOn(frame,RooFit::DrawOption("BX"));
        frame->Draw("BX");
    }
    TCanvas* c_amp_errors = new TCanvas("c_amp_errors","Amplitude errors",1280,800);
    c_amp_errors->Divide(3,2);
    for(int i = 0; i < 6; i++) {
        c_amp_errors->cd(plotMap[i]);
        frame = c.errors[i]->frame(RooFit::AutoRange(*dataset));
        if(frame->GetXaxis()->GetXmin() == frame->GetXaxis()->GetXmax()) {
            frame->GetXaxis()->SetLimits(frame->GetXaxis()->GetXmax()-1,frame->GetXaxis()->GetXmax()+1);
        }
        dataset->plotOn(frame,RooFit::DrawOption("BX"));
        frame->Draw("BX");
    }

    TCanvas* c_phiw = new TCanvas("c_phiw","Weak phase",1280,400);
    c_phiw->Divide(3);
    c_phiw->cd(1);
    frame = c.pull_phiw->frame(RooFit::AutoRange(*dataset));
    if(frame->GetXaxis()->GetXmin() == frame->GetXaxis()->GetXmax()) {
        frame->GetXaxis()->SetLimits(frame->GetXaxis()->GetXmax()-1,frame->GetXaxis()->GetXmax()+1);
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
    if(frame->GetXaxis()->GetXmin() == frame->GetXaxis()->GetXmax()) {
        frame->GetXaxis()->SetLimits(frame->GetXaxis()->GetXmax()-1,frame->GetXaxis()->GetXmax()+1);
    }
    dataset->plotOn(frame,RooFit::DrawOption("BX"));
    frame->Draw();

    c_phiw->cd(3);
    frame = c.phiwe->frame(RooFit::AutoRange(*dataset));
    if(frame->GetXaxis()->GetXmin() == frame->GetXaxis()->GetXmax()) {
        frame->GetXaxis()->SetLimits(frame->GetXaxis()->GetXmax()-1,frame->GetXaxis()->GetXmax()+1);
    }
    dataset->plotOn(frame,RooFit::DrawOption("BX"));
    frame->Draw();
//
    TCanvas* c_rs_pulls = new TCanvas("c_rs_pulls","r and strong phase pulls",1280,800);
    c_rs_pulls->Divide(3,2);
    for(int i = 7; i < 13; i++) {
        c_rs_pulls->cd(i-6);
        frame = c.pulls[i]->frame(RooFit::AutoRange(*dataset));
        if(frame->GetXaxis()->GetXmin() == frame->GetXaxis()->GetXmax()) {
            frame->GetXaxis()->SetLimits(frame->GetXaxis()->GetXmax()-1,frame->GetXaxis()->GetXmax()+1);
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
    for(int i = 7; i < 13; i++) {
        c_rs_residuals->cd(i-6);
        frame = c.residuals[i]->frame(RooFit::AutoRange(*dataset));
        if(frame->GetXaxis()->GetXmin() == frame->GetXaxis()->GetXmax()) {
            frame->GetXaxis()->SetLimits(frame->GetXaxis()->GetXmax()-1,frame->GetXaxis()->GetXmax()+1);
        }
        dataset->plotOn(frame,RooFit::DrawOption("BX"));
        frame->Draw("BX");
    }
    TCanvas* c_rs_errors = new TCanvas("c_rs_errors","r and strong phase errors",1280,800);
    c_rs_errors->Divide(3,2);
    for(int i = 7; i < 13; i++) {
        c_rs_errors->cd(i-6);
        frame = c.errors[i]->frame(RooFit::AutoRange(*dataset));
        if(frame->GetXaxis()->GetXmin() == frame->GetXaxis()->GetXmax()) {
            frame->GetXaxis()->SetLimits(frame->GetXaxis()->GetXmax()-1,frame->GetXaxis()->GetXmax()+1);
        }
        dataset->plotOn(frame,RooFit::DrawOption("BX"));
        frame->Draw("BX");
    }

    c_amp_pulls->SaveAs(".gif");
    c_amp_residuals->SaveAs(".gif");
    c_amp_errors->SaveAs(".gif");
    c_phiw->SaveAs(".gif");
    c_rs_pulls->SaveAs(".gif");
    c_rs_residuals->SaveAs(".gif");
    c_rs_errors->SaveAs(".gif");

    return;
}












