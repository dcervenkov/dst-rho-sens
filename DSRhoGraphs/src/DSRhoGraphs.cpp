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
#include "RooDataHist.h"
#include "RooArgSet.h"
#include "RooArgList.h"
#include "RooGaussian.h"

//#include "ObservablesCollection.h"
#include "ObservablesCollectionCartesian.h"
#include "DSRhoGraphs.h"

#define GRAPHIC

const int num_cmd_options = 1;

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

    if(argc < 3){
        printf("ERROR: Wrong number of arguments.\n");
        printf("Usage: DSRhoGraph MODE FILE-1 [FILE-2]...\n");
        return 85;
    }

    gStyle->SetOptStat(0);
    gStyle->SetPaintTextFormat("3.2g");
    gStyle->SetMarkerSize(2);
    gEnv->SetValue("Canvas.PrintDirectory","plots");

    if(!strcmp(argv[1],"1")){
        if(argc == num_cmd_options + 3)
            CreatePullsAndMaps(argv[2],argv[3]);
        else
            CreatePullsAndMaps(argv[2],nullptr);
    }
    else if(!strcmp(argv[1],"2"))
        CreateErrorProgression(argc,argv);
    else
        printf("ERROR: Unknown mode '%s'.\n",argv[1]);


#ifdef GRAPHIC
    printf("\nProgram execution has finished.\n");
    rootapp->Run();
#endif
    return 0;
}

double CalculateMedian(std::vector<double> vals) {
    double median;
    size_t size = vals.size();
    sort(vals.begin(), vals.end());
    if(size % 2 == 0) {
        median = (vals[size/2 - 1] + vals[size/2]) / 2;
    } else {
        median = vals[size/2];
    }

    return median;
}

void CreateErrorProgression(int argc, char* argv[]){
    if(argc < num_cmd_options + 3){
        printf("ERROR: Error Progression mode requires at least 2 files.");
        return;
    }
    int num_files = argc-1 - num_cmd_options ;
    ObservablesCollection c;
    RooDataSet* dataset;

//    TH1D* h_errors[7];
//    for(int error_no = 0; error_no < 7; error_no++){
//        TString name("h_errors");
//        name += error_no+1;
//        h_errors = new TH1D(name,name,)
//    }
    std::vector<double> mean_errors[num_files];

    for(int file_no = 0; file_no < num_files; file_no++){
        dataset = RooDataSet::read(argv[file_no + num_cmd_options+1], c.CreateArgList());
        c.BindToDataset(dataset);
        //c.AdjustInputSForPeriodicity(dataset);

        const RooArgSet* args = dataset->get();
        RooArgList error_vars(*args->find("phiwe"),*args->find("rpe"),*args->find("r0e"),*args->find("rte"),*args->find("spe"),*args->find("s0e"),*args->find("ste"));


        std::vector<double> errors[error_vars.getSize()];

        for(int entry_no = 0; entry_no < dataset->numEntries(); entry_no++){
            dataset->get(entry_no);
            for(int error_no = 0; error_no < error_vars.getSize();error_no++){
                errors[error_no].push_back(dynamic_cast<RooRealVar*>(&error_vars[error_no])->getValV());
            }
        }

        for(int error_no = 0; error_no < error_vars.getSize(); error_no++){
            mean_errors[file_no].push_back(CalculateMedian(errors[error_no]));
        }
        //delete dataset;
    }

    TH2D* h_errors[mean_errors[0].size()];
    for(int error_no = 0; error_no < mean_errors[0].size(); error_no++){
        TString name("h_errors_");
        name += error_no+1;
        double max_error = 0;
        for(int file_no = 0; file_no < num_files; file_no++){
            if(max_error < mean_errors[file_no][error_no])
                max_error = mean_errors[file_no][error_no];
        }
        h_errors[error_no] = new TH2D(name,name,100,0,num_files+1,100,0,max_error*1.1);
        for(int file_no = 0; file_no < num_files; file_no++){
            h_errors[error_no]->Fill(file_no+1,mean_errors[file_no][error_no]);
        }
    }
    TCanvas* c_errors = new TCanvas("c_errors","Errors",1200,1000);
    c_errors->Divide(3,3);
    for(int error_no = 0; error_no < mean_errors[0].size(); error_no++){
        c_errors->cd(error_no+1);
        h_errors[error_no]->SetMarkerStyle(kPlus);
        h_errors[error_no]->SetMarkerSize(1);
        h_errors[error_no]->Draw();
    }

}


void CreatePullsAndMaps(char* successful_file, char* all_file) {
    ObservablesCollection c;
    RooDataSet* dataset = RooDataSet::read(successful_file,c.CreateArgList());
    dataset = static_cast<RooDataSet*> (dataset->reduce("phiwe < 1"));
    c.BindToDataset(dataset);
    //c.AdjustInputSForPeriodicity(dataset);
    c.CreateResidualsAndPulls(dataset);

    CreateGeneralPlots(dataset,c);

    /// If 2 arguments were supplied, create maps as well
    if(all_file != nullptr){
        ObservablesCollection c_all;
        RooDataSet* dataset_all = RooDataSet::read(all_file,c_all.CreateArgList());
        c_all.BindToDataset(dataset_all);
        //c_all.AdjustInputSForPeriodicity(dataset_all);
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
            gaus[i]->plotOn(frame,RooFit::LineWidth(2));
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
        gaus[6]->plotOn(frame,RooFit::LineWidth(2));
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
            gaus[i]->plotOn(frame,RooFit::LineWidth(2));
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












