#include "TROOT.h"
#include "TFile.h"
#include "TCanvas.h"
#include "RooFit.h"
#include "RooPlot.h"
#include "RooCurve.h"
#include "RooRealVar.h"
#include "RooHist.h"
#include "RooPolynomial.h"

int overlay(char* filename1, char* filename2, char* output_filename){
    TFile *file1 = new TFile(filename1);
    RooPlot* plot1 = (RooPlot*)file1->Get("proj_dt_b");
    TFile *file2 = new TFile(filename2);
    RooPlot* plot2 = (RooPlot*)file2->Get("proj_dt_b");

    TFile* output_file = new TFile(output_filename, "RECREATE");

    TCanvas* c1 = new TCanvas("c1","c1",800,600);

    RooCurve* curve2 = (RooCurve*)plot2->findObject("pdf_a_Int[phit,thb,tht]_Norm[dt,phit,thb,tht]");
    curve2->SetLineColor(kRed);
    RooHist* hist2 = (RooHist*)plot2->findObject("data");
    hist2->SetMarkerColor(kOrange);

    RooHist* hist1 = (RooHist*)plot1->findObject("data");

    RooCurve* curve2clone = (RooCurve*)curve2->Clone();
    RooHist* hist2clone = (RooHist*)hist2->Clone();
    plot1->addPlotable(hist2clone,"P");
    plot1->addPlotable(curve2clone,"L");
    plot1->Draw();
    plot1->Write();

    TCanvas* c2 = new TCanvas("c2","c2",800,1000);
    c2->Divide(1,2);
    c2->cd(1);
    RooHist* histdiff = new RooHist(*hist1,*hist2,1,-1,RooAbsData::SumW2,1);

    RooPlot* plot3 = plot1->emptyClone("plot3");
    plot3->addPlotable(histdiff,"P");
    RooRealVar dt("dt","dt",-3,+3);
    RooPolynomial p("p","p",dt);
    p.plotOn(plot3,RooFit::Normalization(0,RooAbsReal::Raw));
    plot3->Draw();
    plot3->Write();
    
    c2->cd(2);
    RooHist* pull = plot3->pullHist();
    RooPlot* plot4 = plot1->emptyClone("plot4");
    plot4->addPlotable(pull,"P");
    p.plotOn(plot4,RooFit::Normalization(0,RooAbsReal::Raw));
    plot4->Draw();
    plot4->Write();

    output_file->Close();
    file1->Close();
    file2->Close();
}

int overlay(){
    overlay("projections_evt75.root", "projections_gen75.root", "test.root");
}

int main(int argc,char *argv[]){
    //evt_file, gen_file
    overlay(argv[1],argv[2],argv[3]);
    //overlay();
}


