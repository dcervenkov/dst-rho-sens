#include "TROOT.h"
#include "TApplication.h"
#include "TCanvas.h"
#include "TFile.h"
#include "RooFit.h"
#include "RooPlot.h"
#include "RooCurve.h"
#include "RooRealVar.h"
#include "RooHist.h"
#include "RooPolynomial.h"

int sum(char* target, TList* sourceList);

int main( int argc, char **argv ) {

    //TApplication* rootapp = new TApplication("sum",&argc,argv);
    //argc = rootapp->Argc();
    //for(int i = 0; i < argc; i++) {
    //    argv[i] = rootapp->Argv(i);
    //}

    TList* fileList = new TList();

    //if ( argc < 4 ) {
    //    printf("Usage: %s <target> <source1> <source2> ...\n",argv[0]);
    //    printf("supply at least two source files for this to make sense... ;-)\n");
    //    return (-1);
    //}

    printf("Target file: %s\n",argv[1]);

    for ( int i = 2; i < argc; i++ ) {
        printf("Source file %i: %s\n",i-1,argv[i]);
        fileList->Add( TFile::Open( argv[i] ) );
    }
    sum( argv[1], fileList );

    //printf("\nProgram execution has finished.\n");
    //rootapp->Run();
}


int sum(char* target, TList* sourceList){
    TFile* nextSource = (TFile*)sourceList->First();

    TFile* previousSource;

    RooPlot* nextPlot = 0;
    RooPlot* residPlot = 0;
    RooHist* nextHist = 0;
    RooHist* sumHist = 0;
    RooHist* tempHist = 0;

    while (nextSource) {
	nextPlot = (RooPlot*)nextSource->Get("plot3");
	nextHist = (RooHist*)nextPlot->findObject("data");

	if (sumHist == 0){
	    sumHist = nextHist;
	    residPlot = nextPlot->emptyClone("residPlot");
	} else {
	    tempHist = new RooHist(*sumHist,*nextHist,1,1,RooAbsData::SumW2,1);
	    sumHist = tempHist;
	}

	//sumHist->Print("V");
	//nextPlot->Print("V");
	//residPlot->Print("V");

	previousSource = nextSource;
	nextSource = (TFile*)sourceList->After( nextSource );
	if (nextSource)
	    previousSource->Close();
    }

    TFile* outputFile = TFile::Open( target, "RECREATE" );

    TCanvas* c2 = new TCanvas("c2","c2",800,1000);
    c2->Divide(1,2);
    c2->cd(1);

    //gDirectory->cd(target);
    //sumHist->Print("V");
    //residPlot->Print("V");
    residPlot->addPlotable(sumHist,"P");
    
    RooRealVar dt("dt","dt",-3,+3);
    RooPolynomial p("p","p",dt);
    p.plotOn(residPlot,RooFit::Normalization(0,RooAbsReal::Raw));
    residPlot->Draw();
    residPlot->Write();
    
    c2->cd(2);
    RooHist* pullHist = residPlot->pullHist();
    RooPlot* pullPlot = residPlot->emptyClone("pullPlot");
    pullPlot->addPlotable(pullHist,"P");
    p.plotOn(pullPlot,RooFit::Normalization(0,RooAbsReal::Raw));
    pullPlot->Draw();
    pullPlot->Write();

    c2->SaveAs("tpulls.gif");
    c2->SaveAs("tpulls.pdf");
    
    previousSource->Close();
    outputFile->Close();

}

