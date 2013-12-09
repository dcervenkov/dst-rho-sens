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
#include "RooGenericPdf.h"
#include "RooFitResult.h"
#include "RooGaussian.h"

#include "ObservablesCollectionCartesian.h"
#include "DSRhoRecover.h"
#include "MyPDF.h"

//#define GRAPHIC

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

    if(argc != 3){
        printf("ERROR: Wrong number of arguments.\n");
        printf("Usage: DSRhoRecover OUTPUT-FILE INPUT-FILE\n");
        return 85;
    }

    const char* outputFile = argv[1];
    const char* inputFile = argv[2];

//    gStyle->SetOptStat(0);
//    gStyle->SetPaintTextFormat("3.2g");
//    gStyle->SetMarkerSize(1);
    gEnv->SetValue("Canvas.PrintDirectory","plots");

    ObservablesCollection c;
    RooDataSet* dataset = RooDataSet::read(inputFile,c.CreateArgList());
    dataset = static_cast<RooDataSet*> (dataset->reduce("xpe < 0.06"));
    dataset = static_cast<RooDataSet*> (dataset->reduce("xbpe < 0.06"));
    c.BindToDataset(dataset);

    RooRealVar rp("rp","rp",0.01,0,1);
    RooRealVar r0("r0","r0",0.01,0,1);
    RooRealVar rt("rt","rt",0.01,0,1);

    RooRealVar sp("sp","sp",-PI,-2*PI,2*PI);
    RooRealVar s0("s0","s0",-PI,-2*PI,2*PI);
    RooRealVar st("st","st",-PI,-2*PI,2*PI);

    RooRealVar phiw("phiw","phiw",1.79371,0,2*PI);

    MyPDF pdf("pdf","pdf",*c.xp,*c.x0,*c.xt,*c.yp,*c.y0,*c.yt,*c.xbp,*c.xb0,*c.xbt,*c.ybp,*c.yb0,*c.ybt,\
              *c.xpe,*c.x0e,*c.xte,*c.ype,*c.y0e,*c.yte,*c.xbpe,*c.xb0e,*c.xbte,*c.ybpe,*c.yb0e,*c.ybte,\
              *c.rp,*c.r0,*c.rt,*c.sp,*c.s0,*c.st,*c.phiw);

    RooFitResult* result;
    RooDataSet* datasetLine = new RooDataSet("datasetLine","datasetLine",c.CreateArgList());

    FILE* pFile;
    pFile = fopen(outputFile,"w");
    if(pFile == NULL) {
        printf("ERROR: couldn't open file %s for writing!\n",outputFile);
    }

    for(int line = 0; line < dataset->numEntries(); line++) {
        const RooArgSet* row = dataset->get(line);
        datasetLine->add(*row);
        c.CalculateInitialPolarVals();
        c.InitializePolarVals();
        result = pdf.fitTo(*datasetLine,RooFit::Save(),RooFit::Timer(true),RooFit::Minimizer("Minuit2"),\
                           RooFit::Minos(0),RooFit::Hesse(1),RooFit::Strategy(1),RooFit::NumCPU(1));
        datasetLine->reset();
        c.SaveParameters(pFile);
    }

    fclose(pFile);

#ifdef GRAPHIC
    printf("\nProgram execution has finished.\n");
    rootapp->Run();
#endif
    return 0;
}













