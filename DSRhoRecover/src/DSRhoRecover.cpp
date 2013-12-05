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

    if(argc != 2){
        printf("ERROR: Wrong number of arguments.\n");
        printf("Usage: DSRhoRecover FILE\n");
        return 85;
    }

    const char* file = argv[1];

//    gStyle->SetOptStat(0);
//    gStyle->SetPaintTextFormat("3.2g");
//    gStyle->SetMarkerSize(1);
    gEnv->SetValue("Canvas.PrintDirectory","plots");

    ObservablesCollection c;
    RooDataSet* dataset = RooDataSet::read(file,c.CreateArgList());
    dataset = static_cast<RooDataSet*> (dataset->reduce("xpe < 0.06"));
    dataset = static_cast<RooDataSet*> (dataset->reduce("xbpe < 0.06"));
    c.BindToDataset(dataset);

    RooRealVar xp("xp","xp",0.3,-0.1,0.1);
    RooRealVar x0("x0","x0",0.3,-0.1,0.1);
    RooRealVar xt("xt","xt",0.3,-0.1,0.1);
    RooRealVar yp("yp","yp",0.3,-0.1,0.1);
    RooRealVar y0("y0","y0",0.3,-0.1,0.1);
    RooRealVar yt("yt","yt",0.3,-0.1,0.1);

    RooRealVar xbp("xbp","xbp",0.3,-0.1,0.1);
    RooRealVar xb0("xb0","xb0",0.3,-0.1,0.1);
    RooRealVar xbt("xbt","xbt",0.3,-0.1,0.1);
    RooRealVar ybp("ybp","ybp",0.3,-0.1,0.1);
    RooRealVar yb0("yb0","yb0",0.3,-0.1,0.1);
    RooRealVar ybt("ybt","ybt",0.3,-0.1,0.1);

    RooRealVar xpe("xpe","xpe",0.3,0,1);
    RooRealVar x0e("x0e","x0e",0.3,0,1);
    RooRealVar xte("xte","xte",0.3,0,1);
    RooRealVar ype("ype","ype",0.3,0,1);
    RooRealVar y0e("y0e","y0e",0.3,0,1);
    RooRealVar yte("yte","yte",0.3,0,1);

    RooRealVar xbpe("xbpe","xbpe",0.3,0,1);
    RooRealVar xb0e("xb0e","xb0e",0.3,0,1);
    RooRealVar xbte("xbte","xbte",0.3,0,1);
    RooRealVar ybpe("ybpe","ybpe",0.3,0,1);
    RooRealVar yb0e("yb0e","yb0e",0.3,0,1);
    RooRealVar ybte("ybte","ybte",0.3,0,1);

    RooRealVar rp("rp","rp",0.01,0,1);
    RooRealVar r0("r0","r0",0.01,0,1);
    RooRealVar rt("rt","rt",0.01,0,1);

    RooRealVar sp("sp","sp",0.3,-PI,PI);
    RooRealVar s0("s0","s0",0.3,-PI,PI);
    RooRealVar st("st","st",0.3,-PI,PI);

    RooRealVar phiw("phiw","phiw",1.79,0,2*PI);

    MyPDF pdf("pdf","pdf",xp,x0,xt,yp,y0,yt,xbp,xb0,xbt,ybp,yb0,ybt,\
              xpe,x0e,xte,ype,y0e,yte,xbpe,xb0e,xbte,ybpe,yb0e,ybte,\
              rp,r0,rt,sp,s0,st,phiw);

    RooFitResult* result;
    result = pdf.fitTo(*dataset,RooFit::Save(),RooFit::Timer(true),RooFit::Minimizer("Minuit2"),RooFit::Minos(0),RooFit::Hesse(1),RooFit::Strategy(1),RooFit::NumCPU(1));



//
//    RooArgSet cartesianVars(xp,x0,xt,yp,y0,yt);
//    cartesianVars.add(RooArgSet(xbp,xb0,xbt,ybp,yb0,ybt));
//
//    RooDataSet* dataset_new = new RooDataSet("dataset_new","dataset_new",cartesianVars,RooFit::StoreError(cartesianVars));
//
//    for(int entry_no = 0; entry_no < dataset->numEntries(); entry_no++) {
//        dataset->get(entry_no);
//
//        xp.setVal(c.xp->getValV());
//        xp.setError(c.xpe->getValV());
//        x0.setVal(c.x0->getValV());
//        x0.setError(c.x0e->getValV());
//        xt.setVal(c.xt->getValV());
//        xt.setError(c.xte->getValV());
//        yp.setVal(c.yp->getValV());
//        yp.setError(c.ype->getValV());
//        y0.setVal(c.y0->getValV());
//        y0.setError(c.y0e->getValV());
//        yt.setVal(c.yt->getValV());
//        yt.setError(c.yte->getValV());
//
//        xbp.setVal(c.xbp->getValV());
//        xbp.setError(c.xbpe->getValV());
//        xb0.setVal(c.xb0->getValV());
//        xb0.setError(c.xb0e->getValV());
//        xbt.setVal(c.xbt->getValV());
//        xbt.setError(c.xbte->getValV());
//        ybp.setVal(c.ybp->getValV());
//        ybp.setError(c.ybpe->getValV());
//        yb0.setVal(c.yb0->getValV());
//        yb0.setError(c.yb0e->getValV());
//        ybt.setVal(c.ybt->getValV());
//        ybt.setError(c.ybte->getValV());
//
//        dataset_new->add(cartesianVars);
//    }
//
//    TCanvas* myCanvas = new TCanvas("myCanvas","My Canvas",1000,1000);
//    TH1* hh_data = dataset_new->createHistogram("xp,yp",20,20) ;
//    hh_data->Draw("BOX");





#ifdef GRAPHIC
    printf("\nProgram execution has finished.\n");
    rootapp->Run();
#endif
    return 0;
}












