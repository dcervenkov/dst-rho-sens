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

#include "DSRhoGraphs.h"

int main() {
    RooDataSet* dataset = readDataSet("fit5.res");
    const RooArgSet* args = dataset->get();
    args->Print();
    return 0;
}

RooDataSet* readDataSet(const char* filename){
    RooRealVar chi2a ("chi2a","chi2a",0,10);
    RooRealVar chi2ab ("chi2ab","chi2ab",0,10);
    RooRealVar chi2b ("chi2b","chi2b",0,10);
    RooRealVar chi2bb ("chi2bb","chi2bb",0,10);

    RooRealVar ap ("ap","ap",0,1);
    RooRealVar api ("api","api",0,1);
    RooRealVar ape ("ape","ape",0,1);
    RooRealVar apa("apa","apa",0,2*PI);
    RooRealVar apai("apai","apai",0,2*PI);
    RooRealVar apae("apae","apae",0,2*PI);
    RooRealVar a0 ("a0","a0",0,1);
    RooRealVar a0i ("a0i","a0i",0,1);
    RooRealVar a0e ("a0e","a0e",0,1);
    RooRealVar a0a("a0a","a0a",0,2*PI);
    RooRealVar a0ai("a0ai","a0ai",0,2*PI);
    RooRealVar a0ae("a0ae","a0ae",0,2*PI);
    RooRealVar at ("at","at",0,1);
    RooRealVar ati ("ati","ati",0,1);
    RooRealVar ate ("ate","ate",0,1);
    RooRealVar ata("ata","ata",0,2*PI);
    RooRealVar atai("atai","atai",0,2*PI);
    RooRealVar atae("atae","atae",0,2*PI);

    RooRealVar phiw("phiw","phiw",0,2*PI);
    RooRealVar phiwi("phiwi","phiwi",0,2*PI);
    RooRealVar phiwe("phiwe","phiwe",0,2*PI);

    RooRealVar rp ("rp","rp",0,1);
    RooRealVar rpi ("rpi","rpi",0,1);
    RooRealVar rpe ("rpe","rpe",0,1);
    RooRealVar r0 ("r0","r0",0,1);
    RooRealVar r0i ("r0i","r0i",0,1);
    RooRealVar r0e ("r0e","r0e",0,1);
    RooRealVar rt ("rt","rt",0,1);
    RooRealVar rti ("rti","rti",0,1);
    RooRealVar rte ("rte","rte",0,1);

    RooRealVar sp ("sp","sp",-PI,PI);
    RooRealVar spi ("spi","spi",-PI,PI);
    RooRealVar spe ("spe","spe",0,20);
    RooRealVar s0 ("s0","s0",-PI,PI);
    RooRealVar s0i ("s0i","s0i",-PI,PI);
    RooRealVar s0e ("s0e","s0e",0,20);
    RooRealVar st ("st","st",-PI,PI);
    RooRealVar sti ("sti","sti",-PI,PI);
    RooRealVar ste ("ste","ste",0,20);

    RooCategory sep("sep","sep");
    sep.defineType("|",1);
    sep.defineType("||",2);

    RooDataSet* dataset;
    RooArgList argList(chi2a,chi2ab,chi2b,chi2bb,sep,sep);
    argList.add(RooArgList(api,ap,ape,sep,apai,apa,apae,sep));
    argList.add(RooArgList(a0i,a0,a0e,sep,a0ai,a0a,a0ae,sep));
    argList.add(RooArgList(ati,at,ate,sep,atai,ata,atae,sep));
    argList.add(RooArgList(sep,phiwi,phiw,phiwe,sep,sep));
    argList.add(RooArgList(rpi,rp,rpe,sep));
    argList.add(RooArgList(r0i,r0,r0e,sep));
    argList.add(RooArgList(rti,rt,rte,sep,sep));
    argList.add(RooArgList(spi,sp,spe,sep));
    argList.add(RooArgList(s0i,s0,s0e,sep));
    argList.add(RooArgList(sti,st,ste));
    dataset = RooDataSet::read(filename,argList);
    return dataset;
}
