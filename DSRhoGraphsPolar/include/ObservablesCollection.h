#ifndef OBSERVABLESCOLLECTION_H
#define OBSERVABLESCOLLECTION_H

#include "RooFit.h"
#include "RooGlobalFunc.h"
#include "RooRealVar.h"
#include "RooCategory.h"
#include "RooPlot.h"
#include "RooDataSet.h"
#include "RooArgSet.h"

class ObservablesCollection
{
    public:
        ObservablesCollection();
        ~ObservablesCollection();

        RooArgList CreateArgList();
        void BindToDataset(RooDataSet*);
        void CreateResidualsAndPulls(RooDataSet* dataset);
        void AdjustInputSForPeriodicity(RooDataSet*& dataset);


        RooRealVar* chi2a;
        RooRealVar* chi2ab;
        RooRealVar* chi2b;
        RooRealVar* chi2bb;

        RooRealVar* ap;
        RooRealVar* api;
        RooRealVar* ape;
        RooRealVar* apa;
        RooRealVar* apai;
        RooRealVar* apae;
        RooRealVar* a0;
        RooRealVar* a0i;
        RooRealVar* a0e;
        RooRealVar* a0a;
        RooRealVar* a0ai;
        RooRealVar* a0ae;
        RooRealVar* at;
        RooRealVar* ati;
        RooRealVar* ate;
        RooRealVar* ata;
        RooRealVar* atai;
        RooRealVar* atae;

        RooRealVar* phiw;
        RooRealVar* phiwi;
        RooRealVar* phiwe;

        RooRealVar* rp;
        RooRealVar* rpi;
        RooRealVar* rpe;
        RooRealVar* r0;
        RooRealVar* r0i;
        RooRealVar* r0e;
        RooRealVar* rt;
        RooRealVar* rti;
        RooRealVar* rte;

        RooRealVar* sp;
        RooRealVar* spi;
        RooRealVar* spe;
        RooRealVar* s0;
        RooRealVar* s0i;
        RooRealVar* s0e;
        RooRealVar* st;
        RooRealVar* sti;
        RooRealVar* ste;

        RooRealVar* residual_ap;
        RooRealVar* residual_apa;
        RooRealVar* residual_a0;
        RooRealVar* residual_a0a;
        RooRealVar* residual_at;
        RooRealVar* residual_ata;

        RooRealVar* residual_phiw;

        RooRealVar* residual_rp;
        RooRealVar* residual_r0;
        RooRealVar* residual_rt;

        RooRealVar* residual_sp;
        RooRealVar* residual_s0;
        RooRealVar* residual_st;

        RooRealVar* pull_ap;
        RooRealVar* pull_apa;
        RooRealVar* pull_a0;
        RooRealVar* pull_a0a;
        RooRealVar* pull_at;
        RooRealVar* pull_ata;

        RooRealVar* pull_phiw;

        RooRealVar* pull_rp;
        RooRealVar* pull_r0;
        RooRealVar* pull_rt;

        RooRealVar* pull_sp;
        RooRealVar* pull_s0;
        RooRealVar* pull_st;

        RooRealVar* residuals[13];
        RooRealVar* pulls[13];
        RooRealVar* errors[13];


    protected:
    private:
};

#endif // OBSERVABLESCOLLECTION_H
