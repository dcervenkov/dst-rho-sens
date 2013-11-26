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

        RooRealVar* xp;
        RooRealVar* xpi;
        RooRealVar* xpe;
        RooRealVar* x0;
        RooRealVar* x0i;
        RooRealVar* x0e;
        RooRealVar* xt;
        RooRealVar* xti;
        RooRealVar* xte;

        RooRealVar* yp;
        RooRealVar* ypi;
        RooRealVar* ype;
        RooRealVar* y0;
        RooRealVar* y0i;
        RooRealVar* y0e;
        RooRealVar* yt;
        RooRealVar* yti;
        RooRealVar* yte;

        RooRealVar* residual_ap;
        RooRealVar* residual_apa;
        RooRealVar* residual_a0;
        RooRealVar* residual_a0a;
        RooRealVar* residual_at;
        RooRealVar* residual_ata;

        RooRealVar* residual_phiw;

        RooRealVar* residual_xp;
        RooRealVar* residual_x0;
        RooRealVar* residual_xt;

        RooRealVar* residual_yp;
        RooRealVar* residual_y0;
        RooRealVar* residual_yt;

        RooRealVar* pull_ap;
        RooRealVar* pull_apa;
        RooRealVar* pull_a0;
        RooRealVar* pull_a0a;
        RooRealVar* pull_at;
        RooRealVar* pull_ata;

        RooRealVar* pull_phiw;

        RooRealVar* pull_xp;
        RooRealVar* pull_x0;
        RooRealVar* pull_xt;

        RooRealVar* pull_yp;
        RooRealVar* pull_y0;
        RooRealVar* pull_yt;

        RooRealVar* residuals[13];
        RooRealVar* pulls[13];
        RooRealVar* errors[13];


    protected:
    private:
};

#endif // OBSERVABLESCOLLECTION_H
