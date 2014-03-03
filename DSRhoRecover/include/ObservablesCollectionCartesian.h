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
        void AdjustResultsForPeriodicity();
        void CalculateInitialPolarVals();
        void InitializePolarVals();
        void SaveParameters(FILE* pFile);

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

        RooRealVar* xbp;
        RooRealVar* xbpi;
        RooRealVar* xbpe;
        RooRealVar* xb0;
        RooRealVar* xb0i;
        RooRealVar* xb0e;
        RooRealVar* xbt;
        RooRealVar* xbti;
        RooRealVar* xbte;

        RooRealVar* ybp;
        RooRealVar* ybpi;
        RooRealVar* ybpe;
        RooRealVar* yb0;
        RooRealVar* yb0i;
        RooRealVar* yb0e;
        RooRealVar* ybt;
        RooRealVar* ybti;
        RooRealVar* ybte;

        RooRealVar* residual_ap;
        RooRealVar* residual_apa;
        RooRealVar* residual_a0;
        RooRealVar* residual_a0a;
        RooRealVar* residual_at;
        RooRealVar* residual_ata;

        RooRealVar* residual_xp;
        RooRealVar* residual_x0;
        RooRealVar* residual_xt;

        RooRealVar* residual_yp;
        RooRealVar* residual_y0;
        RooRealVar* residual_yt;

        RooRealVar* residual_xbp;
        RooRealVar* residual_xb0;
        RooRealVar* residual_xbt;

        RooRealVar* residual_ybp;
        RooRealVar* residual_yb0;
        RooRealVar* residual_ybt;

        RooRealVar* pull_ap;
        RooRealVar* pull_apa;
        RooRealVar* pull_a0;
        RooRealVar* pull_a0a;
        RooRealVar* pull_at;
        RooRealVar* pull_ata;

        RooRealVar* pull_xp;
        RooRealVar* pull_x0;
        RooRealVar* pull_xt;

        RooRealVar* pull_yp;
        RooRealVar* pull_y0;
        RooRealVar* pull_yt;

        RooRealVar* pull_xbp;
        RooRealVar* pull_xb0;
        RooRealVar* pull_xbt;

        RooRealVar* pull_ybp;
        RooRealVar* pull_yb0;
        RooRealVar* pull_ybt;

        RooRealVar* residuals[18];
        RooRealVar* pulls[18];
        RooRealVar* errors[18];

        RooRealVar* phiwi;
        RooRealVar* phiw;
        RooRealVar* rpi;
        RooRealVar* rp;
        RooRealVar* r0i;
        RooRealVar* r0;
        RooRealVar* rti;
        RooRealVar* rt;
        RooRealVar* spi;
        RooRealVar* sp;
        RooRealVar* s0i;
        RooRealVar* s0;
        RooRealVar* sti;
        RooRealVar* st;


    protected:
    private:
};

#endif // OBSERVABLESCOLLECTION_H
