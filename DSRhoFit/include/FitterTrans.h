#ifndef FITTERTRANS_H
#define FITTERTRANS_H

#include "DSRhoPDF.h"

class FitterTrans
{
    public:
        FitterTrans(RooDataSet* outer_dataSet, Double_t* outer_par_input);
        ~FitterTrans();
        Int_t Fit();
        Int_t ComputeChi2();
        void GetRecoveredParameters(Int_t& numParameters, Double_t** recoveredParameters);
        RooRealVar* GetTht(){return tht;};
        RooRealVar* GetThb(){return thb;};
        RooRealVar* GetPhit(){return phit;};
        RooRealVar* GetDt(){return dt;};
        DSRhoPDF* GetPdf(){return pdf_a;};
        RooDataHist* GetBinnedDataSet();
        void FixAllParameters();
        void FixParameter(const char* par);
        void FreeParameter(const char* par);


    protected:
    private:
        TPluginManager* gPluginMgr;

        void CreateBinnedDataSet();

        RooDataSet* dataSet;
        RooDataHist* dataSet_binned;
        Double_t par_input[11];
        Bool_t doFit;
        RooChi2Var* chi2Var;
        RooFitResult* result;
        RooSimultaneous* simPdf;
        RooArgSet* parameters;

        Int_t tht_bins;
        Int_t thb_bins;
        Int_t phit_bins;
        Int_t dt_bins;

        RooRealVar* tht;
        RooRealVar* thb;
        RooRealVar* phit;
        RooRealVar* dt;
        RooCategory* decType;
        RooRealVar* gamma;

        RooRealVar* ap;
        RooRealVar* apa;
        RooFormulaVar* apr;
        RooFormulaVar* api;
        RooRealVar* a0;
        RooRealVar* a0a;
        RooFormulaVar* a0r;
        RooFormulaVar* a0i;
        RooFormulaVar* at;
        RooRealVar* ata;
        RooFormulaVar* atr;
        RooFormulaVar* ati;

        RooFormulaVar* ap0r;
        RooFormulaVar* a0ti;
        RooFormulaVar* apti;

        /// Time-dep additional vars

        RooRealVar* dm;
        RooRealVar* phiw;

        /// s is strong phase; delta_polarization in BN419
        RooRealVar* st;
        RooRealVar* sp;
        RooRealVar* s0;

        RooRealVar* rt;
        RooRealVar* rp;
        RooRealVar* r0;

        RooFormulaVar* ap0i;
        RooFormulaVar* a0tr;
        RooFormulaVar* aptr;

        RooFormulaVar* At2_a;
        RooFormulaVar* Ap2_a;
        RooFormulaVar* A02_a;

        RooFormulaVar* Ap0r_a;
        RooFormulaVar* A0ti_a;
        RooFormulaVar* Apti_a;

        RooFormulaVar* At2_ab;
        RooFormulaVar* Ap2_ab;
        RooFormulaVar* A02_ab;

        RooFormulaVar* Ap0r_ab;
        RooFormulaVar* A0ti_ab;
        RooFormulaVar* Apti_ab;

        RooFormulaVar* At2_b;
        RooFormulaVar* Ap2_b;
        RooFormulaVar* A02_b;

        RooFormulaVar* Ap0r_b;
        RooFormulaVar* A0ti_b;
        RooFormulaVar* Apti_b;

        RooFormulaVar* At2_bb;
        RooFormulaVar* Ap2_bb;
        RooFormulaVar* A02_bb;

        RooFormulaVar* Ap0r_bb;
        RooFormulaVar* A0ti_bb;
        RooFormulaVar* Apti_bb;

        /// a decays are favored, b corresponding suppressed; a is B0 -> D*- + rho+
        char* formula_a;
        char* formula_ab;
        char* formula_b;
        char* formula_bb;

        RooArgSet* varSet_a;
        RooArgSet* varSet_b;
        RooArgSet* varSet_ab;
        RooArgSet* varSet_bb;

        DSRhoPDF* pdf_a;
        DSRhoPDF* pdf_b;
        DSRhoPDF* pdf_ab;
        DSRhoPDF* pdf_bb;

        RooArgSet* varSet;

        /// numFitParameters holds # of NON-constant fit parameters
        RooArgSet* fitParameters;
        Int_t numFitParameters;

        char* pdfFormula;

        RooGenericPdf* pdf;

        DSRhoPDF* myPdf_a;
        DSRhoPDF* myPdf_b;
};

#endif // FITTERTRANS_H
