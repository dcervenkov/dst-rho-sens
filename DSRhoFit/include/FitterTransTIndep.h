#ifndef FITTERTRANSTINDEP_H
#define FITTERTRANSTINDEP_H

#include "DSRhoPDFTIndep.h"

class FitterTransTIndep
{
    public:
        FitterTransTIndep(RooDataSet* outer_dataSet, Double_t* outer_par_input);
        ~FitterTransTIndep();
        Int_t Fit();
        Int_t ComputeChi2();
        Double_t GetChi2();
        Double_t SaveChi2Maps(const char* type);
        void SaveNllPlot(const char* par);
        void SaveNllPlot(const char* par1, const char* par2);
        void GetRecoveredParameters(Int_t& numParameters, Double_t** recoveredParameters);
        RooRealVar* GetTht(){return tht;};
        RooRealVar* GetThb(){return thb;};
        RooRealVar* GetPhit(){return phit;};
        DSRhoPDFTIndep* GetPdf(){return pdf;};
        RooDataHist* GetBinnedDataSet();
        RooDataSet* GetReducedDataSet();
        RooDataSet* GetDataSet(){return dataSet;};
        void FixAllParameters();
        void FixParameter(const char* par);
        void FreeParameter(const char* par);
        void PrintParameter(const char* par);
        void GetHelParameters(Double_t* params);
        void SaveResiduals();


    protected:
    private:
        void SaveNllPlot(RooRealVar* var);
        void SaveNllPlot(RooRealVar* var1, RooRealVar* var2);

        void CreateBinnedDataSet();
        Double_t GetVPrecise(DSRhoPDFTIndep* pdf);
        Double_t GetVPrecise1D(const int i,DSRhoPDFTIndep* pdf,RooDataSet* loc_dataset);

        RooDataSet* dataSet;
        RooDataSet* dataSet_reduced;
        RooDataHist* dataSet_binned;
        Int_t binnedNumEntries;
        Double_t par_input[11];
        Bool_t doFit;
        RooChi2Var* chi2Var;
        RooFitResult* result;
        RooArgSet* parameters;

        Int_t tht_bins;
        Int_t thb_bins;
        Int_t phit_bins;
        Int_t vars_bins[3];

        RooRealVar* tht;
        RooRealVar* thb;
        RooRealVar* phit;
        RooRealVar* vars[3];

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

        RooArgSet* varSet;

        /// numFitParameters holds # of NON-constant fit parameters
        RooArgSet* fitParameters;
        Int_t numFitParameters;

        DSRhoPDFTIndep* pdf;

};

#endif // FITTERTRANSTINDEP_H
