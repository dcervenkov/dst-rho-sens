#ifndef FITTERTRANS_H
#define FITTERTRANS_H

#include "DSRhoPDF.h"

class FitterTrans
{
    public:
        FitterTrans(Double_t* outer_par_input);
        ~FitterTrans();
        Int_t Fit();
        Int_t ComputeChi2(const char* type);
        Double_t GetChi2(const char* type);
        Double_t SaveChi2Maps(const char* type);
        void SaveNllPlot(const char* par);
        void SaveNllPlot(const char* par1, const char* par2);
        void SaveParameters(char* file);
        RooRealVar* GetTht(){return tht;};
        RooRealVar* GetThb(){return thb;};
        RooRealVar* GetPhit(){return phit;};
        RooRealVar* GetDt(){return dt;};
        DSRhoPDF* GetPdf(){return pdf_a;};
        DSRhoPDF* GetPdfBar(){return pdf_ab;};
        RooDataHist* GetBinnedDataSet();
        RooDataSet* GetReducedDataSet();
        RooDataSet* GetDataSet(){return dataSet;};
        void FixAllParameters();
        void FixParameter(const char* par);
        void FreeParameter(const char* par);
        void PrintParameter(const char* par);
        void GetHelParameters(Double_t* params);
        void SaveResiduals();
        void GenerateDataSet(Int_t numEvents);
        void ReadDataSet(const char* file);

    protected:
    private:
        void SaveNllPlot(RooRealVar* var);
        void SaveNllPlot(RooRealVar* var1, RooRealVar* var2);

        void CreateReducedDataSet(const char* type);
        void CreateBinnedDataSet(const char* type);
        Double_t GetVPrecise(DSRhoPDF* pdf);
        Double_t GetVPrecise1D(const int i,DSRhoPDF* pdf,RooDataSet* loc_dataset);
        Double_t GetVPrecise1D(const int i,RooSimultaneous* spdf,RooDataSet* loc_dataset);

        RooDataSet* dataSet;
        RooDataSet* dataSet_reduced;
        RooDataHist* dataSet_binned;
        Int_t binnedNumEntries;
        Double_t par_input[16];
        Bool_t doFit;
        RooChi2Var* chi2Var;
        RooFitResult* result;
        RooSimultaneous* simPdf;
        RooArgSet* parameters;
        RooArgList* variables;

        Int_t tht_bins;
        Int_t thb_bins;
        Int_t phit_bins;
        Int_t dt_bins;

        RooRealVar* tht;
        RooRealVar* thb;
        RooRealVar* phit;
        RooRealVar* dt;

        Int_t vars_bins[4];
        RooRealVar* vars[4];

        RooCategory* decType;

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

        /// Time-dep additional vars

        RooRealVar* dm;

        RooRealVar* xt;
        RooRealVar* xp;
        RooRealVar* x0;

        RooRealVar* yt;
        RooRealVar* yp;
        RooRealVar* y0;

        RooRealVar* xtb;
        RooRealVar* xpb;
        RooRealVar* x0b;

        RooRealVar* ytb;
        RooRealVar* ypb;
        RooRealVar* y0b;

        DSRhoPDF* pdf_a;
        DSRhoPDF* pdf_b;
        DSRhoPDF* pdf_ab;
        DSRhoPDF* pdf_bb;

        /// numFitParameters holds # of NON-constant fit parameters
        RooArgSet* fitParameters;
        Int_t numFitParameters;

};

#endif // FITTERTRANS_H
