#ifndef FITTERTRANS_H
#define FITTERTRANS_H

// ROOT includes
#include "RooChi2Var.h"
#include "RooSimultaneous.h"
#include "TCanvas.h"
#include "TFile.h"

// Local includes
#include "DSRhoPDF.h"
#include "DSRhoPDFTIndep.h"

class FitterTrans {
   public:
    FitterTrans(Double_t* outer_par_input, bool time_dependent = true);
    ~FitterTrans();
    Int_t Fit();
    Int_t ComputeChi2(const char* type);
    Double_t GetChi2(const char* type);
    Double_t SaveChi2Maps(const char* type);
    void SaveNllPlot(const char* par);
    void SaveNllPlot(const char* par1, const char* par2);
    void SaveParameters(char* file);
    RooRealVar* GetTht() { return tht; };
    RooRealVar* GetThb() { return thb; };
    RooRealVar* GetPhit() { return phit; };
    RooRealVar* GetDt() { return dt; };
    DSRhoPDF* GetPdf() { return pdf_a; };
    RooDataHist* GetBinnedDataSet();
    RooDataSet* GetReducedDataSet();
    RooDataSet* GetDataSet() { return dataSet; };
    void FixAllParameters();
    void FixParameter(const char* par);
    void FreeParameter(const char* par);
    void PrintParameter(const char* par);
    void GetHelParameters(Double_t* params);
    void SaveResiduals();
    void GenerateDataSet(Int_t numEvents);
    void ReadDataSet(const char* file);
    void SetNumCPUs(unsigned int CPUs) { num_CPUs = CPUs; };
    void SaveVarPlot(RooRealVar* var);
    void SaveDtPlots();
    bool IsTimeDependent() { return time_dependent; };

   protected:
   private:
    unsigned int num_CPUs = 1;
    TFile* plot_output_file = NULL;

    void SaveNllPlot(RooRealVar* var);
    void SaveNllPlot(RooRealVar* var1, RooRealVar* var2);

    void CreateReducedDataSet(const char* type);
    void CreateBinnedDataSet(const char* type);
    Double_t GetVPrecise(DSRhoPDF* pdf);
    Double_t GetVPrecise1D(const int i, DSRhoPDF* pdf, RooDataSet* loc_dataset);
    Double_t GetVPrecise1D(const int i, RooSimultaneous* spdf, RooDataSet* loc_dataset);

    void DrawResidualFrame(RooPlot* frame, RooRealVar var, TCanvas* canvas, Int_t padNumber);

    RooDataSet* dataSet;
    RooDataSet* dataSet_reduced;
    RooDataHist* dataSet_binned;
    Int_t binnedNumEntries;
    Double_t par_input[11];
    Bool_t doFit;
    RooChi2Var* chi2Var;
    RooFitResult* result;
    RooArgSet* parameters;
    RooArgList* variables;

    bool time_dependent;

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
    RooRealVar* phiw;

    /// s is strong phase; delta_polarization in BN419
    RooRealVar* st;
    RooRealVar* sp;
    RooRealVar* s0;

    RooRealVar* rt;
    RooRealVar* rp;
    RooRealVar* r0;

    RooAbsPdf* pdf;

    DSRhoPDFTIndep* pdf_tindep;

    DSRhoPDF* pdf_a;
    DSRhoPDF* pdf_b;
    DSRhoPDF* pdf_ab;
    DSRhoPDF* pdf_bb;

    RooSimultaneous* pdf_sim;

    /// numFitParameters holds # of NON-constant fit parameters
    RooArgSet* fitParameters;
    Int_t numFitParameters;
};

#endif  // FITTERTRANS_H
