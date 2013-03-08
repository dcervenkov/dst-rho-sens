#include "DSRhoPDF.h"
#include "DSRhoPDFTIndep.h"

int ProcessTrans(RooDataSet* dataSet, Double_t* par_input, Int_t doFit, Int_t doPlot);
int ToyProcessTrans(RooDataSet* dataSet, Double_t* par_input, Int_t doFit, Int_t doPlot);
Double_t SaveChi2Maps(RooDataHist* data_binned, Int_t numEvents, RooGenericPdf* pdf, RooRealVar var1, RooRealVar var2, RooRealVar var3);
void SavePlots(RooDataSet* dataSet, DSRhoPDF* pdf, const RooRealVar& var1, const RooRealVar& var2, const RooRealVar& var3, RooRealVar& dt);
void SavePlotsTIndep(RooDataSet* dataSet, DSRhoPDFTIndep* pdf, const RooRealVar& var1, const RooRealVar& var2, const RooRealVar& var3);
void ConvertBetweenHelAndTrans(Double_t* par_input);
Double_t Round(Double_t number, Int_t digits);
