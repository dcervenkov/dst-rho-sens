int ProcessHel(RooDataSet* dataSet, RooRealVar& tha, RooRealVar& thb, RooRealVar& chi, RooRealVar& dt, Double_t* par_input, Bool_t doFit, Bool_t doPlot);
int ProcessTrans(RooDataSet* dataSet, Double_t* par_input, Bool_t doFit, Bool_t doPlot);
Double_t SaveChi2Maps(RooDataHist* data_binned, Int_t numEvents, RooGenericPdf* pdf, RooRealVar var1, RooRealVar var2, RooRealVar var3);
Double_t GetChi2(RooDataHist* data_binned, Int_t numEvents, RooGenericPdf* pdf, RooRealVar var1, RooRealVar var2, RooRealVar var3);
void WriteToFile(Int_t num, Double_t* vars, char* file);
void SavePlots(RooDataSet* dataSet, RooGenericPdf* pdf, const RooRealVar& var1, const RooRealVar& var2, const RooRealVar& var3, RooRealVar& dt);
void ConvertTransToHel(Double_t* par_input);
void ConvertHelToTrans(Double_t* par_input);
Double_t Round(Double_t number, Int_t digits);
