int FitHel(RooDataSet* dataSet, RooRealVar& tha, RooRealVar& thb, RooRealVar& chi, Double_t* par_input);
int FitTrans(RooDataSet* dataSet, RooRealVar& tht, RooRealVar& thb, RooRealVar& phit, Double_t* par_input);
Double_t SaveChi2Maps(RooDataHist* data_binned, Int_t numEvents, RooGenericPdf* pdf, RooRealVar var1, RooRealVar var2, RooRealVar var3);
Double_t GetChi2(RooDataHist* data_binned, Int_t numEvents, RooGenericPdf* pdf, RooRealVar var1, RooRealVar var2, RooRealVar var3);
void WriteToFile(Int_t num, Double_t* vars);
void SavePlots(RooDataSet* dataSet, RooGenericPdf* pdf, const RooRealVar& var1, const RooRealVar& var2, const RooRealVar& var3);
