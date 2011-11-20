int FitHel(RooDataSet* dataSet, RooRealVar& tha, RooRealVar& thb, RooRealVar& chi);
int FitTrans(RooDataSet* dataSet, RooRealVar& tht, RooRealVar& thb, RooRealVar& phit);
Double_t GetChi2(RooDataHist* data, RooDataHist* pdf);
Double_t GetChi2Trans(RooDataHist* data, RooDataHist* pdf);
void WriteToFile(Int_t num, Double_t* vars);
void SavePlots(RooDataSet* dataSet, RooGenericPdf* pdf, const RooRealVar& var1, const RooRealVar& var2, const RooRealVar& var3);
