#include "Constants.h"

char* inputFile = 0;
char* outputFile = 0;
const Int_t tha_bins = 20;
const Int_t thb_bins = 20;
const Int_t chi_bins = 40;

const Int_t tht_bins = 20;
// thb_bins already defined
const Int_t phit_bins = 40;

RooRealVar tha("tha","tha",0,PI);
RooRealVar thb("thb","thb",0,PI);
RooRealVar chi("chi","chi",0,2*PI);

RooRealVar tht("tht","tht",0,PI);
// thb already defined
RooRealVar phit("phit","phit",-PI,PI);


int FitHel(RooDataSet* dataSet, RooRealVar& tha, RooRealVar& thb, RooRealVar& chi);
int FitTrans(RooDataSet* dataSet, RooRealVar& tht, RooRealVar& thb, RooRealVar& phit);
Double_t GetChi2(RooDataHist* data, RooDataHist* pdf);
Double_t GetChi2Trans(RooDataHist* data, RooDataHist* pdf);
void WriteToFile(Double_t var1);//, Double_t var2, Double_t var3, Double_t var4);
