#ifndef DSRHOFIT_H_
#define DSRHOFIT_H_

// ROOT includes
#include "TCanvas.h"
#include "RooGenericPdf.h"

// Local includes
#include "DSRhoPDF.h"
#include "DSRhoPDFTIndep.h"
#include "FitterTrans.h"

struct fitter_options {
	int num_CPUs;
	bool num_CPUs_set;
    bool time_independent;
    bool time_independent_set;
	bool make_plots;
	bool make_plots_set;
	bool fit;
	bool fit_set;
	int num_events;
	bool num_events_set;
};

int ProcessCmdLineOptions(const int argc, char* const argv[], char**& optionless_argv, fitter_options& options);
int ProcessTrans(FitterTrans* fitter, Int_t doFit, Int_t doPlot);
Double_t SaveChi2Maps(RooDataHist* data_binned, Int_t numEvents, RooGenericPdf* pdf, RooRealVar var1, RooRealVar var2, RooRealVar var3);
void SavePlots(RooDataSet* dataSet, RooAbsPdf* pdf, const RooRealVar& var1, const RooRealVar& var2, const RooRealVar& var3, RooRealVar& dt);
void ConvertBetweenHelAndTrans(Double_t* par_input);
Double_t Round(Double_t number, Int_t digits);

#endif /* DSRHOFIT_H_ */