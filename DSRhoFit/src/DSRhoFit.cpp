#include "DSRhoFit.h"

// Standard includes
#include <getopt.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>

// ROOT includes
#include "RooAddPdf.h"
#include "RooArgSet.h"
#include "RooCategory.h"
#include "RooChi2Var.h"
#include "RooDataHist.h"
#include "RooDataSet.h"
#include "RooExponential.h"
#include "RooFit.h"
#include "RooFitResult.h"
#include "RooGaussian.h"
#include "RooGenericPdf.h"
#include "RooGlobalFunc.h"
#include "RooHist.h"
#include "RooMinuit.h"
#include "RooNLLVar.h"
#include "RooPlot.h"
#include "RooPolynomial.h"
#include "RooProdPdf.h"
#include "RooRandom.h"
#include "RooRealVar.h"
#include "TApplication.h"
#include "TCanvas.h"
#include "TComplex.h"
#include "TFile.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TLine.h"
#include "TMath.h"
#include "TNtuple.h"
#include "TPad.h"
#include "TROOT.h"
#include "TRandom1.h"
#include "TStopwatch.h"
#include "TStyle.h"
#include "TTree.h"

// Local includes
#include "ASSERT.h"
#include "Constants.h"
#include "FitterTrans.h"

#define DEBUG
//#define VERBOSE

char* inputFile = 0;
char* outputFile = 0;

int main(int argc, char* argv[]) {
    char** optionless_argv = NULL;
    // The {} causes the struct's members to be initialized to 0. Without it
    // they would have unspecified values
    fitter_options options = {};
    const int optionless_argc = ProcessCmdLineOptions(argc, argv, optionless_argv, options);

    int num_mandatory_arguments = options.time_independent_set ? 7 : 14;
    if (optionless_argc != num_mandatory_arguments) {
        printf("ERROR: Wrong number of arguments.\n");
        if (options.time_independent_set) {
            printf("Usage: DSRhoFit inputFile outputFile ap apa a0 ata\n");
        } else {
            printf("Usage: DSRhoFit inputFile outputFile ap apa a0 ata phiw rp r0 rt sp s0 st\n");
        }
        return 85;
    }

    TStopwatch timer;
    timer.Start();

    /// This is so I have to change only the next block if I change the
    /// ordering, etc. of arguments
    inputFile = optionless_argv[1];
    outputFile = optionless_argv[2];
    Int_t numPars = optionless_argc - 3;
    Double_t par_input[numPars];
    for (Int_t i = 0; i < numPars; i++) {
        par_input[i] = atof(optionless_argv[i + 3]);
    }

    Int_t doFit = false;
    Int_t doPlot = false;

    FitterTrans* fitter = new FitterTrans(
        par_input, options.time_independent_set ? (!options.time_independent) : true);

    if (options.num_CPUs_set) fitter->SetNumCPUs(options.num_CPUs);
    if (options.make_plots_set) doPlot = options.make_plots;
    if (options.fit_set) doFit = options.fit;

    if (doFit == 3) {
        ConvertBetweenHelAndTrans(par_input);
    } else if (doFit == 4) {
        fitter->GenerateDataSet(100000);
        fitter->GetDataSet()->write(inputFile);
    } else {
        // ConvertBetweenHelAndTrans(par_input);
        fitter->ReadDataSet(inputFile);
        ProcessTrans(fitter, doFit, doPlot);
    }

    timer.Stop();
    timer.Print();

#ifdef GRAPHIC
    printf("\nProgram execution has finished.\n");
    rootapp->Run();
#endif
    return 0;
}

int ProcessTrans(FitterTrans* fitter, Int_t doFit, Int_t doPlot) {
    if (doFit) {
        fitter->FixAllParameters();
        fitter->FreeParameter("ap");
        fitter->FreeParameter("apa");
        fitter->FreeParameter("a0");
        fitter->FreeParameter("ata");
        // fitter->FreeParameter("phiw");
        // fitter->FreeParameter("rp");
        // fitter->FreeParameter("r0");
        // fitter->FreeParameter("rt");
        // fitter->FreeParameter("sp");
        // fitter->FreeParameter("s0");
        // fitter->FreeParameter("st");
        fitter->Fit();
        fitter->SaveParameters(outputFile);
    }

    if (doPlot == kTRUE) {
        // SaveChi2Maps(fitter->GetBinnedDataSet(),dataSet->numEntries(),fitter->GetPdf(),*(fitter->GetTht()),*(fitter->GetThb()),*(fitter->GetPhit()));
        // Double_t mychi2 = fitter->SaveChi2Maps("a");
        // fitter->SaveResiduals();
        // fitter->SaveNllPlot("r0");
        // printf("mychi2 from SaveChi2Maps = %f\n",mychi2);
        // SavePlots(fitter->GetDataSet(),fitter->GetPdf(),*(fitter->GetTht()),*(fitter->GetThb()),*(fitter->GetPhit()),*(fitter->GetDt()));
        fitter->SaveVarPlot(fitter->GetTht());
        fitter->SaveVarPlot(fitter->GetThb());
        fitter->SaveVarPlot(fitter->GetPhit());
        if (fitter->IsTimeDependent()) {
            fitter->SaveDtPlots();
        }
    }

    return 0;
}

Double_t SaveChi2Maps(RooDataHist* data_binned, Int_t numEvents, RooGenericPdf* pdf,
                      RooRealVar var1, RooRealVar var2, RooRealVar var3) {
    /// Create a histogram from the pdf with the expected number of events with no statistical
    /// fluctuation
    RooDataHist* pdf_binned = pdf->generateBinned(RooArgSet(var1, var2, var3), numEvents, kTRUE);

    const Int_t var1_bins = 30;
    const Int_t var2_bins = 30;
    const Int_t var3_bins = 30;

    Double_t mychi2 = 0;
    Double_t dchi2 = 0;
    Double_t n = 0;
    Double_t v = 0;

    TH1F* h1_chi2 = new TH1F("h1_chi2", "h1_chi2", 100, 0, 100);
    TH2F* h2_chi2_1 = new TH2F("h2_chi2_1", "h2_chi2_1", var1_bins, var1.getMin(), var1.getMax(),
                               var2_bins, var2.getMin(), var2.getMax());
    TH2F* h2_chi2_2 = new TH2F("h2_chi2_2", "h2_chi2_2", var1_bins, var1.getMin(), var1.getMax(),
                               var3_bins, var3.getMin(), var3.getMax());
    TH2F* h2_chi2_3 = new TH2F("h2_chi2_3", "h2_chi2_3", var2_bins, var2.getMin(), var2.getMax(),
                               var3_bins, var3.getMin(), var3.getMax());

    /// Cycle through the centers of all bins
    /// I'm getting width of the first bin, because all bins are of equal width
    for (var1 = var1.getMin() + var1.getBinWidth(0) / 2; var1.getVal() < var1.getMax();
         var1.setVal(var1.getVal() + var1.getBinWidth(0))) {
        for (var2 = var2.getMin() + var2.getBinWidth(0) / 2; var2.getVal() < var2.getMax();
             var2.setVal(var2.getVal() + var2.getBinWidth(0))) {
            for (var3 = var3.getMin() + var3.getBinWidth(0) / 2; var3.getVal() < var3.getMax();
                 var3.setVal(var3.getVal() + var3.getBinWidth(0))) {
                /// Weight is actually the bin content
                n = data_binned->weight(RooArgSet(var1, var2, var3), 0);
                v = pdf_binned->weight(RooArgSet(var1, var2, var3), 0);
                dchi2 = (n - v) * (n - v) / v;
                if (dchi2 > h1_chi2->GetXaxis()->GetXmax() - 1)
                    dchi2 = h1_chi2->GetXaxis()->GetXmax() - 1;
                h1_chi2->Fill(dchi2);
                h2_chi2_1->Fill(var1.getVal(), var2.getVal(), dchi2);
                h2_chi2_2->Fill(var1.getVal(), var3.getVal(), dchi2);
                h2_chi2_3->Fill(var2.getVal(), var3.getVal(), dchi2);
                mychi2 += dchi2;
            }
        }
    }

    TFile* file = new TFile("plots/chi2maps.root", "RECREATE");
    TCanvas* c2 = new TCanvas("c2", "c2", 600, 600);
    TString path;

    c2->SetLogy(kTRUE);
    h1_chi2->GetXaxis()->SetTitle("dchi");
    h1_chi2->GetYaxis()->SetTitle("num bins");
    h1_chi2->Draw();
    h1_chi2->Write();
    c2->SaveAs("plots/chi2_delta.png");

    c2->SetLogy(kFALSE);

    h2_chi2_1->SetOption("colz");
    // hdchi2->SetMinimum(0);
    // hdchi2->SetMaximum(100);
    h2_chi2_1->SetStats(kFALSE);
    h2_chi2_1->GetXaxis()->SetTitle(var1.GetName());
    h2_chi2_1->GetYaxis()->SetTitle(var2.GetName());
    h2_chi2_1->Draw();
    h2_chi2_1->Write();
    path = "plots/chi2map_";
    path += var1.GetName();
    path += "_";
    path += var2.GetName();
    path += ".png";
    c2->SaveAs(path);

    h2_chi2_2->SetOption("colz");
    // hdchi2->SetMinimum(0);
    // hdchi2->SetMaximum(100);
    h2_chi2_2->SetStats(kFALSE);
    h2_chi2_2->GetXaxis()->SetTitle(var1.GetName());
    h2_chi2_2->GetYaxis()->SetTitle(var3.GetName());
    h2_chi2_2->Draw();
    h2_chi2_2->Write();
    path = "plots/chi2map_";
    path += var1.GetName();
    path += "_";
    path += var3.GetName();
    path += ".png";
    c2->SaveAs(path);

    h2_chi2_3->SetOption("colz");
    // hdchi2->SetMinimum(0);
    // hdchi2->SetMaximum(100);
    h2_chi2_3->SetStats(kFALSE);
    h2_chi2_3->GetXaxis()->SetTitle(var2.GetName());
    h2_chi2_3->GetYaxis()->SetTitle(var3.GetName());
    h2_chi2_3->Draw();
    h2_chi2_3->Write();
    path = "plots/chi2map_";
    path += var2.GetName();
    path += "_";
    path += var3.GetName();
    path += ".png";
    c2->SaveAs(path);

    file->Close();

    delete h1_chi2;
    delete h2_chi2_1;
    delete h2_chi2_2;
    delete h2_chi2_3;
    delete c2;

    return mychi2;
}

void ConvertBetweenHelAndTrans(Double_t* par_input) {
    /// The variables are named as if converting parameters from helicity to transversity
    /// but the transformation is symmetric and can be used to convert from transversity
    /// to helicity as well.

    TComplex hp(par_input[0], par_input[1], true);
    TComplex h0(par_input[2], 0, true);
    TComplex hm(sqrt(1 - hp.Rho2() - h0.Rho2()), par_input[3], true);

    TComplex hRhop(par_input[5], par_input[8], true);
    TComplex hRho0(par_input[6], par_input[9], true);
    TComplex hRhom(par_input[7], par_input[10], true);

    TComplex hps = hRhop * hp;
    TComplex h0s = hRho0 * h0;
    TComplex hms = hRhom * hm;

    TComplex ap = (hp + hm) / sqrt(2);
    TComplex a0 = h0;
    TComplex at = (hp - hm) / sqrt(2);

    TComplex aps = (hps + hms) / sqrt(2);
    TComplex a0s = h0s;
    TComplex ats = (hps - hms) / sqrt(2);

    TComplex tRhop = aps / ap;
    TComplex tRho0 = a0s / a0;
    TComplex tRhot = ats / at;

    printf("original hel:\t");
    printf(
        "hp = %.3f\thpa = %.2f\th0 = %.3f\thma = %.2f\tphiw = %.4f\trp = %.3f\tr0 = %.3f\trm = "
        "%.3f\tsp = %.3f\ts0 = %.3f\tsm = %.3f\n",
        par_input[0], par_input[1], par_input[2], par_input[3], par_input[4], par_input[5],
        par_input[6], par_input[7], par_input[8], par_input[9], par_input[10]);

    par_input[0] = ap.Rho();
    par_input[1] = ap.Theta();
    par_input[2] = a0.Rho();
    par_input[3] = at.Theta();
    par_input[5] = tRhop.Rho();
    par_input[6] = tRho0.Rho();
    par_input[7] = tRhot.Rho();
    par_input[8] = tRhop.Theta();
    par_input[9] = tRho0.Theta();
    par_input[10] = tRhot.Theta();

    printf("converted trans:");
    printf(
        "ap = %.3f\tapa = %.2f\ta0 = %.3f\tata = %.2f\tphiw = %.4f\trp = %.3f\tr0 = %.3f\trt = "
        "%.3f\tsp = %.3f\ts0 = %.3f\tst = %.3f\n",
        par_input[0], par_input[1], par_input[2], par_input[3], par_input[4], par_input[5],
        par_input[6], par_input[7], par_input[8], par_input[9], par_input[10]);
}

Double_t Round(Double_t number, Int_t digits) {
    number = number * pow(10, digits);
    if (fmod(number, 1) > 0.5)
        number = ceil(number);
    else
        number = floor(number);

    return number / pow(10, digits);
}

/*
 * Parses command line input and extracts switches and options from it, e.g., -h or --help.
 * Then it acts accordingly, e.g., displaying help or setting variables in an option struct.
 * It also returns optionless_argv and optionless_argc (return value) for easy integration with
 * existing code.
 *
 * @param argc Standard argc
 * @param argv Standard argv
 * @param optionless_argv Pointer where to write the new argv with processed switches removed
 * @param options Struct which holds the variables acted upon by switches
 */
int ProcessCmdLineOptions(const int argc, char* const argv[], char**& optionless_argv,
                          fitter_options& options) {
    int c;
    struct option long_options[] = {{"cpus", required_argument, 0, 'c'},
                                    {"events", required_argument, 0, 'e'},
                                    {"time-independent", no_argument, 0, 'i'},
                                    {"fit", no_argument, 0, 'f'},
                                    {"plot", no_argument, 0, 'p'},
                                    {"help", no_argument, 0, 'h'},
                                    {NULL, no_argument, NULL, 0}};
    int option_index = 0;
    while ((c = getopt_long(argc, argv, "c:e:ifph", long_options, &option_index)) != -1) {
        switch (c) {
            case 0:
                printf("option %s", long_options[option_index].name);
                if (optarg) printf(" with arg %s", optarg);
                printf("\n");
                break;
            case 'c':
                options.num_CPUs = atoi(optarg);
                options.num_CPUs_set = true;
                break;
            case 'e':
                options.num_events = atoi(optarg);
                options.num_events_set = true;
                break;
            case 'i':
                options.time_independent = true;
                options.time_independent_set = true;
                break;
            case 'f':
                options.fit = true;
                options.fit_set = true;
                break;
            case 'p':
                options.make_plots = true;
                options.make_plots_set = true;
                break;
            case 'h':
                printf("Usage: %s [OPTION]... INPUT-FILE OUTPUT_DIR\n\n", argv[0]);
                printf("Mandatory arguments to long options are mandatory for short options too.\n");
                printf("-c, --cpus=NUM_CPUS     number of CPU cores to use for fitting and plotting\n");
                printf("-e, --events=NUM_EVENTS number of events to be imported from the input file\n");
                printf("-h, --help              display this text and exit\n");
                printf("-p, --plot              create angular/dt plots\n");
                printf("-f, --fit               fit\n");
                printf("-i, --time-independent  use time-independent PDF\n");
                exit(0);
                break;
            default:
                printf("?? getopt returned character code 0%o ??\n", c);
        }
    }

    // Create a char** that will become the new argv, with the options removed
    const int optionless_argc = argc - optind + 1;
    optionless_argv = new char*[optionless_argc];
    // We want to keep the program name argument
    optionless_argv[0] = argv[0];
    for (int i = 1; i < optionless_argc; i++) {
        optionless_argv[i] = argv[i - 1 + optind];
    }

    return optionless_argc;
}