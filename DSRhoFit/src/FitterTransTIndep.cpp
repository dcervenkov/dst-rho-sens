#include "RooRealVar.h"
#include "RooFormulaVar.h"
#include "RooCategory.h"
#include "RooGenericPdf.h"
#include "Constants.h"
#include "RooDataSet.h"
#include "RooRandom.h"
#include "RooChi2Var.h"
#include "RooSimultaneous.h"
#include "Minuit2/Minuit2Minimizer.h"
#include "TPluginManager.h"
#include "TMath.h"
#include "TIterator.h"
#include "TLine.h"

#include "DSRhoPDFTIndep.h"
#include "FitterTransTIndep.h"

#include "TH1D.h"
#include "TH2F.h"
#include "TCanvas.h"
#include "TFile.h"

FitterTransTIndep::FitterTransTIndep(RooDataSet* outer_dataSet, Double_t* outer_par_input)
{
    gPluginMgr = new TPluginManager;
    gPluginMgr->AddHandler("ROOT::Math::Minimizer", "Minuit2", "Minuit2Minimizer", "Minuit2", "Minuit2Minimizer(const char *)");

    dataSet = outer_dataSet;
    for(int i = 0; i < 11; i++)
        par_input[i] = outer_par_input[i];

    chi2Var = 0;
    result = 0;

    tht_bins = 60;
    thb_bins = 60;
    phit_bins = 60;

    vars_bins[0] = tht_bins;
    vars_bins[1] = thb_bins;
    vars_bins[2] = phit_bins;

    thb = new RooRealVar("thb","thb",0,PI);
    tht = new RooRealVar("tht","tht",0,PI);
    phit = new RooRealVar("phit","phit",-PI,PI);

    vars[0] = tht;
    vars[1] = thb;
    vars[2] = phit;

    ap = new RooRealVar("ap","ap",par_input[0],0.1,0.4);
    apa = new RooRealVar("apa","apa",par_input[1],0.3,0.7);
    apr = new RooFormulaVar("apr","ap*cos(apa)",RooArgSet(*ap,*apa));
    api = new RooFormulaVar("api","ap*sin(apa)",RooArgSet(*ap,*apa));
    a0 = new RooRealVar("a0","a0",par_input[2],0.8,1);
    a0a = new RooRealVar("a0a","a0a",0);
    a0r = new RooFormulaVar("a0r","a0*cos(a0a)",RooArgSet(*a0,*a0a));
    a0i = new RooFormulaVar("a0i","a0*sin(a0a)",RooArgSet(*a0,*a0a));
    at = new RooFormulaVar("at","sqrt(1-ap*ap-a0*a0)",RooArgSet(*ap,*a0));
    ata = new RooRealVar("ata","ata",par_input[3],2.7,3);
    atr = new RooFormulaVar("atr","at*cos(ata)",RooArgSet(*at,*ata));
    ati = new RooFormulaVar("ati","at*sin(ata)",RooArgSet(*at,*ata));

    ap0r = new RooFormulaVar("ap0r","ap*a0*cos(-apa+a0a)",RooArgSet(*ap,*apa,*a0,*a0a));
    a0ti = new RooFormulaVar("a0ti","a0*at*sin(-a0a+ata)",RooArgSet(*a0,*a0a,*at,*ata));
    apti = new RooFormulaVar("apti","ap*at*sin(-apa+ata)",RooArgSet(*ap,*apa,*at,*ata));


    varSet_a = new RooArgSet(*Ap2_a,*At2_a,*A02_a,*Ap0r_a,*A0ti_a,*Apti_a,*tht,*thb,*phit);
    varSet_a->add(*dt);
    varSet_a->add(*gamma);
    varSet_b = new RooArgSet(*Ap2_b,*At2_b,*A02_b,*Ap0r_b,*A0ti_b,*Apti_b,*tht,*thb,*phit);
    varSet_b->add(*dt);
    varSet_b->add(*gamma);
    varSet_ab = new RooArgSet(*Ap2_ab,*At2_ab,*A02_ab,*Ap0r_ab,*A0ti_ab,*Apti_ab,*tht,*thb,*phit);
    varSet_ab->add(*dt);
    varSet_ab->add(*gamma);
    varSet_bb = new RooArgSet(*Ap2_bb,*At2_bb,*A02_bb,*Ap0r_bb,*A0ti_bb,*Apti_bb,*tht,*thb,*phit);
    varSet_bb->add(*dt);
    varSet_bb->add(*gamma);

    pdf_a = new DSRhoPDF("pdf_a","pdf_a","a",*tht,*thb,*phit,*dt,*ap,*apa,*a0,*ata,*phiw,*rp,*r0,*rt,*sp,*s0,*st);
    pdf_b = new DSRhoPDF("pdf_b","pdf_b","b",*tht,*thb,*phit,*dt,*ap,*apa,*a0,*ata,*phiw,*rp,*r0,*rt,*sp,*s0,*st);
    pdf_ab = new DSRhoPDF("pdf_ab","pdf_ab","ab",*tht,*thb,*phit,*dt,*ap,*apa,*a0,*ata,*phiw,*rp,*r0,*rt,*sp,*s0,*st);
    pdf_bb = new DSRhoPDF("pdf_bb","pdf_bb","bb",*tht,*thb,*phit,*dt,*ap,*apa,*a0,*ata,*phiw,*rp,*r0,*rt,*sp,*s0,*st);

    simPdf = new RooSimultaneous("simPdf","simPdf",*decType);
    simPdf->addPdf(*pdf_a,"a");
    simPdf->addPdf(*pdf_ab,"ab");
    simPdf->addPdf(*pdf_b,"b");
    simPdf->addPdf(*pdf_bb,"bb");

    parameters = new RooArgSet(*ap,*apa,*a0,*ata,*phiw,*rp,*r0,*rt);
    parameters->add(*sp);
    parameters->add(*s0);
    parameters->add(*st);

    /// numFitParameters holds # of NON-constant fit parameters
    numFitParameters = (parameters->selectByAttrib("Constant",kFALSE))->getSize();

//    delete dataSet;
//
//    dataSet = pdf_a->generate(RooArgSet(*vars[0],*vars[1],*vars[2],*vars[3],*decType),80000);
//    dataSet->write("data/dataset_from_pdf_a");
//    delete dataSet;
//
//    dataSet = pdf_ab->generate(RooArgSet(*vars[0],*vars[1],*vars[2],*vars[3]),80000);
//    dataSet->write("data/dataset_from_pdf_ab");
//    delete dataSet;
//
//    dataSet = pdf_b->generate(RooArgSet(*vars[0],*vars[1],*vars[2],*vars[3]),20000);
//    dataSet->write("data/dataset_from_pdf_b");
//    delete dataSet;
//
//    dataSet = pdf_bb->generate(RooArgSet(*vars[0],*vars[1],*vars[2],*vars[3]),20000);
//    dataSet->write("data/dataset_from_pdf_bb");
//    delete dataSet;


////    dataSet = simPdf->generate(RooArgSet(*vars[0],*vars[1],*vars[2],*vars[3],*decType),10000,RooFit::Extended());

}

FitterTransTIndep::~FitterTransTIndep()
{
    delete gPluginMgr;
    delete thb;
    delete tht;
    delete phit;
    delete dt;
    delete decType;
    delete gamma;
    delete ap;
    delete apa;
    delete apr;
    delete api;
    delete a0;
    delete a0a;
    delete a0r;
    delete a0i;
    delete at;
    delete ata;
    delete atr;
    delete ati;
    delete ap0r;
    delete a0ti;
    delete apti;
    delete dm;
    delete phiw;
    delete rp;
    delete r0;
    delete rt;
    delete sp;
    delete s0;
    delete st;
    delete ap0i;
    delete a0tr;
    delete aptr;
    delete At2_a;
    delete Ap2_a;
    delete A02_a;
    delete Ap0r_a;
    delete A0ti_a;
    delete Apti_a;
    delete At2_ab;
    delete Ap2_ab;
    delete A02_ab;
    delete Ap0r_ab;
    delete A0ti_ab;
    delete Apti_ab;
    delete At2_b;
    delete Ap2_b;
    delete A02_b;
    delete Ap0r_b;
    delete A0ti_b;
    delete Apti_b;
    delete At2_bb;
    delete Ap2_bb;
    delete A02_bb;
    delete Ap0r_bb;
    delete A0ti_bb;
    delete Apti_bb;
    delete varSet_a;
    delete varSet_b;
    delete varSet_ab;
    delete varSet_bb;
    delete pdf_a;
    delete pdf_b;
    delete pdf_ab;
    delete pdf_bb;
    delete simPdf;
    delete parameters;
    delete fitParameters;
}

Int_t FitterTransTIndep::Fit()
{
    numFitParameters = (parameters->selectByAttrib("Constant",kFALSE))->getSize();
    //TPluginManager* gPluginMgr = new TPluginManager;
    //gPluginMgr->AddHandler("ROOT::Math::Minimizer", "Minuit2", "Minuit2Minimizer", "Minuit2", "Minuit2Minimizer(const char *)");
    //result = pdf->fitTo(*dataSet,RooFit::Save(),RooFit::Timer(true),RooFit::Minimizer("Minuit2"));//,RooFit::NumCPU(2));
    result = simPdf->fitTo(*dataSet,RooFit::Save(),RooFit::Timer(true),RooFit::Minimizer("Minuit"),RooFit::Minos(),RooFit::Hesse(),RooFit::Strategy(1));//,RooFit::NumCPU(2));
    //result->Print();
}

void FitterTransTIndep::CreateBinnedDataSet(const char* type)
{
    tht->setBins(tht_bins);
    thb->setBins(thb_bins);
    phit->setBins(phit_bins);
    dt->setBins(dt_bins);

    TString cut = "decType==decType::";
    cut += type;

    RooRandom::randomGenerator()->SetSeed(0);
    dataSet_reduced = (RooDataSet*)dataSet->reduce(cut);
    binnedNumEntries = dataSet_reduced->numEntries();

    /// Create a binned dataSet which is needed for chi2 calculation
    dataSet_binned = new RooDataHist("dataSet_binned","dataSet_binned",RooArgSet(*tht,*thb,*phit,*dt),*dataSet_reduced);

    //dataSet_binned = new RooDataHist("dataSet_binned","dataSet_binned",RooArgSet(*tht,*thb,*phit,*dt),*pdf_a);
    //RooDataHist* dataSet_binned = pdf->generateBinned(RooArgSet(var1,var2,var3),dataSet->numEntries(),kFALSE);
}

void FitterTransTIndep::CreateReducedDataSet(const char* type)
{
    TString cut = "decType==decType::";
    cut += type;
    RooRandom::randomGenerator()->SetSeed(0);
    dataSet_reduced = (RooDataSet*)dataSet->reduce(cut);
}

Int_t FitterTransTIndep::ComputeChi2(const char* type)
{
    //if(dataSet_binned == NULL)
        CreateBinnedDataSet(type);

    numFitParameters = (parameters->selectByAttrib("Constant",kFALSE))->getSize();

    if      (strcmp(type,"a") == 0)   chi2Var = new RooChi2Var("chi2Var","chi2Var",*pdf_a,*dataSet_binned);
    else if (strcmp(type,"b") == 0)   chi2Var = new RooChi2Var("chi2Var","chi2Var",*pdf_b,*dataSet_binned);
    else if (strcmp(type,"ab") == 0)  chi2Var = new RooChi2Var("chi2Var","chi2Var",*pdf_ab,*dataSet_binned);
    else if (strcmp(type,"bb") == 0)  chi2Var = new RooChi2Var("chi2Var","chi2Var",*pdf_bb,*dataSet_binned);

	//chi2Var = new RooChi2Var("chi2Var","chi2Var",*simPdf,*dataSet_binned);

	RooRealVar* ndof     = new RooRealVar("ndof","number of degrees of freedom",0);
	RooRealVar* chi2red  = new RooRealVar("chi2red","reduced chi^2",0);
	RooRealVar* prob     = new RooRealVar("prob","prob(chi2,ndof)",0);

	ndof->setVal(dataSet_binned->numEntries()-numFitParameters);
	chi2red->setVal(chi2Var->getVal()/ndof->getVal()) ;
	prob->setVal(TMath::Prob(chi2Var->getVal(),static_cast<int>(ndof->getVal())));

	delete dataSet_binned;

	printf("%s:\tchi2 = %f\n%s:\tndof = %f\n%s:\tchi2red = %f\n%s:\tprob = %f\n\n",type,chi2Var->getVal(),type,ndof->getVal(),type,chi2red->getVal(),type,prob->getVal());

}


Double_t FitterTransTIndep::GetChi2(const char* type)
{
    Double_t mychi2 = 0;
    Double_t n = 0;
    Double_t v = 0;

//    TH1D* h_dchi2 = new TH1D("dchi2","dchi2",100,0,100000);


    //if(dataSet_binned == NULL)
        CreateBinnedDataSet(type);

    DSRhoPDF* pdf = 0;

    if(strcmp(type,"a") == 0)       pdf = pdf_a;
    else if(strcmp(type,"b") == 0)  pdf = pdf_b;
    else if(strcmp(type,"ab") == 0) pdf = pdf_ab;
    else if(strcmp(type,"bb") == 0) pdf = pdf_bb;

    Double_t binVolume = tht->getBinWidth(0)*thb->getBinWidth(0)*phit->getBinWidth(0)*dt->getBinWidth(0);
    Int_t numBins = dataSet_binned->numEntries();

    Int_t numVPrecise = 0;

    /// Cycle through the centers of all bins
    /// I'm getting width of the first bin, because all bins are of equal width
    for(*tht = tht->getMin()+tht->getBinWidth(0)/2; tht->getVal() < tht->getMax(); tht->setVal(tht->getVal()+tht->getBinWidth(0)))
    {
        for(*thb = thb->getMin()+thb->getBinWidth(0)/2; thb->getVal() < thb->getMax(); thb->setVal(thb->getVal()+thb->getBinWidth(0)))
        {
            for(*phit = phit->getMin()+phit->getBinWidth(0)/2; phit->getVal() < phit->getMax(); phit->setVal(phit->getVal()+phit->getBinWidth(0)))
            {
                for(*dt = dt->getMin()+dt->getBinWidth(0)/2; dt->getVal() < dt->getMax(); dt->setVal(dt->getVal()+dt->getBinWidth(0)))
                {
                    /// Weight is actually the bin content
                    n = dataSet_binned->weight(RooArgSet(*tht,*thb,*phit,*dt),0);
                    if(n < 1)
                    {
                        numBins--;
                        continue;
                    }

                    v = pdf->getVal(RooArgSet(*tht,*thb,*phit,*dt))*binVolume*binnedNumEntries;

                    if(((n-v)*(n-v)/v) > 1)
                    {
                        v = GetVPrecise(pdf);
                        numVPrecise++;
                    }

//                    printf("%.10f\t%.10f\n",v,GetVPrecise(pdf));

                    mychi2 += (n-v)*(n-v)/v;
//                    h_dchi2->Fill((n-v)*(n-v)/v);
                }
            }
        }
    }

    printf("# VPrecise called: %i\n",numVPrecise);

//    TCanvas* c1 = new TCanvas("c1","c1",800,600);
//    c1->SetLogy();
//    h_dchi2->Draw();

    delete dataSet_binned;

    //printf("binVolume = %f\n",binVolume);
    printf("%s: numEntries = %i\n",type,binnedNumEntries);
    printf("%s: numBins = %i\n",type,numBins);
    printf("%s: mychi2 = %f\n",type,mychi2);
    printf("%s: mychi2red = %f\n",type,mychi2/numBins);
    printf("%s: prob = %.10f\n\n",type,TMath::Prob(mychi2,numBins));
    return mychi2;
}

Double_t FitterTransTIndep::GetVPrecise(DSRhoPDF* pdf)
{
    /// The pdf seems to be varying quite rapidly at some places, so the approximation of constant pdf in a voxel is sometimes bad.
    /// This function gives a better approximation of a pdf value in a voxel by averaging through multiple points inside it.

    Double_t v = 0;

    Int_t tht_subbins = 3;
    Int_t thb_subbins = 3;
    Int_t phit_subbins = 3;
    Int_t dt_subbins = 3;

    Double_t binVolume = tht->getBinWidth(0)*thb->getBinWidth(0)*phit->getBinWidth(0)*dt->getBinWidth(0);

    Double_t tht_binmin = tht->getVal() - tht->getBinWidth(0) /2;
    Double_t tht_binmax = tht->getVal() + tht->getBinWidth(0) /2;
    Double_t thb_binmin = thb->getVal() - thb->getBinWidth(0) /2;
    Double_t thb_binmax = thb->getVal() + thb->getBinWidth(0) /2;
    Double_t phit_binmin = phit->getVal() - phit->getBinWidth(0) /2;
    Double_t phit_binmax = phit->getVal() + phit->getBinWidth(0) /2;
    Double_t dt_binmin = dt->getVal() - dt->getBinWidth(0) /2;
    Double_t dt_binmax = dt->getVal() + dt->getBinWidth(0) /2;

    Int_t numPasses = 0;
    Int_t pass = 0;

    /// Using *_binmax - 0.001 because when one is at a boundary of e.g. thb, thb->getVal() < thb_binmax is never violated, even though it should be equal.
    for(*tht = tht_binmin+tht->getBinWidth(0)/(2*tht_subbins); tht->getVal() < tht_binmax-0.001; tht->setVal(tht->getVal()+tht->getBinWidth(0)/tht_subbins))
    {
        for(*thb = thb_binmin+thb->getBinWidth(0)/(2*thb_subbins); thb->getVal() < thb_binmax-0.001; thb->setVal(thb->getVal()+thb->getBinWidth(0)/thb_subbins))
        {
            for(*phit = phit_binmin+phit->getBinWidth(0)/(2*phit_subbins); phit->getVal() < phit_binmax-0.001; phit->setVal(phit->getVal()+phit->getBinWidth(0)/phit_subbins))
            {
                for(*dt = dt_binmin+dt->getBinWidth(0)/(2*dt_subbins); dt->getVal() < dt_binmax-0.001; dt->setVal(dt->getVal()+dt->getBinWidth(0)/dt_subbins))
                {
                    v += pdf->getVal(RooArgSet(*tht,*thb,*phit,*dt))*binVolume*binnedNumEntries;
                }
            }
        }
    }

    /// These lines return the variables to their initial states, before calling GetVPrecise()
    tht->setVal(tht_binmin + tht->getBinWidth(0)/2);
    thb->setVal(thb_binmin + thb->getBinWidth(0)/2);
    phit->setVal(phit_binmin + phit->getBinWidth(0)/2);
    dt->setVal(dt_binmin + dt->getBinWidth(0)/2);

    v = v/(tht_subbins*thb_subbins*phit_subbins*dt_subbins);

    return v;
}

Double_t FitterTransTIndep::GetVPrecise1D(const int i,DSRhoPDF* pdf,RooDataSet* loc_dataset)
{
    const Double_t num_subbins = 5;

    RooArgSet intSet;
    for(int j = 0; j < 4; j++) if(j != i) intSet.add(*vars[j]);

    Double_t original_var = vars[i]->getVal();
    Double_t bin_min = original_var - vars[i]->getBinWidth(0)/2;
    Double_t bin_max = original_var + vars[i]->getBinWidth(0)/2;
    Double_t v = 0;
    RooAbsReal* vr;

    for(*vars[i] = bin_min+tht->getBinWidth(0)/(2*num_subbins); vars[i]->getVal() < bin_max-0.001; vars[i]->setVal(vars[i]->getVal()+vars[i]->getBinWidth(0)/num_subbins))
    {
        vr = pdf->createIntegral(intSet,RooArgSet(*vars[0],*vars[1],*vars[2],*vars[3]));
        v += vr->getVal();
    }

    vars[i]->setVal(original_var);
    v *= vars[i]->getBinWidth(0)*loc_dataset->numEntries()/num_subbins;

    return v;
}

Double_t FitterTransTIndep::GetVPrecise1D(const int i,RooSimultaneous* spdf,RooDataSet* loc_dataset)
{
    const Double_t num_subbins = 5;

    RooArgSet intSet;
    for(int j = 0; j < 4; j++) if(j != i) intSet.add(*vars[j]);

    Double_t original_var = vars[i]->getVal();
    Double_t bin_min = original_var - vars[i]->getBinWidth(0)/2;
    Double_t bin_max = original_var + vars[i]->getBinWidth(0)/2;
    Double_t v = 0;
    RooAbsReal* vr;

    for(*vars[i] = bin_min+tht->getBinWidth(0)/(2*num_subbins); vars[i]->getVal() < bin_max-0.001; vars[i]->setVal(vars[i]->getVal()+vars[i]->getBinWidth(0)/num_subbins))
    {
        vr = spdf->createIntegral(intSet,RooArgSet(*vars[0],*vars[1],*vars[2],*vars[3]));
        v += vr->getVal();
    }

    vars[i]->setVal(original_var);
    v *= vars[i]->getBinWidth(0)*loc_dataset->numEntries()/num_subbins;

    return v;
}

void FitterTransTIndep::SaveResiduals()
{

    TFile* file = new TFile("plots/residuals.root","RECREATE");
    TCanvas* c_residuals = new TCanvas("c2","c2",800,600);

    Double_t n = 0;
    Double_t v = 0;

    /// Variables have to be binned to be able to call ->weight to get bin content
    for(int i = 0; i < 4; i++)
        vars[i]->setBins(vars_bins[i]);

    TString name;
    TString path;
    TH2F* h2_residual[7];
    TH1F* h1_residual_bar[7];
    TH1F* h1_pull[7];
    Double_t chi2[7] = {0,0,0,0,0,0,0};
    Int_t ndof[7] = {0,0,0,0,0,0,0};

    /// Loop for tht, thb and phit. Loop for dt_{a,ab,b,bb} follows.
    for(int i = 0; i < 3; i++)
    {
        name = "h2_residual_";
        name += vars[i]->GetName();
        h2_residual[i] = new TH2F(name,name,vars_bins[i],vars[i]->getMin(),vars[i]->getMax(),50,-5,5);
        name = "h1_residual_bar_";
        name += vars[i]->GetName();
        h1_residual_bar[i] = new TH1F(name,name,vars_bins[i],vars[i]->getMin(),vars[i]->getMax());
        name = "h1_pull_";
        name += vars[i]->GetName();
        h1_pull[i] = new TH1F(name,name,40,-7,7);

        /// Binned dataset with *only one* dimension is created from the whole dataset, because tht,thb,phit distributions
        /// are the same for all four decay types.
        RooDataHist* dataSet_binned_1D = new RooDataHist("dataSet_binned_1D","dataSet_binned_1D",RooArgSet(*vars[i]),*dataSet);
        RooAbsReal* vr;

        for(*vars[i] = vars[i]->getMin()+vars[i]->getBinWidth(0)/2; vars[i]->getVal() < vars[i]->getMax(); vars[i]->setVal(vars[i]->getVal()+vars[i]->getBinWidth(0)))
        {
            n = dataSet_binned_1D->weight(RooArgSet(*vars[i]),0);
            if(n <= 1) continue;
            v = GetVPrecise1D(i,pdf_a,dataSet);

            h2_residual[i]->Fill(vars[i]->getVal(),(n-v)/sqrt(n));
            h1_residual_bar[i]->Fill(vars[i]->getVal(),(n-v)/sqrt(n));
            h1_pull[i]->Fill((n-v)/sqrt(n));

            chi2[i] += ((n-v)*(n-v))/v;
            ndof[i]++;
        }

        delete dataSet_binned_1D;

        h2_residual[i]->GetXaxis()->SetTitle(vars[i]->GetName());
        h2_residual[i]->SetMarkerStyle(7);

        /// Prepare 3-sigma lines for the residual histograms
        TLine baseline(vars[i]->getMin(),0,vars[i]->getMax(),0);
        TLine three_sigma_up(vars[i]->getMin(),3,vars[i]->getMax(),3);
        TLine three_sigma_down(vars[i]->getMin(),-3,vars[i]->getMax(),-3);
        three_sigma_up.SetLineColor(2);
        three_sigma_down.SetLineColor(2);

        c_residuals->SetGrid();
        h2_residual[i]->Draw();
        baseline.Draw();
        three_sigma_up.Draw();
        three_sigma_down.Draw();
        h2_residual[i]->Write();
        path = "plots/residual_";
        path += vars[i]->GetName();
        path += ".png";
        c_residuals->SaveAs(path);
        c_residuals->SetGrid(0,0);

        /// These residual_bar histograms are mainly for debugging purposes, may be disabled
        h1_residual_bar[i]->GetXaxis()->SetTitle(vars[i]->GetName());
        c_residuals->SetGrid();
        h1_residual_bar[i]->Draw();
        h1_residual_bar[i]->Write();
        path = "plots/residual_bar_";
        path += vars[i]->GetName();
        path += ".png";
        c_residuals->SaveAs(path);
        c_residuals->SetGrid(0,0);

        h1_pull[i]->Fit("gaus");
        h1_pull[i]->GetXaxis()->SetTitle(vars[i]->GetName());
        h1_pull[i]->Draw();
        h1_pull[i]->Write();
        path = "plots/pull_";
        path += vars[i]->GetName();
        path += ".png";
        c_residuals->SaveAs(path);
    }

    DSRhoPDF* pdf = 0;
    /// Loop for dt_{a,ab,b,bb}
    for(int i = 3; i < 7; i++)
    {
        char* type;

        switch (i)
        {
        case 3:
            type = "a";
            pdf = pdf_a;
            break;

        case 4:
            type = "ab";
            pdf = pdf_ab;
            break;

        case 5:
            type = "b";
            pdf = pdf_b;
            break;

        case 6:
            type = "bb";
            pdf = pdf_bb;
            break;

        default:
            break;
        }

        name = "h2_residual_";
        name += dt->GetName();
        name += "_";
        name += type;
        h2_residual[i] = new TH2F(name,name,dt_bins,dt->getMin(),dt->getMax(),50,-5,5);

        name = "h2_residual_bar_";
        name += dt->GetName();
        name += "_";
        name += type;
        h1_residual_bar[i] = new TH1F(name,name,dt_bins,dt->getMin(),dt->getMax());

        name = "h1_pull_";
        name += dt->GetName();
        name += "_";
        name += type;
        h1_pull[i] = new TH1F(name,name,40,-7,7);

        CreateReducedDataSet(type);
        RooDataHist* dataSet_binned_1D = new RooDataHist("dataSet_binned_1D","dataSet_binned_1D",RooArgSet(*dt),*dataSet_reduced);

        for(*dt = dt->getMin()+dt->getBinWidth(0)/2; dt->getVal() < dt->getMax(); dt->setVal(dt->getVal()+dt->getBinWidth(0)))
        {
            /// It seems the first and last bin collect some under/overflow and are therefore skipped.
            /// TODO: Should be more thoroughly investigated.
            if(dt->getVal() < dt->getMin()+dt->getBinWidth(0)) continue;
            else if(dt->getVal() > dt->getMax()-dt->getBinWidth(0)) continue;

            n = dataSet_binned_1D->weight(RooArgSet(*dt),0);
            if(n <= 10) continue;
            v = GetVPrecise1D(3,pdf,dataSet_reduced);

            h2_residual[i]->Fill(dt->getVal(),(n-v)/sqrt(n));
            h1_residual_bar[i]->Fill(dt->getVal(),(n-v)/sqrt(n));
            h1_pull[i]->Fill((n-v)/sqrt(n));

            chi2[i] += ((n-v)*(n-v))/v;
            ndof[i]++;
        }

        delete dataSet_reduced;
        delete dataSet_binned_1D;

        h2_residual[i]->GetXaxis()->SetTitle(dt->GetName());
        h2_residual[i]->SetMarkerStyle(7);

        /// Prepare 3-sigma lines for the residual histograms
        TLine baseline(dt->getMin(),0,dt->getMax(),0);
        TLine three_sigma_up(dt->getMin(),3,dt->getMax(),3);
        TLine three_sigma_down(dt->getMin(),-3,dt->getMax(),-3);
        three_sigma_up.SetLineColor(2);
        three_sigma_down.SetLineColor(2);

        c_residuals->SetGrid();
        h2_residual[i]->Draw();
        baseline.Draw();
        three_sigma_up.Draw();
        three_sigma_down.Draw();
        h2_residual[i]->Write();
        path = "plots/residual_";
        path += dt->GetName();
        path += "_";
        path += type;
        path += ".png";
        c_residuals->SaveAs(path);
        c_residuals->SetGrid(0,0);

        /// These residual_bar histograms are mainly for debugging purposes, may be disabled
        h1_residual_bar[i]->GetXaxis()->SetTitle(dt->GetName());
        c_residuals->SetGrid();
        h1_residual_bar[i]->Draw();
        h1_residual_bar[i]->Write();
        path = "plots/residual_bar_";
        path += dt->GetName();
        path += "_";
        path += type;
        path += ".png";
        c_residuals->SaveAs(path);
        c_residuals->SetGrid(0,0);

        h1_pull[i]->Fit("gaus");
        h1_pull[i]->GetXaxis()->SetTitle(dt->GetName());
        h1_pull[i]->Draw();
        h1_pull[i]->Write();
        path = "plots/pull_";
        path += dt->GetName();
        path += "_";
        path += type;
        path += ".png";
        c_residuals->SaveAs(path);
        c_residuals->SetGrid(0,0);
    }

    file->Close();
    delete c_residuals;

    ///This is outside of the preceding loop because it would be intersparsed by different messages
    for(int i = 0; i < 3; i++)
        printf("%s\tchi2: %.2f\tndof: %i\tchi2red: %.3f\tprob: %f\n",vars[i]->GetName(),chi2[i],ndof[i],chi2[i]/ndof[i],TMath::Prob(chi2[i],ndof[i]));

    for(int i = 3; i < 7; i++)
        printf("%s_%i\tchi2: %.2f\tndof: %i\tchi2red: %.3f\tprob: %f\n",dt->GetName(),i,chi2[i],ndof[i],chi2[i]/ndof[i],TMath::Prob(chi2[i],ndof[i]));

}

Double_t FitterTransTIndep::SaveChi2Maps(const char* type)
{


//    Double_t mychi2 = 0;
//    Double_t dchi2 = 0;
//    Double_t n = 0;
//    Double_t v = 0;
//
//    RooRealVar* var1 = tht;
//    RooRealVar* var2 = thb;
//    RooRealVar* var3 = phit;
//    RooRealVar* var4 = dt;
//    Int_t var1_bins = tht_bins;
//    Int_t var2_bins = thb_bins;
//    Int_t var3_bins = phit_bins;
//    Int_t var4_bins = dt_bins;
//
//    DSRhoPDF* pdf = 0;
//
//    TFile* file = new TFile("plots/plots.root","RECREATE");
//    TCanvas* c2 = new TCanvas("c2","c2",800,600);
//    TString path;
//
//    if(strcmp(type,"a") == 0)       pdf = pdf_a;
//    else if(strcmp(type,"b") == 0)  pdf = pdf_b;
//    else if(strcmp(type,"ab") == 0) pdf = pdf_ab;
//    else if(strcmp(type,"bb") == 0) pdf = pdf_bb;

//    //if(dataSet_binned == NULL)
//        CreateBinnedDataSet(type);
//
//    Double_t binVolume = tht->getBinWidth(0)*thb->getBinWidth(0)*phit->getBinWidth(0)*dt->getBinWidth(0);
//
//    TH1F* h1_chi2 = new TH1F("h1_chi2","h1_chi2",100,0,10);
//    TH2F* h2_chi2_1 = new TH2F("h2_chi2_1","h2_chi2_1",var1_bins,var1->getMin(),var1->getMax(),var2_bins,var2->getMin(),var2->getMax());
//    TH2F* h2_chi2_2 = new TH2F("h2_chi2_2","h2_chi2_2",var1_bins,var1->getMin(),var1->getMax(),var3_bins,var3->getMin(),var3->getMax());
//    TH2F* h2_chi2_3 = new TH2F("h2_chi2_3","h2_chi2_3",var1_bins,var1->getMin(),var1->getMax(),var4_bins,var4->getMin(),var4->getMax());
//    TH2F* h2_chi2_4 = new TH2F("h2_chi2_4","h2_chi2_4",var2_bins,var2->getMin(),var2->getMax(),var3_bins,var3->getMin(),var3->getMax());
//    TH2F* h2_chi2_5 = new TH2F("h2_chi2_5","h2_chi2_5",var2_bins,var2->getMin(),var2->getMax(),var4_bins,var4->getMin(),var4->getMax());
//    TH2F* h2_chi2_6 = new TH2F("h2_chi2_6","h2_chi2_6",var3_bins,var3->getMin(),var3->getMax(),var4_bins,var4->getMin(),var4->getMax());


    /// Cycle through the centers of all bins
    /// I'm getting width of the first bin, because all bins are of equal width
//    for(*var1 = var1->getMin()+var1->getBinWidth(0)/2; var1->getVal() < var1->getMax(); var1->setVal(var1->getVal()+var1->getBinWidth(0)))
//    {
//        for(*var2 = var2->getMin()+var2->getBinWidth(0)/2; var2->getVal() < var2->getMax(); var2->setVal(var2->getVal()+var2->getBinWidth(0)))
//        {
//            for(*var3 = var3->getMin()+var3->getBinWidth(0)/2; var3->getVal() < var3->getMax(); var3->setVal(var3->getVal()+var3->getBinWidth(0)))
//            {
//                for(*var4 = var4->getMin()+var4->getBinWidth(0)/2; var4->getVal() < var4->getMax(); var4->setVal(var4->getVal()+var4->getBinWidth(0)))
//                {
//
//                    /// Weight is actually the bin content
//                    n = dataSet_binned->weight(RooArgSet(*var1,*var2,*var3,*var4),0);
//                    //if(n == 0) continue;
//                    if(n == 0) continue;
//                    v = pdf->getVal(RooArgSet(*var1,*var2,*var3,*var4))*binVolume*binnedNumEntries;
//
//                    if(((n-v)*(n-v)/v) > 1)
//                    {
//                        v = GetVPrecise(pdf);
//                    }
//
//                    dchi2 = (n-v)*(n-v)/v;
//
//                    //h1_chi2->Fill(dchi2);
//                    h2_chi2_1->Fill(var1->getVal(),var2->getVal(),dchi2);
//                    h2_chi2_2->Fill(var1->getVal(),var3->getVal(),dchi2);
//                    h2_chi2_3->Fill(var1->getVal(),var4->getVal(),dchi2);
//                    h2_chi2_4->Fill(var2->getVal(),var3->getVal(),dchi2);
//                    h2_chi2_5->Fill(var2->getVal(),var4->getVal(),dchi2);
//                    h2_chi2_6->Fill(var3->getVal(),var4->getVal(),dchi2);
//                    mychi2 += dchi2;
//
//                    for(int i = 0; i < 4; i++)
//                        h1_resid[i]->Fill(vars[i]->getVal(),(n-v)/v);
//                }
//            }
//        }
//    }

//    delete dataSet_binned;



    //c2->SetLogy(kTRUE);
//    h1_chi2->GetXaxis()->SetTitle("dchi");
//    h1_chi2->GetYaxis()->SetTitle("num bins");
//    h1_chi2->Draw();
//    h1_chi2->Write();
//    c2->SaveAs("plots/chi2_delta.png");
//
//    c2->SetLogy(kFALSE);

//    h2_chi2_1->SetOption("colz");
//    //hdchi2->SetMinimum(0);
//    //hdchi2->SetMaximum(100);
//    h2_chi2_1->SetStats(kFALSE);
//    h2_chi2_1->GetXaxis()->SetTitle(var1->GetName());
//    h2_chi2_1->GetYaxis()->SetTitle(var2->GetName());
//    h2_chi2_1->Draw();
//    h2_chi2_1->Write();
//    path = "plots/chi2map_";
//    path += var1->GetName();
//    path += "_";
//    path += var2->GetName();
//    path += ".png";
//    c2->SaveAs(path);
//
//    h2_chi2_2->SetOption("colz");
//    //hdchi2->SetMinimum(0);
//    //hdchi2->SetMaximum(100);
//    h2_chi2_2->SetStats(kFALSE);
//    h2_chi2_2->GetXaxis()->SetTitle(var1->GetName());
//    h2_chi2_2->GetYaxis()->SetTitle(var3->GetName());
//    h2_chi2_2->Draw();
//    h2_chi2_2->Write();
//    path = "plots/chi2map_";
//    path += var1->GetName();
//    path += "_";
//    path += var3->GetName();
//    path += ".png";
//    c2->SaveAs(path);
//
//    h2_chi2_3->SetOption("colz");
//    //hdchi2->SetMinimum(0);
//    //hdchi2->SetMaximum(100);
//    h2_chi2_3->SetStats(kFALSE);
//    h2_chi2_3->GetXaxis()->SetTitle(var1->GetName());
//    h2_chi2_3->GetYaxis()->SetTitle(var4->GetName());
//    h2_chi2_3->Draw();
//    h2_chi2_3->Write();
//    path = "plots/chi2map_";
//    path += var1->GetName();
//    path += "_";
//    path += var4->GetName();
//    path += ".png";
//    c2->SaveAs(path);
//
//    h2_chi2_4->SetOption("colz");
//    //hdchi2->SetMinimum(0);
//    //hdchi2->SetMaximum(100);
//    h2_chi2_4->SetStats(kFALSE);
//    h2_chi2_4->GetXaxis()->SetTitle(var2->GetName());
//    h2_chi2_4->GetYaxis()->SetTitle(var3->GetName());
//    h2_chi2_4->Draw();
//    h2_chi2_4->Write();
//    path = "plots/chi2map_";
//    path += var2->GetName();
//    path += "_";
//    path += var3->GetName();
//    path += ".png";
//    c2->SaveAs(path);
//
//    h2_chi2_5->SetOption("colz");
//    //hdchi2->SetMinimum(0);
//    //hdchi2->SetMaximum(100);
//    h2_chi2_5->SetStats(kFALSE);
//    h2_chi2_5->GetXaxis()->SetTitle(var2->GetName());
//    h2_chi2_5->GetYaxis()->SetTitle(var4->GetName());
//    h2_chi2_5->Draw();
//    h2_chi2_5->Write();
//    path = "plots/chi2map_";
//    path += var2->GetName();
//    path += "_";
//    path += var4->GetName();
//    path += ".png";
//    c2->SaveAs(path);
//
//    h2_chi2_6->SetOption("colz");
//    //h2_chi2_6->SetMinimum(0);
//    //h2_chi2_6->SetMaximum(100000);
//    h2_chi2_6->SetStats(kFALSE);
//    h2_chi2_6->GetXaxis()->SetTitle(var3->GetName());
//    h2_chi2_6->GetYaxis()->SetTitle(var4->GetName());
//    h2_chi2_6->Draw();
//    h2_chi2_6->Write();
//    path = "plots/chi2map_";
//    path += var3->GetName();
//    path += "_";
//    path += var4->GetName();
//    path += ".png";
//    c2->SaveAs(path);

//    file->Close();

//    delete h1_chi2;
//    delete h2_chi2_1;
//    delete h2_chi2_2;
//    delete h2_chi2_3;
//    delete h2_chi2_4;
//    delete h2_chi2_5;
//    delete h2_chi2_6;
//    delete c2;
//
//    return mychi2;
}

void FitterTransTIndep::SaveNllPlot(RooRealVar* var)
{
    const int steps = 100;
    Double_t orig_val = var->getVal();
    TCanvas c_nll("c_nll","c_nll",800,600);
    TString name;
    TString basename;
    TString path;

    basename = "nll_";
    basename += var->GetName();
    name = "h1_" + basename;

    TH1F h1_nll(name,name,steps,var->getMin(),var->getMax());
    RooAbsReal* nll;
    nll = simPdf->createNLL(*dataSet,RooFit::NumCPU(2));
    for(Int_t i = 0; i < steps; i++)
    {
        var->setVal(i*var->getMax()/steps+(var->getMax()-var->getMin())/(2*steps));
        printf("Computing %i/%i likelihood function.\n",i+1,steps);
        h1_nll.Fill(var->getVal(),2*nll->getVal());
    }
    delete nll;
    h1_nll.GetXaxis()->SetTitle(var->GetName());
    h1_nll.SetStats(kFALSE);
    h1_nll.Draw();
    h1_nll.Write();
    path = "plots/";
    path += basename;
    path += ".png";
    c_nll.SaveAs(path);
    var->setVal(orig_val);
}

void FitterTransTIndep::SaveNllPlot(RooRealVar* var1, RooRealVar* var2)
{
    const int steps1 = 30;
    const int steps2 = 30;
    Double_t orig_val1 = var1->getVal();
    Double_t orig_val2 = var2->getVal();
    TCanvas c_nll("c_nll","c_nll",600,600);
    TString name;
    TString basename;
    TString path;
    basename = "nll2_";
    basename += var1->GetName();
    basename += "_";
    basename += var2->GetName();
    name = "h2_" + basename;

    TH2F h2_nll(name,name,steps1,var1->getMin(),var1->getMax(),steps2,var2->getMin(),var2->getMax());
    RooAbsReal* nll;
    nll = simPdf->createNLL(*dataSet,RooFit::NumCPU(2));
    for(Int_t i = 0; i < steps1; i++)
    {
        for(Int_t j = 0; j < steps2; j++)
        {
            var1->setVal(i*var1->getMax()/steps1+(var1->getMax()-var1->getMin())/(2*steps1));
            var2->setVal(j*var2->getMax()/steps2+(var2->getMax()-var2->getMin())/(2*steps2));
            printf("Computing %i/%i likelihood function.\n",i*steps2+j+1,steps1*steps2);
            h2_nll.Fill(var1->getVal(),var2->getVal(),2*nll->getVal());
        }
    }

    delete nll;

    h2_nll.GetXaxis()->SetTitle(var1->GetName());
    h2_nll.GetYaxis()->SetTitle(var2->GetName());
    h2_nll.SetStats(kFALSE);
    h2_nll.SetOption("colz");
    h2_nll.Draw();
    h2_nll.Write();
    path = "plots/";
    path += basename;
    path += ".png";
    c_nll.SaveAs(path);
    var1->setVal(orig_val1);
    var2->setVal(orig_val2);
}

void FitterTransTIndep::GetRecoveredParameters(Int_t& numParameters, Double_t** recoveredParameters)
{
    RooRealVar* chi2red  = new RooRealVar("chi2red","reduced chi^2",0);
    chi2red->setVal(chi2Var->getVal()/(dataSet_binned->numEntries()-numFitParameters));

    numParameters = 17;
    Double_t* parameters = new Double_t[numParameters];

    parameters[0] = chi2red->getVal();
    parameters[1] = ap->getVal();
    parameters[2] = ap->getError();
    parameters[3] = apa->getVal();
    parameters[4] = apa->getError();
    parameters[5] = a0->getVal();
    parameters[6] = a0->getError();
    parameters[7] = a0a->getVal();
    parameters[8] = a0a->getError();
    parameters[9] = at->getVal();
    parameters[10] = at->getPropagatedError(*result);
    parameters[11] = ata->getVal();
    parameters[12] = ata->getError();
    parameters[13] = par_input[0];
    parameters[14] = par_input[1];
    parameters[15] = par_input[2];
    parameters[16] = par_input[3];

    *recoveredParameters = parameters;
}

RooDataHist* FitterTransTIndep::GetBinnedDataSet()
{
    if(dataSet_binned == NULL)
        CreateBinnedDataSet("a");

    return dataSet_binned;
}

RooDataSet* FitterTransTIndep::GetReducedDataSet()
{
    if(dataSet_reduced == NULL)
        CreateBinnedDataSet("a");

    return dataSet_reduced;
}

void FitterTransTIndep::FixAllParameters()
{
    RooRealVar* rooPar = 0;
    TIterator* parIter = parameters->createIterator();
    while(rooPar = (RooRealVar*)parIter->Next())
        rooPar->setConstant();

}

void FitterTransTIndep::FixParameter(const char* par)
{
    RooRealVar* rooPar = 0;
    rooPar = (RooRealVar*)parameters->find(par);
    if(rooPar != 0)
        rooPar->setConstant();
}

void FitterTransTIndep::FreeParameter(const char* par)
{
    RooRealVar* rooPar = 0;
    rooPar = (RooRealVar*)parameters->find(par);
    if(rooPar != 0)
        rooPar->setConstant(kFALSE);
}

void FitterTransTIndep::PrintParameter(const char* par)
{
    RooRealVar* rooPar = 0;
    rooPar = (RooRealVar*)parameters->find(par);
    if(rooPar != 0)
        rooPar->Print();

    else if(strcmp(par,at->GetName()) == 0)
        printf("%s = %f +- %f\n",par,at->getVal(),at->getPropagatedError(*result));
    else if(strcmp(par,a0a->GetName()) == 0)
        a0a->Print();
    else if(strcmp(par,"hp") == 0)
    {
        RooFormulaVar hpr("hpr","(apr + atr)/sqrt(2)",RooArgSet(*apr,*atr));
        RooFormulaVar hpi("hpi","(api + ati)/sqrt(2)",RooArgSet(*api,*ati));
        RooFormulaVar hp("hp","sqrt(hpr*hpr+hpi*hpi)",RooArgSet(hpr,hpi));
        RooFormulaVar hpa("hpa","atan2(hpi,hpr)",RooArgSet(hpr,hpi));

        printf("%s = %f +- %f\n","hp",hp.getVal(),hp.getPropagatedError(*result));
        printf("%s = %f +- %f\n","hpa",hpa.getVal(),hpa.getPropagatedError(*result));
    }
    else if(strcmp(par,"hm") == 0)
    {
        RooFormulaVar hmr("hmr","(apr - atr)/sqrt(2)",RooArgSet(*apr,*atr));
        RooFormulaVar hmi("hmi","(api - ati)/sqrt(2)",RooArgSet(*api,*ati));
        RooFormulaVar hm("hm","sqrt(hmr*hmr+hmi*hmi)",RooArgSet(hmr,hmi));
        RooFormulaVar hma("hma","atan2(hmi,hmr)",RooArgSet(hmr,hmi));

        printf("%s = %f +- %f\n","hm",hm.getVal(),hm.getPropagatedError(*result));
        printf("%s = %f +- %f\n","hma",hma.getVal(),hma.getPropagatedError(*result));
    }
    else
        printf("ERROR: Parameter '%s' doesn't exist! Can't print its value.\n",par);
}

void FitterTransTIndep::GetHelParameters(Double_t* params)
{
    RooFormulaVar hpr("hpr","(apr + atr)/sqrt(2)",RooArgSet(*apr,*atr));
    RooFormulaVar hpi("hpi","(api + ati)/sqrt(2)",RooArgSet(*api,*ati));
    RooFormulaVar hp("hp","sqrt(hpr*hpr+hpi*hpi)",RooArgSet(hpr,hpi));
    RooFormulaVar hpa("hpa","atan2(hpi,hpr)",RooArgSet(hpr,hpi));

    RooFormulaVar hmr("hmr","(apr - atr)/sqrt(2)",RooArgSet(*apr,*atr));
    RooFormulaVar hmi("hmi","(api - ati)/sqrt(2)",RooArgSet(*api,*ati));
    RooFormulaVar hm("hm","sqrt(hmr*hmr+hmi*hmi)",RooArgSet(hmr,hmi));
    RooFormulaVar hma("hma","atan2(hmi,hmr)",RooArgSet(hmr,hmi));

    params[0] = hp.getVal();
    params[1] = hp.getPropagatedError(*result);
    params[2] = hpa.getVal();
    params[3] = hpa.getPropagatedError(*result);
    params[4] = a0->getVal();
    params[5] = a0->getError();
    params[6] = hm.getVal();
    params[7] = hm.getPropagatedError(*result);
    params[8] = hma.getVal();
    params[9] = hma.getPropagatedError(*result);
}

void FitterTransTIndep::SaveNllPlot(const char* par)
{
    RooRealVar* rooPar = 0;
    rooPar = (RooRealVar*)parameters->find(par);
    if(rooPar != 0)
        SaveNllPlot(rooPar);
    else
        printf("ERROR: Parameter '%s' doesn't exist! Can't save Nll plot.\n",par);
}

void FitterTransTIndep::SaveNllPlot(const char* par1, const char* par2)
{
    RooRealVar* rooPar1 = 0;
    RooRealVar* rooPar2 = 0;
    rooPar1 = (RooRealVar*)parameters->find(par1);
    rooPar2 = (RooRealVar*)parameters->find(par2);
    if(rooPar1 != 0 && rooPar2 != 0)
        SaveNllPlot(rooPar1,rooPar2);
    if(rooPar1 == 0)
        printf("ERROR: Parameter %s doesn't exist! Can't save Nll plot.\n",par1);
    if(rooPar2 == 0)
        printf("ERROR: Parameter %s doesn't exist! Can't save Nll plot.\n",par2);
}






