#include "RooRealVar.h"
#include "RooFormulaVar.h"
#include "RooCategory.h"
#include "RooGenericPdf.h"
#include "Constants.h"
#include "RooDataSet.h"
#include "RooRandom.h"
#include "RooChi2Var.h"
#include "RooSimultaneous.h"
#include "RooPlot.h"
#include "RooFitResult.h"

#include "TMath.h"
#include "TIterator.h"
#include "TLine.h"

#include "DSRhoPDF.h"
#include "FitterTrans.h"

#include "TH1D.h"
#include "TH2F.h"
#include "TCanvas.h"
#include "TFile.h"

FitterTrans::FitterTrans(Double_t* outer_par_input, bool outer_time_dependent)
{
    chi2Var = 0;
    result = 0;

    tht_bins = 60;
    thb_bins = 60;
    phit_bins = 60;
    dt_bins = 80;

    thb = new RooRealVar("thb","thb",0,PI);
    tht = new RooRealVar("tht","tht",0,PI);
    phit = new RooRealVar("phit","phit",-PI,PI);
    dt = new RooRealVar("dt","dt",-3,3);

    /// Convenient shorthands that can be used in loops
    vars[0] = tht;
    vars[1] = thb;
    vars[2] = phit;
    vars[3] = dt;
    vars_bins[0] = tht_bins;
    vars_bins[1] = thb_bins;
    vars_bins[2] = phit_bins;
    vars_bins[3] = dt_bins;

    decType = new RooCategory("decType","decType");
    decType->defineType("a",1);
    decType->defineType("ab",2);
    decType->defineType("b",3);
    decType->defineType("bb",4);

    time_dependent = outer_time_dependent;
    int num_pars = time_dependent ? 11 : 4;
    for (int i = 0; i < num_pars; i++)
        par_input[i] = outer_par_input[i];

    /// apr and other seemingly unnecessary vars are used when calculating/printing e.g. "hp"
    ap = new RooRealVar("ap","ap",par_input[0],0,0.5);
    apa = new RooRealVar("apa","apa",par_input[1],0,1);
    apr = new RooFormulaVar("apr","ap*cos(apa)",RooArgSet(*ap,*apa));
    api = new RooFormulaVar("api","ap*sin(apa)",RooArgSet(*ap,*apa));
    a0 = new RooRealVar("a0","a0",par_input[2],0.8,1);
    a0a = new RooRealVar("a0a","a0a",0);
    a0r = new RooFormulaVar("a0r","a0*cos(a0a)",RooArgSet(*a0,*a0a));
    a0i = new RooFormulaVar("a0i","a0*sin(a0a)",RooArgSet(*a0,*a0a));
    at = new RooFormulaVar("at","sqrt(1-ap*ap-a0*a0)",RooArgSet(*ap,*a0));
    ata = new RooRealVar("ata","ata",par_input[3],2,4);
    atr = new RooFormulaVar("atr","at*cos(ata)",RooArgSet(*at,*ata));
    ati = new RooFormulaVar("ati","at*sin(ata)",RooArgSet(*at,*ata));

    if (time_dependent) {
        /// Time-dep additional vars

        dm = new RooRealVar("dm","dm",0.507e12);
        phiw = new RooRealVar("phiw","phiw",par_input[4],0,2*PI);

        rp = new RooRealVar("rp","rp",par_input[5],0.001,0.2);
        r0 = new RooRealVar("r0","r0",par_input[6],0.001,0.2);
        rt = new RooRealVar("rt","rt",par_input[7],0.001,0.2); /// eq. (100) in BN419 approximates this

        /// s is strong phase; delta_polarization in BN419
        sp = new RooRealVar("sp","sp",par_input[8],-3*PI,3*PI);
        s0 = new RooRealVar("s0","s0",par_input[9],-3*PI,3*PI);
        st = new RooRealVar("st","st",par_input[10],-3*PI,3*PI);

        pdf_a = new DSRhoPDF("pdf_a","pdf_a","a",*tht,*thb,*phit,*dt,*ap,*apa,*a0,*ata,*phiw,*rp,*r0,*rt,*sp,*s0,*st);
        pdf_b = new DSRhoPDF("pdf_b","pdf_b","b",*tht,*thb,*phit,*dt,*ap,*apa,*a0,*ata,*phiw,*rp,*r0,*rt,*sp,*s0,*st);
        pdf_ab = new DSRhoPDF("pdf_ab","pdf_ab","ab",*tht,*thb,*phit,*dt,*ap,*apa,*a0,*ata,*phiw,*rp,*r0,*rt,*sp,*s0,*st);
        pdf_bb = new DSRhoPDF("pdf_bb","pdf_bb","bb",*tht,*thb,*phit,*dt,*ap,*apa,*a0,*ata,*phiw,*rp,*r0,*rt,*sp,*s0,*st);

        pdf_sim = new RooSimultaneous("pdf_sim","pdf_sim",*decType);
        pdf_sim->addPdf(*pdf_a,"a");
        pdf_sim->addPdf(*pdf_ab,"ab");
        pdf_sim->addPdf(*pdf_b,"b");
        pdf_sim->addPdf(*pdf_bb,"bb");

        parameters = new RooArgSet(*ap,*apa,*a0,*ata,*phiw,*rp,*r0,*rt);
        parameters->add(*sp);
        parameters->add(*s0);
        parameters->add(*st);

        variables = new RooArgList(*tht,*thb,*phit,*dt,*decType);
    } else {
        pdf_tindep = new DSRhoPDFTIndep("pdf_tindep","pdf_tindep",*tht,*thb,*phit,*ap,*apa,*a0,*ata);
        parameters = new RooArgSet(*ap,*apa,*a0,*ata);
        variables = new RooArgList(*tht,*thb,*phit);
    }

    /// numFitParameters holds # of NON-constant fit parameters
    numFitParameters = (parameters->selectByAttrib("Constant",kFALSE))->getSize();

}

FitterTrans::~FitterTrans()
{
    delete thb;
    delete tht;
    delete phit;
    delete dt;
    delete decType;

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

    delete dm;
    delete phiw;
    delete rp;
    delete r0;
    delete rt;
    delete sp;
    delete s0;
    delete st;
    delete pdf_a;
    delete pdf_b;
    delete pdf_ab;
    delete pdf_bb;
    delete pdf_tindep;
    delete pdf_sim;
    delete parameters;
    delete fitParameters;
}

Int_t FitterTrans::Fit()
{
    numFitParameters = (parameters->selectByAttrib("Constant",kFALSE))->getSize();

    RooAbsPdf* pdf;
    if (time_dependent) {
        pdf = pdf_sim;
    } else {
        pdf = pdf_tindep;
    }
    printf("*********** cpus = %i *********\n", num_CPUs);
    result = pdf->fitTo(*dataSet,RooFit::Save(),RooFit::Timer(true),\
                           RooFit::Minos(0),RooFit::Hesse(1),RooFit::Strategy(1),RooFit::NumCPU(num_CPUs));

    const TMatrixDSym& cor = result->correlationMatrix();
    result->Print();
    cor.Print();
    //TCanvas c1;
    //result->correlationHist()->Draw("colz");
    //c1.SaveAs("corr.gif");
}

void FitterTrans::ReadDataSet(const char* file) {
    dataSet = RooDataSet::read(file,*variables);
}

void FitterTrans::GenerateDataSet(Int_t numEvents)
{
    RooRandom::randomGenerator()->SetSeed(0);

    /// TODO: For some reason this is the generated fraction of favored decays, have to check why
    const Double_t fracFav = 0.812;
    const Double_t fracSup = 1 - fracFav;

    if(dataSet != 0)
        delete dataSet;

    RooDataSet* temp_dataSet;
    temp_dataSet = pdf_a->generate(RooArgSet(*tht,*thb,*phit,*dt),TMath::Nint(numEvents*fracFav/2));
    decType->setLabel("a");
    temp_dataSet->addColumn(*decType);

    dataSet = (RooDataSet*)temp_dataSet->Clone();
    delete temp_dataSet;

    temp_dataSet = pdf_ab->generate(RooArgSet(*tht,*thb,*phit,*dt),TMath::Nint(numEvents*fracFav/2));
    decType->setLabel("ab");
    temp_dataSet->addColumn(*decType);
    dataSet->append(*temp_dataSet);
    delete temp_dataSet;

    temp_dataSet = pdf_b->generate(RooArgSet(*tht,*thb,*phit,*dt),TMath::Nint(numEvents*fracSup/2));
    decType->setLabel("b");
    temp_dataSet->addColumn(*decType);
    dataSet->append(*temp_dataSet);
    delete temp_dataSet;

    temp_dataSet = pdf_bb->generate(RooArgSet(*tht,*thb,*phit,*dt),TMath::Nint(numEvents*fracSup/2));
    decType->setLabel("bb");
    temp_dataSet->addColumn(*decType);
    dataSet->append(*temp_dataSet);
    delete temp_dataSet;
}

void FitterTrans::CreateBinnedDataSet(const char* type)
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

void FitterTrans::CreateReducedDataSet(const char* type)
{
    TString cut = "decType==decType::";
    cut += type;
    RooRandom::randomGenerator()->SetSeed(0);
    dataSet_reduced = (RooDataSet*)dataSet->reduce(cut);
}

Int_t FitterTrans::ComputeChi2(const char* type)
{
    //if(dataSet_binned == NULL)
        CreateBinnedDataSet(type);

    numFitParameters = (parameters->selectByAttrib("Constant",kFALSE))->getSize();

    if      (strcmp(type,"a") == 0)   chi2Var = new RooChi2Var("chi2Var","chi2Var",*pdf_a,*dataSet_binned);
    else if (strcmp(type,"b") == 0)   chi2Var = new RooChi2Var("chi2Var","chi2Var",*pdf_b,*dataSet_binned);
    else if (strcmp(type,"ab") == 0)  chi2Var = new RooChi2Var("chi2Var","chi2Var",*pdf_ab,*dataSet_binned);
    else if (strcmp(type,"bb") == 0)  chi2Var = new RooChi2Var("chi2Var","chi2Var",*pdf_bb,*dataSet_binned);

	//chi2Var = new RooChi2Var("chi2Var","chi2Var",*pdf_sim,*dataSet_binned);

	RooRealVar* ndof     = new RooRealVar("ndof","number of degrees of freedom",0);
	RooRealVar* chi2red  = new RooRealVar("chi2red","reduced chi^2",0);
	RooRealVar* prob     = new RooRealVar("prob","prob(chi2,ndof)",0);

	ndof->setVal(dataSet_binned->numEntries()-numFitParameters);
	chi2red->setVal(chi2Var->getVal()/ndof->getVal()) ;
	prob->setVal(TMath::Prob(chi2Var->getVal(),static_cast<int>(ndof->getVal())));

	delete dataSet_binned;

	printf("%s:\tchi2 = %f\n%s:\tndof = %f\n%s:\tchi2red = %f\n%s:\tprob = %f\n\n",type,chi2Var->getVal(),type,ndof->getVal(),type,chi2red->getVal(),type,prob->getVal());

}


Double_t FitterTrans::GetChi2(const char* type)
{
    Double_t mychi2 = 0;
    Double_t n = 0;
    Double_t v = 0;

//    TH1D* h_dchi2 = new TH1D("dchi2","dchi2",100,0,100000);


    //if(dataSet_binned == NULL)
        CreateBinnedDataSet(type);

    DSRhoPDF* pdf_temp = 0;

    if(strcmp(type,"a") == 0)       pdf_temp = pdf_a;
    else if(strcmp(type,"b") == 0)  pdf_temp = pdf_b;
    else if(strcmp(type,"ab") == 0) pdf_temp = pdf_ab;
    else if(strcmp(type,"bb") == 0) pdf_temp = pdf_bb;

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

                    v = pdf_temp->getVal(RooArgSet(*tht,*thb,*phit,*dt))*binVolume*binnedNumEntries;

                    if(((n-v)*(n-v)/v) > 1)
                    {
                        v = GetVPrecise(pdf_temp);
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

Double_t FitterTrans::GetVPrecise(DSRhoPDF* pdf)
{
    /// The pdf seems to be varying quite rapidly at some places, so the approximation of constant pdf_temp in a voxel is sometimes bad.
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

Double_t FitterTrans::GetVPrecise1D(const int i,DSRhoPDF* pdf,RooDataSet* loc_dataset)
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

Double_t FitterTrans::GetVPrecise1D(const int i,RooSimultaneous* spdf,RooDataSet* loc_dataset)
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

void FitterTrans::SaveResiduals()
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

    DSRhoPDF* pdf_temp = 0;
    /// Loop for dt_{a,ab,b,bb}
    for(int i = 3; i < 7; i++)
    {
        char* type;

        switch (i)
        {
        case 3:
            type = "a";
            pdf_temp = pdf_a;
            break;

        case 4:
            type = "ab";
            pdf_temp = pdf_ab;
            break;

        case 5:
            type = "b";
            pdf_temp = pdf_b;
            break;

        case 6:
            type = "bb";
            pdf_temp = pdf_bb;
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
            v = GetVPrecise1D(3,pdf_temp,dataSet_reduced);

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

Double_t FitterTrans::SaveChi2Maps(const char* type) {

    Double_t mychi2 = 0;
    Double_t dchi2 = 0;
    Double_t n = 0;
    Double_t v = 0;

    RooRealVar* var1 = tht;
    RooRealVar* var2 = thb;
    RooRealVar* var3 = phit;
    RooRealVar* var4 = dt;
    Int_t var1_bins = tht_bins;
    Int_t var2_bins = thb_bins;
    Int_t var3_bins = phit_bins;
    Int_t var4_bins = dt_bins;

    DSRhoPDF* pdf_temp = 0;

    TFile* file = new TFile("plots/plots.root","RECREATE");
    TCanvas* c2 = new TCanvas("c2","c2",800,600);
    TString path;

    if(strcmp(type,"a") == 0)       pdf_temp = pdf_a;
    else if(strcmp(type,"b") == 0)  pdf_temp = pdf_b;
    else if(strcmp(type,"ab") == 0) pdf_temp = pdf_ab;
    else if(strcmp(type,"bb") == 0) pdf_temp = pdf_bb;

    //if(dataSet_binned == NULL)
    CreateBinnedDataSet(type);

    Double_t binVolume = tht->getBinWidth(0)*thb->getBinWidth(0)*phit->getBinWidth(0)*dt->getBinWidth(0);

    TH1F* h1_chi2 = new TH1F("h1_chi2","h1_chi2",100,0,1000000);
    TH2F* h2_chi2_1 = new TH2F("h2_chi2_1","h2_chi2_1",var1_bins,var1->getMin(),var1->getMax(),var2_bins,var2->getMin(),var2->getMax());
    TH2F* h2_chi2_2 = new TH2F("h2_chi2_2","h2_chi2_2",var1_bins,var1->getMin(),var1->getMax(),var3_bins,var3->getMin(),var3->getMax());
    TH2F* h2_chi2_3 = new TH2F("h2_chi2_3","h2_chi2_3",var1_bins,var1->getMin(),var1->getMax(),var4_bins,var4->getMin(),var4->getMax());
    TH2F* h2_chi2_4 = new TH2F("h2_chi2_4","h2_chi2_4",var2_bins,var2->getMin(),var2->getMax(),var3_bins,var3->getMin(),var3->getMax());
    TH2F* h2_chi2_5 = new TH2F("h2_chi2_5","h2_chi2_5",var2_bins,var2->getMin(),var2->getMax(),var4_bins,var4->getMin(),var4->getMax());
    TH2F* h2_chi2_6 = new TH2F("h2_chi2_6","h2_chi2_6",var3_bins,var3->getMin(),var3->getMax(),var4_bins,var4->getMin(),var4->getMax());


    /// Cycle through the centers of all bins
    /// I'm getting width of the first bin, because all bins are of equal width
    for(*var1 = var1->getMin()+var1->getBinWidth(0)/2; var1->getVal() < var1->getMax(); var1->setVal(var1->getVal()+var1->getBinWidth(0))) {
        for(*var2 = var2->getMin()+var2->getBinWidth(0)/2; var2->getVal() < var2->getMax(); var2->setVal(var2->getVal()+var2->getBinWidth(0))) {
            for(*var3 = var3->getMin()+var3->getBinWidth(0)/2; var3->getVal() < var3->getMax(); var3->setVal(var3->getVal()+var3->getBinWidth(0))) {
                for(*var4 = var4->getMin()+var4->getBinWidth(0)/2; var4->getVal() < var4->getMax(); var4->setVal(var4->getVal()+var4->getBinWidth(0))) {

                    /// Weight is actually the bin content
                    n = dataSet_binned->weight(RooArgSet(*var1,*var2,*var3,*var4),0);
                    //if(n == 0) continue;
                    if(n == 0) continue;
                    v = pdf_temp->getVal(RooArgSet(*var1,*var2,*var3,*var4))*binVolume*binnedNumEntries;

                    if(((n-v)*(n-v)/v) > 1) {
                        v = GetVPrecise(pdf_temp);
                    }

                    dchi2 = (n-v)*(n-v)/v;

                    h1_chi2->Fill(dchi2);
                    h2_chi2_1->Fill(var1->getVal(),var2->getVal(),dchi2);
                    h2_chi2_2->Fill(var1->getVal(),var3->getVal(),dchi2);
                    h2_chi2_3->Fill(var1->getVal(),var4->getVal(),dchi2);
                    h2_chi2_4->Fill(var2->getVal(),var3->getVal(),dchi2);
                    h2_chi2_5->Fill(var2->getVal(),var4->getVal(),dchi2);
                    h2_chi2_6->Fill(var3->getVal(),var4->getVal(),dchi2);
                    mychi2 += dchi2;

//                    for(int i = 0; i < 4; i++)
//                        h1_resid[i]->Fill(vars[i]->getVal(),(n-v)/v);
                }
            }
        }
    }

    delete dataSet_binned;

    c2->SetLogy(kTRUE);
    h1_chi2->GetXaxis()->SetTitle("dchi");
    h1_chi2->GetYaxis()->SetTitle("num bins");
    h1_chi2->Draw();
    h1_chi2->Write();
    c2->SaveAs("plots/chi2_delta.png");

    c2->SetLogy(kFALSE);

    h2_chi2_1->SetOption("colz");
    //hdchi2->SetMinimum(0);
    //hdchi2->SetMaximum(100);
    h2_chi2_1->SetStats(kFALSE);
    h2_chi2_1->GetXaxis()->SetTitle(var1->GetName());
    h2_chi2_1->GetYaxis()->SetTitle(var2->GetName());
    h2_chi2_1->Draw();
    h2_chi2_1->Write();
    path = "plots/chi2map_";
    path += var1->GetName();
    path += "_";
    path += var2->GetName();
    path += ".png";
    c2->SaveAs(path);

    h2_chi2_2->SetOption("colz");
    //hdchi2->SetMinimum(0);
    //hdchi2->SetMaximum(100);
    h2_chi2_2->SetStats(kFALSE);
    h2_chi2_2->GetXaxis()->SetTitle(var1->GetName());
    h2_chi2_2->GetYaxis()->SetTitle(var3->GetName());
    h2_chi2_2->Draw();
    h2_chi2_2->Write();
    path = "plots/chi2map_";
    path += var1->GetName();
    path += "_";
    path += var3->GetName();
    path += ".png";
    c2->SaveAs(path);

    h2_chi2_3->SetOption("colz");
    //hdchi2->SetMinimum(0);
    //hdchi2->SetMaximum(100);
    h2_chi2_3->SetStats(kFALSE);
    h2_chi2_3->GetXaxis()->SetTitle(var1->GetName());
    h2_chi2_3->GetYaxis()->SetTitle(var4->GetName());
    h2_chi2_3->Draw();
    h2_chi2_3->Write();
    path = "plots/chi2map_";
    path += var1->GetName();
    path += "_";
    path += var4->GetName();
    path += ".png";
    c2->SaveAs(path);

    h2_chi2_4->SetOption("colz");
    //hdchi2->SetMinimum(0);
    //hdchi2->SetMaximum(100);
    h2_chi2_4->SetStats(kFALSE);
    h2_chi2_4->GetXaxis()->SetTitle(var2->GetName());
    h2_chi2_4->GetYaxis()->SetTitle(var3->GetName());
    h2_chi2_4->Draw();
    h2_chi2_4->Write();
    path = "plots/chi2map_";
    path += var2->GetName();
    path += "_";
    path += var3->GetName();
    path += ".png";
    c2->SaveAs(path);

    h2_chi2_5->SetOption("colz");
    //hdchi2->SetMinimum(0);
    //hdchi2->SetMaximum(100);
    h2_chi2_5->SetStats(kFALSE);
    h2_chi2_5->GetXaxis()->SetTitle(var2->GetName());
    h2_chi2_5->GetYaxis()->SetTitle(var4->GetName());
    h2_chi2_5->Draw();
    h2_chi2_5->Write();
    path = "plots/chi2map_";
    path += var2->GetName();
    path += "_";
    path += var4->GetName();
    path += ".png";
    c2->SaveAs(path);

    h2_chi2_6->SetOption("colz");
    //h2_chi2_6->SetMinimum(0);
    //h2_chi2_6->SetMaximum(100000);
    h2_chi2_6->SetStats(kFALSE);
    h2_chi2_6->GetXaxis()->SetTitle(var3->GetName());
    h2_chi2_6->GetYaxis()->SetTitle(var4->GetName());
    h2_chi2_6->Draw();
    h2_chi2_6->Write();
    path = "plots/chi2map_";
    path += var3->GetName();
    path += "_";
    path += var4->GetName();
    path += ".png";
    c2->SaveAs(path);

    file->Close();

//    delete h1_chi2;
//    delete h2_chi2_1;
//    delete h2_chi2_2;
//    delete h2_chi2_3;
//    delete h2_chi2_4;
//    delete h2_chi2_5;
//    delete h2_chi2_6;
//    delete c2;

    return mychi2;
}

void FitterTrans::SaveNllPlot(RooRealVar* var)
{
    const int steps = 100;
    Double_t stepsize = (var->getMax()-var->getMin())/steps;
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
    nll = pdf_sim->createNLL(*dataSet,RooFit::NumCPU(4));

    for(Int_t i = 0; i < steps; i++)
    {
        var->setVal(i*stepsize + var->getMin() + stepsize/2);
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

void FitterTrans::SaveNllPlot(RooRealVar* var1, RooRealVar* var2)
{
    const int steps1 = 30;
    const int steps2 = 30;
    Double_t stepsize1 = (var1->getMax()-var1->getMin())/steps1;
    Double_t stepsize2 = (var2->getMax()-var2->getMin())/steps2;
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
    nll = pdf_sim->createNLL(*dataSet,RooFit::NumCPU(2));
    for(Int_t i = 0; i < steps1; i++)
    {
        for(Int_t j = 0; j < steps2; j++)
        {
            var1->setVal(i*stepsize1 + var1->getMin() + stepsize1/2);
            var2->setVal(j*stepsize2 + var2->getMin() + stepsize2/2);
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

void FitterTrans::SaveParameters(char* file)
{
    const RooArgSet* args = dataSet->get();
    RooPlot* frame = 0;

    const Int_t numParameters = time_dependent ? 39 : 18;
    Double_t* parameters = new Double_t[numParameters];

    // if (time_dependent) {
    //     /// The next 2 lines enable getting category items' names and therefore reduced datasets in a loop
    //     const RooCategory* cat = (RooCategory*)args->find("decType");
    //     RooDataSet* datacut;
    //     const Int_t org_pdf_type = pdf_a->getType();

    //     /// Getting 1D chi^2 for all 4 decay types
    //     for(int i = 1; i <= 4; i++)
    //     {
    //         frame = dt->frame();
    //         TString type = (char*)cat->lookupType(i)->GetName();
    //         TString cut = "decType==decType::" + type;
    //         datacut = (RooDataSet*)dataSet->reduce(*dt,cut);
    //         datacut->plotOn(frame,RooFit::Name("data"));

    //         pdf_a->setType(i);
    //         pdf_a->plotOn(frame,RooFit::Project(RooArgSet(*tht,*thb,*phit)));

    //         parameters[i-1] = frame->chiSquare(11);

    //         delete frame;
    //     }

    //     pdf_a->setType(org_pdf_type);
    // }

    FILE* pFile;
    pFile = fopen (file,"w");
    if (pFile == NULL)
    {
        printf("ERROR: couldn't open file %s for writing!\n",file);
        delete[] parameters;
        return;
    }

    parameters[0] = par_input[0];
    parameters[1] = ap->getVal();
    parameters[2] = ap->getError();
    parameters[3] = par_input[1];
    parameters[4] = apa->getVal();
    parameters[5] = apa->getError();
    parameters[6] = par_input[2];
    parameters[7] = a0->getVal();
    parameters[8] = a0->getError();
    parameters[9] = 0;
    parameters[10] = a0a->getVal();
    parameters[11] = a0a->getError();
    parameters[12] = sqrt(1-par_input[0]*par_input[0]-par_input[2]*par_input[2]);
    parameters[13] = at->getVal();
    if (result == 0){
        parameters[14] = 0;
    }else{
        parameters[14] = at->getPropagatedError(*result);
    }
    parameters[15] = par_input[3];
    parameters[16] = ata->getVal();
    parameters[17] = ata->getError();

    if (time_dependent) {
        parameters[18] = par_input[4];
        parameters[19] = phiw->getVal();
        parameters[20] = phiw->getError();
        parameters[21] = par_input[5];
        parameters[22] = rp->getVal();
        parameters[23] = rp->getError();
        parameters[24] = par_input[6];
        parameters[25] = r0->getVal();
        parameters[26] = r0->getError();
        parameters[27] = par_input[7];
        parameters[28] = rt->getVal();
        parameters[29] = rt->getError();
        parameters[30] = par_input[8];
        /// To prevent the fitter from becoming stuck at the limit of strong phase
        /// range, the range has been extended from the [-PI,PI] interval.
        /// This is kosher as the strong phase is of course 2*PI periodic.
        /// The following code ensures the final result is in the [-PI,PI] interval.
        if(sp->getVal() > PI){
            sp->setVal(sp->getVal()-2*PI);
        }else if(sp->getVal() < -PI){
            sp->setVal(sp->getVal()+2*PI);
        }
        if(s0->getVal() > PI){
            s0->setVal(s0->getVal()-2*PI);
        }else if(s0->getVal() < -PI){
            s0->setVal(s0->getVal()+2*PI);
        }
        if(st->getVal() > PI){
            st->setVal(st->getVal()-2*PI);
        }else if(st->getVal() < -PI){
            st->setVal(st->getVal()+2*PI);
        }
        parameters[31] = sp->getVal();
        parameters[32] = sp->getError();
        parameters[33] = par_input[9];
        parameters[34] = s0->getVal();
        parameters[35] = s0->getError();
        parameters[36] = par_input[10];
        parameters[37] = st->getVal();
        parameters[38] = st->getError();
    }

    Int_t separators[] = {0,0,1, 0,0,1, 0,0,1, 0,0,1, 0,0,1, 0,0,2, 0,0,2, 0,0,1, 0,0,1, 0,0,2, 0,0,1, 0,0,1, 0,0,2};

    for(Int_t i = 0; i < numParameters; i++)
    {
        /// These parameters can be both positive and negative therefore the
        /// added + in front of positive numbers keeps the columns aligned
        if(i==30||i==31||i==33||i==34||i==36||i==37){
            fprintf(pFile,"%+.5f ",parameters[i]);
        } else {
            fprintf(pFile,"%.5f ",parameters[i]);
        }
        if (i == numParameters - 1) continue;
        if (separators[i] == 1)
            fprintf(pFile,"| ");
        else if (separators[i] == 2)
            fprintf(pFile,"|| ");
    }
    fprintf(pFile,"\n");
    fclose (pFile);
    delete[] parameters;
}

RooDataHist* FitterTrans::GetBinnedDataSet()
{
    if(dataSet_binned == NULL)
        CreateBinnedDataSet("a");

    return dataSet_binned;
}

RooDataSet* FitterTrans::GetReducedDataSet()
{
    if(dataSet_reduced == NULL)
        CreateBinnedDataSet("a");

    return dataSet_reduced;
}

void FitterTrans::FixAllParameters()
{
    RooRealVar* rooPar = 0;
    TIterator* parIter = parameters->createIterator();
    while(rooPar = (RooRealVar*)parIter->Next())
        rooPar->setConstant();

}

void FitterTrans::FixParameter(const char* par)
{
    RooRealVar* rooPar = 0;
    rooPar = (RooRealVar*)parameters->find(par);
    if(rooPar != 0)
        rooPar->setConstant();
}

void FitterTrans::FreeParameter(const char* par)
{
    RooRealVar* rooPar = 0;
    rooPar = (RooRealVar*)parameters->find(par);
    if(rooPar != 0)
        rooPar->setConstant(kFALSE);
}

void FitterTrans::PrintParameter(const char* par)
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

void FitterTrans::GetHelParameters(Double_t* params)
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

void FitterTrans::SaveNllPlot(const char* par)
{
    RooRealVar* rooPar = 0;
    rooPar = (RooRealVar*)parameters->find(par);
    if(rooPar != 0)
        SaveNllPlot(rooPar);
    else
        printf("ERROR: Parameter '%s' doesn't exist! Can't save Nll plot.\n",par);
}

void FitterTrans::SaveNllPlot(const char* par1, const char* par2)
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






