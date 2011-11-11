#include <stdio.h>
#include <vector>
#include <stdlib.h>

#include "TROOT.h"
#include "TStyle.h"
#include "TApplication.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TNtuple.h"
#include "TTree.h"
#include "TRandom1.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TRotation.h"
#include "TLorentzRotation.h"
#include "TLorentzVector.h"
#include "TVector3.h"
#include "TStopwatch.h"
#include "TMath.h"

#include "RooFit.h"
#include "RooRealVar.h"
#include "RooPlot.h"
#include "RooDataSet.h"
#include "RooDataHist.h"
#include "RooArgSet.h"
#include "RooGenericPdf.h"
#include "RooFitResult.h"
#include "RooChi2Var.h"
#include "RooGaussian.h"
#include "RooProdPdf.h"
#include "RooNLLVar.h"
#include "RooMinuit.h"
#include "RooPolynomial.h"
#include "RooAddPdf.h"


#include "DSRhoSens.h"
#include "Particle.h"
#include "Constants.h"
#include "ASSERT.h"

#define DEBUG
//#define VERBOSE
//#define GRAPHIC
#define HELICITY
//#define TRANSVERSITY

#ifdef HELICITY
TH1D *hChi = new TH1D("hChi", "Chi distribution", 50, 0, 2*PI);
TH1D *hThA = new TH1D("hThA", "Theta A distribution", 50, 0, PI);
TH1D *hThB = new TH1D("hThB", "Theta B distribution", 50, 0, PI);
TH2D *hG = new TH2D("hG","Gamma distribution", 30, 0, PI, 30, 0, PI);
#endif

#ifdef TRANSVERSITY
TH1D *hPhiT = new TH1D("hPhiT", "Phi T distribution", 50, -PI, PI);
TH1D *hThT = new TH1D("hThT", "Theta T distribution", 50, 0, PI);
TH1D *hThB = new TH1D("hThB", "Theta B distribution", 50, 0, PI);
#endif


int main(int argc, char* argv[])
{
    #ifdef GRAPHIC
    TApplication* rootapp = new TApplication("DSRhoSens",&argc,argv); //For graphic-output apps only
    /**
     * When using TApplication it removes arguments it "handles" from
     * the argument array. E.g. -b, -x, -q, --help, <dir>, <file>, etc.
     * For more info read TApplication's GetOptions function help.
     * The solution is to use rootapp->Argc() and rootapp->Argv(i).
     * The next few lines are for compatibility of GRAPHIC vs. non-GRAPHIC -
     * they recreate the original argc and argv when using GRAPHIC.
     **/
    argc = rootapp->Argc();
    for(int i = 0; i < argc; i++)
    {
        argv[i] = rootapp->Argv(i);
    }
    #endif

    TStopwatch timer;
    timer.Start();

    #ifdef HELICITY
    RooRealVar tha("tha","tha",0,PI);
    RooRealVar thb("thb","thb",0,PI);
    RooRealVar chi("chi","chi",0,2*PI);

    RooDataSet* dataSet = new RooDataSet("data","data",RooArgList(tha,thb,chi));
    #endif

    #ifdef TRANSVERSITY
    RooRealVar phit("phit","phit",-PI,PI);
    RooRealVar tht("tht","tht",0,PI);
    RooRealVar thb("thb","thb",0,PI);
    RooRealVar ap2("ap2","ap2",0,1);
    RooDataSet* dataSet = new RooDataSet("data","data",RooArgSet(phit,tht,thb));
    RooGenericPdf* pdf = new RooGenericPdf("pdf","Generic PDF","ap2*sin(phit)*sin(phit)*sin(tht)*sin(tht)*sin(tht)*sin(thb)*sin(thb)*sin(thb)",RooArgSet(phit,tht,thb,ap2));
    #endif

    for(int i = 2; i < argc; i++)
    {
        ReadEvents(argv[i]);
        Analyze(dataSet);
    }

    //dataSet->write(argv[1]);
    //dataSet = RooDataSet::read("dataset",RooArgList(tha,thb,chi));
    Fit(dataSet,tha,thb,chi);

    timer.Stop();
    timer.Print();

    #ifdef GRAPHIC
//    TCanvas* cc = new TCanvas("cc","Canvas",800,600);

        #ifdef HELICITY
//    	RooPlot* frame1 = chi.frame();
//        RooPlot* frame2 = tha.frame();
//        RooPlot* frame3 = thb.frame();
//        binnedDataSet->plotOn(frame1);
//        pdf->plotOn(frame1);
//        binnedDataSet->plotOn(frame2);
//        pdf->plotOn(frame2);
//        binnedDataSet->plotOn(frame3);
//        pdf->plotOn(frame3);
//        cc->Divide(2,2);
//        cc->cd(1);
//        frame1->Draw();
//        cc->cd(2);
//    	frame2->Draw();
//    	cc->cd(3);
//    	frame3->Draw();
//    	printf("chi1red = %f\nchi2red = %f\nchi3red = %f\n",frame1->chiSquare(5),frame2->chiSquare(5),frame3->chiSquare(5));
        #endif

        #ifdef TRANSVERSITY
        cc->Divide(2,2);
        cc->cd(1);
        hPhiT->Draw();
        cc->cd(2);
        hThT->Draw();
        cc->cd(3);
        hThB->Draw();
        cc->cd(4);

        TCanvas* c1 = new TCanvas("c1","c1");
        pdf->fitTo(*dataSet);
        RooPlot* xframe = phit.frame();
        dataSet->plotOn(xframe);
        pdf->plotOn(xframe);
        xframe->Draw();

        TCanvas* c2 = new TCanvas("c2","c2");
        RooPlot* yframe = tht.frame();
        dataSet->plotOn(yframe);
        pdf->plotOn(yframe);
        yframe->Draw();

        TCanvas* c3 = new TCanvas("c3","c3");
        RooPlot* zframe = thb.frame();
        dataSet->plotOn(zframe);
        pdf->plotOn(zframe);
        zframe->Draw();
        #endif

    printf("\nProgram execution has finished.\n");
    rootapp->Run(); //For graphic-output apps only
    #endif
    return 0;
}

void Fit(RooDataSet* dataSet, RooRealVar& tha, RooRealVar& thb, RooRealVar& chi)
{

    /**
     * The f one gets from fitting with the following PDF is actually (|H+|^2 + |H-|^2)
     * viz formula (35) in BN419. Equvivalentely (1-f) = |H0|^2.
     **/
//    RooRealVar f("f","fraction",0.5,0,1);
//    RooRealVar c1("c1","c1 var",0.05,0,0.5);
//    RooRealVar c2("c2","c2 var",0.05,0,0.5);

    //RooRealVar hp("hp","hp",0.152,0,0.3);
    RooRealVar hp("hp","hp",0.152);
    //RooRealVar hpa("hpa","hpa",1.47,0,2*PI);
    RooRealVar hpa("hpa","hpa",1.47);
    RooFormulaVar hpr("hpr","hp*cos(hpa)",RooArgSet(hp,hpa));
    RooFormulaVar hpi("hpi","hp*sin(hpa)",RooArgSet(hp,hpa));
    //RooRealVar h0("h0","h0",0.936,0.8,1);
    RooRealVar h0("h0","h0",0.936);
    //RooRealVar h0a("h0a","h0a",0,0,2*PI);
    RooRealVar h0a("h0a","h0a",0);
    RooFormulaVar h0r("h0r","h0*cos(h0a)",RooArgSet(h0,h0a));
    RooFormulaVar h0i("h0i","h0*sin(h0a)",RooArgSet(h0,h0a));
    RooFormulaVar hm("hm","sqrt(1-hp*hp-h0*h0)",RooArgSet(hp,h0));
    //RooRealVar hma("hma","hma",0.19,0,2*PI);
    RooRealVar hma("hma","hma",0.19);
    RooFormulaVar hmr("hmr","hm*cos(hma)",RooArgSet(hm,hma));
    RooFormulaVar hmi("hmi","hm*sin(hma)",RooArgSet(hm,hma));

    RooFormulaVar hptr("hptr","hp*hm*cos(hpa-hma)",RooArgSet(hp,hpa,hm,hma));
    RooFormulaVar hpti("hpti","hp*hm*sin(hpa-hma)",RooArgSet(hp,hpa,hm,hma));

    /**
     * hat: h - helicity, a - addition, t - transverse
     * hst: s - subtraction
     */
    RooFormulaVar hatr("hatr","hp*cos(hpa) + hm*cos(hma)",RooArgSet(hp,hpa,hm,hma));
    RooFormulaVar hati("hati","hp*sin(hpa) + hm*sin(hma)",RooArgSet(hp,hpa,hm,hma));
    RooFormulaVar hat("hat","sqrt(hatr*hatr + hati*hati)",RooArgSet(hatr,hati));
    RooFormulaVar hata("hata","atan2(hati,hatr)",RooArgSet(hati,hatr));

    RooFormulaVar hstr("hstr","hp*cos(hpa) - hm*cos(hma)",RooArgSet(hp,hpa,hm,hma));
    RooFormulaVar hsti("hsti","hp*sin(hpa) - hm*sin(hma)",RooArgSet(hp,hpa,hm,hma));
    RooFormulaVar hst("hst","sqrt(hstr*hstr + hsti*hsti)",RooArgSet(hstr,hsti));
    RooFormulaVar hsta("hsta","atan2(hsti,hstr)",RooArgSet(hsti,hstr));

    RooArgSet varSet(tha,thb,hp,hpa,hpr,hpi,h0,h0a,h0r);
    varSet.add(h0i);
    varSet.add(hm);
    varSet.add(hma);
    varSet.add(hmr);
    varSet.add(hmi);
    varSet.add(hat);
    varSet.add(hata);
    varSet.add(hst);
    varSet.add(hsta);
    varSet.add(chi);

    const char* pdfFormula =   "(hp*hp+hm*hm)*sin(tha)*sin(tha)*sin(tha)*sin(thb)*sin(thb)*sin(thb)+\
                                4*h0*h0*cos(tha)*cos(tha)*sin(tha)*cos(thb)*cos(thb)*sin(thb)+\
                                2*(hp*hm*cos(hpa-hma)*cos(2*chi)-hp*hm*sin(hpa-hma)*sin(2*chi))*sin(tha)*sin(tha)*sin(tha)*sin(thb)*sin(thb)*sin(thb)+\
                                (hat*h0*cos(hata-h0a)*cos(chi)-hst*h0*sin(hsta-h0a)*sin(chi))*sin(2*tha)*sin(tha)*sin(2*thb)*sin(thb)";

    RooGenericPdf* pdf = new RooGenericPdf("pdf","Generic PDF",pdfFormula,varSet);

    //RooGenericPdf* pdf_int_chi = new RooGenericPdf("pdf_int_chi","pdf_int_chi","f*sin(tha)*sin(tha)*sin(tha)*sin(thb)*sin(thb)*sin(thb)+(1-f)*4*cos(tha)*cos(tha)*sin(tha)*cos(thb)*cos(thb)*sin(thb)",RooArgSet(tha,thb,f));
    //RooGenericPdf* pdf_int_chi_thb = new RooGenericPdf("pdf_int_chi","pdf_int_chi","f*(4/3)*sin(tha)*sin(tha)*sin(tha)+(1-f)*(2/3)*4*cos(tha)*cos(tha)*sin(tha)",RooArgSet(tha,f));
    //RooGenericPdf* pdf_int_tha_thb = new RooGenericPdf("pdf_int_tha_thb","pdf_int_tha_thb","1+2*(c1*cos(2*chi)-c2*sin(2*chi))",RooArgSet(c1,c2,chi));


    //RooDataSet* data = pdf->generate(RooArgSet(tha,thb),10000);

//    RooFitResult* result = pdf_int_chi->fitTo(*dataSet,RooFit::Save(),RooFit::Timer(true));//,RooFit::NumCPU(2));
//    result->Print();
//    RooPlot* frame1 = tha.frame(RooFit::Bins(100));
//    dataSet->plotOn(frame1);
//    pdf_int_chi->plotOn(frame1);
//    pdf_int_chi->paramOn(frame1);
//    frame1->Draw();
//    printf("(int_chi) chi^2 = %f\n",frame1->chiSquare(1));

//    RooFormulaVar h0fit("h0fit","sqrt(1-f)",f);
//    h0 = h0fit;
//    h0fit.Print();
//    h0.Print();
//    h0.setConstant(true);
    //c1 = 0.0138;
    //c2 = 0.0462;
    //c1.setConstant(true);
    //c2.setConstant(true);

//    RooFitResult* result2 = pdf_int_tha_thb->fitTo(*dataSet,RooFit::Save(),RooFit::Timer(true));
//    result2->Print();
//    c1.Print();
//    c2.Print();
//    h0.Print();
    //printf("h0 = %f\n",h0.getVal());
//    c1.setConstant(true);
//    c2.setConstant(true);
//
//	RooPlot* frame2 = chi.frame(RooFit::Bins(100));
//	dataSet->plotOn(frame2);
//	pdf_int_tha_thb->plotOn(frame2);
//	pdf_int_tha_thb->paramOn(frame2);
//	frame2->Draw();
//	printf("(int_tha_thb) chi^2 = %f\n",frame2->chiSquare(2));

	//RooRealVar c1err("c1err","c1err",c1.getError());
	//RooRealVar c2err("c2err","c2err",c2.getError());
	//RooGaussian* c1_constraint = new RooGaussian("c1_constraint","c1_constraint",hptr,c1,c1err);
	//RooGaussian* c2_constraint = new RooGaussian("c2_constraint","c2_constraint",hpti,c2,c2err);

    //RooFitResult* result3 = pdf->fitTo(*dataSet,RooFit::ExternalConstraints(RooArgSet(*c1_constraint,*c2_constraint)),RooFit::Save(),RooFit::Timer(true));
    //RooFitResult* result3 = pdf->fitTo(*dataSet,RooFit::Save(),RooFit::Timer(true));
    //result3->Print();

    hp.Print();
    hpa.Print();
    h0.Print();
    h0a.Print();
    hm.Print();
    hma.Print();
//    c1.Print();
//    c2.Print();
//    hptr.Print();
//    hpti.Print();

    chi.setBins(40);
	tha.setBins(20);
	thb.setBins(20);

	RooDataHist* binnedDataSet = new RooDataHist("binnedDataSet","binnedDataSet",RooArgSet(chi,tha,thb),*dataSet);
	RooChi2Var chi2Var("chi2Var","chi2Var",*pdf,*binnedDataSet);

	RooRealVar* _chi2     = new RooRealVar("chi2","chi^2",0) ;
	RooRealVar* _ndof     = new RooRealVar("ndof","number of degrees of freedom",0) ;
	RooRealVar* _chi2red  = new RooRealVar("chi2red","reduced chi^2",0) ;
	RooRealVar* _prob     = new RooRealVar("prob","prob(chi2,ndof)",0) ;

	_chi2->setVal(chi2Var.getVal()) ;
	//RooArgSet* floatPars = (RooArgSet*) fitParams()->selectByAttrib("Constant",kFALSE);
	// Should use the above line instead of 5
	_ndof->setVal(binnedDataSet->numEntries()-5) ;
	_chi2red->setVal(_chi2->getVal()/_ndof->getVal()) ;
	_prob->setVal(TMath::Prob(_chi2->getVal(),static_cast<int>(_ndof->getVal())));

	printf("chi2 = %f\nndof = %f\nchi2red = %f\nprob = %f\n",_chi2->getVal(),_ndof->getVal(),_chi2red->getVal(),_prob->getVal());

	TFile* chi2File = new TFile("chi2File.root","UPDATE");
	TNtuple* ntuple = (TNtuple*)chi2File->Get("ntuple");
	if (!ntuple)
	{
		printf("ntuple not found, creating it.\n");
		ntuple = new TNtuple("ntuple","chi2 ntuple","chi2");
	}

	ntuple->Fill((double)chi2Var.getVal());
	ntuple->Write("",TObject::kOverwrite);
	chi2File->Close();
	delete chi2File;

}



void Analyze(RooDataSet* dataSet)
{
    printf("# Events = %i\n",events->size());

    for(int i = 0; i < events->size(); i++)
    //for(int i = 0; i < 3; i++)
    {
        //printf("\n\nEvent %i # particles = %i\n",i+1,(*events)[i].size());
        //PrintEvent(i);
        Particle* B0 = 0;
        Particle* DS = 0;
        Particle* DSD0 = 0;
        Particle* DSPi = 0;
        Particle* Rho = 0;
        Particle* RhoPi0 = 0;
        Particle* RhoPi = 0;

        #ifdef HELICITY
        double dChi = 0;
        double theta_a = 0;
        double theta_b = 0;
        #endif

        #ifdef TRANSVERSITY
        double phi_t = 0;
        double theta_t = 0;
        double theta_b = 0;
        #endif

        /// I want pointers to the particles returned, therefore I have to
        /// pass the function a pointer address, so the function can write in it
        if(GetRelevantParticles(i,&B0, &DS, &DSD0, &DSPi, &Rho, &RhoPi0, &RhoPi))
        {
            #ifdef VERBOSE
            printf("All relevant particles found\n");
            #endif

            #ifdef HELICITY
            TransformHel(B0,DS,DSD0,DSPi,Rho,RhoPi0,RhoPi);
            GetAnglesHel(DSD0,RhoPi,dChi,theta_a,theta_b);
            hChi->Fill(dChi);
            hThA->Fill(theta_a);
            hThB->Fill(theta_b);
            hG->Fill(theta_a,theta_b);

            RooRealVar tha("tha","tha",0,PI);
            RooRealVar thb("thb","thb",0,PI);
            RooRealVar chi("chi","chi",0,2*PI);
            tha = theta_a;
            thb = theta_b;
            chi = dChi;
            dataSet->add(RooArgSet(tha,thb,chi));
            //printf("chi: %0.4f\tth_a: %0.4f\tth_b: %0.4f\n",chi,theta_a,theta_b);
            #endif

            #ifdef TRANSVERSITY
            TransformTrans(B0,DS,DSD0,DSPi,Rho,RhoPi0,RhoPi);
            GetAnglesTrans(DSD0,RhoPi,theta_t,phi_t,theta_b);
            hPhiT->Fill(phi_t);
            hThT->Fill(theta_t);
            hThB->Fill(theta_b);

            RooRealVar phit("phit","phit",-PI,PI);
            RooRealVar tht("tht","tht",0,PI);
            RooRealVar thb("thb","thb",0,PI);
            phit = phi_t;
            tht = theta_t;
            thb = theta_b;
            dataSet->add(RooArgSet(phit,tht,thb));
            //printf("th_t: %0.4f\tphi_t: %0.4f\tth_b: %0.4f\n",theta_t,phi_t,theta_b);
            #endif

        }
        else
            printf("WARNING: Could not find all relevant particles in event #%i\n",i+1);
    }
}

int ReadEvents(char fileName[])
{
    /// This function reads a ROOT file and creates an object
    /// with the following structure:
    ///
    /// events - event1 - particle1
    ///                 - particle2
    ///                 -    ...
    ///
    ///        - event2 - particle1
    ///                 - particle2
    ///                 -    ...
    ///
    ///        -  ...

    printf("Opening file %s...\n",fileName);
    TFile* file = new TFile(fileName);
    if(file->IsOpen() == kTRUE)
        printf("File %s opened.\n",fileName);
    TTree* tree = (TTree*)file->Get("Tree");

    events->clear();
    events->resize(tree->GetEntries());

    /// Cycle through every event in the file
    for(int eventNo = 0; eventNo < tree->GetEntries(); eventNo++)
    {
        /// This code block extracts the number of particles in the event
        /// so arrays of appropriate size can be created in the next block
        int numParticles;
        tree->SetBranchStatus("*",0);
        tree->SetBranchStatus("numParticles",1);
        tree->SetBranchAddress("numParticles",&numParticles);
        tree->GetEntry(eventNo);

        int id[numParticles];
        int idhep[numParticles];
        int mother[numParticles];
        int da1[numParticles];
        int da2[numParticles];
        double p[numParticles][5];
        double v[numParticles][4];

        /// Activate readout of all branches except numParticles
        /// which is no longer necessary
        tree->SetBranchStatus("*",1);
        tree->SetBranchStatus("numParticles",0);

        tree->SetBranchAddress("id",id);
        tree->SetBranchAddress("idhep",idhep);
        tree->SetBranchAddress("mother",mother);
        tree->SetBranchAddress("da1",da1);
        tree->SetBranchAddress("da2",da2);
        tree->SetBranchAddress("p",p);
        tree->SetBranchAddress("v",v);

        tree->GetEntry(eventNo);

        //if(eventNo == 2044)
        //    PrintEventOrig(eventNo,numParticles,id,idhep,mother,da1,da2,p,v);

        /// Resize the vector to its final size, so it doesn't need to
        /// resize itself many times (explained in std::vector documentation)
        (*events)[eventNo].resize(numParticles);

        /// This loop populates the event vector with all particles of the
        /// event and sets non-relation variables: idhep, momentum, position
        for(int i = 0; i < numParticles; i++)
        {
            Particle part;
            part.SetIdhep(idhep[i]);

            for(int j = 0; j < 4; j++)
            {
                part.SetP(j,p[i][j]);
                part.SetV(j,v[i][j]);
            }

            (*events)[eventNo][i] = part;
        }

        /// This loop creates relations (mother, daughters) between the particles
        /// created in the previous loop
        for(int i = 0; i < numParticles; i++)
        {
            if((mother[i]-1) >= 0)
                (*events)[eventNo][i].SetMother((*events)[eventNo][mother[i]-1]);

            if(da1[i] != 0)
            {
                for(int j = 0; j <= (da2[i] - da1[i]); j++)
                    (*events)[eventNo][i].SetDaughter(j,(*events)[eventNo][da1[i]-1+j]);
            }
        }
    }

    file->Close();
    delete file;
    file = 0;
    return 0;
}

void PrintEvent(int evtNo)
{
    std::vector<Particle> particles = (*events)[evtNo];
    printf("Event %i\n",evtNo);
    printf("ID\tIDHEP\tP1\tP2\tP3\tE\tM\n");
    for(int i = 0; i < particles.size(); i++)
    {
        printf("%i\t%i\t%.2f\t%.2f\t%.2f\t%.2f\t%.3f\n",i+1,particles[i].GetIdhep(),\
               particles[i].GetP(0),particles[i].GetP(1),particles[i].GetP(2),particles[i].GetP(3),particles[i].GetM());
    }
    printf("\n");
}

void PrintEventOrig(int evtNo, int numParticles, int* id, int* idhep, int* mother, \
                    int* da1, int* da2, double (*p)[5], double (*v)[4])
{
    printf("Event %i\n",evtNo);
    printf("ID\tIDHEP\tMOTHER\tDA1\tDA2\tP1\tP2\tP3\tE\tM\n");
    for(int i = 0; i < numParticles; i++)
    {
        printf("%i\t%i\t%i\t%i\t%i\t%.2f\t%.2f\t%.2f\t%.2f\t%.3f\n",id[i],idhep[i],\
            mother[i],da1[i],da2[i],p[i][0],p[i][1],p[i][2],p[i][3],p[i][4]);
    }
    printf("\n");
}

void PrintRelevantParticles(const Particle* DS,const Particle* DSD0,const Particle* DSPi, \
                            const Particle* Rho,const Particle* RhoPi0,const Particle* RhoPi)
{
    printf("Part\tPx\tPy\tPz\tE\tM\n");
    printf("D*\t%0.2f\t%0.2f\t%0.2f\t%0.2f\t%0.2f\t\n", DS->GetP(0), DS->GetP(1), DS->GetP(2), DS->GetP(3), DS->GetM());
    printf("D0\t%0.2f\t%0.2f\t%0.2f\t%0.2f\t%0.2f\t\n", DSD0->GetP(0), DSD0->GetP(1), DSD0->GetP(2), DSD0->GetP(3), DSD0->GetM());
    printf("Pi\t%0.2f\t%0.2f\t%0.2f\t%0.2f\t%0.2f\t\n\n", DSPi->GetP(0), DSPi->GetP(1), DSPi->GetP(2), DSPi->GetP(3), DSPi->GetM());

    printf("Rho\t%0.2f\t%0.2f\t%0.2f\t%0.2f\t%0.2f\t\n", Rho->GetP(0), Rho->GetP(1), Rho->GetP(2), Rho->GetP(3), Rho->GetM());
    printf("Pi0\t%0.2f\t%0.2f\t%0.2f\t%0.2f\t%0.2f\t\n", RhoPi0->GetP(0), RhoPi0->GetP(1), RhoPi0->GetP(2), RhoPi0->GetP(3), RhoPi0->GetM());
    printf("Pi\t%0.2f\t%0.2f\t%0.2f\t%0.2f\t%0.2f\t\n\n\n", RhoPi->GetP(0), RhoPi->GetP(1), RhoPi->GetP(2), RhoPi->GetP(3), RhoPi->GetM());

}

/// When you pass a pointer to a function, a new pointer in the function's
/// scope is created and initialized to point at the same address as the
/// original pointer. Therefore you cannot reassign the original pointer
/// inside the function; you can only change the actual object the
/// original pointer points at. To change the address the pointer points
/// at you must pass the pointer by reference - you pass an address
/// of a pointer, in other words a pointer to a pointer. Therefore, the
/// following function takes pointer to a pointer to a particle (Particle**)
/// as argument(s).
bool GetRelevantParticles(int eventNo, Particle** B0 ,Particle** DS, Particle** DSD0, Particle** DSPi, \
                          Particle** Rho, Particle** RhoPi0, Particle** RhoPi)
{
    /// This cycles through all particles in an event
    for(int i = 0; i < (*events)[eventNo].size(); i++)
    {
        /// You don't want to have e.g. DSD0 = 0, because that would set
        /// the POINTER TO A POINTER that points to a particle to null.
        /// What you want to do is set the POINTER TO A PARTICLE to null,
        /// because you don't yet know which particle it should point to.
        /// Therefore you must use one dereference.
        *B0 = 0;
        *DS = 0;
        *Rho = 0;
        *DSD0 = 0;
        *DSPi = 0;
        *RhoPi0 = 0;
        *RhoPi = 0;

        if(abs((*events)[eventNo][i].GetIdhep()) == B0_IDHEP)
        {
            #ifdef VERBOSE
            printf("Found a B0\n");
            #endif
            *B0 = &(*events)[eventNo][i];
            if(GetDSRhoFromB0(*B0,DS,Rho))
            {
                if(GetD0PiFromDS(*DS,DSD0,DSPi) && GetPi0PiFromRho(*Rho,RhoPi0,RhoPi))
                    return 1;
            }
            #ifdef VERBOSE
            printf("This branch doesn't match\n");
            #endif
        } // If block searching for B0
    } // For cycle going through all particles in an event

    return 0;
}

bool GetDSRhoFromB0(const Particle* const B0, Particle** DS, Particle** Rho)
{
    bool foundDS = 0;
    bool foundRho = 0;
    /// Viz GetRelevantParticles()
    *DS = 0;
    *Rho = 0;

    for(int i = 0; i < B0->GetNumDaughters(); i++)
    {
        switch(abs(B0->GetDaughter(i)->GetIdhep()))
        {
        case DS_IDHEP:
            #ifdef VERBOSE
            printf("Found a D*\n");
            #endif
            *DS = B0->GetDaughter(i);
            foundDS = 1;
            break;

        case RHO_IDHEP:
            #ifdef VERBOSE
            printf("Found a Rho\n");
            #endif
            *Rho = B0->GetDaughter(i);
            foundRho = 1;
            break;

        case PHOTON_IDHEP:
            #ifdef VERBOSE
            printf("Found a gamma\n");
            #endif
            break;

        default:
            return 0;
        }
    }

    if(foundDS && foundRho)
        return 1;
    else
    {
        *DS = 0;
        *Rho = 0;
        return 0;
    }
}

bool GetD0PiFromDS(const Particle* const DS, Particle** D0, Particle** Pi)
{
    bool foundD0 = 0;
    bool foundPi = 0;
    /// GetRelevantParticles()
    *D0 = 0;
    *Pi = 0;

    for(int i = 0; i < DS->GetNumDaughters(); i++)
    {
        switch(abs(DS->GetDaughter(i)->GetIdhep()))
        {
        case D0_IDHEP:
            #ifdef VERBOSE
            printf("Found a D0 from D*\n");
            #endif
            *D0 = DS->GetDaughter(i);
            foundD0 = 1;
            break;

        case PI_IDHEP:
            #ifdef VERBOSE
            printf("Found a Pi from D*\n");
            #endif
            *Pi = DS->GetDaughter(i);;
            foundPi = 1;
            break;

        case PHOTON_IDHEP:
            #ifdef VERBOSE
            printf("Found a gamma from D*\n");
            #endif
            break;

        default:
            return 0;
        }
    }

    if(foundD0 && foundPi)
        return 1;
    else
    {
        *D0 = 0;
        *Pi = 0;
        return 0;
    }
}

bool GetPi0PiFromRho(const Particle* const Rho, Particle** Pi0, Particle** Pi)
{
    bool foundPi0 = 0;
    bool foundPi = 0;
    /// Viz GetRelevantParticles
    *Pi0 = 0;
    *Pi = 0;

    for(int i = 0; i < Rho->GetNumDaughters(); i++)
    {
        switch(abs(Rho->GetDaughter(i)->GetIdhep()))
        {
        case PI0_IDHEP:
            #ifdef VERBOSE
            printf("Found a Pi0 from Rho\n");
            #endif
            *Pi0 = Rho->GetDaughter(i);
            foundPi0 = 1;
            break;

        case PI_IDHEP:
            #ifdef VERBOSE
            printf("Found a Pi from Rho\n");
            #endif
            *Pi = Rho->GetDaughter(i);
            foundPi = 1;
            break;

        case PHOTON_IDHEP:
            #ifdef VERBOSE
            printf("Found a gamma from Rho\n");
            #endif
            break;

        default:
            return 0;
        }
    }

    if(foundPi0 && foundPi)
        return 1;
    else
    {
        *Pi0 = 0;
        *Pi = 0;
        return 0;
    }
}

/// This function returns a rotation needed to transform the momentum of the
/// given particle so that it is wholly in the direction of the Z axis
TRotation GetRotationToZ(const Particle* const part)
{
    TVector3 tMom; // tMom stands for threeMomentum
    TRotation rot;

    tMom.SetXYZ(part->GetP(0),part->GetP(1),part->GetP(2));

    double phi = acos(tMom.Z()/(sqrt((tMom.Y())*(tMom.Y())+(tMom.Z())*(tMom.Z()))));
    if(tMom.Y() < 0)
        rot.RotateX(-phi);
    else
        rot.RotateX(phi);

    tMom = rot*tMom;
    double theta = acos(tMom.Z()/(sqrt((tMom.Z())*(tMom.Z())+(tMom.X())*(tMom.X()))));
    if(tMom.X() < 0)
        rot.RotateY(theta);
    else
        rot.RotateY(-theta);

    return rot;
}

TRotation GetZRotationToX(const Particle* const part)
{
    TVector3 tMom; // tMom stands for threeMomentum
    TRotation rot;

    tMom.SetXYZ(part->GetP(0),part->GetP(1),part->GetP(2));

    double phi = acos(tMom.X()/(sqrt((tMom.X())*(tMom.X())+(tMom.Y())*(tMom.Y()))));
    if(tMom.Y() < 0)
        rot.RotateZ(phi);
    else
        rot.RotateZ(-phi);

    return rot;
}

void TransformTrans(Particle* B0, Particle* DS, Particle* DSD0, Particle* DSPi, \
                  Particle* Rho, Particle* RhoPi0, Particle* RhoPi)
{
    //PrintRelevantParticles(DS, DSD0, DSPi, Rho, RhoPi0, RhoPi);

    TVector3 B0Beta = B0->GetBoost();
    DS->BoostP(-B0Beta);
    DSD0->BoostP(-B0Beta);
    DSPi->BoostP(-B0Beta);
    Rho->BoostP(-B0Beta);
    RhoPi->BoostP(-B0Beta);
    RhoPi0->BoostP(-B0Beta);

    //PrintRelevantParticles(DS, DSD0, DSPi, Rho, RhoPi0, RhoPi);

    TRotation rot = GetRotationToZ(DS);
    DS->RotateP(rot);
    DSD0->RotateP(rot);
    DSPi->RotateP(rot);
    Rho->RotateP(rot);
    RhoPi->RotateP(rot);
    RhoPi0->RotateP(rot);

    //PrintRelevantParticles(DS, DSD0, DSPi, Rho, RhoPi0, RhoPi);

    TRotation rot2 = GetZRotationToX(RhoPi);
    DS->RotateP(rot2);
    DSD0->RotateP(rot2);
    DSPi->RotateP(rot2);
    Rho->RotateP(rot2);
    RhoPi->RotateP(rot2);
    RhoPi0->RotateP(rot2);

    //PrintRelevantParticles(DS, DSD0, DSPi, Rho, RhoPi0, RhoPi);

    TRotation rot3;
    rot3.RotateX(PI/2);
    rot3.RotateZ(PI/2);
    DS->RotateP(rot3);
    DSD0->RotateP(rot3);
    DSPi->RotateP(rot3);
    Rho->RotateP(rot3);
    RhoPi->RotateP(rot3);
    RhoPi0->RotateP(rot3);

    //PrintRelevantParticles(DS, DSD0, DSPi, Rho, RhoPi0, RhoPi);

    TVector3 DSBeta = DS->GetBoost();
    TVector3 RhoBeta = Rho->GetBoost();
    DS->BoostP(-DSBeta);
    DSD0->BoostP(-DSBeta);
    DSPi->BoostP(-DSBeta);
    Rho->BoostP(-RhoBeta);
    RhoPi->BoostP(-RhoBeta);
    RhoPi0->BoostP(-RhoBeta);

    //PrintRelevantParticles(DS, DSD0, DSPi, Rho, RhoPi0, RhoPi);
}

void TransformHel(Particle* B0, Particle* DS, Particle* DSD0, Particle* DSPi, \
                  Particle* Rho, Particle* RhoPi0, Particle* RhoPi)
{
    //PrintRelevantParticles(DS, DSD0, DSPi, Rho, RhoPi0, RhoPi);

    TVector3 B0Beta = B0->GetBoost();
    DS->BoostP(-B0Beta);
    DSD0->BoostP(-B0Beta);
    DSPi->BoostP(-B0Beta);
    Rho->BoostP(-B0Beta);
    RhoPi->BoostP(-B0Beta);
    RhoPi0->BoostP(-B0Beta);

    //PrintRelevantParticles(DS, DSD0, DSPi, Rho, RhoPi0, RhoPi);

    TRotation rot = GetRotationToZ(DS);
    DS->RotateP(rot);
    DSD0->RotateP(rot);
    DSPi->RotateP(rot);
    Rho->RotateP(rot);
    RhoPi->RotateP(rot);
    RhoPi0->RotateP(rot);

    //PrintRelevantParticles(DS, DSD0, DSPi, Rho, RhoPi0, RhoPi);

    TVector3 DSBeta = DS->GetBoost();
    TVector3 RhoBeta = Rho->GetBoost();
    DS->BoostP(-DSBeta);
    DSD0->BoostP(-DSBeta);
    DSPi->BoostP(-DSBeta);
    Rho->BoostP(-RhoBeta);
    RhoPi->BoostP(-RhoBeta);
    RhoPi0->BoostP(-RhoBeta);

    //PrintRelevantParticles(DS, DSD0, DSPi, Rho, RhoPi0, RhoPi);
}

void GetAnglesHel(Particle* a, Particle* b, double& chi, double& theta_a, double& theta_b)
{
    double a_phi = a->GetPhi();
    double b_phi = b->GetPhi();

    /// GetPhi() returns phi from -PI to PI. The formula for chi works only when phi is
    /// from 0 to 2PI and that's the reason for these two ifs.
    if(a_phi < 0)
        a_phi = 2*PI + a_phi;

    if(b_phi < 0)
        b_phi = 2*PI + b_phi;

    if(a_phi >= b_phi)
        chi = a_phi - b_phi;
    else
        chi = 2*PI - (b_phi - a_phi);

    theta_a = a->GetTheta();
    theta_b = PI - b->GetTheta();
}

void GetAnglesTrans(Particle* a, Particle* b, double& theta_t, double& phi_t, double& theta_b)
{
    theta_t = a->GetTheta();
    phi_t = a->GetPhi();
    theta_b = PI - b->GetPhi();
}
