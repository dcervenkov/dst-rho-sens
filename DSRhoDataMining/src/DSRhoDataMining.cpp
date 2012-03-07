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
#include "RooCategory.h"
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


#include "DSRhoDataMining.h"
#include "Particle.h"
#include "Constants.h"
#include "ASSERT.h"

#define DEBUG
//#define VERBOSE
//#define HELICITY
#define TRANSVERSITY

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
    if (argc < 3)
    {
        printf("ERROR: Too few arguments!\n");
        printf("Usage: DSRhoDataMining output_file input_file(s)\n");
        return -1;
    }

    TStopwatch timer;
    timer.Start();

    RooRealVar dt("dt","dt",-10,10);
    RooCategory decType("decType","decType");

    #ifdef HELICITY
    RooRealVar tha("tha","tha",0,PI);
    RooRealVar thb("thb","thb",0,PI);
    RooRealVar chi("chi","chi",0,2*PI);
    RooDataSet* dataSet = new RooDataSet("data","data",RooArgList(tha,thb,chi,dt,decType));
    #endif

    #ifdef TRANSVERSITY
    RooRealVar tht("tht","tht",0,PI);
    RooRealVar thb("thb","thb",0,PI);
    RooRealVar phit("phit","phit",-PI,PI);
    RooDataSet* dataSet = new RooDataSet("data","data",RooArgSet(tht,thb,phit,dt,decType));
    #endif

    /// The argv[0] is path and name of the program itself, argv[1] is the output filename
    /// so input files start at argv[2].
    for(int i = 2; i < argc; i++)
    {
        ReadEvents(argv[i]);
        Analyze(dataSet);
    }

    dataSet->write(argv[1]);

    timer.Stop();
    timer.Print();

    return 0;
}


void Analyze(RooDataSet* dataSet)
{
    printf("# Events = %i\n",events->size());

    int not_found = 0;

    for(int i = 0; i < events->size(); i++)
    //for(int i = 4; i < 6; i++)
    {
        //printf("\n\nEvent %i # particles = %i\n",i,(*events)[i].size());
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

        double delta_t = 0;
        int dec_type = 0;

        /// I want pointers to the particles returned, therefore I have to
        /// pass the function a pointer address, so the function can write in it
        if(GetRelevantParticles(i,&B0, &DS, &DSD0, &DSPi, &Rho, &RhoPi0, &RhoPi, delta_t, dec_type))
        {
            #ifdef VERBOSE
            printf("All relevant particles found\n");
            #endif

            RooRealVar dt("dt","dt",-10,10);
            dt = delta_t;

            RooCategory decType("decType","decType");
            decType.defineType("a",1);
            decType.defineType("ab",2);
            decType.defineType("b",3);
            decType.defineType("bb",4);
            decType.setIndex(dec_type);

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
            dataSet->add(RooArgSet(tha,thb,chi,dt,decType));
            //printf("chi: %0.4f\tth_a: %0.4f\tth_b: %0.4f\n",chi,theta_a,theta_b);
            #endif

            #ifdef TRANSVERSITY
            TransformTrans(B0,DS,DSD0,DSPi,Rho,RhoPi0,RhoPi);
            GetAnglesTrans(DSD0,RhoPi,theta_t,phi_t,theta_b);
            hPhiT->Fill(phi_t);
            hThT->Fill(theta_t);
            hThB->Fill(theta_b);

            RooRealVar tht("tht","tht",0,PI);
            RooRealVar thb("thb","thb",0,PI);
            RooRealVar phit("phit","phit",-PI,PI);
            tht = theta_t;
            thb = theta_b;
            phit = phi_t;
            dataSet->add(RooArgSet(tht,thb,phit,dt,decType));
            //printf("th_t: %0.4f\tphi_t: %0.4f\tth_b: %0.4f\n",theta_t,phi_t,theta_b);
            #endif

        }
        else
        {
            //printf("WARNING: Could not find all relevant particles in event #%i\n",i+1);
            not_found++;
        }
    }
    printf("WARNING: Could not find all relevant particles in %i events.\n",not_found);
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
                          Particle** Rho, Particle** RhoPi0, Particle** RhoPi, double& delta_t, int& dec_type)
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

    Particle* tmpB0 = 0;
    Particle* tmpDS = 0;
    Particle* tmpRho = 0;
    Particle* tmpDSD0 = 0;
    Particle* tmpDSPi = 0;
    Particle* tmpRhoPi0 = 0;
    Particle* tmpRhoPi = 0;

    int found_sig = 0;
    int found_tag = 0;
    int B0s = 0;
    int B0Bars = 0;

    double t_tag = 0;
    double t_sig = 0;

    /// This cycles through all particles in an event
    for(int i = 0; i < (*events)[eventNo].size(); i++)
    {
        bool found_sig_this_iter = 0;

        if(abs((*events)[eventNo][i].GetIdhep()) == B0_IDHEP)
        {
            #ifdef VERBOSE
            printf("Found a B0\n");
            #endif

            tmpB0 = &(*events)[eventNo][i];

            if(tmpB0->GetIdhep() > 0)
                B0s++;
            else
                B0Bars++;

            if(GetDSRhoFromB0(tmpB0,&tmpDS,&tmpRho))
            {
                if(GetD0PiFromDS(tmpDS,&tmpDSD0,&tmpDSPi) && GetPi0PiFromRho(tmpRho,&tmpRhoPi0,&tmpRhoPi))
                {
                    t_sig = tmpDS->GetV(3);
                    found_sig += 1;
                    found_sig_this_iter = 1;

                    *B0 = tmpB0;
                    *DS = tmpDS;
                    *Rho = tmpRho;
                    *DSD0 = tmpDSD0;
                    *DSPi = tmpDSPi;
                    *RhoPi0 = tmpRhoPi0;
                    *RhoPi = tmpRhoPi;

                }
            }

            if(!found_sig_this_iter)
            {
                /// This is a safeguard against decays, where B0tag doesn't decay for
                /// some reason
                if(tmpB0->GetNumDaughters() == 0)
                    return 0;

                t_tag = tmpB0->GetDaughter(0)->GetV(3);
                found_tag += 1;

                #ifdef VERBOSE
                printf("This branch doesn't match\n");
                #endif
            }
        } // If block searching for B0



        /// What should happen if 2 signal events are found? Now they are discarded
        if((found_sig == 1) && (found_tag == 1))
        {
            //printf("Now printing in foundsig && foundtag\n");
            //PrintRelevantParticles(*DS, *DSD0, *DSPi, *Rho, *RhoPi0, *RhoPi);

            if((B0s == 1) && (B0Bars == 1))
            {
                if((*DS)->GetIdhep() > 0)
                    dec_type = 2;
                else
                    dec_type = 1;
            }
            /// EvtGen simulates mixing by decaying Ups(4S) into two
            /// mesons of the same kind with the correct time distribution.
            else if((B0s == 2) && (B0Bars == 0))
                dec_type = 4;
            else if((B0s == 0) && (B0Bars == 2))
                dec_type = 3;

            /// dec_type:   1:  B0      -> D*- + rho+   (a  - favored)
            ///             2:  B0Bar   -> D*+ + rho-   (ab - favored)
            ///             3:  B0      -> D*+ + rho-   (b  - suppressed)
            ///             4:  B0Bar   -> D*- + rho+   (bb - suppressed)

            delta_t = t_sig - t_tag;
            return 1;
        }
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
