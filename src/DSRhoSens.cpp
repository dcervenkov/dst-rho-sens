#include <stdio.h>
#include <vector>
#include <stdlib.h>

#include "TApplication.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TTree.h"
#include "TRandom1.h"
#include "TH1F.h"
#include "TRotation.h"
#include "TLorentzRotation.h"
#include "TLorentzVector.h"
#include "TVector3.h"

#include "DSRhoSens.h"
#include "Particle.h"
#include "Constants.h"
#include "ASSERT.h"

//#define VERBOSE
#define GRAPHIC

TH1D *hChi = new TH1D("hChi", "Chi distribution", 100, 0, 4*PI);
TH1D *hThA = new TH1D("hThA", "Theta A distribution", 100, 0, PI);
TH1D *hThB = new TH1D("hThB", "Theta B distribution", 100, 0, PI);

int main(int argc, char* argv[])
{
    #ifdef GRAPHIC
    TApplication* rootapp = new TApplication("DSRhoSens",&argc,argv); //For graphic-output apps only
    /// When using TApplication it removes arguments it "handles" from
    /// the argument array. E.g. -b, -x, -q, --help, <dir>, <file>, etc.
    /// For more info read TApplication's GetOptions function help.
    /// The solution is to use rootapp->Argc() and rootapp->Argv(i).
    /// The next few lines are for compatibility of GRAPHIC vs. non-GRAPHIC.
    argc = rootapp->Argc();
    for(int i = 0; i < argc; i++)
    {
        argv[i] = rootapp->Argv(i);
    }
    #endif


    for(int i = 1; i < argc; i++)
    {
        ReadEvents(argv[i]);
        Analyze();
    }

    #ifdef GRAPHIC
    TCanvas* c1 = new TCanvas("c1","Canvas",800,600);
    c1->Divide(2,2);
    c1->cd(1);
    hChi->Draw();
    c1->cd(2);
    hThA->Draw();
    c1->cd(3);
    hThB->Draw();

    printf("\nProgram execution has finished.\n");
    rootapp->Run(); //For graphic-output apps only
    #endif

    return 0;
}

void Analyze()
{
    printf("# Events = %i\n",events->size());

    for(int i = 0; i < events->size(); i++)
    {
        //if(i != 2043)
        //    continue;
        //printf("\n\nEvent %i # particles = %i\n",i+1,(*events)[i].size());
        //PrintEvent(i);
        Particle* B0 = 0;
        Particle* DS = 0;
        Particle* DSD0 = 0;
        Particle* DSPi = 0;
        Particle* Rho = 0;
        Particle* RhoPi0 = 0;
        Particle* RhoPi = 0;
        double chi = 0;
        double theta_a = 0;
        double theta_b = 0;

        /// I want pointers to the particles returned, therefore I have to
        /// pass the function a pointer address, so the function can write in it
        if(GetRelevantParticles(i,&B0, &DS, &DSD0, &DSPi, &Rho, &RhoPi0, &RhoPi))
        {
            #ifdef VERBOSE
            printf("All relevant particles found\n");
            #endif
            TransformHel(B0,DS,DSD0,DSPi,Rho,RhoPi0,RhoPi);
            GetAngles(DSD0,RhoPi,chi,theta_a,theta_b);
            //printf("chi: %0.4f\tth_a: %0.4f\tth_b: %0.4f\n",chi,theta_a,theta_b);
            hChi->Fill(chi);
            hThA->Fill(theta_a);
            hThB->Fill(theta_b);
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
        /// event and sets a few variables: idhep, mass, etc.
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
                {
                    (*events)[eventNo][i].SetDaughter(j,(*events)[eventNo][da1[i]-1+j]);
                    //if(j > 9)
                        //printf("eventno: %i\n",eventNo);
                }
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

void GetAngles(Particle* a, Particle* b, double& chi, double& theta_a, double& theta_b)
{
    chi = 2*PI + a->GetPhi() - b->GetPhi();
    theta_a = a->GetTheta();
    theta_b = PI - b->GetTheta();
}



