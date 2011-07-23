#include <stdio.h>
#include <vector>
#include <stdlib.h>

#include "TApplication.h"
#include "TFile.h"
#include "TTree.h"
#include "TRandom1.h"
#include "TH1F.h"

#include "DSRhoSens.h"
#include "Particle.h"
#include "Constants.h"
#include "ASSERT.h"

int main(int argc, char* argv[])
{
    //TApplication* rootapp = new TApplication("DSRhoSens",&argc,argv); //For graphic-output apps only

    ReadEvents("data/DSRho_exp07-0.root");
    printf("# Events = %i\n",events->size());
    for(int i = 0; i < 3; i++)
    {
        printf("\n\nEvent %i # particles = %i\n",i+1,(*events)[i].size());
        PrintEvent(i);
        Particle *DSD0, *DSPi, *RhoPi0, *RhoPi;
        /// I want pointers to the particles returned, therefore I have to
        /// pass the function a pointer address, so the function can write in it
        if(GetRelevantParticles(i, &DSD0, &DSPi, &RhoPi0, &RhoPi))
        {
            printf("All relevant particles found\n");
            //printf("D0 IDHEP: %i\tE: %f GeV\n", DSD0->GetIdhep(),DSD0->GetP(0));
            //printf("Pi0 IDHEP: %i\tE: %f GeV\n", RhoPi0->GetIdhep(),RhoPi0->GetP(0));
        }

    }

    //rootapp->Run(); //For graphic-output apps only
    return 0;
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

        //printf("\n\nID\tIDHEP\tMOTHER\tDA1\tDA2\tP1\tP2\tP3\tE\tM\n");

        /// Resize the vector to its final size, so it doesn't need to
        /// resize itself many times (explained in std::vector documentation)
        (*events)[eventNo].resize(numParticles);

        /// This loop populates the event vector with all particles of the
        /// event and sets a few variables: idhep, mass, etc.
        for(int i = 0; i < numParticles; i++)
        {
            //printf("%i\t%i\t%i\t%i\t%i\t%.2f\t%.2f\t%.2f\t%.2f\t%.3f\n",id[i],idhep[i],\
            mother[i],da1[i],da2[i],p[i][0],p[i][1],p[i][2],p[i][3],p[i][4]);

            Particle part;
            part.SetIdhep(idhep[i]);
            part.SetM(p[i][4]);

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
                    (*events)[eventNo][i].SetDaughter(j,(*events)[eventNo][da1[i]-1+j]);;
            }
        }
    }

    file->Close();
    delete file;
    return 0;
}

void PrintEvent(int evtNo)
{
    std::vector<Particle> particles = (*events)[evtNo];
    printf("ID\tIDHEP\tP1\tP2\tP3\tE\tM\n");
    for(int i = 0; i < particles.size(); i++)
    {
        printf("%i\t%i\t%.2f\t%.2f\t%.2f\t%.2f\t%.3f\n",i+1,particles[i].GetIdhep(),\
               particles[i].GetP(0),particles[i].GetP(1),particles[i].GetP(2),particles[i].GetP(3),particles[i].GetM());
    }
    printf("\n");
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
bool GetRelevantParticles(int eventNo, Particle** DSD0, Particle** DSPi, Particle** RhoPi0, Particle** RhoPi)
{
    /// Setting pointers to null address, so a check that all particles
    /// were found and pointed to can be done at the end of this function
    Particle* B0 = 0;
    Particle* DS = 0;
    Particle* Rho = 0;

    /// This cycles through all particles in an event
    for(int i = 0; i < (*events)[eventNo].size(); i++)
    {
        /// All pointers are initialized to null for safety
        B0 = 0;
        DS = 0;
        Rho = 0;
        /// You don't want to have e.g. DSD0 = 0, because that would set
        /// the POINTER TO A POINTER that points to a particle to null.
        /// What you want to do is set the POINTER TO A PARTICLE to null,
        /// because you don't yet know which particle it should point to.
        /// Therefore you must use one dereference.
        *DSD0 = 0;
        *DSPi = 0;
        *RhoPi0 = 0;
        *RhoPi = 0;

        if(abs((*events)[eventNo][i].GetIdhep()) == B0_IDHEP)
        {
            printf("Found a B0\n");
            B0 = &(*events)[eventNo][i];
            if(GetDSRhoFromB0(B0,&DS,&Rho))
            {
                if(GetD0PiFromDS(DS,DSD0,DSPi) && GetPi0PiFromRho(Rho,RhoPi0,RhoPi))
                    return 1;
            }
            printf("This branch doesn't match\n");
        } // If block searching for B0
    } // For cycle going through all particles in an event

    return 0;
}

bool GetDSRhoFromB0(Particle* B0, Particle** DS, Particle** Rho)
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
            printf("Found a D*\n");
            *DS = B0->GetDaughter(i);
            foundDS = 1;
            break;

        case RHO_IDHEP:
            printf("Found a Rho\n");
            *Rho = B0->GetDaughter(i);
            foundRho = 1;
            break;

        case PHOTON_IDHEP:
            printf("Found a gamma\n");
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

bool GetD0PiFromDS(Particle* DS, Particle** D0, Particle** Pi)
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
            printf("Found a D0 from D*\n");
            *D0 = DS->GetDaughter(i);
            foundD0 = 1;
            break;

        case PI_IDHEP:
            printf("Found a Pi from D*\n");
            *Pi = DS->GetDaughter(i);;
            foundPi = 1;
            break;

        case PHOTON_IDHEP:
            printf("Found a gamma from D*\n");
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

bool GetPi0PiFromRho(Particle* Rho, Particle** Pi0, Particle** Pi)
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
            printf("Found a Pi0 from Rho\n");
            *Pi0 = Rho->GetDaughter(i);
            foundPi0 = 1;
            break;

        case PI_IDHEP:
            printf("Found a Pi from Rho\n");
            *Pi = Rho->GetDaughter(i);
            foundPi = 1;
            break;

        case PHOTON_IDHEP:
            printf("Found a gamma from Rho\n");
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

