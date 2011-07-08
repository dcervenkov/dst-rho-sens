#include <stdio.h>
#include <vector>

#include "TApplication.h"
#include "TFile.h"
#include "TTree.h"
#include "TRandom1.h"
#include "TH1F.h"

#include "DSRhoSens.h"
#include "Particle.h"
#include "ASSERT.h"

int CreateEvents(std::vector< std::vector<Particle> > &events, char fileName[])
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
    TFile *file = new TFile(fileName);
    if (file->IsOpen() == kTRUE)
        printf("File %s opened.\n",fileName);
    TTree *tree = (TTree*)file->Get("Tree");

    events.resize(tree->GetEntries());

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

        //printf("TEMP\tID\tIDHEP\tMOTHER\tDA1\tDA2\tP1\tP2\tP3\tE\n");

        /// Resize the vector to its final size, so it doesn't need to
        /// resize itself many times (explained in std::vector documentation)
        events[eventNo].resize(numParticles);

        /// This loop populates the event vector with all particles of the
        /// event and sets a few variables: idhep, mass, etc.
        for (int i = 0; i < numParticles; i++)
        {
            //printf("%i\t%i\t%i\t%i\t%i\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\n",id[i],idhep[i],\
                   mother[i],da1[i],da2[i],p[i][0],p[i][1],p[i][2],p[i][3],p[i][4]);

            Particle part;
            part.SetIdhep(idhep[i]);
            part.SetM(p[i][4]);
            for (int j = 0; j < 4; j++)
            {
                part.SetP(j,p[i][j]);
                part.SetV(j,v[i][j]);
            }
            events[eventNo][i] = part;
        }

        /// This loop creates relations (mother, daughters) between the particles
        /// created in the previous loop
        for (int i = 0; i < numParticles; i++)
        {
            if((mother[i]-1) >= 0)
                events[eventNo][i].SetMother(events[eventNo][mother[i]-1]);

            if(da1[i] != 0)
            {
                for (int j = 0; j <= (da2[i] - da1[i]); j++)
                    events[eventNo][i].SetDaughter(j,events[eventNo][da1[i]-1+j]);;
            }
        }
        //printf("\n\n\n");
    }

    file->Close();
    delete file;
    return 0;
}

int main(int argc, char* argv[])
{
    //TApplication* rootapp = new TApplication("DSRhoSens",&argc,argv); //For graphic-output apps only

    std::vector< std::vector<Particle> >* events = new std::vector< std::vector<Particle> >;
    CreateEvents(*events,"data/3events.root");

    printf("# Events = %i\n",events->size());
    for (int i = 0; i < events->size(); i++)
        printf("Event %i # particles = %i\n",i+1,(*events)[i].size());

    //rootapp->Run(); //For graphic-output apps only
    return 0;
}