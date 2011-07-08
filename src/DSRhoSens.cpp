#include "DSRhoSens.h"

int treer()
{
    printf("Opening file %s...\n",fileName);
    TFile *file = new TFile(fileName);
    if (file->IsOpen() == kTRUE)
        printf("File %s opened.\n",fileName);
    TTree *tree = (TTree*)file->Get("Tree");
    for(int eventNo = 0; eventNo < tree->GetEntries(); eventNo++)
    {
        //This code block extracts the number of particles in the event
        //so arrays of appropriate size can be created
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

        printf("ID\tIDHEP\tMOTHER\tDA1\tDA2\tP1\tP2\tP3\tE\tM\n");

        for (int i = 0; i < numParticles; i++)
        {
            printf("%i\t%i\t%i\t%i\t%i\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\n",id[i],idhep[i],\
                   mother[i],da1[i],da2[i],p[i][0],p[i][1],p[i][2],p[i][3],p[i][4]);
        }
        printf("\n\n\n");
    }

    file->Close();
    delete file;
    return 0;
}

int main(int argc, char* argv[])
{
    //TApplication* rootapp = new TApplication("DSRhoSens",&argc,argv); //For graphic-output apps only
    treer();
    //rootapp->Run(); //For graphic-output apps only
    return 0;
}