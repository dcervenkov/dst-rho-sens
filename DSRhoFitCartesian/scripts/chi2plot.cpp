#include "Riostream.h"
void chi2plot() {
//  Read data from an ascii file and create a root file with an histogram and an ntuple.
//   see a variant of this macro in basic2.C
//Author: Rene Brun

   ifstream in;
   in.open("../data/org_many/chi/chi");

   Double_t chi2;
   Int_t nlines = 0;
   Float_t low = 15000;
   Float_t high = 17500;
   Float_t ndof = 15995;
   TFile *f = new TFile("chi2.root","RECREATE");
   TH1F *h1 = new TH1F("h1","chi2 distribution with 15995 d.o.f.",50,low,high);
   TH1F *h2 = new TH1F("h2","chi2red distribution with 15995 d.o.f.",50,low/ndof,high/ndof);

   while (1) {
      in >> chi2;
      if (!in.good()) break;
      h1->Fill(chi2);
      h2->Fill(chi2/ndof);
      nlines++;
   }
   printf(" found %d points\n",nlines);

   in.close();

   TCanvas* cc = new TCanvas("cc","Chi2 distributions",1300,600);
   cc->Divide(2);
   cc->cd(1);
   h1->Draw();
   cc->cd(2);
   h2->Draw();

   f->Write();
}
