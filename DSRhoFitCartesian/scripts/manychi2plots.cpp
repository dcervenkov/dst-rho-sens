#include "Riostream.h"
void chi2plot() {
//  Read data from an ascii file and create a root file with an histogram and an ntuple.
//   see a variant of this macro in basic2.C
//Author: Rene Brun

   ifstream in;
   in.open("../data/org_many/chi/chi");

   Double_t chi2;
   Double_t chi2tha;
   Double_t chi2thb;
   Double_t chi2chi;
   Int_t nlines = 0;
   Float_t low = 15000;
   Float_t high = 19000;
   Float_t ndof = 15995;
   TFile *f = new TFile("chi2.root","RECREATE");
   //TH1F *h1 = new TH1F("h1","chi2 distribution with 1995 d.o.f.",50,low,high);
   TH1F *h2 = new TH1F("h2","chi2red distribution with 15995 d.o.f.",50,low/ndof,high/ndof);
   TH1F *h3 = new TH1F("h3","chi2tha distribution with 15 d.o.f.",50,0,3);
   TH1F *h4 = new TH1F("h4","chi2thb distribution with 15 d.o.f.",50,0,3);
   TH1F *h5 = new TH1F("h5","chi2chi distribution with 15 d.o.f.",50,0,3);

   while (1) {
      in >> chi2 >> chi2tha >> chi2thb >> chi2chi;
      if (!in.good()) break;
      //h1->Fill(chi2);
      h2->Fill(chi2/ndof);
      h3->Fill(chi2tha);
      h4->Fill(chi2thb);
      h5->Fill(chi2chi);
      nlines++;
   }
   printf(" found %d points\n",nlines);

   in.close();

   TCanvas* cc = new TCanvas("cc","Chi2 distributions",1300,600);
   cc->Divide(2,2);
   cc->cd(1);
   h2->Draw();
   cc->cd(2);
   h3->Draw();
   cc->cd(3);
   h4->Draw();
   cc->cd(4);
   h5->Draw();

   f->Write();
}
