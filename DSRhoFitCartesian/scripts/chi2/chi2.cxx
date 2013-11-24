void draw_chi2(TF1* chi2){
    const Int_t ndof = 89;
    TH1F *htemp = (TH1F*)gPad->GetPrimitive("htemp");
    Double_t range_min = htemp->GetXaxis()->GetXmin();
    Double_t range_max = htemp->GetXaxis()->GetXmax();
    Double_t bin_width = htemp->GetXaxis()->GetBinWidth(1);
    Double_t num_events = htemp->GetEntries();
    chi2->SetParameters(ndof,1);
    Double_t integral = chi2->Integral(range_min,range_max);
    //printf("range_min = %f, range_max = %f, int = %f, bin_width = %f, num_events = %f\n",range_min,range_max,integral,bin_width,num_events);
    chi2->SetParameters(ndof,num_events*bin_width/integral);
    chi2->Draw("same");
}

void chi2(const bool drawchi2) {
    TFile *f = new TFile("chi2.root","RECREATE");
    TTree *T_fit = new TTree("T_fit","T_fit");
    TTree *T_nofit = new TTree("T_nofit","T_nofit");
    T_fit->ReadFile("chi2_fit","a:ab:b:bb");
    T_nofit->ReadFile("chi2_nofit","a:ab:b:bb");

    TF1* chi2 = new TF1("chi2","ROOT::Math::chisquared_pdf(x*[0],[0])*[1]",0,10);

    TCanvas* c_fit = new TCanvas("c_fit","Canvas fit",900,800);
    c_fit->Divide(2,2);
    c_fit->cd(1);
    T_fit->Draw("a");
    if (drawchi2) draw_chi2(chi2);
    c_fit->cd(2);
    T_fit->Draw("ab");
    if (drawchi2) draw_chi2(chi2);
    c_fit->cd(3);
    T_fit->Draw("b");
    if (drawchi2) draw_chi2(chi2);
    c_fit->cd(4);
    T_fit->Draw("bb");
    if (drawchi2) draw_chi2(chi2);
    T_fit->Write();

    TCanvas* c_nofit = new TCanvas("c_nofit","Canvas nofit",900,800);
    c_nofit->Divide(2,2);
    c_nofit->cd(1);
    T_nofit->Draw("a");
    if (drawchi2) draw_chi2(chi2);
    c_nofit->cd(2);
    T_nofit->Draw("ab");
    if (drawchi2) draw_chi2(chi2);
    c_nofit->cd(3);
    T_nofit->Draw("b");
    if (drawchi2) draw_chi2(chi2);
    c_nofit->cd(4);
    T_nofit->Draw("bb");
    if (drawchi2) draw_chi2(chi2);
    T_nofit->Write();
}
