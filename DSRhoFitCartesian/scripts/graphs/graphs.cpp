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

void graphs() {
    const bool drawchi2 = 1;

    TTree* tree = new TTree("tree","tree");
    tree->ReadFile("fit5.res","chi2a/D:chi2ab:chi2b:chi2bb:sep1/C:api/D:ap:ape:sep2/C:apai/D:apa:apae:sep3/C:a0i/D:a0:a0e:sep4/C:a0ai/D:a0a:a0ae:sep5/C:ati/D:at:ate:sep6/C:atai/D:ata:atae:sep7/C:phiwi/D:phiw:phiwe:sep8/C:rpi/D:rp:rpe:sep9/C:r0i/D:r0:r0e:sep10/C:rti/D:rt:rte:sep11/C:spi/D:sp:spe:sep12/C:s0i/D:s0:s0e:sep13/C:sti/D:st:ste");

    TCanvas* c_a_residuals = new TCanvas("c_a_residuals","Transversity amplitudes residuals",1200,800);
    c_a_residuals->Divide(3,2);
    c_a_residuals->cd(1);
    tree->Draw("(ap-api)");
    c_a_residuals->cd(2);
    tree->Draw("(a0-a0i)");
    c_a_residuals->cd(3);
    tree->Draw("(at-ati)");
    c_a_residuals->cd(4);
    tree->Draw("(apa-apai)");
    c_a_residuals->cd(5);
    tree->Draw("(a0a-a0ai)");
    c_a_residuals->cd(6);
    tree->Draw("(ata-atai)");

    TH1F* htemp;

    gStyle->SetOptFit(111);
    gStyle->SetOptStat(1100);
    TCanvas* c_a_pulls = new TCanvas("c_a_pulls","Transversity amplitudes pulls",1200,800);
    c_a_pulls->Divide(3,2);
    c_a_pulls->cd(1);
    tree->Draw("(ap-api)/ape");
    htemp = (TH1F*)gPad->GetPrimitive("htemp");    
    htemp->Fit("gaus");
    c_a_pulls->cd(2);
    tree->Draw("(a0-a0i)/a0e");
    htemp = (TH1F*)gPad->GetPrimitive("htemp");    
    htemp->Fit("gaus");
    c_a_pulls->cd(3);
    tree->Draw("(at-ati)/ate");
    htemp = (TH1F*)gPad->GetPrimitive("htemp");    
    htemp->Fit("gaus");
    c_a_pulls->cd(4);
    tree->Draw("(apa-apai)/apae");
    htemp = (TH1F*)gPad->GetPrimitive("htemp");    
    htemp->Fit("gaus");
    c_a_pulls->cd(5);
    tree->Draw("(a0a-a0ai)/a0ae");
    htemp = (TH1F*)gPad->GetPrimitive("htemp");    
    htemp->Fit("gaus");
    c_a_pulls->cd(6);
    tree->Draw("(ata-atai)/atae");
    htemp = (TH1F*)gPad->GetPrimitive("htemp");    
    htemp->Fit("gaus");

    TCanvas* c_a_errors = new TCanvas("c_a_errors","Transversity amplitudes errors",1200,800);
    c_a_errors->Divide(3,2);
    c_a_errors->cd(1);
    tree->Draw("ape");
    c_a_errors->cd(2);
    tree->Draw("a0e");
    c_a_errors->cd(3);
    tree->Draw("ate");
    c_a_errors->cd(4);
    tree->Draw("apae");
    c_a_errors->cd(5);
    tree->Draw("a0ae");
    c_a_errors->cd(6);
    tree->Draw("atae");

    TCanvas* c_rs_residuals = new TCanvas("c_rs_residuals","r and strong phase residuals",1200,800);
    c_rs_residuals->Divide(3,2);
    c_rs_residuals->cd(1);
    tree->Draw("(rp-rpi)");
    c_rs_residuals->cd(2);
    tree->Draw("(r0-r0i)");
    c_rs_residuals->cd(3);
    tree->Draw("(rt-rti)");
    c_rs_residuals->cd(4);
    tree->Draw("(sp-spi)");
    c_rs_residuals->cd(5);
    tree->Draw("(s0-s0i)");
    c_rs_residuals->cd(6);
    tree->Draw("(st-sti)");

    TCanvas* c_rs_pulls = new TCanvas("c_rs_pulls","r and strong phase pulls",1200,800);
    c_rs_pulls->Divide(3,2);
    c_rs_pulls->cd(1);
    tree->Draw("(rp-rpi)/rpe");
    htemp = (TH1F*)gPad->GetPrimitive("htemp");    
    htemp->Fit("gaus");
    c_rs_pulls->cd(2);
    tree->Draw("(r0-r0i)/r0e");
    htemp = (TH1F*)gPad->GetPrimitive("htemp");    
    htemp->Fit("gaus");
    c_rs_pulls->cd(3);
    tree->Draw("(rt-rti)/rte");
    htemp = (TH1F*)gPad->GetPrimitive("htemp");    
    htemp->Fit("gaus");
    c_rs_pulls->cd(4);
    tree->Draw("(sp-spi)/spe");
    htemp = (TH1F*)gPad->GetPrimitive("htemp");    
    htemp->Fit("gaus");
    c_rs_pulls->cd(5);
    tree->Draw("(s0-s0i)/s0e");
    htemp = (TH1F*)gPad->GetPrimitive("htemp");    
    htemp->Fit("gaus");
    c_rs_pulls->cd(6);
    tree->Draw("(st-sti)/ste");
    htemp = (TH1F*)gPad->GetPrimitive("htemp");    
    htemp->Fit("gaus");

    TCanvas* c_rs_errors = new TCanvas("c_rs_errors","r and strong phase errors",1200,800);
    c_rs_errors->Divide(3,2);
    c_rs_errors->cd(1);
    tree->Draw("rpe");
    c_rs_errors->cd(2);
    tree->Draw("r0e");
    c_rs_errors->cd(3);
    tree->Draw("rte");
    c_rs_errors->cd(4);
    tree->Draw("spe");
    c_rs_errors->cd(5);
    tree->Draw("s0e");
    c_rs_errors->cd(6);
    tree->Draw("ste");

    TCanvas* c_phiw = new TCanvas("c_phiw","phiw",1200,400);
    c_phiw->Divide(3,1);
    c_phiw->cd(1);
    tree->Draw("(phiw-phiwi)");
    c_phiw->cd(2);
    tree->Draw("(phiw-phiwi)/phiwe");
    htemp = (TH1F*)gPad->GetPrimitive("htemp");    
    htemp->Fit("gaus");
    c_phiw->cd(3);
    tree->Draw("phiwe");

    TF1* chi2 = new TF1("chi2","ROOT::Math::chisquared_pdf(x*[0],[0])*[1]",0,10);

    TCanvas* c_chi2 = new TCanvas("c_chi2","chi2",900,800);
    c_chi2->Divide(2,2);
    c_chi2->cd(1);
    tree->Draw("chi2a");
    if (drawchi2) draw_chi2(chi2);
    c_chi2->cd(2);
    tree->Draw("chi2ab");
    if (drawchi2) draw_chi2(chi2);
    c_chi2->cd(3);
    tree->Draw("chi2b");
    if (drawchi2) draw_chi2(chi2);
    c_chi2->cd(4);
    tree->Draw("chi2bb");
    if (drawchi2) draw_chi2(chi2);
}
