const Double_t PI = 3.1415;

void rgraphs() {
    RooRealVar chi2a ("chi2a","chi2a",0,10);
    RooRealVar chi2ab ("chi2ab","chi2ab",0,10);
    RooRealVar chi2b ("chi2b","chi2b",0,10);
    RooRealVar chi2bb ("chi2bb","chi2bb",0,10);

    RooRealVar ap ("ap","ap",0,1);
    RooRealVar apa("apa","apa",0,2*PI);
    RooRealVar a0 ("a0","a0",0,1);
    RooRealVar a0a("a0a","a0a",0,2*PI);
    RooRealVar at ("at","at",0,1);
    RooRealVar ata("ata","ata",0,2*PI);

    RooRealVar phi("phiw","phiw",0,2*PI);

    RooRealVar rp ("rp","rp",0,0.2);
    RooRealVar r0 ("r0","r0",0,0.2);
    RooRealVar rt ("rt","rt",0,0.2);

    RooRealVar sp ("sp","sp",-PI,PI);
    RooRealVar s0 ("s0","s0",-PI,PI);
    RooRealVar st ("st","st",-PI,PI);

    RooCategory sep("sep","sep");
    sep->defineType("|",1);
    sep->defineType("||",2);

    RooDataSet* dataset;
    dataset = RooDataSet::read("fit.res",RooArgList(chi2a,chi2ab,chi2b,chi2bb,sep,sep));

    RooPlot* frame = chi2a.frame();
    dataset->plotOn(frame);
    frame->Draw();

}
