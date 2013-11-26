#include "ObservablesCollectionCartesian.h"
#include "DSRhoGraphs.h"

ObservablesCollection::ObservablesCollection()
{
    chi2a = new RooRealVar("chi2a","chi2a",0,10);
    chi2ab = new RooRealVar("chi2ab","chi2ab",0,10);
    chi2b = new RooRealVar("chi2b","chi2b",0,10);
    chi2bb = new RooRealVar("chi2bb","chi2bb",0,10);

    ap = new RooRealVar("ap","ap",0,1);
    api = new RooRealVar("api","api",0,1);
    ape = new RooRealVar("ape","ape",0,1);
    apa= new RooRealVar("apa","apa",0,2*PI);
    apai= new RooRealVar("apai","apai",0,2*PI);
    apae= new RooRealVar("apae","apae",0,2*PI);
    a0 = new RooRealVar("a0","a0",0,1);
    a0i = new RooRealVar("a0i","a0i",0,1);
    a0e = new RooRealVar("a0e","a0e",0,1);
    a0a= new RooRealVar("a0a","a0a",0,2*PI);
    a0ai= new RooRealVar("a0ai","a0ai",0,2*PI);
    a0ae= new RooRealVar("a0ae","a0ae",0,2*PI);
    at = new RooRealVar("at","at",0,1);
    ati = new RooRealVar("ati","ati",0,1);
    ate = new RooRealVar("ate","ate",0,1);
    ata= new RooRealVar("ata","ata",0,2*PI);
    atai= new RooRealVar("atai","atai",0,2*PI);
    atae= new RooRealVar("atae","atae",0,2*PI);

    phiw= new RooRealVar("phiw","phiw",0,2*PI);
    phiwi= new RooRealVar("phiwi","phiwi",0,2*PI);
    phiwe= new RooRealVar("phiwe","phiwe",0,2*PI);

    xp = new RooRealVar("xp","xp",-1,1);
    xpi = new RooRealVar("xpi","xpi",-1,1);
    xpe = new RooRealVar("xpe","xpe",0,1);
    x0 = new RooRealVar("x0","x0",-1,1);
    x0i = new RooRealVar("x0i","x0i",-1,1);
    x0e = new RooRealVar("x0e","x0e",0,1);
    xt = new RooRealVar("xt","xt",-1,1);
    xti = new RooRealVar("xti","xti",-1,1);
    xte = new RooRealVar("xte","xte",0,1);

    yp = new RooRealVar("yp","yp",-PI,PI);
    ypi = new RooRealVar("ypi","ypi",-1*PI,1*PI);
    ype = new RooRealVar("ype","ype",0,20);
    y0 = new RooRealVar("y0","y0",-PI,PI);
    y0i = new RooRealVar("y0i","y0i",-1*PI,1*PI);
    y0e = new RooRealVar("y0e","y0e",0,20);
    yt = new RooRealVar("yt","yt",-PI,PI);
    yti = new RooRealVar("yti","yti",-1*PI,1*PI);
    yte = new RooRealVar("yte","yte",0,20);
}

ObservablesCollection::~ObservablesCollection()
{
    //dtor
}

RooArgList ObservablesCollection::CreateArgList()
{
    RooCategory* sep = new RooCategory("sep","sep");
    sep->defineType("|",1);
    sep->defineType("||",2);

    RooArgList argList(*chi2a,*chi2ab,*chi2b,*chi2bb,*sep,*sep);
    argList.add(RooArgList(*api,*ap,*ape,*sep,*apai,*apa,*apae,*sep));
    argList.add(RooArgList(*a0i,*a0,*a0e,*sep,*a0ai,*a0a,*a0ae,*sep));
    argList.add(RooArgList(*ati,*at,*ate,*sep,*atai,*ata,*atae,*sep));
    argList.add(RooArgList(*sep,*phiwi,*phiw,*phiwe,*sep,*sep));
    argList.add(RooArgList(*xpi,*xp,*xpe,*sep));
    argList.add(RooArgList(*x0i,*x0,*x0e,*sep));
    argList.add(RooArgList(*xti,*xt,*xte,*sep,*sep));
    argList.add(RooArgList(*ypi,*yp,*ype,*sep));
    argList.add(RooArgList(*y0i,*y0,*y0e,*sep));
    argList.add(RooArgList(*yti,*yt,*yte));
    return argList;
}

void ObservablesCollection::BindToDataset(RooDataSet* dataset){
    const RooArgSet* vars = dataset->get();
    chi2a = static_cast<RooRealVar*>(vars->find(chi2a->GetName()));
    chi2ab = static_cast<RooRealVar*>(vars->find(chi2ab->GetName()));
    chi2b = static_cast<RooRealVar*>(vars->find(chi2b->GetName()));
    chi2bb = static_cast<RooRealVar*>(vars->find(chi2bb->GetName()));

    ap = static_cast<RooRealVar*>(vars->find(ap->GetName()));
    api = static_cast<RooRealVar*>(vars->find(api->GetName()));
    ape = static_cast<RooRealVar*>(vars->find(ape->GetName()));
    apa= static_cast<RooRealVar*>(vars->find(apa->GetName()));
    apai= static_cast<RooRealVar*>(vars->find(apai->GetName()));
    apae= static_cast<RooRealVar*>(vars->find(apae->GetName()));
    a0 = static_cast<RooRealVar*>(vars->find(a0->GetName()));
    a0i = static_cast<RooRealVar*>(vars->find(a0i->GetName()));
    a0e = static_cast<RooRealVar*>(vars->find(a0e->GetName()));
    a0a= static_cast<RooRealVar*>(vars->find(a0a->GetName()));
    a0ai= static_cast<RooRealVar*>(vars->find(a0ai->GetName()));
    a0ae= static_cast<RooRealVar*>(vars->find(a0ae->GetName()));
    at = static_cast<RooRealVar*>(vars->find(at->GetName()));
    ati = static_cast<RooRealVar*>(vars->find(ati->GetName()));
    ate = static_cast<RooRealVar*>(vars->find(ate->GetName()));
    ata= static_cast<RooRealVar*>(vars->find(ata->GetName()));
    atai= static_cast<RooRealVar*>(vars->find(atai->GetName()));
    atae= static_cast<RooRealVar*>(vars->find(atae->GetName()));

    phiw= static_cast<RooRealVar*>(vars->find(phiw->GetName()));
    phiwi= static_cast<RooRealVar*>(vars->find(phiwi->GetName()));
    phiwe= static_cast<RooRealVar*>(vars->find(phiwe->GetName()));

    xp = static_cast<RooRealVar*>(vars->find(xp->GetName()));
    xpi = static_cast<RooRealVar*>(vars->find(xpi->GetName()));
    xpe = static_cast<RooRealVar*>(vars->find(xpe->GetName()));
    x0 = static_cast<RooRealVar*>(vars->find(x0->GetName()));
    x0i = static_cast<RooRealVar*>(vars->find(x0i->GetName()));
    x0e = static_cast<RooRealVar*>(vars->find(x0e->GetName()));
    xt = static_cast<RooRealVar*>(vars->find(xt->GetName()));
    xti = static_cast<RooRealVar*>(vars->find(xti->GetName()));
    xte = static_cast<RooRealVar*>(vars->find(xte->GetName()));

    yp = static_cast<RooRealVar*>(vars->find(yp->GetName()));
    ypi = static_cast<RooRealVar*>(vars->find(ypi->GetName()));
    ype = static_cast<RooRealVar*>(vars->find(ype->GetName()));
    y0 = static_cast<RooRealVar*>(vars->find(y0->GetName()));
    y0i = static_cast<RooRealVar*>(vars->find(y0i->GetName()));
    y0e = static_cast<RooRealVar*>(vars->find(y0e->GetName()));
    yt = static_cast<RooRealVar*>(vars->find(yt->GetName()));
    yti = static_cast<RooRealVar*>(vars->find(yti->GetName()));
    yte = static_cast<RooRealVar*>(vars->find(yte->GetName()));

    errors[0] = ape;
    errors[1] = apae;
    errors[2] = a0e;
    errors[3] = a0ae;
    errors[4] = ate;
    errors[5] = atae;

    errors[6] = phiwe;

    errors[7] = xpe;
    errors[8] = x0e;
    errors[9] = xte;
    errors[10] = ype;
    errors[11] = y0e;
    errors[12] = yte;

}

void ObservablesCollection::CreateResidualsAndPulls(RooDataSet* dataset){
    RooFormulaVar f_residual_ap("residual_ap","(ap-api)",RooArgSet(*ap,*api,*ape));
    RooFormulaVar f_residual_apa("residual_apa","(apa-apai)",RooArgSet(*apa,*apai,*apae));
    RooFormulaVar f_residual_a0("residual_a0","(a0-a0i)",RooArgSet(*a0,*a0i,*a0e));
    RooFormulaVar f_residual_a0a("residual_a0a","(a0a-a0ai)",RooArgSet(*a0a,*a0ai,*a0ae));
    RooFormulaVar f_residual_at("residual_at","(at-ati)",RooArgSet(*at,*ati,*ate));
    RooFormulaVar f_residual_ata("residual_ata","(ata-atai)",RooArgSet(*ata,*atai,*atae));

    RooFormulaVar f_residual_phiw("residual_phiw","(phiw-phiwi)",RooArgSet(*phiw,*phiwi,*phiwe));

    RooFormulaVar f_residual_xp("residual_xp","(xp-xpi)",RooArgSet(*xp,*xpi,*xpe));
    RooFormulaVar f_residual_x0("residual_x0","(x0-x0i)",RooArgSet(*x0,*x0i,*x0e));
    RooFormulaVar f_residual_xt("residual_xt","(xt-xti)",RooArgSet(*xt,*xti,*xte));

    RooFormulaVar f_residual_yp("residual_yp","(yp-ypi)",RooArgSet(*yp,*ypi,*ype));
    RooFormulaVar f_residual_y0("residual_y0","(y0-y0i)",RooArgSet(*y0,*y0i,*y0e));
    RooFormulaVar f_residual_yt("residual_yt","(yt-yti)",RooArgSet(*yt,*yti,*yte));

    residual_ap = static_cast<RooRealVar*>(dataset->addColumn(f_residual_ap));
    residual_apa = static_cast<RooRealVar*>(dataset->addColumn(f_residual_apa));
    residual_a0 = static_cast<RooRealVar*>(dataset->addColumn(f_residual_a0));
    residual_a0a = static_cast<RooRealVar*>(dataset->addColumn(f_residual_a0a));
    residual_at = static_cast<RooRealVar*>(dataset->addColumn(f_residual_at));
    residual_ata = static_cast<RooRealVar*>(dataset->addColumn(f_residual_ata));

    residual_phiw = static_cast<RooRealVar*>(dataset->addColumn(f_residual_phiw));

    residual_xp = static_cast<RooRealVar*>(dataset->addColumn(f_residual_xp));
    residual_x0 = static_cast<RooRealVar*>(dataset->addColumn(f_residual_x0));
    residual_xt = static_cast<RooRealVar*>(dataset->addColumn(f_residual_xt));

    residual_yp = static_cast<RooRealVar*>(dataset->addColumn(f_residual_yp));
    residual_y0 = static_cast<RooRealVar*>(dataset->addColumn(f_residual_y0));
    residual_yt = static_cast<RooRealVar*>(dataset->addColumn(f_residual_yt));

    residuals[0] = residual_ap;
    residuals[1] = residual_apa;
    residuals[2] = residual_a0;
    residuals[3] = residual_a0a;
    residuals[4] = residual_at;
    residuals[5] = residual_ata;

    residuals[6] = residual_phiw;

    residuals[7] = residual_xp;
    residuals[8] = residual_x0;
    residuals[9] = residual_xt;
    residuals[10] = residual_yp;
    residuals[11] = residual_y0;
    residuals[12] = residual_yt;


    RooFormulaVar f_pull_ap("pull_ap","(ap-api)/ape",RooArgSet(*ap,*api,*ape));
    RooFormulaVar f_pull_apa("pull_apa","(apa-apai)/apae",RooArgSet(*apa,*apai,*apae));
    RooFormulaVar f_pull_a0("pull_a0","(a0-a0i)/a0e",RooArgSet(*a0,*a0i,*a0e));
    RooFormulaVar f_pull_a0a("pull_a0a","(a0a-a0ai)/a0ae",RooArgSet(*a0a,*a0ai,*a0ae));
    RooFormulaVar f_pull_at("pull_at","(at-ati)/ate",RooArgSet(*at,*ati,*ate));
    RooFormulaVar f_pull_ata("pull_ata","(ata-atai)/atae",RooArgSet(*ata,*atai,*atae));

    RooFormulaVar f_pull_phiw("pull_phiw","(phiw-phiwi)/phiwe",RooArgSet(*phiw,*phiwi,*phiwe));

    RooFormulaVar f_pull_xp("pull_xp","(xp-xpi)/xpe",RooArgSet(*xp,*xpi,*xpe));
    RooFormulaVar f_pull_x0("pull_x0","(x0-x0i)/x0e",RooArgSet(*x0,*x0i,*x0e));
    RooFormulaVar f_pull_xt("pull_xt","(xt-xti)/xte",RooArgSet(*xt,*xti,*xte));

    RooFormulaVar f_pull_yp("pull_yp","(yp-ypi)/ype",RooArgSet(*yp,*ypi,*ype));
    RooFormulaVar f_pull_y0("pull_y0","(y0-y0i)/y0e",RooArgSet(*y0,*y0i,*y0e));
    RooFormulaVar f_pull_yt("pull_yt","(yt-yti)/yte",RooArgSet(*yt,*yti,*yte));

    pull_ap = static_cast<RooRealVar*>(dataset->addColumn(f_pull_ap));
    pull_apa = static_cast<RooRealVar*>(dataset->addColumn(f_pull_apa));
    pull_a0 = static_cast<RooRealVar*>(dataset->addColumn(f_pull_a0));
    pull_a0a = static_cast<RooRealVar*>(dataset->addColumn(f_pull_a0a));
    pull_at = static_cast<RooRealVar*>(dataset->addColumn(f_pull_at));
    pull_ata = static_cast<RooRealVar*>(dataset->addColumn(f_pull_ata));

    pull_phiw = static_cast<RooRealVar*>(dataset->addColumn(f_pull_phiw));

    pull_xp = static_cast<RooRealVar*>(dataset->addColumn(f_pull_xp));
    pull_x0 = static_cast<RooRealVar*>(dataset->addColumn(f_pull_x0));
    pull_xt = static_cast<RooRealVar*>(dataset->addColumn(f_pull_xt));

    pull_yp = static_cast<RooRealVar*>(dataset->addColumn(f_pull_yp));
    pull_y0 = static_cast<RooRealVar*>(dataset->addColumn(f_pull_y0));
    pull_yt = static_cast<RooRealVar*>(dataset->addColumn(f_pull_yt));

    pulls[0] = pull_ap;
    pulls[1] = pull_apa;
    pulls[2] = pull_a0;
    pulls[3] = pull_a0a;
    pulls[4] = pull_at;
    pulls[5] = pull_ata;

    pulls[6] = pull_phiw;

    pulls[7] = pull_xp;
    pulls[8] = pull_x0;
    pulls[9] = pull_xt;
    pulls[10] = pull_yp;
    pulls[11] = pull_y0;
    pulls[12] = pull_yt;
}
