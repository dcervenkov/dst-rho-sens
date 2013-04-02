#include "ObservablesCollection.h"
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

    rp = new RooRealVar("rp","rp",0,1);
    rpi = new RooRealVar("rpi","rpi",0,1);
    rpe = new RooRealVar("rpe","rpe",0,1);
    r0 = new RooRealVar("r0","r0",0,1);
    r0i = new RooRealVar("r0i","r0i",0,1);
    r0e = new RooRealVar("r0e","r0e",0,1);
    rt = new RooRealVar("rt","rt",0,1);
    rti = new RooRealVar("rti","rti",0,1);
    rte = new RooRealVar("rte","rte",0,1);

    sp = new RooRealVar("sp","sp",-PI,PI);
    spi = new RooRealVar("spi","spi",-PI,PI);
    spe = new RooRealVar("spe","spe",0,20);
    s0 = new RooRealVar("s0","s0",-PI,PI);
    s0i = new RooRealVar("s0i","s0i",-PI,PI);
    s0e = new RooRealVar("s0e","s0e",0,20);
    st = new RooRealVar("st","st",-PI,PI);
    sti = new RooRealVar("sti","sti",-PI,PI);
    ste = new RooRealVar("ste","ste",0,20);
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

    RooDataSet* dataset;
    RooArgList argList(*chi2a,*chi2ab,*chi2b,*chi2bb,*sep,*sep);
    argList.add(RooArgList(*api,*ap,*ape,*sep,*apai,*apa,*apae,*sep));
    argList.add(RooArgList(*a0i,*a0,*a0e,*sep,*a0ai,*a0a,*a0ae,*sep));
    argList.add(RooArgList(*ati,*at,*ate,*sep,*atai,*ata,*atae,*sep));
    argList.add(RooArgList(*sep,*phiwi,*phiw,*phiwe,*sep,*sep));
    argList.add(RooArgList(*rpi,*rp,*rpe,*sep));
    argList.add(RooArgList(*r0i,*r0,*r0e,*sep));
    argList.add(RooArgList(*rti,*rt,*rte,*sep,*sep));
    argList.add(RooArgList(*spi,*sp,*spe,*sep));
    argList.add(RooArgList(*s0i,*s0,*s0e,*sep));
    argList.add(RooArgList(*sti,*st,*ste));
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

    rp = static_cast<RooRealVar*>(vars->find(rp->GetName()));
    rpi = static_cast<RooRealVar*>(vars->find(rpi->GetName()));
    rpe = static_cast<RooRealVar*>(vars->find(rpe->GetName()));
    r0 = static_cast<RooRealVar*>(vars->find(r0->GetName()));
    r0i = static_cast<RooRealVar*>(vars->find(r0i->GetName()));
    r0e = static_cast<RooRealVar*>(vars->find(r0e->GetName()));
    rt = static_cast<RooRealVar*>(vars->find(rt->GetName()));
    rti = static_cast<RooRealVar*>(vars->find(rti->GetName()));
    rte = static_cast<RooRealVar*>(vars->find(rte->GetName()));

    sp = static_cast<RooRealVar*>(vars->find(sp->GetName()));
    spi = static_cast<RooRealVar*>(vars->find(spi->GetName()));
    spe = static_cast<RooRealVar*>(vars->find(spe->GetName()));
    s0 = static_cast<RooRealVar*>(vars->find(s0->GetName()));
    s0i = static_cast<RooRealVar*>(vars->find(s0i->GetName()));
    s0e = static_cast<RooRealVar*>(vars->find(s0e->GetName()));
    st = static_cast<RooRealVar*>(vars->find(st->GetName()));
    sti = static_cast<RooRealVar*>(vars->find(sti->GetName()));
    ste = static_cast<RooRealVar*>(vars->find(ste->GetName()));

    errors[0] = ape;
    errors[1] = apae;
    errors[2] = a0e;
    errors[3] = a0ae;
    errors[4] = ate;
    errors[5] = atae;

    errors[6] = phiwe;

    errors[7] = rpe;
    errors[8] = r0e;
    errors[9] = rte;
    errors[10] = spe;
    errors[11] = s0e;
    errors[12] = ste;

}

void ObservablesCollection::CreateResidualsAndPulls(RooDataSet* dataset){
    RooFormulaVar f_residual_ap("residual_ap","(ap-api)",RooArgSet(*ap,*api,*ape));
    RooFormulaVar f_residual_apa("residual_apa","(apa-apai)",RooArgSet(*apa,*apai,*apae));
    RooFormulaVar f_residual_a0("residual_a0","(a0-a0i)",RooArgSet(*a0,*a0i,*a0e));
    RooFormulaVar f_residual_a0a("residual_a0a","(a0a-a0ai)",RooArgSet(*a0a,*a0ai,*a0ae));
    RooFormulaVar f_residual_at("residual_at","(at-ati)",RooArgSet(*at,*ati,*ate));
    RooFormulaVar f_residual_ata("residual_ata","(ata-atai)",RooArgSet(*ata,*atai,*atae));

    RooFormulaVar f_residual_phiw("residual_phiw","(phiw-phiwi)",RooArgSet(*phiw,*phiwi,*phiwe));

    RooFormulaVar f_residual_rp("residual_rp","(rp-rpi)",RooArgSet(*rp,*rpi,*rpe));
    RooFormulaVar f_residual_r0("residual_r0","(r0-r0i)",RooArgSet(*r0,*r0i,*r0e));
    RooFormulaVar f_residual_rt("residual_rt","(rt-rti)",RooArgSet(*rt,*rti,*rte));

    RooFormulaVar f_residual_sp("residual_sp","(sp-spi)",RooArgSet(*sp,*spi,*spe));
    RooFormulaVar f_residual_s0("residual_s0","(s0-s0i)",RooArgSet(*s0,*s0i,*s0e));
    RooFormulaVar f_residual_st("residual_st","(st-sti)",RooArgSet(*st,*sti,*ste));

    residual_ap = static_cast<RooRealVar*>(dataset->addColumn(f_residual_ap));
    residual_apa = static_cast<RooRealVar*>(dataset->addColumn(f_residual_apa));
    residual_a0 = static_cast<RooRealVar*>(dataset->addColumn(f_residual_a0));
    residual_a0a = static_cast<RooRealVar*>(dataset->addColumn(f_residual_a0a));
    residual_at = static_cast<RooRealVar*>(dataset->addColumn(f_residual_at));
    residual_ata = static_cast<RooRealVar*>(dataset->addColumn(f_residual_ata));

    residual_phiw = static_cast<RooRealVar*>(dataset->addColumn(f_residual_phiw));

    residual_rp = static_cast<RooRealVar*>(dataset->addColumn(f_residual_rp));
    residual_r0 = static_cast<RooRealVar*>(dataset->addColumn(f_residual_r0));
    residual_rt = static_cast<RooRealVar*>(dataset->addColumn(f_residual_rt));

    residual_sp = static_cast<RooRealVar*>(dataset->addColumn(f_residual_sp));
    residual_s0 = static_cast<RooRealVar*>(dataset->addColumn(f_residual_s0));
    residual_st = static_cast<RooRealVar*>(dataset->addColumn(f_residual_st));

    residuals[0] = residual_ap;
    residuals[1] = residual_apa;
    residuals[2] = residual_a0;
    residuals[3] = residual_a0a;
    residuals[4] = residual_at;
    residuals[5] = residual_ata;

    residuals[6] = residual_phiw;

    residuals[7] = residual_rp;
    residuals[8] = residual_r0;
    residuals[9] = residual_rt;
    residuals[10] = residual_sp;
    residuals[11] = residual_s0;
    residuals[12] = residual_st;


    RooFormulaVar f_pull_ap("pull_ap","(ap-api)/ape",RooArgSet(*ap,*api,*ape));
    RooFormulaVar f_pull_apa("pull_apa","(apa-apai)/apae",RooArgSet(*apa,*apai,*apae));
    RooFormulaVar f_pull_a0("pull_a0","(a0-a0i)/a0e",RooArgSet(*a0,*a0i,*a0e));
    RooFormulaVar f_pull_a0a("pull_a0a","(a0a-a0ai)/a0ae",RooArgSet(*a0a,*a0ai,*a0ae));
    RooFormulaVar f_pull_at("pull_at","(at-ati)/ate",RooArgSet(*at,*ati,*ate));
    RooFormulaVar f_pull_ata("pull_ata","(ata-atai)/atae",RooArgSet(*ata,*atai,*atae));

    RooFormulaVar f_pull_phiw("pull_phiw","(phiw-phiwi)/phiwe",RooArgSet(*phiw,*phiwi,*phiwe));

    RooFormulaVar f_pull_rp("pull_rp","(rp-rpi)/rpe",RooArgSet(*rp,*rpi,*rpe));
    RooFormulaVar f_pull_r0("pull_r0","(r0-r0i)/r0e",RooArgSet(*r0,*r0i,*r0e));
    RooFormulaVar f_pull_rt("pull_rt","(rt-rti)/rte",RooArgSet(*rt,*rti,*rte));

    RooFormulaVar f_pull_sp("pull_sp","(sp-spi)/spe",RooArgSet(*sp,*spi,*spe));
    RooFormulaVar f_pull_s0("pull_s0","(s0-s0i)/s0e",RooArgSet(*s0,*s0i,*s0e));
    RooFormulaVar f_pull_st("pull_st","(st-sti)/ste",RooArgSet(*st,*sti,*ste));

    pull_ap = static_cast<RooRealVar*>(dataset->addColumn(f_pull_ap));
    pull_apa = static_cast<RooRealVar*>(dataset->addColumn(f_pull_apa));
    pull_a0 = static_cast<RooRealVar*>(dataset->addColumn(f_pull_a0));
    pull_a0a = static_cast<RooRealVar*>(dataset->addColumn(f_pull_a0a));
    pull_at = static_cast<RooRealVar*>(dataset->addColumn(f_pull_at));
    pull_ata = static_cast<RooRealVar*>(dataset->addColumn(f_pull_ata));

    pull_phiw = static_cast<RooRealVar*>(dataset->addColumn(f_pull_phiw));

    pull_rp = static_cast<RooRealVar*>(dataset->addColumn(f_pull_rp));
    pull_r0 = static_cast<RooRealVar*>(dataset->addColumn(f_pull_r0));
    pull_rt = static_cast<RooRealVar*>(dataset->addColumn(f_pull_rt));

    pull_sp = static_cast<RooRealVar*>(dataset->addColumn(f_pull_sp));
    pull_s0 = static_cast<RooRealVar*>(dataset->addColumn(f_pull_s0));
    pull_st = static_cast<RooRealVar*>(dataset->addColumn(f_pull_st));

    pulls[0] = pull_ap;
    pulls[1] = pull_apa;
    pulls[2] = pull_a0;
    pulls[3] = pull_a0a;
    pulls[4] = pull_at;
    pulls[5] = pull_ata;

    pulls[6] = pull_phiw;

    pulls[7] = pull_rp;
    pulls[8] = pull_r0;
    pulls[9] = pull_rt;
    pulls[10] = pull_sp;
    pulls[11] = pull_s0;
    pulls[12] = pull_st;
}

