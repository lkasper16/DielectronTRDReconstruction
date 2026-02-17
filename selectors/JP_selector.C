#define JP_selector_cxx
// The class definition in JP_selector.h has been generated automatically
// by the ROOT utility TTree::MakeSelector(). This class is derived
// from the ROOT class TSelector. For more information on the TSelector
// framework see $ROOTSYS/README/README.SELECTOR or the ROOT User Manual.

// The following methods are defined in this file:
//    Begin():        called every time a loop on the tree starts,
//                    a convenient place to create your histograms.
//    SlaveBegin():   called after Begin(), when on PROOF called only on the
//                    slave servers.
//    Process():      called for each event, in this function you decide what
//                    to read and fill your histograms.
//    SlaveTerminate: called at the end of the loop on the tree, when on PROOF
//                    called only on the slave servers.
//    Terminate():    called at the end of the loop on the tree,
//                    a convenient place to draw/fit your histograms.
//
// To use this file, try the following session on your Tree T:
//
// Root > T->Process("JP_selector.C")
// Root > T->Process("JP_selector.C","some options")
// Root > T->Process("JP_selector.C+")
//

#include "JP_selector.h"
#include <TH2.h>
#include <TStyle.h>

#define PI 3.14159265
#define p_mass2 0.880354
#define p_mass 0.938272
#define pi_mass2 0.01948
#define pi_mass 0.13957
TFile *fout;
TFile *hflux;
TH1F *myflux;
TAxis *fax;

// output for event-by-event fits
TNtuple *ntuple;
float prompt,fweight,tweight;

TFile *bcal_bg;
TH1F *bhbg_factorb;
TAxis *xbhb;
TFile *fcal_bg;
TH1F *bhbg_factorf;
TAxis *xbhf;
TF1 *sige;
TF1 *dsige;
TF1 *fcalbg;
TF1 *bcalbg;
TF1 *mcxsec;
TF1 *bpove_cut;
TF1 *bpove_cuts;
TF1 *fpove_cut;
TF1 *fpove_cuts;
float smgl; //mass of the gluons(2 or 3) squared
float tslope; //used in jpsi simulations
int ncount;
int ncombos;
int acombos;
double JP_selector::ecorr(double pp, double thp, double ebeam)
{
#define p_mass2 0.880354
#define p_mass 0.938272
double minv0=sqrt(-2.*(ebeam+p_mass)*(sqrt(pp*pp+p_mass2)-p_mass)+2.*ebeam*pp*cos(thp));
double minv1=sqrt(-2.*(ebeam*1.0025+p_mass)*(sqrt(pp*pp+p_mass2)-p_mass)+2.*ebeam*1.0025*pp*cos(thp));
return minv1-minv0;
//return minv1;
}
double JP_selector::ftmin(double jmass, double ebeam)
{
float pmass=0.9383;
//float jmass=3.097;
float s=pmass*pmass+2.*pmass*ebeam;
float ss=sqrt(s);
float ecm1=(s-pmass*pmass)/2./ss;
float ecm2=(s+jmass*jmass-pmass*pmass)/2./ss;
float pcm1=ecm1;
float pcm2=sqrt(ecm2*ecm2-jmass*jmass);
double tmn=2.*(ecm1*ecm2-pcm1*pcm2)-jmass*jmass;
return tmn;
}

double JP_selector::ftmax(double jmass, double ebeam)  {
float pmass=0.9383;
//float jmass=3.097;
float s=pmass*pmass+2.*pmass*ebeam;
float ss=sqrt(s);
float ecm1=(s-pmass*pmass)/2./ss;
float ecm2=(s+jmass*jmass-pmass*pmass)/2./ss;
float pcm1=ecm1;
float pcm2=sqrt(ecm2*ecm2-jmass*jmass);
double tmx=2.*(ecm1*ecm2+pcm1*pcm2)-jmass*jmass;
return tmx;
}


void JP_selector::Begin(TTree * /*tree*/)
{
   // The Begin() function is called at the start of the query.
   // When running with PROOF Begin() is only called on the client.
   // The tree argument is deprecated (on PROOF 0 is passed).

//cokcok
 nrsigmas = 3.0; //number of sigmas to cut E/p
// nrsigmas = 6.0; //number of sigmas to cut E/p
 nlsigmas = 3.0; //number of sigmas to cut E/p
// nlsigmas = 6.0; //number of sigmas to cut E/p
 ndesigmas= 6.0; //number of sigmas for dE/dx cut

   TString option = GetOption();


smgl=2.2;
tslope=1.;
ncount=0;
ncombos=0;
acombos=0;


dsige = new TF1("dsige","0.900-0.412*x+0.038*x*x",8.2,12.);
sige = new TF1("sige","-0.00902003+0.00111854*x",8.2,12.);

mcxsec = new TF1("mcxsec","pol6",8.2,11.8);
mcxsec->SetParameters(-899.577,501.601,-115.956,14.2234,-0.976859,0.0356609,-0.000540931);

hflux = new TFile("flux_80bins.root");
myflux = (TH1F*)hflux->Get("myflux");
//myflux->Scale(1./myflux->Integral());
//coldcold cokcok myflux->Scale(1./4.);
fax = (TAxis*)myflux->GetXaxis();

fcal_bg = new TFile("mPovE_FCAL.root");
bcal_bg = new TFile("mPovE_BCAL.root");
   char foutname[80];
   sprintf(foutname,"s%dde%d.root",(int)nrsigmas,(int)ndesigmas);
   cout<<foutname<<endl;
fout = new TFile(foutname,"recreate");
bhbg_factorf = (TH1F*)fcal_bg->Get("bhbg_factorf");
xbhf = (TAxis*)bhbg_factorf->GetXaxis();
bhbg_factorb = (TH1F*)bcal_bg->Get("bhbg_factorb");
xbhb = (TAxis*)bhbg_factorb->GetXaxis();

fcalbg = new TF1("fcalbg","0.225530+0.220108*x-0.0881692*x*x",1.2,3.5);
bcalbg = new TF1("bcalbg","1.79669-2.21354*x+1.09814*x*x-0.173206*x*x*x",1.2,3.5);


bcalbg = new TF1("bcalbg","-0.674666+1.62783*x-0.833527*x*x+0.130829*x*x*x",1.2,3.5);
fcalbg = new TF1("fcalbg","-0.0406324+0.451861*x-0.133072*x*x",1.2,3.5);
  
//last
fcalbg = new TF1("fcalbg","0.00576061+0.456857*x-0.1281580*x*x",1.2,3.5);
bcalbg = new TF1("bcalbg","0.16108+0.198688*x-0.058957*x*x",1.2,3.5);

bpove_cut = new TF1("bpove_cut","1.2989-0.272492*x+0.0879449*x*x-0.00949335*x*x*x",1.2,3.2);
bpove_cuts = new TF1("bpove_cuts","0.158354-0.107121*x+0.0413965*x*x-0.00552854*x*x*x",1.2,3.2);
fpove_cut = new TF1("fpove_cut","1.1256-0.0620255*x+0.0271978*x*x-0.00404517*x*x*x",1.2,3.2);
fpove_cuts = new TF1("fpove_cuts","0.0719702-0.018369*x+0.00628064*x*x-0.000741237*x*x*x",1.2,3.2);

ntuple = new TNtuple("ntuple","data for event-by-event fit","Minv:prompt:t:ebeam:fweight:tweight");

}

void JP_selector::SlaveBegin(TTree * /*tree*/)
{
   // The SlaveBegin() function is called after the Begin() function.
   // When running with PROOF SlaveBegin() is called on each slave server.
   // The tree argument is deprecated (on PROOF 0 is passed).

   TString option = GetOption();

// defining the cuts:

//cokcok
theta_cut=2.*PI/180.;
//theta_cut=0.*PI/180.;

//cokcok
chi2_cut=5000.;
//chi2_cut=500.;
//chi2_cut=12.;
//chi2_cut=5000.e+10;

//cokcok
p_cut=0.4;
//p_cut=0.0;

delta_cut=0.0;

dedx_cut=3.127e-6-ndesigmas*0.5e-6; //dEdx cut

pbcal_cut=0.03;
//cokcok
//pbcal_cut=0.03;
//pbcal_cut=0.00;

theta_cm_cut = 90.-0.; //in deg

bh1=1.2;
bh2=2.5;

// p/E cuts:
float bm_data0=1.04616;
float bs_data0=6.82571e-02;


float fm_data0=1.07380;
float fs_data0=3.91639e-02;
//2018
bm_data0=1.0521;
bs_data0=6.490e-02;
fm_data0= 1.0780;
fs_data0=5.230e-02;

//final JP
bm_data0=1.0245;
bs_data0=4.960e-02;
fm_data0= 1.0748;
fs_data0=5.300e-02;


//cokcok
//no KF
//fm_data0= 1.0840;
//fs_data0=13.230e-02;


bcut_data0=bm_data0+bs_data0*nrsigmas;
fcut_data0=fm_data0+fs_data0*nrsigmas;
bcut_data1=bm_data0-bs_data0*nlsigmas;
fcut_data1=fm_data0-fs_data0*nlsigmas;

cout<<" bcut0, fcut0 ="<<bcut_data0<<" "<<fcut_data0<<" 1/ "<<1./bcut_data0<<" "<<1./fcut_data0<<endl;
cout<<" bcut0, fcut0 ="<<bcut_data1<<" "<<fcut_data1<<" 1/ "<<1./bcut_data1<<" "<<1./fcut_data1<<endl;



// MC scale
// int events_bh=2248653.;
int events_bh=2084890;
events_bh=9393506;

int events_jp=930428+3756692;
events_jp=394076+1627477; //t2.9
events_jp=388401+1410042; //t1.4
events_jp=879306; //t1.4
//int events_phi=196719.; //old
//cold int events_phi=160302.+145921.; //2017LIrnd+2017HIrnd >8.2
int events_phi=160302.+145921.; //2017LIrnd+2017HIrnd >8.2
//int events_phi=154861.+141336.;
events_phi=183925+1628250;


// lumi
//float lumi=(10590.+20520.+16080.); // nb^-1
//float lumi=(10260.); // nb^-1
//float lumi=(20366.); // nb^-1
//correct but need to match jpsi_eff_v3.C with old flux float lumi=(11041.8 + 49553.1); // nb^-1
//cold float lumi=(11068.9 + 39128.7); // nb^-1 used for ANroo !!!!!!!!!!!!!!!!!!!!
//float lumi=(14874. + 53725.3); // nb^-1
lumi=320255.; // nb^-1 
lumi16=13827.300;
lumi17=53387.400;
lumi18s=154630.00;
lumi18f=98409.800;
lumi=lumi16+lumi17+lumi18s+lumi18f;

//cokcok very important - update when redoing MC !!!!! 
events_jp16=415019.00;
events_jp17=454295.00;
events_jp18s=456299.00;
events_jp18f=454525.00;
events_jp=events_jp16+events_jp17+events_jp18s+events_jp18f;
//alex
//events_jp18s=899052.00;


events_bh16=1.85676e+06;
events_bh17=1.95588e+06;
events_bh18s=1.97545e+06;
events_bh18f=1.95498e+06;
events_bh=events_bh16+events_bh17+events_bh18s+events_bh18f;


//cross-sections and branching ratios:
//Mike float xsec_bh=1000. ; // just conversion from pb to nb, the BH x-section is in the weight
//HallB float xsec_bh=0.001/(2.*PI) ; // just conversion from pb to nb, the BH x-section is in the weight
xsec_bh=0.001/(2.*PI) ; // just conversion from pb to nb, the BH x-section is in the weight AND correction for the 2*pi problem 
xsec_jp=0.750/1.36; // nb //adjusted to match jpsi MC and data
br_jp=0.06; // jp->e+e- BR
float xsec_phi=550.; // nb
float br_phi=0.0003; // phi->e+e- BR

scale_bh = lumi*xsec_bh/events_bh;
float scale_jpsi18s = lumi18s*xsec_jp*br_jp/events_jp18s;
cout<<" ---------- scale_jpsi18s ==============="<<scale_jpsi18s<<endl;
float scale_jpsi16 = lumi16*xsec_jp*br_jp/events_jp16;
cout<<" ---------- scale_jpsi16 ==============="<<scale_jpsi16<<endl;

scale_jpsi = lumi*xsec_jp*br_jp/events_jp;
scale_phi = lumi*xsec_phi*br_phi/events_phi;

cout<<" sacle_bh "<<scale_bh<<endl;
cout<<" sacle_jpsi "<<scale_jpsi<<endl;
cout<<" sacle_phi "<<scale_phi<<endl;

cout<<" lumi, events_jp"<<lumi<<" "<<events_jp<<endl;


   run_old=0; event_old=0; 
 revent.open("revent_jp.dat");
   nbins=600;
   //nbins=300;
   mbeg=0.500;
   mend=3.500;
//nbins=16;
//mbeg=1.2;
//mend=8*0.8/3+mbeg;

   //mbeg=0.2500; //Minv^2
   //mend=12.2500; //Minv^2
   //nbins=8;
   //mbeg=1.2;
   //mend=3.33333333333;

   nbins2=500;
   mbeg2=-5.000;
   mend2=5.000;

   hweight = new TH2F("hweight","",100,8.2,12.,100,0.,100.);
   htint = new TH2F("htint","",100,8.2,12.,100,0.,1.0);
   minv_nokf= new TH1F("minv_nokf","MC Minv no KF",nbins,mbeg,mend);
   minvr= new TH1F("minvr","MC Mrec",nbins,mbeg,mend);

   minva = new TH1F("minva","MC Minv all (J/#psi)",nbins,mbeg,mend);
   minvb = new TH1F("minvb","MC Minv accidentials (J/#psi)",nbins,mbeg,mend);
   minvc = new TH1F("minvc","MC J/#psi peak",nbins,mbeg,mend);
   minv = new TH1F("minv","MC J/#psi",nbins,mbeg,mend);
   minvap = new TH1F("minvap","MC Minv all (#phi)",nbins,mbeg,mend);
   minvbp = new TH1F("minvbp","MC Minv accidentials (#phi)",nbins,mbeg,mend);
   minvcp = new TH1F("minvcp","MC #phi peak",nbins,mbeg,mend);
   minvp = new TH1F("minvp","MC #phi",nbins,mbeg,mend);

   minva->Sumw2();
   minvb->Sumw2();
   minvc->Sumw2();
   minv->Sumw2();
   minvap->Sumw2();
   minvbp->Sumw2();
   minvcp->Sumw2();
   minvp->Sumw2();


   dminv_nokf= new TH1F("dminv_nokf","Minv no KF",nbins,mbeg,mend);
   dminvr= new TH1F("dminvr","Mrec",nbins,mbeg,mend);
   dminvra= new TH1F("dminvra","Mrec",nbins,mbeg,mend);
   dminvrb= new TH1F("dminvrb","Mrec",nbins,mbeg,mend);
   dminvrc= new TH1F("dminvrc","Mrec",nbins,mbeg,mend);

   dminva = new TH1F("dminva","Minv all (J/#psi)",nbins,mbeg,mend);
   dminvb = new TH1F("dminvb","Minv accidentials (J/#psi)",nbins,mbeg,mend);
   dminvc = new TH1F("dminvc","J/#psi peak",nbins,mbeg,mend);
   dminv = new TH1F("dminv","J/#psi",nbins,mbeg,mend);
   dminvap = new TH1F("dminvap","Minv all (#phi)",nbins,mbeg,mend);
   dminvbp = new TH1F("dminvbp","Minv accidentials (#phi)",nbins,mbeg,mend);
   dminvcp = new TH1F("dminvcp","#phi peak",nbins,mbeg,mend);
   dminvp = new TH1F("dminvp","#phi",nbins,mbeg,mend);

       dminva_ecomb = new TH2F("dminva_ecomb","Minv vs ecomb in time",nbins,mbeg,mend,nbins2,mbeg2,mend2);
       dminvb_ecomb = new TH2F("dminvb_ecomb","Minv vs ecomb out of time",nbins,mbeg,mend,nbins2,mbeg2,mend2);
       dminva_tcomb = new TH2F("dminva_tcomb","Minv vs tcomb in time",nbins,mbeg,mend,nbins2,mbeg2,mend2);
       dminvb_tcomb = new TH2F("dminvb_tcomb","Minv vs tcomb out of time",nbins,mbeg,mend,nbins2,mbeg2,mend2);

       minva_ecomb = new TH2F("minva_ecomb","Minv vs ecomb in time",nbins,mbeg,mend,nbins2,mbeg2,mend2);
       minvb_ecomb = new TH2F("minvb_ecomb","Minv vs ecomb out of time",nbins,mbeg,mend,nbins2,mbeg2,mend2);
       minva_tcomb = new TH2F("minva_tcomb","Minv vs tcomb in time",nbins,mbeg,mend,nbins2,mbeg2,mend2);
       minvb_tcomb = new TH2F("minvb_tcomb","Minv vs tcomb out of time",nbins,mbeg,mend,nbins2,mbeg2,mend2);

       dninva_ecomb = new TH1F("dninva_ecomb","Minv vs ecomb in time",nbins,mbeg,mend);
       dninvb_ecomb = new TH1F("dninvb_ecomb","Minv vs ecomb out of time",nbins,mbeg,mend);
       dninva_tcomb = new TH1F("dninva_tcomb","Minv vs tcomb in time",nbins,mbeg,mend);
       dninvb_tcomb = new TH1F("dninvb_tcomb","Minv vs tcomb out of time",nbins,mbeg,mend);

       ninva_ecomb = new TH1F("ninva_ecomb","Minv vs ecomb in time",nbins,mbeg,mend);
       ninvb_ecomb = new TH1F("ninvb_ecomb","Minv vs ecomb out of time",nbins,mbeg,mend);
       ninva_tcomb = new TH1F("ninva_tcomb","Minv vs tcomb in time",nbins,mbeg,mend);
       ninvb_tcomb = new TH1F("ninvb_tcomb","Minv vs tcomb out of time",nbins,mbeg,mend);


   dminva->Sumw2();
   dminvb->Sumw2();
   dminvc->Sumw2();
   dminv->Sumw2();
   dminvap->Sumw2();
   dminvbp->Sumw2();
   dminvcp->Sumw2();
   dminvp->Sumw2();

  char hnam[128];
  char htit[128];

  nebins=18;
  ebeg=8.2;
  eend=11.44;
  estep=(eend-ebeg)/nebins;


for (int ebin=0;ebin<nebins;ebin++){
       float e1=ebeg+ebin*estep;
       float e2=e1+estep;

// MC

       sprintf(hnam,"eminva%d",ebin);  
       sprintf(htit,"MC Minv for %d < E <%d in time (J/#psi)",e1,e2);
       eminva[ebin]= new TH1F(hnam,htit,nbins,mbeg,mend);
       eminva[ebin]->Sumw2();

       sprintf(hnam,"eminvb%d",ebin);  
       sprintf(htit,"MC Minv for %f < E <%f out of time (J/#psi)",e1,e2);
       eminvb[ebin]= new TH1F(hnam,htit,nbins,mbeg,mend);
       eminvb[ebin]->Sumw2();

       sprintf(hnam,"eminvc%d",ebin);  
       sprintf(htit,"MC J/#psi Minv for %f < E <%f )",e1,e2);
       eminvc[ebin]= new TH1F(hnam,htit,nbins,mbeg,mend);
       eminvc[ebin]->Sumw2();

       sprintf(hnam,"eminv%d",ebin);  
       sprintf(htit,"MC J/#psi fit of Minv for %f < E <%f accidential corrected",e1,e2);
       eminv[ebin]= new TH1F(hnam,htit,nbins,mbeg,mend);
       eminv[ebin]->Sumw2();

       sprintf(hnam,"eminvap%d",ebin);  
       sprintf(htit,"MC Minv for %f < E <%f in time (#phi)",e1,e2);
       eminvap[ebin]= new TH1F(hnam,htit,nbins,mbeg,mend);
       eminvap[ebin]->Sumw2();

       sprintf(hnam,"eminvbp%d",ebin);  
       sprintf(htit,"MC Minv for %f < E <%f out of time (#phi)",e1,e2);
       eminvbp[ebin]= new TH1F(hnam,htit,nbins,mbeg,mend);
       eminvbp[ebin]->Sumw2();

       sprintf(hnam,"eminvcp%d",ebin);  
       sprintf(htit,"MC #phi Minv for %f < E <%f ",e1,e2);
       eminvcp[ebin]= new TH1F(hnam,htit,nbins,mbeg,mend);
       eminvcp[ebin]->Sumw2();

       sprintf(hnam,"eminvp%d",ebin);  
       sprintf(htit,"MC #phi fit of Minv for %f < E <%f ",e1,e2);
       eminvp[ebin]= new TH1F(hnam,htit,nbins,mbeg,mend);
       eminvp[ebin]->Sumw2();

       sprintf(hnam,"eminva_ecomb%d",ebin);
       sprintf(htit,"MC Minv vs ecomb for %d < E <%d in time (J/#psi)",e1,e2);
       eminva_ecomb[ebin]= new TH2F(hnam,htit,nbins,mbeg,mend,10,-0.5,9.5);

       sprintf(hnam,"eminva_tcomb%d",ebin);
       sprintf(htit,"MC Minv vs tcomb for %d < E <%d in time (J/#psi)",e1,e2);
       eminva_tcomb[ebin]= new TH2F(hnam,htit,nbins,mbeg,mend,10,-0.5,9.5);

       sprintf(hnam,"eminvb_ecomb%d",ebin);
       sprintf(htit,"MC Minv vs ecomb for %d < E <%d out of time (J/#psi)",e1,e2);
       eminvb_ecomb[ebin]= new TH2F(hnam,htit,nbins,mbeg,mend,10,-0.5,9.5);

       sprintf(hnam,"eminvb_tcomb%d",ebin);
       sprintf(htit,"MC Minv vs tcomb for %d < E <%d out of time (J/#psi)",e1,e2);
       eminvb_tcomb[ebin]= new TH2F(hnam,htit,nbins,mbeg,mend,10,-0.5,9.5);


// data

// combinations 
       sprintf(hnam,"deminva_ecomb%d",ebin);
       sprintf(htit,"Minv vs ecomb for %d < E <%d in time (J/#psi)",e1,e2);
       deminva_ecomb[ebin]= new TH2F(hnam,htit,nbins,mbeg,mend,10,-0.5,9.5);

       sprintf(hnam,"deminva_tcomb%d",ebin);
       sprintf(htit,"Minv vs tcomb for %d < E <%d in time (J/#psi)",e1,e2);
       deminva_tcomb[ebin]= new TH2F(hnam,htit,nbins,mbeg,mend,10,-0.5,9.5);

       sprintf(hnam,"deminvb_ecomb%d",ebin);
       sprintf(htit,"Minv vs ecomb for %d < E <%d out of time (J/#psi)",e1,e2);
       deminvb_ecomb[ebin]= new TH2F(hnam,htit,nbins,mbeg,mend,10,-0.5,9.5);

       sprintf(hnam,"deminvb_tcomb%d",ebin);
       sprintf(htit,"Minv vs tcomb for %d < E <%d out of time (J/#psi)",e1,e2);
       deminvb_tcomb[ebin]= new TH2F(hnam,htit,nbins,mbeg,mend,10,-0.5,9.5);

//

       sprintf(hnam,"deminva%d",ebin);
       sprintf(htit,"%6.2f < E <%6.2f",e1,e2);
       deminva[ebin]= new TH1F(hnam,htit,nbins,mbeg,mend);
       deminva[ebin]->Sumw2();

       sprintf(hnam,"deminvb%d",ebin);
       sprintf(htit,"%6.2f < E <%6.2f",e1,e2);
       deminvb[ebin]= new TH1F(hnam,htit,nbins,mbeg,mend);
       deminvb[ebin]->Sumw2();

       sprintf(hnam,"deminvc%d",ebin);
       sprintf(htit,"%6.2f < E <%6.2f )",e1,e2);
       deminvc[ebin]= new TH1F(hnam,htit,nbins,mbeg,mend);
       deminvc[ebin]->Sumw2();

       sprintf(hnam,"deminv%d",ebin);
       sprintf(htit,"%6.2f < E <%6.2f",e1,e2);
       deminv[ebin]= new TH1F(hnam,htit,nbins,mbeg,mend);
       deminv[ebin]->Sumw2();

       sprintf(hnam,"deminvap%d",ebin);
       sprintf(htit,"Minv for %f < E <%f in time (#phi)",e1,e2);
       deminvap[ebin]= new TH1F(hnam,htit,nbins,mbeg,mend);
       deminvap[ebin]->Sumw2();

       sprintf(hnam,"deminvbp%d",ebin);
       sprintf(htit,"Minv for %f < E <%f out of time (#phi)",e1,e2);
       deminvbp[ebin]= new TH1F(hnam,htit,nbins,mbeg,mend);
       deminvbp[ebin]->Sumw2();

       sprintf(hnam,"deminvcp%d",ebin);
       sprintf(htit,"#phi Minv for %f < E <%f ",e1,e2);
       deminvcp[ebin]= new TH1F(hnam,htit,nbins,mbeg,mend);
       deminvcp[ebin]->Sumw2();

       sprintf(hnam,"deminvp%d",ebin);
       sprintf(htit,"#phi fit of Minv for %f < E <%f ",e1,e2);
       deminvp[ebin]= new TH1F(hnam,htit,nbins,mbeg,mend);
       deminvp[ebin]->Sumw2();

}

  //cok ntbins=7;
  ntbins=15;
  ntbins=1;
  tbeg=0.00;
  //coktend=1.05;
  tend=2.25;
  
  tstep=(tend-tbeg)/ntbins;

  temin=ebeg;
  temax=eend;
  temin=10.;

  tstep=(tend-tbeg)/ntbins;

  //cokcok etmin=0.;
  //cokcok etmax=10.;
  etmin=0.;
  etmax=10.;

  ethist=new TH2D("ethist","",nebins,ebeg,eend,ntbins,tbeg,tend);


for (int tbin=0;tbin<ntbins;tbin++){
       float t1=tbeg+tbin*tstep;
       float t2=t1+tstep;

// MC

       sprintf(hnam,"tminva%d",tbin);  
       sprintf(htit,"MC Minv for %d < t <%d in time (J/#psi)",t1,t2);
       tminva[tbin]= new TH1F(hnam,htit,nbins,mbeg,mend);
       tminva[tbin]->Sumw2();

       sprintf(hnam,"tminvb%d",tbin);  
       sprintf(htit,"MC Minv for %f < t <%f out of time (J/#psi)",t1,t2);
       tminvb[tbin]= new TH1F(hnam,htit,nbins,mbeg,mend);
       tminvb[tbin]->Sumw2();

       sprintf(hnam,"tminvc%d",tbin);  
       sprintf(htit,"MC J/#psi Minv for %f < t <%f )",t1,t2);
       tminvc[tbin]= new TH1F(hnam,htit,nbins,mbeg,mend);
       tminvc[tbin]->Sumw2();

       sprintf(hnam,"tminv%d",tbin);  
       sprintf(htit,"MC J/#psi fit of Minv for %f < t <%f accidential corrected",t1,t2);
       tminv[tbin]= new TH1F(hnam,htit,nbins,mbeg,mend);
       tminv[tbin]->Sumw2();

       sprintf(hnam,"tminvap%d",tbin);  
       sprintf(htit,"MC Minv for %d < t <%d in time (phi)",t1,t2);
       tminvap[tbin]= new TH1F(hnam,htit,nbins,mbeg,mend);
       tminvap[tbin]->Sumw2();

       sprintf(hnam,"tminvbp%d",tbin);  
       sprintf(htit,"MC Minv for %f < t <%f out of time (phi)",t1,t2);
       tminvbp[tbin]= new TH1F(hnam,htit,nbins,mbeg,mend);
       tminvbp[tbin]->Sumw2();

       sprintf(hnam,"tminvcp%d",tbin);  
       sprintf(htit,"MC phi Minv for %f < t <%f )",t1,t2);
       tminvcp[tbin]= new TH1F(hnam,htit,nbins,mbeg,mend);
       tminvcp[tbin]->Sumw2();

       sprintf(hnam,"tminvp%d",tbin);  
       sprintf(htit,"MC phi fit of Minv for %f < t <%f accidential corrected",t1,t2);
       tminvp[tbin]= new TH1F(hnam,htit,nbins,mbeg,mend);
       tminvp[tbin]->Sumw2();

// data

       sprintf(hnam,"dtminva%d",tbin);
       sprintf(htit,"Minv for %d < t <%d in time (J/#psi)",t1,t2);
       dtminva[tbin]= new TH1F(hnam,htit,nbins,mbeg,mend);
       dtminva[tbin]->Sumw2();

       sprintf(hnam,"dtminvb%d",tbin);
       sprintf(htit,"Minv for %f < t <%f out of time (J/#psi)",t1,t2);
       dtminvb[tbin]= new TH1F(hnam,htit,nbins,mbeg,mend);
       dtminvb[tbin]->Sumw2();

       sprintf(hnam,"dtminvc%d",tbin);
       sprintf(htit,"J/#psi Minv for %f < t <%f )",t1,t2);
       dtminvc[tbin]= new TH1F(hnam,htit,nbins,mbeg,mend);
       dtminvc[tbin]->Sumw2();

       sprintf(hnam,"dtminv%d",tbin);
       sprintf(htit,"J/#psi fit of Minv for %f < t <%f accidential corrected",t1,t2);
       dtminv[tbin]= new TH1F(hnam,htit,nbins,mbeg,mend);
       dtminv[tbin]->Sumw2();

       sprintf(hnam,"dtminvap%d",tbin);
       sprintf(htit,"Minv for %d < t <%d in time (phi)",t1,t2);
       dtminvap[tbin]= new TH1F(hnam,htit,nbins,mbeg,mend);
       dtminvap[tbin]->Sumw2();

       sprintf(hnam,"dtminvbp%d",tbin);
       sprintf(htit,"Minv for %f < t <%f out of time (phi)",t1,t2);
       dtminvbp[tbin]= new TH1F(hnam,htit,nbins,mbeg,mend);
       dtminvbp[tbin]->Sumw2();

       sprintf(hnam,"dtminvcp%d",tbin);
       sprintf(htit,"phi Minv for %f < t <%f )",t1,t2);
       dtminvcp[tbin]= new TH1F(hnam,htit,nbins,mbeg,mend);
       dtminvcp[tbin]->Sumw2();

       sprintf(hnam,"dtminvp%d",tbin);
       sprintf(htit,"phi fit of Minv for %f < t <%f accidential corrected",t1,t2);
       dtminvp[tbin]= new TH1F(hnam,htit,nbins,mbeg,mend);
       dtminvp[tbin]->Sumw2();

// combinations
       sprintf(hnam,"dtminva_ecomb%d",tbin);
       sprintf(htit,"Minv vs ecomb for %d < t <%d in time (J/#psi)",t1,t2);
       dtminva_ecomb[tbin]= new TH2F(hnam,htit,nbins,mbeg,mend,10,-0.5,9.5);

       sprintf(hnam,"dtminva_tcomb%d",tbin);
       sprintf(htit,"Minv vs tcomb for %d < t <%d in time (J/#psi)",t1,t2);
       dtminva_tcomb[tbin]= new TH2F(hnam,htit,nbins,mbeg,mend,10,-0.5,9.5);

       sprintf(hnam,"dtminvb_ecomb%d",tbin);
       sprintf(htit,"Minv vs ecomb for %d < t <%d out of time (J/#psi)",t1,t2);
       dtminvb_ecomb[tbin]= new TH2F(hnam,htit,nbins,mbeg,mend,10,-0.5,9.5);

       sprintf(hnam,"dtminvb_tcomb%d",tbin);
       sprintf(htit,"Minv vs tcomb for %d < t <%d out of time (J/#psi)",t1,t2);
       dtminvb_tcomb[tbin]= new TH2F(hnam,htit,nbins,mbeg,mend,10,-0.5,9.5);

}

      ncpe = new TH1F("ncpe","Number of combos per event",10,-0.5,9.5);
      iminv=0;
      iebeam=0;
      anecomb=0;
      bnecomb=0;
      antcomb=0;
      bntcomb=0;

      for (int i=0;i<100;i++){
        minv_old[i]=0.;
        dtrf_old[i]=0.;
        tbeam_old[i]=0.;
        trf_old[i]=0.;
        chi2_old[i]=0.;
        ebeam_old[i]=0.;
        t_old[i]=0.;
      }



}

Bool_t JP_selector::Process(Long64_t entry)
{
   // The Process() function is called for each entry in the tree (or possibly
   // keyed object in the case of PROOF) to be processed. The entry argument
   // specifies which entry in the currently loaded tree is to be processed.
   // It can be passed to either JP_selector::GetEntry() or TBranch::GetEntry()
   // to read either all or the required parts of the data. When processing
   // keyed objects with PROOF, the object is already loaded and is available
   // via the fObject pointer.
   //
   // This function should contain the "body" of the analysis. It can contain
   // simple or elaborate selection criteria, run algorithms on the data
   // of the event and typically fill histograms.
   //
   // The processing can be stopped by calling Abort().
   //
   // Use fStatus to set the return value of TTree::Process().
   //
   // The return value is currently not used.

    weight=-100.;
    bool data=true;
    bool new_event=false;
    GetEntry(entry);
    if(weight!=-100.){
       data=false;
    } 

//cokcok
//     Minv = Mrec_m;
//    pp=pp_m;
//    pep=pep_m;
//    pem=pem_m;
//    thp=thp_m;
//    thep=thep_m;
//    them=them_m;
//check
//    if(run_no>11553)return 0;
    
    if(run_no!=run_old){
      run_old=run_no;
      new_event=true;
      trf1=(2.*(run_no<11569||run_no>11663)+1.*(run_no>11568&&run_no<11664));
      trf2=trf1+(4.*(run_no<11569||run_no>11663)+2.*(run_no>11568&&run_no<11664))*0.;
      trf3=trf1+(4.*(run_no<11569||run_no>11663)+2.*(run_no>11568&&run_no<11664))*3.;
      //float trf3=trf1+(4.*(run_no<11569||run_no>11663)+2.*(run_no>11568&&run_no<11664))*1.;
      //final2017LI if(data&&run_no>30000)trf3=trf1+(4.*(run_no<11569||run_no>11663)+2.*(run_no>11568&&run_no<11664))*1.;
      //if(data&&run_no>30000)trf3=trf1+(4.*(run_no<11569||run_no>11663)+2.*(run_no>11568&&run_no<11664))*1.;
      //v5 float trfd=(1./2.*(run_no<11336)+1./6.*(run_no>11335));
      if(!data)fdedxem=1.;
      trfd=trf1/(trf3-trf2)*fdedxem;
      //cold float  trfd=trf1/(trf3-trf2);

      cout<<" run number ="<<run_no<<endl;
      if(run_no<19999){
       scale_bh = lumi16*xsec_bh/events_bh16;
       scale_jpsi = lumi16*xsec_jp*br_jp/events_jp16;
       //cout<<" sacle_bh16 "<<scale_bh<<endl;
       //cout<<" sacle_jpsi16 "<<scale_jpsi<<endl;
      } else if(run_no<39999) {
       scale_bh = lumi17*xsec_bh/events_bh17;
       scale_jpsi = lumi17*xsec_jp*br_jp/events_jp17;
       //cout<<" sacle_bh17 "<<scale_bh<<endl;
       //cout<<" sacle_jpsi17 "<<scale_jpsi<<endl;
      } else if(run_no<49999) {
       scale_bh = lumi18s*xsec_bh/events_bh18s;
       scale_jpsi = lumi18s*xsec_jp*br_jp/events_jp18s;
       //cout<<" sacle_bh18s "<<scale_bh<<endl;
       //cout<<" sacle_jpsi18s "<<scale_jpsi<<endl;
      } else if(run_no<59999) {
       scale_bh = lumi18f*xsec_bh/events_bh18f;
       scale_jpsi = lumi18f*xsec_jp*br_jp/events_jp18f;
       //cout<<" sacle_bh18f "<<scale_bh<<endl;
       //cout<<" sacle_jpsi18f "<<scale_jpsi<<endl;
      }
//cokcok
//scale_jpsi*=100.;
//scale_bh*=100.;
//scale_jpsi=1.;
//scale_bh=1.;

    } else if(event_no!=event_old){
      event_old=event_no;
      new_event=true;
    }
       if(new_event){
       ncount++;
       float www=weight_old;
       int naind[100];
       int nbind[100];
       antcomb=0;
       bntcomb=0;
       anecomb=1;
       bnecomb=1;
       if(iminv>0){
       for (int i=1;i<iminv;i++){
       if(abs(dtrf_old[i])<trf1&&abs(dtrf_old[i])<trf1){
       int j=0;
       bool brk=false;
       for (j=0;j<i;j++){
           if(abs(ebeam_old[i]-ebeam_old[j])<1.e-7&&abs(tbeam_old[i]-tbeam_old[j])<1.e-7){
             brk=true;
           }
       }
           if(brk){
           naind[antcomb]=i;
           antcomb++;
           } else {
            anecomb++;
           }
       } //end in-time
       if(abs(dtrf_old[i])>=trf2&&abs(dtrf_old[i])<trf3){
       int j=0;
       bool brk=false;
       for (j=0;j<i;j++){
           if(abs(ebeam_old[i]-ebeam_old[j])<1.e-7&&abs(tbeam_old[i]-tbeam_old[j])<1.e-7){
             brk=true;
           }
       }
           if(brk){
           nbind[bntcomb]=i;
           bntcomb++;
           } else {
            bnecomb++;
           }
       } //end out-time
       } //end i-loop
       for (int ik=0;ik<antcomb;ik++){
           if(data){
             dninva_tcomb->Fill(minv_old[naind[ik]],www);
           } else {
             ninva_tcomb->Fill(minv_old[naind[ik]],www);
           }
       }
       for (int ik=0;ik<bntcomb;ik++){
           if(data){
             dninvb_tcomb->Fill(minv_old[nbind[ik]],www);
           } else {
             ninvb_tcomb->Fill(minv_old[nbind[ik]],www);
           }
       }

      for (int i=0;i<iminv;i++){
        minv_old[i]=0.;
        ebeam_old[i]=0.;
        t_old[i]=0.;
        dtrf_old[i]=0.;
        tbeam_old[i]=0.;
        trf_old[i]=0.;
        chi2_old[i]=0.;
      }
      } //end iminv>0
      iminv=0;
      iebeam=0;
      anecomb=0;
      bnecomb=0;
      antcomb=0;
      bntcomb=0;

    } //end new event
  

      if(!data){
       if (weight>0.) {
       } else {
       }
      } else {
       if (run_no<11664){
       cdedxep*=3.12/2.52;
       cdedxem*=3.12/2.52;
       } else {
       }
       cdedxep*=1./1.74*(1+1.32/(1+abs(cos(thep))));
       cdedxem*=1./1.74*(1+1.32/(1+abs(cos(them))));
      }
      float Missep=sqrt(p_mass2+pi_mass2+2.*sqrt(pp_m*pp_m+p_mass2)*sqrt(pep_m*pep_m+pi_mass2)-2.*pp_m*pep_m*cos(thp_m+thep_m));
      float Missem=sqrt(p_mass2+pi_mass2+2.*sqrt(pp_m*pp_m+p_mass2)*sqrt(pem_m*pem_m+pi_mass2)-2.*pp_m*pem_m*cos(thp_m+them_m));

///*
//cold    bcut_data0=bpove_cut->Eval(Minv)+bpove_cuts->Eval(Minv)*nrsigmas;
//cold    bcut_data1=bpove_cut->Eval(Minv)-bpove_cuts->Eval(Minv)*nlsigmas;


// important cuts are here for TRD
    if(
       ((pep/bcalep<bcut_data0&&bcalep>0&&pbcalep*sin(thep)>pbcal_cut
       &&pep/bcalep>bcut_data1)
       ||(pep/fcalep<fcut_data0&&fcalep>0
        &&pep/fcalep>fcut_data1))
       &&((pem/bcalem<bcut_data0&&bcalem>0&&pbcalem*sin(them)>pbcal_cut
       &&pem/bcalem>bcut_data1)
       ||(pem/fcalem<fcut_data0&&fcalem>0
        &&pem/fcalem>fcut_data1))
     //&&(fcalep>0||fcalem>0)
     //&&(fcalep<=0&&fcalem<=0)
     //cokcok
      &&chi2>0&&chi2<chi2_cut
     //  &&chi2>0
      &&ebeam>8.2&&thep>theta_cut&&them>theta_cut
     // cokcok
     // &&ebeam>8.2&&thep<30.*PI/180.&&them<30.*PI/180.
     //&&ebeam>9.
     //&&ebeam<11.
  //   &&Mang>0.85
     &&abs(Missep-1.232)>delta_cut&&abs(Missem-1.232)>delta_cut&&pep>p_cut&&pem>p_cut
     &&pp>p_cut
    &&pep>p_cut&&pem>p_cut
  //cold  &&((cdedxep>dedx_cut&&thep*57.3>13.)||(thep*57.3<=13.))
  //cold   &&((cdedxem>dedx_cut&&them*57.3>13.)||(them*57.3<=13.))
     //&&((!data)||(cdedxep/1.74*(1+1.32/(1+abs(cos(thep))))>dedx_cut&&thep*57.3>13.)||(thep*57.3<=13.))
     //&&((!data)||(cdedxem/1.74*(1+1.32/(1+abs(cos(thep))))>dedx_cut&&thep*57.3>13.)||(thep*57.3<=13.))
   //  &&(abs(ebeam+0.938-sqrt(pp_m*pp_m+0.88)-pep_m-pem_m)<3.)
      &&((abs(trf-tbeam)<trf1)||(abs(trf-tbeam)>=trf2&&abs(trf-tbeam)<trf3))
     // &&(abs(t-tmin)<0.6)
     //  &&(pol==-1)
    // &&(abs(Theta-90.)<theta_cm_cut)
     //&&((data)||(weight<0.)||(weight<1000.))
   //  &&(ebeam>8.25&&ebeam<11.85)
     //&&(Minv>3.033&&Minv<3.153)
     //&&run_no>71727&&run_no<71856
     //&&run_no>71862&&run_no<71943
     //&&run_no>71942&&run_no<72068
     //&&run_no>72361&&run_no<72436
    // &&run_no>51203&&run_no<51384
    // &&run_no>51496&&run_no<51769
    // &&(run_no<11071||run_no>11346) //CDC hole excluded
    //   &&thep*57.3<40.&&them*57.3<40.
    //cokcok
      // &&abs(t-tmin)<0.6
      ) {
//*/
//cokcok 
// { 
      float www=1.;



    bool new_ecomb=true;
    bool new_comb=true;
    bool new_tcomb=false;
    for (int i=0;i<iminv;i++){
      if(abs(ebeam-ebeam_old[i])<1.e-7)new_ecomb=false; 
    }
    if(iminv==0)new_ecomb=false;
    for (int i=0;i<iminv;i++){
      //if(abs(Minv-minv_old[i])<1.e-7||(
      if(abs(ebeam-ebeam_old[i])<1.e-7
                                  //     (abs(ebeam-ebeam_old[i])>1.e-7)&&
                                  //     ((abs(trf-tbeam)<trf1&&abs(dtrf_old[i])<trf1)
                                  //     ||(abs(trf-tbeam)>=trf2&&abs(trf-tbeam)<trf3&&abs(dtrf_old[i])>=trf2&&abs(dtrf_old[i])<trf3)
                                  //    )
                                  //    )
        )new_tcomb=true; 
    }
    for (int i=0;i<iminv;i++){
      //if(abs(Minv-minv_old[i])<1.e-7&&abs(ebeam-ebeam_old[i])<1.e-7)new_comb=false; 
      if(abs(Minv-minv_old[i])<1.e-7)new_comb=false; 
    }
    
    new_comb=true;
    if(new_comb){                 //same combinations excluded
        minv_old[iminv]=Minv;
        dtrf_old[iminv]=trf-tbeam;
        tbeam_old[iminv]=tbeam;
        chi2_old[iminv]=chi2;
        trf_old[iminv]=trf;
        ebeam_old[iminv]=ebeam;
        t_old[iminv]=-(t-tmin);
        event_no_old[iminv]=event_no;
        iminv++;
    }
     
    if(!new_comb){                 //same combinations excluded
       /* cok
       cout<<" b i n g o o o o   ev_no="<<event_no<<endl;
       cout<<" Minv, ebeam="<<Minv<<" "<<ebeam<<endl;
       cout<<" old_minv=";
       for (int i=0;i<iminv;i++){
       cout<<" minv_old="<<minv_old[i]-Minv;
       }
       cout<<endl;
       cout<<" old_minv=";
       for (int i=0;i<iminv;i++){
       cout<<" ebeam_old="<<ebeam_old[i]-ebeam;
       }
       cout<<endl;
       int inext;
       //cin>>inext;
       */
     } else {                       //all fills done for new combinations only:
       ncombos++;
        
/*
        if(new_event){
        if(abs(trf-tbeam)<trf1){
          anecomb=1;
          antcomb=1;
        } else if (abs(trf-tbeam)>=trf2&&abs(trf-tbeam)<trf3) {
          bnecomb=1;
          bntcomb=1;
        }
        }
*/
        
/*
        if(abs(trf-tbeam)<trf1){
          iebeam++;
          if(new_ecomb)anecomb++;
          if(new_tcomb)antcomb++;
        } else if (abs(trf-tbeam)>=trf2&&abs(trf-tbeam)<trf3) {
          iebeam++;
          if(new_ecomb)bnecomb++;
          if(new_tcomb)bntcomb++;
        }
*/
//if(event_no==5288556&&run_no==10392){
//cout<<"new_ecomb,new_tcomb"<<new_ecomb<<" "<<new_tcomb<<endl;

          //cout<<" anecomb,antcomb,bnecomb,bntcomb="<< anecomb<<" "<<antcomb<<" "<<bnecomb<<" "<<bntcomb<<endl;
          //for (int i=0;i<iminv;i++){
          //cout<<" minv,ebeam,dtrf,trf="<<minv_old[i]<<" "<<ebeam_old[i]<<" "<<dtrf_old[i]<<" "<<trf_old[i]<<endl;
          //}
          //int inext;
          //cin>>inext;

//}

        


        //Minv=Minv*Minv; //Minv^2
        if(!data){
        
        //trf1=2.;
        //trf2=2.;
        //trf3=14.;
        //trfd=trf1/(trf3-trf2);
        if(weight==-2.){
// phi case
          www=scale_phi;
       
        } else if (weight==-3.) {
// J/psi case
          www=scale_jpsi;
          //cokcok -\/ uncomment just to evaluate the J/psi MC flux
          //www=scale_jpsi/mcxsec->Eval(ebeam);
          //??? www*=scale_jpsi*exp(1.66*abs(t));
        } else if (weight>0.) {
// BH case
          www=weight*scale_bh;
          float bhbg=1.;
          if(Minv>1.2&&Minv<3.333){
          bhbg=0.;
          //int mbin=xbhb->FindBin(Minv);
          //if(Minv>3.126)mbin=7;
          //if(bcalep/pep>0.5)bhbg+=1./bhbg_factorb->GetBinContent(mbin);
          //if(bcalem/pem>0.5)bhbg+=1./bhbg_factorb->GetBinContent(mbin);
          //if(fcalep/pep>0.5)bhbg+=1./bhbg_factorf->GetBinContent(mbin);
          //if(fcalem/pem>0.5)bhbg+=1./bhbg_factorf->GetBinContent(mbin);

          if(bcalep/pep>0.5)bhbg+=1./bcalbg->Eval(Minv);
          if(bcalem/pem>0.5)bhbg+=1./bcalbg->Eval(Minv);
          if(fcalep/pep>0.5)bhbg+=1./fcalbg->Eval(Minv);
          if(fcalem/pem>0.5)bhbg+=1./fcalbg->Eval(Minv);

          //if(fcalep/pep>0.5)bhbg+=1./bcalbg->Eval(Minv);
          //if(fcalem/pem>0.5)bhbg+=1./bcalbg->Eval(Minv);
          //if(Minv>3.126)bhbg*=2.;
          bhbg/=2.;
          }
          //cokcok
          //bhbg=1.;
          www=www*bhbg;
        } else {
         www=0.;
        }

        } else {
          www=1.;
          
        }
        double tmn=abs(ftmin(Minv,ebeam));
        double tmx=abs(ftmax(Minv,ebeam));
        if(tmn<etmin)tmn=etmin;
        if(tmx>etmax)tmx=etmax;
        //cokcok
        //www=www/abs(tmx-tmn);
        //www=www/myflux->GetBinContent(fax->FindBin(ebeam));
        tweight=1./abs(tmx-tmn);
        fweight=1./myflux->GetBinContent(fax->FindBin(ebeam));
        prompt=0.;
        if(abs(trf-tbeam)<trf1)prompt=1.;
        //if(-(t)>etmin&&-(t)<etmax&&ebeam>8.56&&ebeam<8.92)
        if(-(t)>etmin&&-(t)<etmax)
        ntuple->Fill(Minv,prompt,t,ebeam,fweight,tweight); 

        double deweight,dtint,eweight,tint,dtweight,tweight;
        if(data){
        deweight=1.;
        dtint=smgl/3.*(1./pow(1.+tmn/smgl,3)-1./pow(1+tmx/smgl,3));
        //deweight/=dtint;
        dtweight=deweight;
        } else {
        eweight=1.;
        //tint=1./tslope*exp(tmn*tslope)*(exp(-tmn*tslope)-exp(-tmx*tslope));
        tint=1./tslope*(exp(-tmn*tslope)-exp(-tmx*tslope));
        tweight=eweight;
        //eweight/=tint;
        }


    if(abs(trf-tbeam)<2.&&(data)&&Minv>3.063&&Minv<3.123)ethist->Fill(ebeam,-t);

        if(Minv>3.&&Minv<3.2&&(!data)){
        hweight->Fill(ebeam,eweight);
        htint->Fill(ebeam,tint); 
        }
        if(abs(trf-tbeam)<trf1){
           if(data){
           dminvra->Fill(Mrec_m,www);
           //cok dminva->Fill(Minv,www);
           dminva->Fill(Minv+ecorr(pp_m,thp_m,ebeam),www);
           dminvap->Fill(Minv,www);
           //dminva->Fill(Mrec_m,www);
           //dminvap->Fill(Mang,www);
           } else {
           minvr->Fill(Mrec_m,www);
           minva->Fill(Minv,www);
           minvap->Fill(Minv,www);
           }
           for (int ebin=0;ebin<nebins;ebin++){
             float e1=ebeg+ebin*estep;
             float e2=e1+estep;
             //cok if(ebeam>e1&&ebeam<=e2&&-(t-tmin)>etmin&&-(t-tmin)<etmax)
             if(ebeam>e1&&ebeam<=e2&&-(t-tmin)>etmin&&-(t-tmin)<etmax) {
             //if(ebeam>e1&&ebeam<=e2&&-(t)>etmin&&-(t)<etmax){
               if(data){
               deminva[ebin]->Fill(Minv,www/deweight);
               deminvap[ebin]->Fill(Minv,www);
               } else {
               eminva[ebin]->Fill(Minv,www/eweight);
               eminvap[ebin]->Fill(Minv,www);
               }
             }
           }
           for (int tbin=0;tbin<ntbins;tbin++){
             float t1=tbeg+tbin*tstep;
             float t2=t1+tstep;
             //cok if(-(t-tmin)>t1&&-(t-tmin)<=t2&&ebeam>temin&&ebeam<temax)
             if(-(t)>t1&&-(t)<=t2&&ebeam>temin&&ebeam<temax){
             //test if(thp*57.3>t1&&thp*57.3<=t2)
               if(data){
               dtminva[tbin]->Fill(Minv,www);
               dtminvap[tbin]->Fill(Minv,www);

               //dtminva[tbin]->Fill(Minv,www/dtweight);
               //dtminvap[tbin]->Fill(Minv,www);
               } else {
               tminva[tbin]->Fill(Minv,www);
               tminvap[tbin]->Fill(Minv,www);
               //tminva[tbin]->Fill(Minv,www/tweight);
               //tminvap[tbin]->Fill(Minv,www);
               }
             }
           }

        } else if (abs(trf-tbeam)>=trf2&&abs(trf-tbeam)<trf3) {
           if(data){
           dminvrb->Fill(Mrec_m,trfd*www);
           dminvb->Fill(Minv+ecorr(pp_m,thp_m,ebeam),trfd*www);
           //cok dminvb->Fill(Minv,trfd*www);
           dminvbp->Fill(Minv,trfd*www);
           } else {
           minvb->Fill(Minv,trfd*www);
           minvbp->Fill(Minv,trfd*www);
           }
           for (int ebin=0;ebin<nebins;ebin++){
             float e1=ebeg+ebin*estep;
             float e2=e1+estep;
             //cok if(ebeam>e1&&ebeam<=e2&&-(t-tmin)>etmin&&-(t-tmin)<etmax)
             if(ebeam>e1&&ebeam<=e2&&-(t-tmin)>etmin&&-(t-tmin)<etmax) {
             //if(ebeam>e1&&ebeam<=e2&&-(t)>etmin&&-(t)<etmax){
               if(data){
               deminvb[ebin]->Fill(Minv,trfd*www/deweight);
               deminvbp[ebin]->Fill(Minv,trfd*www);
               } else {
               eminvb[ebin]->Fill(Minv,trfd*www/eweight);
               eminvbp[ebin]->Fill(Minv,trfd*www);
               }
             }
           }
           for (int tbin=0;tbin<ntbins;tbin++){
             float t1=tbeg+tbin*tstep;
             float t2=t1+tstep;
             //cok if(-(t-tmin)>t1&&-(t-tmin)<=t2&&ebeam>temin&&ebeam<temax)
             if(-(t)>t1&&-(t)<=t2&&ebeam>temin&&ebeam<temax){
             //test if(thp*57.3>t1&&thp*57.3<=t2)
               if(data){
               dtminvb[tbin]->Fill(Minv,trfd*www);
               dtminvbp[tbin]->Fill(Minv,trfd*www);
               //dtminvb[tbin]->Fill(Minv,trfd*www/dtweight);
               //dtminvbp[tbin]->Fill(Minv,trfd*www);
               } else {
               tminvb[tbin]->Fill(Minv,trfd*www);
               tminvbp[tbin]->Fill(Minv,trfd*www);
               //tminvb[tbin]->Fill(Minv,trfd*www/tweight);
               //tminvbp[tbin]->Fill(Minv,trfd*www);
               }
             }
           }

        } //end else if


    } //end old Minv old ebeam

    } // E N D   B I G    I F


//    if(((bcalep/pep>0.5&&pbcalep/pep>(1-bcalep/pep*1.1))||fcalep/pep>0.8)
//     &&((bcalem/pem>0.5&&pbcalem/pem>(1-bcalem/pem*1.1))||fcalem/pem>0.8)
//     &&ebeam>8.2&&chi2>0&&chi2<50&&Minv>0.5&&abs(tbeam-trf)<2.&&Minv<3.5) {

    if((bcalep/pep_m>0.7||fcalep/pep_m>0.7)
     &&(bcalem/pem_m>0.7||fcalem/pem_m>0.7)
     &&((bcalep/pep>0.7&&pbcalep/pep>(1-bcalep/pep))||fcalep/pep>0.7)
     &&((bcalem/pem>0.7&&pbcalem/pem>(1-bcalem/pem))||fcalem/pem>0.7)
     &&(abs(thep-thep_m)<0.4&&abs((phep-phep_m)*sin(thep))<0.2)
     &&(abs(them-them_m)<0.4&&abs((phem-phem_m)*sin(them))<0.2)
     &&(abs(pep-pep_m)<5.&&abs(pem-pem_m)<5.)
     &&(Minv>(0.857*Mrec_m-0.0855))
     &&(Minv>0.8&&Mrec_m>0.8&&ebeam>8.2&&chi2>-10&&abs(tbeam-trf)<2.)) {
     if(data){
     dminv_nokf->Fill(Minv);
     } else {
     minv_nokf->Fill(Minv);
     }
    } // E N D   B I G    I F


   return kTRUE;
}

void JP_selector::SlaveTerminate()
{
   // The SlaveTerminate() function is called after all entries or objects
   // have been processed. When running with PROOF SlaveTerminate() is called
   // on each slave server.


}

void JP_selector::Terminate()
{
   // The Terminate() function is the last function to be called during
   // a query. It always runs on the client, it can be used to present
   // the results graphically or save the results to file.

cout<<" total number of events,combos,acombos ====== "<<ncount<<" "<<ncombos<<" "<<acombos<<endl;

bool fitbg; // itrue - fit background hist and subtract, false - just subtract the histogram
double sum,sigma,bg;
float estep=(eend-ebeg)/nebins;
TH1F *debh = new TH1F("debh","",nebins,ebeg,eend);
TH1F *dejpsi = new TH1F("dejpsi","",nebins,ebeg,eend);
TH1F *dephi = new TH1F("dephi","",nebins,ebeg,eend);
TH1F *dtbh = new TH1F("dtbh","",ntbins,tbeg,tend);
TH1F *dtjpsi = new TH1F("dtjpsi","",ntbins,tbeg,tend);
TH1F *dtphi = new TH1F("dtphi","",ntbins,tbeg,tend);

TH1F *ebh = new TH1F("ebh","",nebins,ebeg,eend);
TH1F *ejpsi = new TH1F("ejpsi","",nebins,ebeg,eend);
TH1F *ephi = new TH1F("ephi","",nebins,ebeg,eend);
TH1F *tbh = new TH1F("tbh","",ntbins,tbeg,tend);
TH1F *tjpsi = new TH1F("tjpsi","",ntbins,tbeg,tend);
TH1F *tphi = new TH1F("tphi","",ntbins,tbeg,tend);
int mcflag;

// MC Minv widths
// 10 bins
double jpsigma=1.53822e-02;
//float mcwidth[10]={0.0056276429,0.0097527122,0.010375901,0.011995159,0.013918157,0.014618722,0.014765249,0.016105106,0.017058151,0.014553332};
float mcwidth[20]={0.00361951,0.00681605,0.00938171,0.00997647,0.0101224,0.0103699,0.0113894,0.0123979,0.0130037,0.0143896,0.0144498,0.0140419,0.0145275,0.0145317,0.0155685,0.01609,0.0164124,0.0165902,0.0139572,0.0149109};

TAxis *xaxis = dminv->GetXaxis();
// data all
//cok fitbg=false; //do not fit the background, just subtract it, integrate the Gaussian to estimate events
fitbg=true; //fit the background and subtract the fit, sum events within 3 sigma (sigma fixed now to 10 MeV)
double herror;
cout<<"[][][][][][][][] data all BH: [][][][][][][][] "<<endl;
dminv->Add(dminva,dminvb,1.,-1.);
float bhintall=dminv->IntegralAndError(xaxis->FindBin(bh1),xaxis->FindBin(bh2),herror,"");
cout<<"BH 1.5 - 2.5 GeV = "<<bhintall<<" +/- "<<herror<<endl;
cout<<endl;
cout<<"[][][][][][][][] Rec. data all J/psi: [][][][][][][][] "<<endl;
mcflag=0;
//FitJpsi(dminvra,dminvrb,dminvrc,dminvr,sum,sigma,bg,fitbg,mcflag,jpsigma);
cout<<"[][][][][][][][] data all J/psi: [][][][][][][][] "<<endl;
//FitJpsi(dminva,dminvb,dminvc,dminv,sum,sigma,bg,fitbg,mcflag,jpsigma);
cout<<"[][][][][][][][] data all phi: [][][][][][][][] "<<endl;
//FitPhi(dminvap,dminvbp,dminvcp,dminvp,sum,sigma,bg,fitbg);

// MC all
cout<<"[][][][][][][][] MC all BH: [][][][][][][][] "<<endl;
minv->Add(minva,minvb,1.,-1.);
bhintall=minv->IntegralAndError(xaxis->FindBin(bh1),xaxis->FindBin(bh2),herror,"");
cout<<"BH 1.5 - 2.5 GeV = "<<bhintall<<" +/- "<<herror<<endl;
cout<<endl;
cout<<"[][][][][][][][] MC all J/psi: [][][][][][][][] "<<endl;
mcflag=1;
//FitJpsi(minva,minvb,minvc,minv,sum,sigma,bg,fitbg,mcflag,jpsigma);
cout<<"[][][][][][][][] MC all phi: [][][][][][][][] "<<endl;
//FitPhi(minvap,minvbp,minvcp,minvp,sum,sigma,bg,fitbg);

TH1F *rminv = (TH1F*) dminv->Clone("rminv");
rminv->Divide(minv);


cout<<endl;
cout<<endl;
fitbg=true; //fit the background and subtract the fit, sum events within 3 sigma (sigma fixed now to 10 MeV)
//cokcok   for (int ebin=0;ebin<nebins;ebin++){
   for (int ebin=0;ebin<nebins;ebin++){
cout<<"[][][][][][][][] data E-slices bin "<<ebin<<" BH [][][][][][][][] "<<endl;
deminv[ebin]->Add(deminva[ebin],deminvb[ebin],1.,-1.);
float bhint=deminv[ebin]->IntegralAndError(xaxis->FindBin(bh1),xaxis->FindBin(bh2),herror,"");
cout<<"BH 1.5 - 2.5 GeV = "<<bhint<<" +/- "<<herror<<endl;
debh->SetBinContent(ebin+1,bhint);
debh->SetBinError(ebin+1,herror);
cout<<endl;
cout<<"[][][][][][][][] data E-slices bin "<<ebin<<" J/psi [][][][][][][][] "<<endl;
mcflag=0;
jpsigma=mcwidth[ebin];
//FitJpsi(deminva[ebin],deminvb[ebin],deminvc[ebin],deminv[ebin],sum,sigma,bg,fitbg,mcflag,jpsigma);
dejpsi->SetBinContent(ebin+1,abs(sum));
dejpsi->SetBinError(ebin+1,sigma);
   } //end ebin loop

cout<<endl;
//fitbg=false; //do not fit the background, just subtract it, integrate the Gaussian to estimate events
   //cok for (int ebin=0;ebin<nebins;ebin++){
   for (int ebin=nebins-1;ebin>=0;ebin--){
cout<<"[][][][][][][][] data E-slices bin "<<ebin<<" phi [][][][][][][][] "<<endl;
//FitPhi(deminvap[ebin],deminvbp[ebin],deminvcp[ebin],deminvp[ebin],sum,sigma,bg,fitbg);
dephi->SetBinContent(ebin+1,abs(sum));
dephi->SetBinError(ebin+1,sigma);
   } //end ebin loop


cout<<endl;
fitbg=true; //fit the background and subtract the fit, sum events within 3 sigma (sigma fixed now to 10 MeV)
   for (int tbin=0;tbin<ntbins;tbin++){
cout<<"[][][][][][][][] data t-slices bin "<<tbin<<" J/psi [][][][][][][][] "<<endl;
dtminv[tbin]->Add(dtminva[tbin],dtminvb[tbin],1.,-1.);
float bhint=dtminv[tbin]->IntegralAndError(xaxis->FindBin(bh1),xaxis->FindBin(bh2),herror,"");
dtbh->SetBinContent(tbin+1,bhint);
dtbh->SetBinError(tbin+1,herror);
cout<<"BH 1.5 - 2.5 GeV = "<<bhint<<" +/- "<<herror<<endl;
mcflag=0;
//FitJpsi(dtminva[tbin],dtminvb[tbin],dtminvc[tbin],dtminv[tbin],sum,sigma,bg,fitbg,mcflag,jpsigma);
dtjpsi->SetBinContent(tbin+1,abs(sum));
dtjpsi->SetBinError(tbin+1,sigma);
   } //end tbin loop

cout<<endl;
//fitbg=false; //do not fit the background, just subtract it, integrate the Gaussian to estimate events
   for (int tbin=0;tbin<ntbins;tbin++){
cout<<"[][][][][][][][] data t-slices bin "<<tbin<<" phi [][][][][][][][] "<<endl;
//FitPhi(dtminvap[tbin],dtminvbp[tbin],dtminvcp[tbin],dtminvp[tbin],sum,sigma,bg,fitbg);
dtphi->SetBinContent(tbin+1,abs(sum));
dtphi->SetBinError(tbin+1,sigma);
   } //end tbin loop

cout<<endl;
cout<<endl;
// MC
fitbg=true; //fit the background and subtract the fit, sum events within 3 sigma (sigma fixed now to 10 MeV)
   for (int ebin=0;ebin<nebins;ebin++){
cout<<"[][][][][][][][] MC E-slices bin "<<ebin<<" BH [][][][][][][][] "<<endl;
eminv[ebin]->Add(eminva[ebin],eminvb[ebin],1.,-1.);
double bhint=eminv[ebin]->IntegralAndError(xaxis->FindBin(bh1),xaxis->FindBin(bh2),herror,"");
cout<<"BH 1.5 - 2.5 GeV = "<<bhint<<" +/- "<<herror<<endl;
ebh->SetBinContent(ebin+1,bhint);
ebh->SetBinError(ebin+1,herror);
cout<<endl;
cout<<"[][][][][][][][] MC E-slices bin "<<ebin<<" J/psi [][][][][][][][] "<<endl;
mcflag=1;
jpsigma=mcwidth[ebin];
//FitJpsi(eminva[ebin],eminvb[ebin],eminvc[ebin],eminv[ebin],sum,sigma,bg,fitbg,mcflag,jpsigma);
ejpsi->SetBinContent(ebin+1,abs(sum));
ejpsi->SetBinError(ebin+1,sigma);
   } //end ebin loop

cout<<endl;
//fitbg=false; //do not fit the background, just subtract it, integrate the Gaussian to estimate events
   for (int ebin=0;ebin<nebins;ebin++){
cout<<"[][][][][][][][] MC E-slices bin "<<ebin<<" phi [][][][][][][][] "<<endl;
//FitPhi(eminvap[ebin],eminvbp[ebin],eminvcp[ebin],eminvp[ebin],sum,sigma,bg,fitbg);
ephi->SetBinContent(ebin+1,abs(sum));
ephi->SetBinError(ebin+1,sigma);
   } //end ebin loop

cout<<endl;
fitbg=true; //fit the background and subtract the fit, sum events within 3 sigma (sigma fixed now to 10 MeV)
   for (int tbin=0;tbin<ntbins;tbin++){
cout<<"[][][][][][][][] MC t-slices bin "<<tbin<<" J/psi [][][][][][][][] "<<endl;
tminv[tbin]->Add(tminva[tbin],tminvb[tbin],1.,-1.);
float bhint=tminv[tbin]->IntegralAndError(xaxis->FindBin(bh1),xaxis->FindBin(bh2),herror,"");
tbh->SetBinContent(tbin+1,bhint);
tbh->SetBinError(tbin+1,herror);
cout<<"BH 1.5 - 2.5 GeV = "<<bhint<<" +/- "<<herror<<endl;
mcflag=1;
//FitJpsi(tminva[tbin],tminvb[tbin],tminvc[tbin],tminv[tbin],sum,sigma,bg,fitbg,mcflag,jpsigma);
tjpsi->SetBinContent(tbin+1,abs(sum));
tjpsi->SetBinError(tbin+1,sigma);
   } //end tbin loop

cout<<endl;
//fitbg=false; //do not fit the background, just subtract it, integrate the Gaussian to estimate events
   for (int tbin=0;tbin<ntbins;tbin++){
cout<<"[][][][][][][][] MC t-slices bin "<<tbin<<" phi [][][][][][][][] "<<endl;
//FitPhi(tminvap[tbin],tminvbp[tbin],tminvcp[tbin],tminvp[tbin],sum,sigma,bg,fitbg);
tphi->SetBinContent(tbin+1,abs(sum));
tphi->SetBinError(tbin+1,sigma);
   } //end tbin loop

ntuple->Write();
fout->Write();

}

void JP_selector::FitPhi(TH1F*aa, TH1F*bb, TH1F*cc, TH1F*vv, double& ssum, double& ssigma, double& bbg, bool fitbg)
{

double bin=3./nbins;

double m1=0.95;
m1=0.9;
double m2=1.25;
//m2=1.30;
double mnv=1.020;
double sigma=0.0106;
double p0,p1,p2,p3,p4,p5;
double err;
expp = new TF1("expp","exp([0]+[1]*x)",m1,m2);
//cold expp = new TF1("expp","[0]+[1]*x",m1,m2);
ffit = new TF1("ffit","exp([0]+[1]*x)+[2]*exp(-0.5*(x-[3])*(x-[3])/([4]*[4]))",m1,m2);
//cold ffit = new TF1("ffit","([0]+[1]*x)+[2]*exp(-0.5*(x-[3])*(x-[3])/([4]*[4]))",m1,m2);

       
   // phi fit

bb->Fit("expp","qr");
for (int i=1;i<nbins+1;i++){
  double ma=aa->GetBinContent(i);
  double ema=aa->GetBinError(i);
  double bin_cent=aa->GetBinCenter(i);
  double mb=expp->Eval(bin_cent);
  vv->SetBinContent(i,ma-mb);
  vv->SetBinError(i,ema);
}
//cold if(!fitbg)vv->Add(aa,bb,1.,-1.);
vv->Add(aa,bb,1.,-1.);
vv->Fit("expp","qr");
//gMinuit->GetParameter(0,p0,err);
//gMinuit->GetParameter(1,p1,err);
p0=expp->GetParameter(0);
p1=expp->GetParameter(1);
ffit->SetParameters(p0,p1,10.,mnv,sigma);
//ffit->FixParameter(4,sigma);
vv->Fit("ffit","qr");
//gMinuit->GetParameter(0,p0,err);
//gMinuit->GetParameter(1,p1,err);
p0=ffit->GetParameter(0);
p1=ffit->GetParameter(1);
expp->SetParameters(p0,p1);
for (int i=1;i<nbins+1;i++){
  double ma=vv->GetBinContent(i);
  double ema=vv->GetBinError(i);
  double bin_cent=vv->GetBinCenter(i);
  double mb=expp->Eval(bin_cent);
  cc->SetBinContent(i,ma-mb);
  cc->SetBinError(i,ema);
}

double m[5][5];
double a;
double s;
double e;
double ae; double se; double sae;

//gMinuit->mnemat(&m[0][0],5);
//ae=m[2][2];
//se=m[4][4];
//sae=m[4][2];
TVirtualFitter *fitter = TVirtualFitter::GetFitter();
TMatrixD matrix(5,5,fitter->GetCovarianceMatrix());
ae=matrix[2][2];
se=matrix[4][4];
sae=matrix[4][2];
//gMinuit->GetParameter(2,a,e);
//gMinuit->GetParameter(4,s,e);
a=ffit->GetParameter(2);
s=ffit->GetParameter(4);
double summ=sqrt(2.*PI)/bin*s*a;
err=sqrt(2.*PI)/bin*sqrt(s*s*ae+a*a*se+2.*a*s*sae);


double cexp; double pexp; double mean;
//gMinuit->GetParameter(0,cexp,e);
//gMinuit->GetParameter(1,pexp,e);
//gMinuit->GetParameter(3,mean,e);
cexp=ffit->GetParameter(0);
pexp=ffit->GetParameter(1);
mean=ffit->GetParameter(3);

//if(!fitbg)cout<<" "<<summ<<" +/- "<<err<<endl;
double alt_err;
if(fitbg){
TAxis *xaxis = cc->GetXaxis();
s=0.0102;
mean=1.016;
int bin1=xaxis->FindBin(mean-3.*s);
int bin2=xaxis->FindBin(mean+3.*s);
//cout<<" --bin1--bin2-- "<<bin1<<" "<<bin2<<endl;
//cold summ=cc->Integral(bin1,bin2);
//cold summ=cc->IntegralAndError(bin1,bin2,alt_err,"");
}

//if(!fitbg){
cout<<" "<<summ<<" +/- "<<err<<endl;
cout<<" background= "<<exp(cexp+pexp*mean)<<endl;
//}
ssum=summ;
ssigma=err;
bbg=exp(cexp+pexp*mean);

}

void JP_selector::FitJpsi(TH1F*aa, TH1F*bb, TH1F*cc, TH1F*vv,double &ssum, double &ssigma, double &bbg, bool fitbg, int mcflag, double &jpsigma)
{

double bin=3./nbins;
cout<<" biiiiiiiiinnnnnnnnnnnnnnnnnnnnnnnnnnnn="<<bin<<endl;
double m1=2.95;
//cok m1=2.9;
m1=2.95;
double m2=3.2;
//cok m2=3.2;
m2=3.2;
//double m2=3.15;
//m2=3.3;
double mnv=3.096;
double sigma=0.0123;
sigma=jpsigma;
//if(mcflag==1){
//mnv=3.096;
//sigma=0.0137;
//}

double p0,p1,p2,p3,p4,p5;
double summ,err;
double m[5][5];
double a;
double s;
double e;
double ae; double se; double sae;
double cexp; double pexp; double mean;

expp = new TF1("expp","[0]+[1]*x",m1,m2);
ffit = new TF1("ffit","[0]+[1]*x+[2]*exp(-0.5*(x-[3])*(x-[3])/([4]*[4]))",m1,m2);


       
bb->Fit("expp","qr");
for (int i=1;i<nbins+1;i++){
  double ma=aa->GetBinContent(i);
  double ema=aa->GetBinError(i);
  double bin_cent=aa->GetBinCenter(i);
  double mb=expp->Eval(bin_cent);
  // important change -> method IV -> Fit the total (not acc. subtracted) distribution
  vv->SetBinContent(i,ma-mb);
  //cok for method IV  vv->SetBinContent(i,ma);
  //cokcok vv->SetBinContent(i,ma);
  vv->SetBinContent(i,ma);
  vv->SetBinError(i,ema);
}
//cold if(!fitbg)vv->Add(aa,bb,1.,-1.);
//
//check this:
//cokcok vv->Add(aa,bb,1.,-1.);
vv->Add(aa,bb,1.,-1.);

if(fitbg){
vv->Fit("expp","qr");
} else {
vv->Fit("expp","qr");
}
//gMinuit->GetParameter(0,p0,err);
//gMinuit->GetParameter(1,p1,err);
p0=expp->GetParameter(0);
p1=expp->GetParameter(1);
ffit->SetParameters(p0,p1,10.,mnv,sigma);
//if(fitbg)ffit->FixParameter(4,sigma);
//cokcok
//ffit->FixParameter(4,sigma);
if(mcflag==1){
vv->Fit("ffit","r");
} else {
//cokcok vv->Fit("ffit","rl");
vv->Fit("ffit","rl");
}
p0=ffit->GetParameter(0);
p1=ffit->GetParameter(1);
//check p0=expp->GetParameter(0);
//check p1=expp->GetParameter(1);
expp->SetParameters(p0,p1);
for (int i=1;i<nbins+1;i++){
  double ma=vv->GetBinContent(i);
  double ema=vv->GetBinError(i);
  double bin_cent=vv->GetBinCenter(i);
  double mb=expp->Eval(bin_cent);
  cc->SetBinContent(i,ma-mb);
  cc->SetBinError(i,ema);
}


//gMinuit->mnemat(&m[0][0],5);
//ae=m[2][2];
//se=m[4][4];
//sae=m[4][2];
TVirtualFitter *fitter = TVirtualFitter::GetFitter();
TMatrixD matrix(5,5,fitter->GetCovarianceMatrix());
ae=matrix[2][2];
se=matrix[4][4];
sae=matrix[4][2];
//gMinuit->GetParameter(2,a,e);
//gMinuit->GetParameter(4,s,e);
a=ffit->GetParameter(2);
mean=ffit->GetParameter(3);
s=ffit->GetParameter(4);
cout<<" biiiiiiiiinnnnnnnnnnnnnnnnnnnnnnnnnnnn="<<bin<<endl;
summ=sqrt(2.*PI)/bin*s*a;
err=sqrt(2.*PI)/bin*sqrt(s*s*ae+a*a*se+2.*a*s*sae);

//gMinuit->GetParameter(0,cexp,e);
//gMinuit->GetParameter(1,pexp,e);
//gMinuit->GetParameter(3,mean,e);
cexp=ffit->GetParameter(0);
pexp=ffit->GetParameter(1);
mean=ffit->GetParameter(3);
//if(!fitbg)cout<<" "<<summ<<" +/- "<<err<<endl;
if(fitbg){
TAxis *xaxis = cc->GetXaxis();
//cok? fixed s=sigma;
//cok? fixed mean=mnv;
s=sigma;
//cold mean=mnv;
int bin1=xaxis->FindBin(mean-3.*s);
int bin2=xaxis->FindBin(mean+3.*s);
cout<<" >>>>>>>>>>>>>>>>>>>>>>>> mean sigma --bin1--bin2-- "<<mean<<" "<<sigma<<" "<<bin1<<" "<<bin2<<endl;
double alt_err;
//cok summ=cc->IntegralAndError(bin1,bin2,alt_err,"");
summ=cc->IntegralAndError(bin1,bin2,alt_err,"");
//cok??? err=sqrt(summ);
err=sqrt(summ);

cout<<" "<<summ<<" +/- "<<err<<endl;
cout<<" "<<summ<<" +/- "<<alt_err<<endl;
}
cout<<" background= "<<(cexp+pexp*mean)<<endl;
ssum=summ;
ssigma=err;
//cok ssigma=sqrt(ssum);
//ssigma=ssum/10.;
bbg=cexp+pexp*mean;


}
