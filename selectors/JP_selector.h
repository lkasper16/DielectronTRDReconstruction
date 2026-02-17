//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Tue Aug 22 10:55:44 2017 by ROOT version 5.34/34
// from TTree JP/results from JP real data
// found on file: ../jpsi_mc/jpsi_lp_out.root
//////////////////////////////////////////////////////////

#ifndef JP_selector_h
#define JP_selector_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TSelector.h>
#include <TH2.h>
#include <TH1.h>
#include <TF1.h>
#include <iostream>
#include <fstream>
#include <TVirtualFitter.h>
#include <TMatrixD.h>
#include <TNtuple.h>


// Header file for the classes stored in the TTree if any.

// Fixed size dimensions of array or collections stored in the TTree if any.

class JP_selector : public TSelector {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
ofstream revent;

   double ftmin(double, double );
   double ftmax(double, double );
   double ecorr(double, double ,double);

   int run_old,event_old;
    float trf1;
    float trf2;
    float trf3;
    float trfd;


  TH2D *ethist;


   TH1F *ncpe;
   float minv_old[100];
   float dtrf_old[100];
   float trf_old[100];
   float tbeam_old[100];
   float chi2_old[100];
   float ebeam_old[100];
   float t_old[100];
   int event_no_old[100];
   int anecomb,bnecomb,antcomb,bntcomb;
   int iminv,iebeam;
   float weight_old;

// MC hists
   TH2F* hweight;
   TH2F* htint;
   TH1F* minv;
   TH1F* minvp;

   TH1F* minvr;
   TH1F* minv_nokf;

   TH1F* minva;
   TH1F* minvb;
   TH1F* minvc;

   TH1F* minvap;
   TH1F* minvbp;
   TH1F* minvcp;

   TH1F* eminva[200];
   TH1F* eminvb[200];
   TH1F* eminvc[200];
   TH1F* eminv[200];
   TH1F* eminvap[200];
   TH1F* eminvbp[200];
   TH1F* eminvcp[200];
   TH1F* eminvp[200];
//combo counting
   TH1F* dninva_ecomb;
   TH1F* dninvb_ecomb;
   TH1F* dninva_tcomb;
   TH1F* dninvb_tcomb;

   TH1F* ninva_ecomb;
   TH1F* ninvb_ecomb;
   TH1F* ninva_tcomb;
   TH1F* ninvb_tcomb;

   TH2F* dminva_ecomb;
   TH2F* dminvb_ecomb;
   TH2F* dminva_tcomb;
   TH2F* dminvb_tcomb;

   TH2F* minva_ecomb;
   TH2F* minvb_ecomb;
   TH2F* minva_tcomb;
   TH2F* minvb_tcomb;

   TH2F* deminva_ecomb[200];
   TH2F* deminvb_ecomb[200];
   TH2F* deminva_tcomb[200];
   TH2F* deminvb_tcomb[200];

   TH2F* dtminva_ecomb[200];
   TH2F* dtminvb_ecomb[200];
   TH2F* dtminva_tcomb[200];
   TH2F* dtminvb_tcomb[200];

   TH2F* eminva_ecomb[200];
   TH2F* eminvb_ecomb[200];
   TH2F* eminva_tcomb[200];
   TH2F* eminvb_tcomb[200];

   TH2F* tminva_ecomb[200];
   TH2F* tminvb_ecomb[200];
   TH2F* tminva_tcomb[200];
   TH2F* tminvb_tcomb[200];


//
   TH1F* tminva[200];
   TH1F* tminvb[200];
   TH1F* tminvc[200];
   TH1F* tminv[200];
   TH1F* tminvap[200];
   TH1F* tminvbp[200];
   TH1F* tminvcp[200];
   TH1F* tminvp[200];

// data hists
   TH1F* dminv;
   TH1F* dminvp;

   TH1F* dminvr;
   TH1F* dminvra;
   TH1F* dminvrb;
   TH1F* dminvrc;
   TH1F* dminv_nokf;

   TH1F* dminva;
   TH1F* dminvb;
   TH1F* dminvc;

   TH1F* dminvap;
   TH1F* dminvbp;
   TH1F* dminvcp;

   TH1F* deminva[200];
   TH1F* deminvb[200];
   TH1F* deminvc[200];
   TH1F* deminv[200];
   TH1F* deminvap[200];
   TH1F* deminvbp[200];
   TH1F* deminvcp[200];
   TH1F* deminvp[200];

   TH1F* dtminva[200];
   TH1F* dtminvb[200];
   TH1F* dtminvc[200];
   TH1F* dtminv[200];
   TH1F* dtminvap[200];
   TH1F* dtminvbp[200];
   TH1F* dtminvcp[200];
   TH1F* dtminvp[200];

// fiting functions
   TF1 *expp;
   TF1 *ffit;

   int nbins;
   float mbeg,mend;
   int nbins2;
   float mbeg2,mend2;
   float bcut_data0,bcut_data1,fcut_data0,fcut_data1,theta_cut,chi2_cut,dedx_cut,pbcal_cut;
   float bcut_data0s,bcut_data1s,fcut_data0s,fcut_data1s,theta_cuts,chi2_cuts;
   float p_cut,delta_cut, theta_cm_cut;
   float p_cuts;

   int nebins;
   float ebeg,eend,estep;
   float etmin,etmax;

   int ntbins;
   float tbeg,tend,tstep;
   float temin,temax;

   float scale_jpsi, scale_phi, scale_bh;
   float lumi,lumi16,lumi17,lumi18s,lumi18f;
   int events_jp,events_jp16,events_jp17,events_jp18s,events_jp18f;
   int events_bh,events_bh16,events_bh17,events_bh18s,events_bh18f;
   float xsec_jp,br_jp,xsec_bh;

float nrsigmas; //number of sigmas to cut E/p
float nlsigmas; //number of sigmas to cut E/p
float ndesigmas; //number of sigmas for dE/dx cut


   float bh1, bh2;

   virtual void    FitPhi(TH1F*,TH1F*,TH1F*,TH1F*,double&,double&,double&,bool);
   virtual void    FitJpsi(TH1F*,TH1F*,TH1F*,TH1F*,double&,double&,double&,bool,int, double&);
   TFile *fout;

   // Declaration of leaf types
   Float_t         ebeam;
   Int_t           pol;
   Float_t         pphi;
   Float_t         Phi;
   Float_t         Theta;
   Float_t         tbeam;
   Float_t         trf;
   Float_t         chi2;
   Float_t         ndf;
   Float_t         pp;
   Float_t         thp;
   Float_t         php;
   Float_t         fcalp;
   Float_t         bcalp;
   Float_t         tp;
   Float_t         pep;
   Float_t         thep;
   Float_t         phep;
   Float_t         fcalep;
   Float_t         bcalep;
   Float_t         pbcalep;
   Float_t         tep;
   Float_t         pem;
   Float_t         them;
   Float_t         phem;
   Float_t         fcalem;
   Float_t         bcalem;
   Float_t         pbcalem;
   Float_t         tem;
   Float_t         Minv;
   Float_t         t;
   Float_t         tmin;
   Float_t         mmiss2;
   Float_t         ptmiss;
   Float_t         fdedxp;
   Float_t         cdedxp;
   Float_t         fdedxep;
   Float_t         cdedxep;
   Float_t         fdedxem;
   Float_t         cdedxem;
   Float_t         missem;
   Float_t         missep;
   Float_t         pp_m;
   Float_t         thp_m;
   Float_t         php_m;
   Float_t         tp_m;
   Float_t         pep_m;
   Float_t         thep_m;
   Float_t         phep_m;
   Float_t         tep_m;
   Float_t         pem_m;
   Float_t         them_m;
   Float_t         phem_m;
   Float_t         tem_m;
   Float_t         weight;
   Float_t         pp_t;
   Float_t         thp_t;
   Float_t         php_t;
   Float_t         tp_t;
   Float_t         pep_t;
   Float_t         thep_t;
   Float_t         phep_t;
   Float_t         tep_t;
   Float_t         pem_t;
   Float_t         them_t;
   Float_t         phem_t;
   Float_t         tem_t;
   Float_t         Minv_t;
   Float_t         t_t;
   Float_t         tmin_t;
   Float_t         ebeam_t;
   Float_t         tbeam_t;
   Float_t         Minv_m;
   Float_t         t_m;
   Float_t         cthcm;
   Float_t         cthcm_m;
   Float_t         cthcm_jp;
   Float_t         cthcm_jpm;
   Float_t         Mrec;
   Float_t         Mrec_m;
   Float_t         Mrec_a;
   Float_t         Mang;
   Float_t         Mpt;
   Float_t         pxp;
   Float_t         pxep;
   Float_t         pxem;
   Float_t         ptep;
   Float_t         ptem;
   Float_t         xp;
   Float_t         yp;
   Float_t         zp;
   Float_t         xep;
   Float_t         yep;
   Float_t         zep;
   Float_t         xem;
   Float_t         yem;
   Float_t         zem;
   Float_t         ndfp;
   Float_t         chi2p;
   Float_t         ndfep;
   Float_t         chi2ep;
   Float_t         ndfem;
   Float_t         chi2em;
   Float_t         cthcm_ep;
   Float_t         cthcm_em;
   Float_t         esum;
   Float_t         pxsum;
   Float_t         pysum;
   Float_t         pzsum;
   Int_t           run_no;
   Int_t           event_no;
   Int_t           no_tracks;
   Int_t           no_hypos;
   Int_t           no_energies;
   Int_t           no_combos;

   // List of branches
   TBranch        *b_ebeam;   //!
   TBranch        *b_pol;   //!
   TBranch        *b_pphi;   //!
   TBranch        *b_Phi;   //!
   TBranch        *b_Theta;   //!
   TBranch        *b_tbeam;   //!
   TBranch        *b_trf;   //!
   TBranch        *b_chi2;   //!
   TBranch        *b_ndf;   //!
   TBranch        *b_pp;   //!
   TBranch        *b_thp;   //!
   TBranch        *b_php;   //!
   TBranch        *b_fcalp;   //!
   TBranch        *b_bcalp;   //!
   TBranch        *b_tp;   //!
   TBranch        *b_pep;   //!
   TBranch        *b_thep;   //!
   TBranch        *b_phep;   //!
   TBranch        *b_fcalep;   //!
   TBranch        *b_bcalep;   //!
   TBranch        *b_pbcalep;   //!
   TBranch        *b_tep;   //!
   TBranch        *b_pem;   //!
   TBranch        *b_them;   //!
   TBranch        *b_phem;   //!
   TBranch        *b_fcalem;   //!
   TBranch        *b_bcalem;   //!
   TBranch        *b_pbcalem;   //!
   TBranch        *b_tem;   //!
   TBranch        *b_Minv;   //!
   TBranch        *b_t;   //!
   TBranch        *b_tmin;   //!
   TBranch        *b_mmiss2;   //!
   TBranch        *b_ptmiss;   //!
   TBranch        *b_fdedxp;   //!
   TBranch        *b_cdedxp;   //!
   TBranch        *b_fdedxep;   //!
   TBranch        *b_cdedxep;   //!
   TBranch        *b_fdedxem;   //!
   TBranch        *b_cdedxem;   //!
   TBranch        *b_missem;   //!
   TBranch        *b_missep;   //!
   TBranch        *b_pp_m;   //!
   TBranch        *b_thp_m;   //!
   TBranch        *b_php_m;   //!
   TBranch        *b_tp_m;   //!
   TBranch        *b_pep_m;   //!
   TBranch        *b_thep_m;   //!
   TBranch        *b_phep_m;   //!
   TBranch        *b_tep_m;   //!
   TBranch        *b_pem_m;   //!
   TBranch        *b_them_m;   //!
   TBranch        *b_phem_m;   //!
   TBranch        *b_tem_m;   //!
   TBranch        *b_weight;   //!
   TBranch        *b_pp_t;   //!
   TBranch        *b_thp_t;   //!
   TBranch        *b_php_t;   //!
   TBranch        *b_tp_t;   //!
   TBranch        *b_pep_t;   //!
   TBranch        *b_thep_t;   //!
   TBranch        *b_phep_t;   //!
   TBranch        *b_tep_t;   //!
   TBranch        *b_pem_t;   //!
   TBranch        *b_them_t;   //!
   TBranch        *b_phem_t;   //!
   TBranch        *b_tem_t;   //!
   TBranch        *b_Minv_t;   //!
   TBranch        *b_t_t;   //!
   TBranch        *b_tmin_t;   //!
   TBranch        *b_ebeam_t;   //!
   TBranch        *b_tbeam_t;   //!
   TBranch        *b_Minv_m;   //!
   TBranch        *b_t_m;   //!
   TBranch        *b_cthcm;   //!
   TBranch        *b_cthcm_m;   //!
   TBranch        *b_cthcm_jp;   //!
   TBranch        *b_cthcm_jpm;   //!
   TBranch        *b_Mrec;   //!
   TBranch        *b_Mrec_m;   //!
   TBranch        *b_Mrec_a;   //!
   TBranch        *b_Mang;   //!
   TBranch        *b_Mpt;   //!
   TBranch        *b_pxp;   //!
   TBranch        *b_pxep;   //!
   TBranch        *b_pxem;   //!
   TBranch        *b_ptep;   //!
   TBranch        *b_ptem;   //!
   TBranch        *b_xp;   //!
   TBranch        *b_yp;   //!
   TBranch        *b_zp;   //!
   TBranch        *b_xep;   //!
   TBranch        *b_yep;   //!
   TBranch        *b_zep;   //!
   TBranch        *b_xem;   //!
   TBranch        *b_yem;   //!
   TBranch        *b_zem;   //!
   TBranch        *b_ndfp;   //!
   TBranch        *b_chi2p;   //!
   TBranch        *b_ndfep;   //!
   TBranch        *b_chi2ep;   //!
   TBranch        *b_ndfem;   //!
   TBranch        *b_chi2em;   //!
   TBranch        *b_cthcm_ep;   //!
   TBranch        *b_cthcm_em;   //!
   TBranch        *b_esum;   //!
   TBranch        *b_pxsum;   //!
   TBranch        *b_pysum;   //!
   TBranch        *b_pzsum;   //!
   TBranch        *b_run_no;   //!
   TBranch        *b_event_no;   //!
   TBranch        *b_no_tracks;   //!
   TBranch        *b_no_hypos;   //!
   TBranch        *b_no_energies;   //!
   TBranch        *b_no_combos;   //!

   JP_selector(TTree * /*tree*/ =0) : fChain(0) { }
   virtual ~JP_selector() { }
   virtual Int_t   Version() const { return 2; }
   virtual void    Begin(TTree *tree);
   virtual void    SlaveBegin(TTree *tree);
   virtual void    Init(TTree *tree);
   virtual Bool_t  Notify();
   virtual Bool_t  Process(Long64_t entry);
   virtual Int_t   GetEntry(Long64_t entry, Int_t getall = 0) { return fChain ? fChain->GetTree()->GetEntry(entry, getall) : 0; }
   virtual void    SetOption(const char *option) { fOption = option; }
   virtual void    SetObject(TObject *obj) { fObject = obj; }
   virtual void    SetInputList(TList *input) { fInput = input; }
   virtual TList  *GetOutputList() const { return fOutput; }
   virtual void    SlaveTerminate();
   virtual void    Terminate();

   ClassDef(JP_selector,0);
};

#endif

#ifdef JP_selector_cxx
void JP_selector::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("ebeam", &ebeam, &b_ebeam);
   fChain->SetBranchAddress("pol", &pol, &b_pol);
   fChain->SetBranchAddress("pphi", &pphi, &b_pphi);
   fChain->SetBranchAddress("Phi", &Phi, &b_Phi);
   fChain->SetBranchAddress("Theta", &Theta, &b_Theta);
   fChain->SetBranchAddress("tbeam", &tbeam, &b_tbeam);
   fChain->SetBranchAddress("trf", &trf, &b_trf);
   fChain->SetBranchAddress("chi2", &chi2, &b_chi2);
   fChain->SetBranchAddress("ndf", &ndf, &b_ndf);
   fChain->SetBranchAddress("pp", &pp, &b_pp);
   fChain->SetBranchAddress("thp", &thp, &b_thp);
   fChain->SetBranchAddress("php", &php, &b_php);
   fChain->SetBranchAddress("fcalp", &fcalp, &b_fcalp);
   fChain->SetBranchAddress("bcalp", &bcalp, &b_bcalp);
   fChain->SetBranchAddress("tp", &tp, &b_tp);
   fChain->SetBranchAddress("pep", &pep, &b_pep);
   fChain->SetBranchAddress("thep", &thep, &b_thep);
   fChain->SetBranchAddress("phep", &phep, &b_phep);
   fChain->SetBranchAddress("fcalep", &fcalep, &b_fcalep);
   fChain->SetBranchAddress("bcalep", &bcalep, &b_bcalep);
   fChain->SetBranchAddress("pbcalep", &pbcalep, &b_pbcalep);
   fChain->SetBranchAddress("tep", &tep, &b_tep);
   fChain->SetBranchAddress("pem", &pem, &b_pem);
   fChain->SetBranchAddress("them", &them, &b_them);
   fChain->SetBranchAddress("phem", &phem, &b_phem);
   fChain->SetBranchAddress("fcalem", &fcalem, &b_fcalem);
   fChain->SetBranchAddress("bcalem", &bcalem, &b_bcalem);
   fChain->SetBranchAddress("pbcalem", &pbcalem, &b_pbcalem);
   fChain->SetBranchAddress("tem", &tem, &b_tem);
   fChain->SetBranchAddress("Minv", &Minv, &b_Minv);
   fChain->SetBranchAddress("t", &t, &b_t);
   fChain->SetBranchAddress("tmin", &tmin, &b_tmin);
   fChain->SetBranchAddress("mmiss2", &mmiss2, &b_mmiss2);
   fChain->SetBranchAddress("ptmiss", &ptmiss, &b_ptmiss);
   fChain->SetBranchAddress("fdedxp", &fdedxp, &b_fdedxp);
   fChain->SetBranchAddress("cdedxp", &cdedxp, &b_cdedxp);
   fChain->SetBranchAddress("fdedxep", &fdedxep, &b_fdedxep);
   fChain->SetBranchAddress("cdedxep", &cdedxep, &b_cdedxep);
   fChain->SetBranchAddress("fdedxem", &fdedxem, &b_fdedxem);
   fChain->SetBranchAddress("cdedxem", &cdedxem, &b_cdedxem);
   fChain->SetBranchAddress("missem", &missem, &b_missem);
   fChain->SetBranchAddress("missep", &missep, &b_missep);
   fChain->SetBranchAddress("pp_m", &pp_m, &b_pp_m);
   fChain->SetBranchAddress("thp_m", &thp_m, &b_thp_m);
   fChain->SetBranchAddress("php_m", &php_m, &b_php_m);
   fChain->SetBranchAddress("tp_m", &tp_m, &b_tp_m);
   fChain->SetBranchAddress("pep_m", &pep_m, &b_pep_m);
   fChain->SetBranchAddress("thep_m", &thep_m, &b_thep_m);
   fChain->SetBranchAddress("phep_m", &phep_m, &b_phep_m);
   fChain->SetBranchAddress("tep_m", &tep_m, &b_tep_m);
   fChain->SetBranchAddress("pem_m", &pem_m, &b_pem_m);
   fChain->SetBranchAddress("them_m", &them_m, &b_them_m);
   fChain->SetBranchAddress("phem_m", &phem_m, &b_phem_m);
   fChain->SetBranchAddress("tem_m", &tem_m, &b_tem_m);
   fChain->SetBranchAddress("weight", &weight, &b_weight);
   fChain->SetBranchAddress("pp_t", &pp_t, &b_pp_t);
   fChain->SetBranchAddress("thp_t", &thp_t, &b_thp_t);
   fChain->SetBranchAddress("php_t", &php_t, &b_php_t);
   fChain->SetBranchAddress("tp_t", &tp_t, &b_tp_t);
   fChain->SetBranchAddress("pep_t", &pep_t, &b_pep_t);
   fChain->SetBranchAddress("thep_t", &thep_t, &b_thep_t);
   fChain->SetBranchAddress("phep_t", &phep_t, &b_phep_t);
   fChain->SetBranchAddress("tep_t", &tep_t, &b_tep_t);
   fChain->SetBranchAddress("pem_t", &pem_t, &b_pem_t);
   fChain->SetBranchAddress("them_t", &them_t, &b_them_t);
   fChain->SetBranchAddress("phem_t", &phem_t, &b_phem_t);
   fChain->SetBranchAddress("tem_t", &tem_t, &b_tem_t);
   fChain->SetBranchAddress("Minv_t", &Minv_t, &b_Minv_t);
   fChain->SetBranchAddress("t_t", &t_t, &b_t_t);
   fChain->SetBranchAddress("tmin_t", &tmin_t, &b_tmin_t);
   fChain->SetBranchAddress("ebeam_t", &ebeam_t, &b_ebeam_t);
   fChain->SetBranchAddress("tbeam_t", &tbeam_t, &b_tbeam_t);
   fChain->SetBranchAddress("Minv_m", &Minv_m, &b_Minv_m);
   fChain->SetBranchAddress("t_m", &t_m, &b_t_m);
   fChain->SetBranchAddress("cthcm", &cthcm, &b_cthcm);
   fChain->SetBranchAddress("cthcm_m", &cthcm_m, &b_cthcm_m);
   fChain->SetBranchAddress("cthcm_jp", &cthcm_jp, &b_cthcm_jp);
   fChain->SetBranchAddress("cthcm_jpm", &cthcm_jpm, &b_cthcm_jpm);
   fChain->SetBranchAddress("Mrec", &Mrec, &b_Mrec);
   fChain->SetBranchAddress("Mrec_m", &Mrec_m, &b_Mrec_m);
   fChain->SetBranchAddress("Mrec_a", &Mrec_a, &b_Mrec_a);
   fChain->SetBranchAddress("Mang", &Mang, &b_Mang);
   fChain->SetBranchAddress("Mpt", &Mpt, &b_Mpt);
   fChain->SetBranchAddress("pxp", &pxp, &b_pxp);
   fChain->SetBranchAddress("pxep", &pxep, &b_pxep);
   fChain->SetBranchAddress("pxem", &pxem, &b_pxem);
   fChain->SetBranchAddress("ptep", &ptep, &b_ptep);
   fChain->SetBranchAddress("ptem", &ptem, &b_ptem);
   fChain->SetBranchAddress("xp", &xp, &b_xp);
   fChain->SetBranchAddress("yp", &yp, &b_yp);
   fChain->SetBranchAddress("zp", &zp, &b_zp);
   fChain->SetBranchAddress("xep", &xep, &b_xep);
   fChain->SetBranchAddress("yep", &yep, &b_yep);
   fChain->SetBranchAddress("zep", &zep, &b_zep);
   fChain->SetBranchAddress("xem", &xem, &b_xem);
   fChain->SetBranchAddress("yem", &yem, &b_yem);
   fChain->SetBranchAddress("zem", &zem, &b_zem);
   fChain->SetBranchAddress("ndfp", &ndfp, &b_ndfp);
   fChain->SetBranchAddress("chi2p", &chi2p, &b_chi2p);
   fChain->SetBranchAddress("ndfep", &ndfep, &b_ndfep);
   fChain->SetBranchAddress("chi2ep", &chi2ep, &b_chi2ep);
   fChain->SetBranchAddress("ndfem", &ndfem, &b_ndfem);
   fChain->SetBranchAddress("chi2em", &chi2em, &b_chi2em);
   fChain->SetBranchAddress("cthcm_ep", &cthcm_ep, &b_cthcm_ep);
   fChain->SetBranchAddress("cthcm_em", &cthcm_em, &b_cthcm_em);
   fChain->SetBranchAddress("esum", &esum, &b_esum);
   fChain->SetBranchAddress("pxsum", &pxsum, &b_pxsum);
   fChain->SetBranchAddress("pysum", &pysum, &b_pysum);
   fChain->SetBranchAddress("pzsum", &pzsum, &b_pzsum);
   fChain->SetBranchAddress("run_no", &run_no, &b_run_no);
   fChain->SetBranchAddress("event_no", &event_no, &b_event_no);
   fChain->SetBranchAddress("no_tracks", &no_tracks, &b_no_tracks);
   fChain->SetBranchAddress("no_hypos", &no_hypos, &b_no_hypos);
   fChain->SetBranchAddress("no_energies", &no_energies, &b_no_energies);
   fChain->SetBranchAddress("no_combos", &no_combos, &b_no_combos);
}

Bool_t JP_selector::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

#endif // #ifdef JP_selector_cxx
