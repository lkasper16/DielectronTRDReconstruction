//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Mon Oct 24 16:55:17 2016 by ROOT version 5.34/34
// from TTree jpsi_lp_Tree/jpsi_lp_Tree
// found on file: phi_timecut.root
//////////////////////////////////////////////////////////

#ifndef jpsi_lp_selector_h
#define jpsi_lp_selector_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TSelector.h>

#include <TH1.h>
#include <TH2.h>
#include <iostream>
#include <fstream>



// Header file for the classes stored in the TTree if any.
#include <TLorentzVector.h>
#include <TClonesArray.h>

// Fixed size dimensions of array or collections stored in the TTree if any.

//ofstream fpolar;
int run_polar[10000];
float polarization[10000];
int rall;

//ofstream fscale;
int run_scale[10000];
float scale[10000];
float escale[10000];
int rall_scale;

class jpsi_lp_selector : public TSelector {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain

    TFile* fOut; //output file
    TTree* JP; //output tree
   float ebeam,pp,thp,pep,pem,thep,them,bcalp,fcalp,bcalep,fcalep,fcalem,bcalem,Minv,t,tmin;
   float pbcalep,pbcalem,chi2,ndf;
   float fdedxp,cdedxp,fdedxep,cdedxep,fdedxem,cdedxem;
   float ebeam_max;
   float tbeam,trf;
   int event_no,run_no;
//cok temp event selection
   int run[1000];
   int event[1000];
   int event2[1000];
   int event3[1000];

   int NumTracks;
   float pp_m,thp_m,pep_m,pem_m,thep_m,them_m,Minv_m,t_m;
   float tp,tep,tem,tp_m,tep_m,tem_m;
   float Theta,php,phep,phem,php_m,phep_m,phem_m;
   float cthcm,cthcm_m,cthcm_jp,cthcm_jpm,Mrec,Mrec_m,Mrec_a,Mang,Mpt;
   float pxp,pxep,pxem,ptep,ptem;
   float xp,yp,zp,xep,yep,zep,xem,yem,zem,xkf,ykf,zkf;
   float ndfp,chi2p,ndfep,chi2ep,ndfem,chi2em;
   float cthcm_ep,cthcm_em;
   float esum,pxsum,pysum,pzsum;
   float jpmass;
   float mmiss2,ptmiss;
   float missem,missep;
   int pol;
   float pphi,Phi;

    int old_rn=0;
    int file_no=0;
    int event_no_max=0;

    double phi_cm,theta_cm;

   // Declaration of leaf types

   UInt_t          RunNumber;
   ULong64_t       EventNumber;
   UInt_t          L1TriggerBits;
   TLorentzVector  *X4_Production;
   UInt_t          NumBeam;
   Int_t           Beam__PID[90000];   //[NumBeam]
   TClonesArray    *Beam__X4_Measured;
   TClonesArray    *Beam__P4_Measured;
   UInt_t          NumNeutralHypos;
   Int_t           NeutralHypo__NeutralID[100000];   //[NumNeutralHypos]
   Int_t           NeutralHypo__PID[100000];   //[NumNeutralHypos]
   TClonesArray    *NeutralHypo__X4_Measured;
   TClonesArray    *NeutralHypo__P4_Measured;
   Float_t         NeutralHypo__Beta_Timing[100000];   //[NumNeutralHypos]
   Float_t         NeutralHypo__ChiSq_Timing[100000];   //[NumNeutralHypos]
   UInt_t          NeutralHypo__NDF_Timing[100000];   //[NumNeutralHypos]
   TClonesArray    *NeutralHypo__X4_Shower;
   Float_t         NeutralHypo__Energy_BCAL[100000];   //[NumNeutralHypos]
   Float_t         NeutralHypo__Energy_BCALPreshower[100000];   //[NumNeutralHypos]
   Float_t         NeutralHypo__Energy_FCAL[100000];   //[NumNeutralHypos]
   Float_t         NeutralHypo__TrackBCAL_DeltaPhi[100000];   //[NumNeutralHypos]
   Float_t         NeutralHypo__TrackBCAL_DeltaZ[100000];   //[NumNeutralHypos]
   Float_t         NeutralHypo__TrackFCAL_DOCA[100000];   //[NumNeutralHypos]
   Float_t         NeutralHypo__PhotonRFDeltaTVar[100000];   //[NumNeutralHypos]
   UInt_t          NumChargedHypos;
   Int_t           ChargedHypo__TrackID[1400000];   //[NumChargedHypos]
   Int_t           ChargedHypo__PID[1400000];   //[NumChargedHypos]
   TClonesArray    *ChargedHypo__X4_Measured;
   TClonesArray    *ChargedHypo__P4_Measured;
   UInt_t          ChargedHypo__NDF_Tracking[1400000];   //[NumChargedHypos]
   Float_t         ChargedHypo__ChiSq_Tracking[1400000];   //[NumChargedHypos]
   UInt_t          ChargedHypo__NDF_DCdEdx[1400000];   //[NumChargedHypos]
   Float_t         ChargedHypo__ChiSq_DCdEdx[1400000];   //[NumChargedHypos]
   Float_t         ChargedHypo__dEdx_CDC[1400000];   //[NumChargedHypos]
   Float_t         ChargedHypo__dEdx_FDC[1400000];   //[NumChargedHypos]
   Float_t         ChargedHypo__HitTime[1400000];   //[NumChargedHypos]
   Float_t         ChargedHypo__RFDeltaTVar[1400000];   //[NumChargedHypos]
   Float_t         ChargedHypo__Beta_Timing[1400000];   //[NumChargedHypos]
   Float_t         ChargedHypo__ChiSq_Timing[1400000];   //[NumChargedHypos]
   UInt_t          ChargedHypo__NDF_Timing[1400000];   //[NumChargedHypos]
   Float_t         ChargedHypo__dEdx_TOF[1400000];   //[NumChargedHypos]
   Float_t         ChargedHypo__dEdx_ST[1400000];   //[NumChargedHypos]
   Float_t         ChargedHypo__Energy_BCAL[1400000];   //[NumChargedHypos]
   Float_t         ChargedHypo__Energy_BCALPreshower[1400000];   //[NumChargedHypos]
   Float_t         ChargedHypo__Energy_FCAL[1400000];   //[NumChargedHypos]
   Float_t         ChargedHypo__TrackBCAL_DeltaPhi[1400000];   //[NumChargedHypos]
   Float_t         ChargedHypo__TrackBCAL_DeltaZ[14000000];   //[NumChargedHypos]
   Float_t         ChargedHypo__TrackFCAL_DOCA[1400000];   //[NumChargedHypos]
   UInt_t          NumCombos;
   Bool_t          IsComboCut[4300000];   //[NumCombos]
   Float_t         RFTime_Measured[4300000];   //[NumCombos]
   Float_t         ChiSq_KinFit[4300000];   //[NumCombos]
   UInt_t          NDF_KinFit[4300000];   //[NumCombos]
   Int_t           ComboBeam__BeamIndex[4300000];   //[NumCombos]
   TClonesArray    *ComboBeam__P4_KinFit;
   TClonesArray    *ComboBeam__X4_KinFit;
   Int_t           Proton__ChargedIndex[4300000];   //[NumCombos]
   TClonesArray    *Proton__P4_KinFit;
   TClonesArray    *Proton__X4_KinFit;
   Float_t         Proton__Beta_Timing_KinFit[4300000];   //[NumCombos]
   Float_t         Proton__ChiSq_Timing_KinFit[4300000];   //[NumCombos]
//   TClonesArray    *DecayingJpsi__P4_KinFit;
   Int_t           Electron__ChargedIndex[4300000];   //[NumCombos]
   TClonesArray    *Electron__P4_KinFit;
   TClonesArray    *Electron__X4_KinFit;
   Float_t         Electron__Beta_Timing_KinFit[4300000];   //[NumCombos]
   Float_t         Electron__ChiSq_Timing_KinFit[4300000];   //[NumCombos]
   Int_t           Positron__ChargedIndex[4300000];   //[NumCombos]
   TClonesArray    *Positron__P4_KinFit;
   TClonesArray    *Positron__X4_KinFit;
   Float_t         Positron__Beta_Timing_KinFit[4300000];   //[NumCombos]
   Float_t         Positron__ChiSq_Timing_KinFit[4300000];   //[NumCombos]

   // List of branches

   TBranch        *b_RunNumber;   //!
   TBranch        *b_EventNumber;   //!
   TBranch        *b_L1TriggerBits;   //!
   TBranch        *b_X4_Production;   //!
   TBranch        *b_NumBeam;   //!
   TBranch        *b_Beam__PID;   //!
   TBranch        *b_Beam__X4_Measured;   //!
   TBranch        *b_Beam__P4_Measured;   //!
   TBranch        *b_NumNeutralHypos;   //!
   TBranch        *b_NeutralHypo__NeutralID;   //!
   TBranch        *b_NeutralHypo__PID;   //!
   TBranch        *b_NeutralHypo__X4_Measured;   //!
   TBranch        *b_NeutralHypo__P4_Measured;   //!
   TBranch        *b_NeutralHypo__Beta_Timing;   //!
   TBranch        *b_NeutralHypo__ChiSq_Timing;   //!
   TBranch        *b_NeutralHypo__NDF_Timing;   //!
   TBranch        *b_NeutralHypo__X4_Shower;   //!
   TBranch        *b_NeutralHypo__Energy_BCAL;   //!
   TBranch        *b_NeutralHypo__Energy_BCALPreshower;   //!
   TBranch        *b_NeutralHypo__Energy_FCAL;   //!
   TBranch        *b_NeutralHypo__TrackBCAL_DeltaPhi;   //!
   TBranch        *b_NeutralHypo__TrackBCAL_DeltaZ;   //!
   TBranch        *b_NeutralHypo__TrackFCAL_DOCA;   //!
   TBranch        *b_NeutralHypo__PhotonRFDeltaTVar;   //!
   TBranch        *b_NumChargedHypos;   //!
   TBranch        *b_ChargedHypo__TrackID;   //!
   TBranch        *b_ChargedHypo__PID;   //!
   TBranch        *b_ChargedHypo__X4_Measured;   //!
   TBranch        *b_ChargedHypo__P4_Measured;   //!
   TBranch        *b_ChargedHypo__NDF_Tracking;   //!
   TBranch        *b_ChargedHypo__ChiSq_Tracking;   //!
   TBranch        *b_ChargedHypo__NDF_DCdEdx;   //!
   TBranch        *b_ChargedHypo__ChiSq_DCdEdx;   //!
   TBranch        *b_ChargedHypo__dEdx_CDC;   //!
   TBranch        *b_ChargedHypo__dEdx_FDC;   //!
   TBranch        *b_ChargedHypo__HitTime;   //!
   TBranch        *b_ChargedHypo__RFDeltaTVar;   //!
   TBranch        *b_ChargedHypo__Beta_Timing;   //!
   TBranch        *b_ChargedHypo__ChiSq_Timing;   //!
   TBranch        *b_ChargedHypo__NDF_Timing;   //!
   TBranch        *b_ChargedHypo__dEdx_TOF;   //!
   TBranch        *b_ChargedHypo__dEdx_ST;   //!
   TBranch        *b_ChargedHypo__Energy_BCAL;   //!
   TBranch        *b_ChargedHypo__Energy_BCALPreshower;   //!
   TBranch        *b_ChargedHypo__Energy_FCAL;   //!
   TBranch        *b_ChargedHypo__TrackBCAL_DeltaPhi;   //!
   TBranch        *b_ChargedHypo__TrackBCAL_DeltaZ;   //!
   TBranch        *b_ChargedHypo__TrackFCAL_DOCA;   //!
   TBranch        *b_NumCombos;   //!
   TBranch        *b_IsComboCut;   //!
   TBranch        *b_RFTime_Measured;   //!
   TBranch        *b_ChiSq_KinFit;   //!
   TBranch        *b_NDF_KinFit;   //!
   TBranch        *b_ComboBeam__BeamIndex;   //!
   TBranch        *b_ComboBeam__P4_KinFit;   //!
   TBranch        *b_ComboBeam__X4_KinFit;   //!
   TBranch        *b_Proton__ChargedIndex;   //!
   TBranch        *b_Proton__P4_KinFit;   //!
   TBranch        *b_Proton__X4_KinFit;   //!
   TBranch        *b_Proton__Beta_Timing_KinFit;   //!
   TBranch        *b_Proton__ChiSq_Timing_KinFit;   //!
//   TBranch        *b_DecayingJpsi__P4_KinFit;   //!
   TBranch        *b_Electron__ChargedIndex;   //!
   TBranch        *b_Electron__P4_KinFit;   //!
   TBranch        *b_Electron__X4_KinFit;   //!
   TBranch        *b_Electron__Beta_Timing_KinFit;   //!
   TBranch        *b_Electron__ChiSq_Timing_KinFit;   //!
   TBranch        *b_Positron__ChargedIndex;   //!
   TBranch        *b_Positron__P4_KinFit;   //!
   TBranch        *b_Positron__X4_KinFit;   //!
   TBranch        *b_Positron__Beta_Timing_KinFit;   //!
   TBranch        *b_Positron__ChiSq_Timing_KinFit;   //!

   jpsi_lp_selector(TTree * /*tree*/ =0) : fChain(0) { }
   virtual ~jpsi_lp_selector() { }
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

  void Kine(double, TLorentzVector , TLorentzVector , TLorentzVector );

   ClassDef(jpsi_lp_selector,0);
};

#endif

#ifdef jpsi_lp_selector_cxx
void jpsi_lp_selector::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set object pointer

   X4_Production = 0;
   Beam__X4_Measured = 0;
   Beam__P4_Measured = 0;
   NeutralHypo__X4_Measured = 0;
   NeutralHypo__P4_Measured = 0;
   NeutralHypo__X4_Shower = 0;
   ChargedHypo__X4_Measured = 0;
   ChargedHypo__P4_Measured = 0;
   ComboBeam__P4_KinFit = 0;
   ComboBeam__X4_KinFit = 0;
   Proton__P4_KinFit = 0;
   Proton__X4_KinFit = 0;
//   DecayingJpsi__P4_KinFit = 0;
   Electron__P4_KinFit = 0;
   Electron__X4_KinFit = 0;
   Positron__P4_KinFit = 0;
   Positron__X4_KinFit = 0;
   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("RunNumber", &RunNumber, &b_RunNumber);
   fChain->SetBranchAddress("EventNumber", &EventNumber, &b_EventNumber);
   fChain->SetBranchAddress("L1TriggerBits", &L1TriggerBits, &b_L1TriggerBits);
   fChain->SetBranchAddress("X4_Production", &X4_Production, &b_X4_Production);
   fChain->SetBranchAddress("NumBeam", &NumBeam, &b_NumBeam);
   fChain->SetBranchAddress("Beam__PID", Beam__PID, &b_Beam__PID);
   fChain->SetBranchAddress("Beam__X4_Measured", &Beam__X4_Measured, &b_Beam__X4_Measured);
   fChain->SetBranchAddress("Beam__P4_Measured", &Beam__P4_Measured, &b_Beam__P4_Measured);
   fChain->SetBranchAddress("NumNeutralHypos", &NumNeutralHypos, &b_NumNeutralHypos);
   fChain->SetBranchAddress("NeutralHypo__NeutralID", &NeutralHypo__NeutralID, &b_NeutralHypo__NeutralID);
   fChain->SetBranchAddress("NeutralHypo__PID", &NeutralHypo__PID, &b_NeutralHypo__PID);
   fChain->SetBranchAddress("NeutralHypo__X4_Measured", &NeutralHypo__X4_Measured, &b_NeutralHypo__X4_Measured);
   fChain->SetBranchAddress("NeutralHypo__P4_Measured", &NeutralHypo__P4_Measured, &b_NeutralHypo__P4_Measured);
   fChain->SetBranchAddress("NeutralHypo__Beta_Timing", &NeutralHypo__Beta_Timing, &b_NeutralHypo__Beta_Timing);
   fChain->SetBranchAddress("NeutralHypo__ChiSq_Timing", &NeutralHypo__ChiSq_Timing, &b_NeutralHypo__ChiSq_Timing);
   fChain->SetBranchAddress("NeutralHypo__NDF_Timing", &NeutralHypo__NDF_Timing, &b_NeutralHypo__NDF_Timing);
   fChain->SetBranchAddress("NeutralHypo__X4_Shower", &NeutralHypo__X4_Shower, &b_NeutralHypo__X4_Shower);
   fChain->SetBranchAddress("NeutralHypo__Energy_BCAL", &NeutralHypo__Energy_BCAL, &b_NeutralHypo__Energy_BCAL);
   fChain->SetBranchAddress("NeutralHypo__Energy_BCALPreshower", &NeutralHypo__Energy_BCALPreshower, &b_NeutralHypo__Energy_BCALPreshower);
   fChain->SetBranchAddress("NeutralHypo__Energy_FCAL", &NeutralHypo__Energy_FCAL, &b_NeutralHypo__Energy_FCAL);
   fChain->SetBranchAddress("NeutralHypo__TrackBCAL_DeltaPhi", &NeutralHypo__TrackBCAL_DeltaPhi, &b_NeutralHypo__TrackBCAL_DeltaPhi);
   fChain->SetBranchAddress("NeutralHypo__TrackBCAL_DeltaZ", &NeutralHypo__TrackBCAL_DeltaZ, &b_NeutralHypo__TrackBCAL_DeltaZ);
   fChain->SetBranchAddress("NeutralHypo__TrackFCAL_DOCA", &NeutralHypo__TrackFCAL_DOCA, &b_NeutralHypo__TrackFCAL_DOCA);
   fChain->SetBranchAddress("NeutralHypo__PhotonRFDeltaTVar", &NeutralHypo__PhotonRFDeltaTVar, &b_NeutralHypo__PhotonRFDeltaTVar);
   fChain->SetBranchAddress("NumChargedHypos", &NumChargedHypos, &b_NumChargedHypos);
   fChain->SetBranchAddress("ChargedHypo__TrackID", ChargedHypo__TrackID, &b_ChargedHypo__TrackID);
   fChain->SetBranchAddress("ChargedHypo__PID", ChargedHypo__PID, &b_ChargedHypo__PID);
   fChain->SetBranchAddress("ChargedHypo__X4_Measured", &ChargedHypo__X4_Measured, &b_ChargedHypo__X4_Measured);
   fChain->SetBranchAddress("ChargedHypo__P4_Measured", &ChargedHypo__P4_Measured, &b_ChargedHypo__P4_Measured);
   fChain->SetBranchAddress("ChargedHypo__NDF_Tracking", ChargedHypo__NDF_Tracking, &b_ChargedHypo__NDF_Tracking);
   fChain->SetBranchAddress("ChargedHypo__ChiSq_Tracking", ChargedHypo__ChiSq_Tracking, &b_ChargedHypo__ChiSq_Tracking);
   fChain->SetBranchAddress("ChargedHypo__NDF_DCdEdx", ChargedHypo__NDF_DCdEdx, &b_ChargedHypo__NDF_DCdEdx);
   fChain->SetBranchAddress("ChargedHypo__ChiSq_DCdEdx", ChargedHypo__ChiSq_DCdEdx, &b_ChargedHypo__ChiSq_DCdEdx);
   fChain->SetBranchAddress("ChargedHypo__dEdx_CDC", ChargedHypo__dEdx_CDC, &b_ChargedHypo__dEdx_CDC);
   fChain->SetBranchAddress("ChargedHypo__dEdx_FDC", ChargedHypo__dEdx_FDC, &b_ChargedHypo__dEdx_FDC);
   fChain->SetBranchAddress("ChargedHypo__HitTime", ChargedHypo__HitTime, &b_ChargedHypo__HitTime);
   fChain->SetBranchAddress("ChargedHypo__RFDeltaTVar", ChargedHypo__RFDeltaTVar, &b_ChargedHypo__RFDeltaTVar);
   fChain->SetBranchAddress("ChargedHypo__Beta_Timing", ChargedHypo__Beta_Timing, &b_ChargedHypo__Beta_Timing);
   fChain->SetBranchAddress("ChargedHypo__ChiSq_Timing", ChargedHypo__ChiSq_Timing, &b_ChargedHypo__ChiSq_Timing);
   fChain->SetBranchAddress("ChargedHypo__NDF_Timing", ChargedHypo__NDF_Timing, &b_ChargedHypo__NDF_Timing);
   fChain->SetBranchAddress("ChargedHypo__dEdx_TOF", ChargedHypo__dEdx_TOF, &b_ChargedHypo__dEdx_TOF);
   fChain->SetBranchAddress("ChargedHypo__dEdx_ST", ChargedHypo__dEdx_ST, &b_ChargedHypo__dEdx_ST);
   fChain->SetBranchAddress("ChargedHypo__Energy_BCAL", ChargedHypo__Energy_BCAL, &b_ChargedHypo__Energy_BCAL);
   fChain->SetBranchAddress("ChargedHypo__Energy_BCALPreshower", ChargedHypo__Energy_BCALPreshower, &b_ChargedHypo__Energy_BCALPreshower);
   fChain->SetBranchAddress("ChargedHypo__Energy_FCAL", ChargedHypo__Energy_FCAL, &b_ChargedHypo__Energy_FCAL);
   fChain->SetBranchAddress("ChargedHypo__TrackBCAL_DeltaPhi", ChargedHypo__TrackBCAL_DeltaPhi, &b_ChargedHypo__TrackBCAL_DeltaPhi);
   fChain->SetBranchAddress("ChargedHypo__TrackBCAL_DeltaZ", ChargedHypo__TrackBCAL_DeltaZ, &b_ChargedHypo__TrackBCAL_DeltaZ);
   fChain->SetBranchAddress("ChargedHypo__TrackFCAL_DOCA", ChargedHypo__TrackFCAL_DOCA, &b_ChargedHypo__TrackFCAL_DOCA);
   fChain->SetBranchAddress("NumCombos", &NumCombos, &b_NumCombos);
   fChain->SetBranchAddress("IsComboCut", IsComboCut, &b_IsComboCut);
   fChain->SetBranchAddress("RFTime_Measured", RFTime_Measured, &b_RFTime_Measured);
   fChain->SetBranchAddress("ChiSq_KinFit", ChiSq_KinFit, &b_ChiSq_KinFit);
   fChain->SetBranchAddress("NDF_KinFit", NDF_KinFit, &b_NDF_KinFit);
   fChain->SetBranchAddress("ComboBeam__BeamIndex", ComboBeam__BeamIndex, &b_ComboBeam__BeamIndex);
   fChain->SetBranchAddress("ComboBeam__P4_KinFit", &ComboBeam__P4_KinFit, &b_ComboBeam__P4_KinFit);
   fChain->SetBranchAddress("ComboBeam__X4_KinFit", &ComboBeam__X4_KinFit, &b_ComboBeam__X4_KinFit);
   fChain->SetBranchAddress("Proton__ChargedIndex", Proton__ChargedIndex, &b_Proton__ChargedIndex);
   fChain->SetBranchAddress("Proton__P4_KinFit", &Proton__P4_KinFit, &b_Proton__P4_KinFit);
   fChain->SetBranchAddress("Proton__X4_KinFit", &Proton__X4_KinFit, &b_Proton__X4_KinFit);
   fChain->SetBranchAddress("Proton__Beta_Timing_KinFit", Proton__Beta_Timing_KinFit, &b_Proton__Beta_Timing_KinFit);
   fChain->SetBranchAddress("Proton__ChiSq_Timing_KinFit", Proton__ChiSq_Timing_KinFit, &b_Proton__ChiSq_Timing_KinFit);
 //  fChain->SetBranchAddress("DecayingJpsi__P4_KinFit", &DecayingJpsi__P4_KinFit, &b_DecayingJpsi__P4_KinFit);
   fChain->SetBranchAddress("Electron__ChargedIndex", Electron__ChargedIndex, &b_Electron__ChargedIndex);
   fChain->SetBranchAddress("Electron__P4_KinFit", &Electron__P4_KinFit, &b_Electron__P4_KinFit);
   fChain->SetBranchAddress("Electron__X4_KinFit", &Electron__X4_KinFit, &b_Electron__X4_KinFit);
   fChain->SetBranchAddress("Electron__Beta_Timing_KinFit", Electron__Beta_Timing_KinFit, &b_Electron__Beta_Timing_KinFit);
   fChain->SetBranchAddress("Electron__ChiSq_Timing_KinFit", Electron__ChiSq_Timing_KinFit, &b_Electron__ChiSq_Timing_KinFit);
   fChain->SetBranchAddress("Positron__ChargedIndex", Positron__ChargedIndex, &b_Positron__ChargedIndex);
   fChain->SetBranchAddress("Positron__P4_KinFit", &Positron__P4_KinFit, &b_Positron__P4_KinFit);
   fChain->SetBranchAddress("Positron__X4_KinFit", &Positron__X4_KinFit, &b_Positron__X4_KinFit);
   fChain->SetBranchAddress("Positron__Beta_Timing_KinFit", Positron__Beta_Timing_KinFit, &b_Positron__Beta_Timing_KinFit);
   fChain->SetBranchAddress("Positron__ChiSq_Timing_KinFit", Positron__ChiSq_Timing_KinFit, &b_Positron__ChiSq_Timing_KinFit);

}

Bool_t jpsi_lp_selector::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

#endif // #ifdef jpsi_lp_selector_cxx
