#define jpsi_lp_selector_cxx
// The class definition in jpsi_lp_selector.h has been generated automatically
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
// Root > T->Process("jpsi_lp_selector.C")
// Root > T->Process("jpsi_lp_selector.C","some options")
// Root > T->Process("jpsi_lp_selector.C+")
//

#include "jpsi_lp_selector.h"
#include <TH2.h>
#include <TStyle.h>
#include <TH1.h>
#include "TROOT.h"
#include <iostream>
#include "TMath.h"
#include "TMatrixD.h"
#include "TMatrixDLazy.h"
#include "TVectorD.h"
#include "TDecompLU.h"
#include "TDecompSVD.h"
#include <sstream>
#include <string>

#define mass 0.938272 //proton mass
#define PI 3.14159265
#define radian 57.29577958
//#define jpmass 3.097 // jp mass


void jpsi_lp_selector::Begin(TTree * /*tree*/)
{
   // The Begin() function is called at the start of the query.
   // When running with PROOF Begin() is only called on the client.
   // The tree argument is deprecated (on PROOF 0 is passed).

jpmass=3.097;

TString option = GetOption();
cout<<" begin "<<endl;

//ifstream fpolar("para_perp.dat");
ifstream fpolar("rcdb_pol.dat");
int nrun=0;
while (fpolar.good()){
	string changle;
	fpolar>>run_polar[nrun]>>polarization[nrun]>>changle;
	nrun++;
	//if(nrun!=0)cout<<"   "<<nrun<<" "<<run_polar[nrun-1]<<" "<<changle<<" "<<polarization[nrun-1]<<endl;
}
rall=nrun;
ifstream fscale("ScalingFactors_Microscope.dat");
int nrun_scale=0;
while (fscale.good()){
	fscale>>run_scale[nrun_scale]>>scale[nrun_scale]>>escale[nrun_scale];
	nrun_scale++;
	if(nrun_scale!=0)cout<<"   "<<nrun_scale<<" "<<run_scale[nrun_scale-1]<<" "<<scale[nrun_scale-1]<<" "<<escale[nrun_scale-1]<<endl;
}
rall_scale=nrun_scale;

}

void jpsi_lp_selector::SlaveBegin(TTree * /*tree*/)
{
   // The SlaveBegin() function is called after the Begin() function.
   // When running with PROOF SlaveBegin() is called on each slave server.
   // The tree argument is deprecated (on PROOF 0 is passed).

   cout<<" begin slave 1 "<<endl;
   TString option = GetOption();
    fOut = new TFile("jpsi_lp_out2.root", "RECREATE");
    JP = new TTree("JP", "results from JP real data");

    JP->Branch("ebeam", &ebeam, "ebeam/F");
    JP->Branch("pol", &pol, "pol/I");
    JP->Branch("pphi", &pphi, "pphi/F");
    JP->Branch("Phi", &Phi, "Phi/F");
    JP->Branch("Theta", &Theta, "Theta/F");
    JP->Branch("tbeam", &tbeam, "tbeam/F");
    JP->Branch("trf", &trf, "trf/F");
    JP->Branch("chi2", &chi2, "chi2/F");
    JP->Branch("ndf", &ndf, "ndf/F");

    JP->Branch("pp", &pp, "pp/F");
    JP->Branch("thp", &thp, "thp/F");
    JP->Branch("php", &php, "php/F");
    JP->Branch("fcalp", &fcalp, "fcalp/F");
    JP->Branch("bcalp", &bcalp, "bcalp/F");
    JP->Branch("tp", &tp, "tp/F");

    JP->Branch("pep", &pep, "pep/F");
    JP->Branch("thep", &thep, "thep/F");
    JP->Branch("phep", &phep, "phep/F");
    JP->Branch("fcalep", &fcalep, "fcalep/F");
    JP->Branch("bcalep", &bcalep, "bcalep/F");
    JP->Branch("pbcalep", &pbcalep, "pbcalep/F");
    JP->Branch("tep", &tep, "tep/F");

    JP->Branch("pem", &pem, "pem/F");
    JP->Branch("them", &them, "them/F");
    JP->Branch("phem", &phem, "phem/F");
    JP->Branch("fcalem", &fcalem, "fcalem/F");
    JP->Branch("bcalem", &bcalem, "bcalem/F");
    JP->Branch("pbcalem", &pbcalem, "pbcalem/F");
    JP->Branch("tem", &tem, "tem/F");

    JP->Branch("Minv", &Minv, "Minv/F");
    JP->Branch("t", &t, "t/F");
    JP->Branch("tmin", &tmin, "tmin/F");
    JP->Branch("mmiss2", &mmiss2, "mmiss2/F");
    JP->Branch("ptmiss", &ptmiss, "ptmiss/F");

    JP->Branch("fdedxp", &fdedxp, "fdedxp/F");
    JP->Branch("cdedxp", &cdedxp, "cdedxp/F");
    JP->Branch("fdedxep", &fdedxep, "fdedxep/F");
    JP->Branch("cdedxep", &cdedxep, "cdedxep/F");
    JP->Branch("fdedxem", &fdedxem, "fdedxem/F");
    JP->Branch("cdedxem", &cdedxem, "cdedxem/F");

    JP->Branch("missem", &missem, "missem/F");
    JP->Branch("missep", &missep, "missep/F");

    JP->Branch("pp_m", &pp_m, "pp_m/F");
    JP->Branch("thp_m", &thp_m, "thp_m/F");
    JP->Branch("php_m", &php_m, "php_m/F");
    JP->Branch("tp_m", &tp_m, "tp_m/F");

    JP->Branch("pep_m", &pep_m, "pep_m/F");
    JP->Branch("thep_m", &thep_m, "thep_m/F");
    JP->Branch("phep_m", &phep_m, "phep_m/F");
    JP->Branch("tep_m", &tep_m, "tep_m/F");

    JP->Branch("pem_m", &pem_m, "pem_m/F");
    JP->Branch("them_m", &them_m, "them_m/F");
    JP->Branch("phem_m", &phem_m, "phem_m/F");
    JP->Branch("tem_m", &tem_m, "tem_m/F");

    JP->Branch("Minv_m", &Minv_m, "Minv_m/F");
    JP->Branch("t_m", &t_m, "t_m/F");

    JP->Branch("cthcm", &cthcm, "cthcm/F");
    JP->Branch("cthcm_m", &cthcm_m, "cthcm_m/F");

    JP->Branch("cthcm_jp", &cthcm_jp, "cthcm_jp/F");
    JP->Branch("cthcm_jpm", &cthcm_jpm, "cthcm_jpm/F");

    JP->Branch("Mrec", &Mrec, "Mrec/F");
    JP->Branch("Mrec_m", &Mrec_m, "Mrec_m/F");
    JP->Branch("Mrec_a", &Mrec_a, "Mrec_a/F");
    JP->Branch("Mang", &Mang, "Mang/F");
    JP->Branch("Mpt", &Mpt, "Mpt/F");

    JP->Branch("pxp", &pxp, "pxp/F");
    JP->Branch("pxep", &pxep, "pxep/F");
    JP->Branch("pxem", &pxem, "pxem/F");
    JP->Branch("ptep", &ptep, "ptep/F");
    JP->Branch("ptem", &ptem, "ptem/F");

    JP->Branch("xp", &xp, "xp/F");
    JP->Branch("yp", &yp, "yp/F");
    JP->Branch("zp", &zp, "zp/F");
    JP->Branch("xep", &xep, "xep/F");
    JP->Branch("yep", &yep, "yep/F");
    JP->Branch("zep", &zep, "zep/F");
    JP->Branch("xem", &xem, "xem/F");
    JP->Branch("yem", &yem, "yem/F");
    JP->Branch("zem", &zem, "zem/F");
    JP->Branch("xkf", &xkf, "xkf/F");
    JP->Branch("ykf", &ykf, "ykf/F");
    JP->Branch("zkf", &zkf, "zkf/F");

    JP->Branch("ndfp", &ndfp, "ndfp/F");
    JP->Branch("chi2p", &chi2p, "chi2p/F");
    JP->Branch("ndfep", &ndfep, "ndfep/F");
    JP->Branch("chi2ep", &chi2ep, "chi2ep/F");
    JP->Branch("ndfem", &ndfem, "ndfem/F");
    JP->Branch("chi2em", &chi2em, "chi2em/F");

    JP->Branch("cthcm_ep", &cthcm_ep, "cthcm_ep/F");
    JP->Branch("cthcm_em", &cthcm_em, "cthcm_em/F");

    JP->Branch("esum", &esum, "esum/F");
    JP->Branch("pxsum", &pxsum, "pxsum/F");
    JP->Branch("pysum", &pysum, "pysum/F");
    JP->Branch("pzsum", &pzsum, "pzsum/F");

    JP->Branch("run_no", &run_no, "run_no/I");
    JP->Branch("event_no", &event_no, "event_no/I");

    JP->Branch("no_tracks", &NumTracks, "no_tracks/I");
    JP->Branch("no_hypos", &NumChargedHypos, "no_hypos/I");
    JP->Branch("no_energies", &NumBeam, "no_energies/I");
    JP->Branch("no_combos", &NumCombos, "no_combos/I");

   
   cout<<" begin slave 2 "<<endl;


///*
for(int i=0;i<1000;i++){
	run[i]=-100;
	event[i]=-100;
}
//v11miss    ifstream fmiss("jpsi_v11miss.dat");
//    ifstream fmiss("jpsi_v11andv6.dat");
//cok cold    ifstream fmiss("jpsi_v6miss.dat");
ifstream fmiss("20f2017all.list");
int ind=0;
while (fmiss.good()){
	fmiss>>run[ind]>>event[ind];
	//fmiss>>run[ind]>>event[ind]>>event2[ind]>>event3[ind];
	cout<<" run,event="<<run[ind]<<" "<<event[ind]<<endl;
	//cout<<" run,event,2,3="<<run[ind]<<" "<<event[ind]<<" "<<event2[ind]<<" "<<event3[ind]<<endl;
	ind++;
}
//*/

old_rn=0;
event_no_max=0;

}

Bool_t jpsi_lp_selector::Process(Long64_t entry)
{
   // The Process() function is called for each entry in the tree (or possibly
   // keyed object in the case of PROOF) to be processed. The entry argument
   // specifies which entry in the currently loaded tree is to be processed.
   // It can be passed to either jpsi_lp_selector::GetEntry() or TBranch::GetEntry()
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

   if(entry==0){
     //cout<<"begin reading new file"<<endl;
   }
   GetEntry(entry);
	
   
   if(entry==0){
     float rscale=1;
     if(RunNumber!=old_rn){
       //cout<<" run, event_no_max="<<old_rn<<" "<<event_no_max<<endl;
       event_no_max=0;
       file_no=0;
       old_rn=RunNumber;
       for (int rnum=0;rnum<rall;rnum++){
          if(RunNumber==run_polar[rnum]){
             pol=polarization[rnum];
             break;
          } 
       }
       for (int rnum=0;rnum<rall_scale;rnum++){
          if(RunNumber==run_scale[rnum]){
             rscale=scale[rnum];
             break;
          } 
       }
     }
     cout<<" run,pol,rscale,file="<<RunNumber<<" "<<pol<<" "<<rscale<<" "<<file_no<<endl;
     fdedxem=rscale;
     file_no++;
   }

   if(event_no>event_no_max)event_no_max=event_no;
//cok test event selection
    bool select_ev=false;
    for (int i=0;i<1000;i++){
      if(RunNumber==run[i]&&EventNumber==event[i])select_ev=true;
/*
      if(RunNumber==run[i]){
        if(event[i]!=event2[i]){
           if(EventNumber<=event[i]||EventNumber>=event2[i])select_ev=true;
        } else {
           if(EventNumber<=event[i])select_ev=true;
        }
        if(EventNumber>event3[i])select_ev=false;
      }
*/
    }
    if(select_ev)
   cout<<" RunNumber, FileNumber, EventNumber = "<<RunNumber<<" "<<file_no-1<<" "<<EventNumber<<endl;
   //return 0;
    
// find the number of actual tracks (i.e. physical segments ignoring PID)
     float tracks[40];
     float energies[40];
     for (int i=0;i<40;i++){
       tracks[i]=-100.;
     } 

     NumTracks=0.;
     bool old_track;
     for(int ih=0;ih<NumChargedHypos;ih++){
      float track_no=ChargedHypo__TrackID[ih];
      if(NumTracks==0){
        old_track=false;
      } else {
      old_track=false;
      for (int i=0;i<40;i++){
       if(tracks[i]==track_no){
         old_track=true;
         break;
       }
      } 
      }
      if(!old_track){
        tracks[NumTracks]=track_no;
        NumTracks++;
      }

     } //end loop over ChargedHypos

  double Combo_Invmass[80] = {0};
  double Combo_ElectronE[80] = {0}; //REPLACE 6->20 when running multiple files
  double Combo_PositronE[80] = {0}; // Combo variables.


   ebeam_max=0.;
for(UInt_t loc_i = 0; loc_i < NumCombos; ++loc_i) {


     Int_t MyPositron_index = Positron__ChargedIndex[loc_i];
     TLorentzVector& MyPositron = *((TLorentzVector*)ChargedHypo__P4_Measured->At(MyPositron_index));
     TLorentzVector& MyPositron_KF = *((TLorentzVector*)Positron__P4_KinFit->At(loc_i));
     TLorentzVector& MyPositronX = *((TLorentzVector*)ChargedHypo__X4_Measured->At(MyPositron_index));
     TLorentzVector& MyPositronX_KF = *((TLorentzVector*)Positron__X4_KinFit->At(loc_i));

     Int_t MyElectron_index = Electron__ChargedIndex[loc_i];
     TLorentzVector& MyElectron = *((TLorentzVector*)ChargedHypo__P4_Measured->At(MyElectron_index));
     TLorentzVector& MyElectron_KF = *((TLorentzVector*)Electron__P4_KinFit->At(loc_i));
     TLorentzVector& MyElectronX = *((TLorentzVector*)ChargedHypo__X4_Measured->At(MyElectron_index));
     TLorentzVector& MyElectronX_KF = *((TLorentzVector*)Electron__X4_KinFit->At(loc_i));

     Int_t MyProton_index = Proton__ChargedIndex[loc_i];
     TLorentzVector& MyProton = *((TLorentzVector*)ChargedHypo__P4_Measured->At(MyProton_index));
     TLorentzVector& MyProton_KF = *((TLorentzVector*)Proton__P4_KinFit->At(loc_i));
     TLorentzVector& MyProtonX = *((TLorentzVector*)ChargedHypo__X4_Measured->At(MyProton_index));
     TLorentzVector& MyProtonX_KF = *((TLorentzVector*)Proton__X4_KinFit->At(loc_i));

     TLorentzVector Beam_Target;
     Beam_Target.SetE(0.938272);
     Int_t locBeamIndex = ComboBeam__BeamIndex[loc_i];
     TLorentzVector& Beam_photon = *((TLorentzVector*)Beam__P4_Measured->At(locBeamIndex));
     TLorentzVector& Beam_photonX = *((TLorentzVector*)Beam__X4_Measured->At(locBeamIndex));

     TLorentzVector jpsi_IM = MyPositron + MyElectron;
     TLorentzVector jpsi_IM_KF = MyPositron_KF + MyElectron_KF;
     TLorentzVector My_initstate = Beam_Target + Beam_photon;
     TLorentzVector My_finalstate = MyPositron + MyElectron + MyProton;
     TLorentzVector My_finalstate_KF = MyPositron_KF + MyElectron_KF + MyProton_KF;

     TLorentzVector My_MissM = My_initstate - My_finalstate;
     TLorentzVector My_MissM_KF = My_initstate - My_finalstate_KF;

     TLorentzVector eminus_p = MyElectron + MyProton;
     TLorentzVector eminus_p_KF = MyElectron_KF + MyProton_KF;

     TLorentzVector eplus_p = MyPositron + MyProton;
     TLorentzVector eplus_p_KF = MyPositron_KF + MyProton_KF;

     TLorentzVector eminus_eplus = MyElectron + MyPositron;
     TLorentzVector eminus_eplus_KF = MyElectron_KF + MyPositron_KF;

     ndfp=ChargedHypo__NDF_Tracking[MyProton_index];
     chi2p=ChargedHypo__ChiSq_Tracking[MyProton_index];
     ndfep=ChargedHypo__NDF_Tracking[MyPositron_index];
     chi2ep=ChargedHypo__ChiSq_Tracking[MyPositron_index];
     ndfem=ChargedHypo__NDF_Tracking[MyElectron_index];
     chi2em=ChargedHypo__ChiSq_Tracking[MyElectron_index];

     Float_t ElectronBCAL = ChargedHypo__Energy_BCAL[MyElectron_index];
     Float_t ElectronPBCAL = ChargedHypo__Energy_BCALPreshower[MyElectron_index];
     Float_t ElectronFCAL = ChargedHypo__Energy_FCAL[MyElectron_index];
     Float_t PositronBCAL = ChargedHypo__Energy_BCAL[MyPositron_index];
     Float_t PositronPBCAL = ChargedHypo__Energy_BCALPreshower[MyPositron_index];
     Float_t PositronFCAL = ChargedHypo__Energy_FCAL[MyPositron_index];
     Float_t ProtonBCAL = ChargedHypo__Energy_BCAL[MyProton_index];
     Float_t ProtonFCAL = ChargedHypo__Energy_FCAL[MyProton_index];

     double Electron_BCAL_EP = 0;
     double Electron_FCAL_EP = 0;
     double Positron_BCAL_EP = 0;
     double Positron_FCAL_EP = 0;
     double Proton_BCAL_EP = 0;
     double Proton_FCAL_EP = 0;

     if (MyElectron.P() > 0) {
       Electron_BCAL_EP = ElectronBCAL / MyElectron.P();
       Electron_FCAL_EP = ElectronFCAL / MyElectron.P();
     }
     if (MyPositron.P() > 0) {
       Positron_BCAL_EP = PositronBCAL / MyPositron.P();
       Positron_FCAL_EP = PositronFCAL / MyPositron.P();
     }
     if (MyProton.P() > 0) {
       Proton_BCAL_EP = ProtonBCAL / MyProton.P();
       Proton_FCAL_EP = ProtonFCAL / MyProton.P();
     }


     double FLAG_skip = 0;

     
    //   for(UInt_t loc_kk = 0; loc_kk < loc_i; ++loc_kk){
    //       if(Combo_Invmass[loc_kk] == jpsi_IM.M()){
    //          FLAG_skip = 1; // Stops double counting
    //       }
    //    }
    
     Combo_Invmass[loc_i] = jpsi_IM.M();

     float ebeam_test=Beam_photon.E();
     //cold !!!! if (ebeam_test>=ebeam_max){
     chi2=ChiSq_KinFit[loc_i];
     ndf=NDF_KinFit[loc_i];
     ebeam=Beam_photon.E();
     ebeam_max=ebeam;
     tbeam=Beam_photonX.T();

     trf=RFTime_Measured[loc_i];

     pp=MyProton_KF.P();
     thp=MyProton_KF.Theta();
     php=MyProton_KF.Phi();
     tp=MyProtonX_KF.T();
     tp_m=MyProtonX.T();

     pep=MyPositron_KF.P();
     thep=MyPositron_KF.Theta();
     phep=MyPositron_KF.Phi();
     tep=MyPositronX_KF.T();
     tep_m=MyPositronX.T();

     pem=MyElectron_KF.P();
     them=MyElectron_KF.Theta();
     phem=MyElectron_KF.Phi();
     tem=MyElectronX_KF.T();
     tem_m=MyElectronX.T();

     bcalp=ProtonBCAL;
     fcalp=ProtonFCAL;
     bcalep=PositronBCAL;
     pbcalep=PositronPBCAL;
     fcalep=PositronFCAL;
     bcalem=ElectronBCAL;
     pbcalem=ElectronPBCAL;
     fcalem=ElectronFCAL;

     fdedxp=ChargedHypo__dEdx_FDC[MyProton_index];
     cdedxp=ChargedHypo__dEdx_CDC[MyProton_index];
     fdedxep=ChargedHypo__dEdx_FDC[MyPositron_index];
     cdedxep=ChargedHypo__dEdx_CDC[MyPositron_index];
     //cold fdedxem=ChargedHypo__dEdx_FDC[MyElectron_index];
     cdedxem=ChargedHypo__dEdx_CDC[MyElectron_index];

     Kine(ebeam,MyElectron_KF,MyPositron_KF,MyProton_KF);
     Minv=jpsi_IM_KF.M();
     TVector3 boost = - jpsi_IM_KF.BoostVector();
     TLorentzVector Electron_rf=MyElectron_KF;
     Electron_rf.Boost(boost);
     pphi=Electron_rf.Phi();
     Phi=phi_cm;
     Theta=theta_cm;
     t=(Beam_photon-jpsi_IM_KF)*(Beam_photon-jpsi_IM_KF);

     pp_m=MyProton.P();
     thp_m=MyProton.Theta();
     php_m=MyProton.Phi();
     pep_m=MyPositron.P();
     pem_m=MyElectron.P();
     thep_m=MyPositron.Theta();
     phep_m=MyPositron.Phi();
     them_m=MyElectron.Theta();
     phem_m=MyElectron.Phi();
     Minv_m=jpsi_IM.M();
     t_m=(Beam_photon-jpsi_IM)*(Beam_photon-jpsi_IM);

     xp=MyProtonX.X();
     yp=MyProtonX.Y();
     zp=MyProtonX.Z();
     xep=MyPositronX.X();
     yep=MyPositronX.Y();
     zep=MyPositronX.Z();
     xem=MyElectronX.X();
     yem=MyElectronX.Y();
     zem=MyElectronX.Z();
     xkf=MyProtonX_KF.X();
     ykf=MyProtonX_KF.Y();
     zkf=MyProtonX_KF.Z();

     esum=(MyProton.E()+MyPositron.E()+MyElectron.E())-(ebeam+mass);
     pxsum=(MyProton.Px()+MyPositron.Px()+MyElectron.Px());
     pysum=(MyProton.Py()+MyPositron.Py()+MyElectron.Py());
     pzsum=(MyProton.Pz()+MyPositron.Pz()+MyElectron.Pz())-ebeam;

     Mrec=(MyProton+MyPositron+MyElectron).M();

     Mrec_m=(My_initstate-MyProton).M();

// begin kin. rec 
//  angles only

     double theta=MyProton.Theta();
     double phi=MyProton.Phi();
     double thetap=MyPositron.Theta();
     double phip=MyPositron.Phi();
     double thetam=MyElectron.Theta();
     double phim=MyElectron.Phi();
     double ct=cos(theta); double st=sin(theta); double cf=cos(phi); double sf=sin(phi);
     double ctp=cos(thetap); double stp=sin(thetap); double cfp=cos(phip); double sfp=sin(phip);
     double ctm=cos(thetam); double stm=sin(thetam); double cfm=cos(phim); double sfm=sin(phim);

     TMatrixD a(3,3); TVectorD b(3);
     a(0,0)=ct; a(0,1)=ctp; a(0,2)=ctm;
     a(1,0)=st*cf; a(1,1)=stp*cfp; a(1,2)=stm*cfm;
     a(2,0)=st*sf; a(2,1)=stp*sfp; a(2,2)=stm*sfm;
     bool ok;
     b(0)=ebeam; b(1)=0.; b(2)=0.;
     if(MyPositron_index!=MyElectron_index&&MyPositron_index!=MyProton_index&&MyProton_index!=MyElectron_index&&ct!=ctp&&ct!=ctm&&ctp!=ctm){
     TDecompLU lu(a);
     TVectorD xx=lu.Solve(b, ok);
     pxp=xx(0); pxep=xx(1); pxem=xx(2);

     Mang=Minv_m/sqrt(pep_m*pem_m)*sqrt(abs(xx(1)*xx(2)));
     Mrec_a=sqrt(2.*(ebeam*mass+mass*mass-(ebeam+mass)*sqrt(xx(0)*xx(0)+mass*mass)+ebeam*xx(0)*ct));
     } else {
     pxp=0.; pxep=0.; pxem=0.; Mang=0.; Mrec_a=0.;

     }


// end angles only
//
// begin pt fit

     TMatrixD aa(2,2); TVectorD bb(2);
     aa(0,0)=stp*cfp; aa(0,1)=stm*cfm;
     aa(1,0)=stp*sfp; aa(1,1)=stm*sfm;
     bb(0)=-st*cf*pp_m; bb(1)=-st*sf*pp_m;
     if(MyPositron_index!=MyElectron_index&&MyPositron_index!=MyProton_index&&MyProton_index!=MyElectron_index&&stp!=stm){
     TDecompLU luu(aa);
     TVectorD xxx=luu.Solve(bb, ok);
     //cout<<" pt+,pt-="<<xxx(0)<<" "<<xxx(1)<<endl;
     ptep=xxx(0); ptem=xxx(1);

     Mpt=Minv_m/sqrt(pep_m*pem_m)*sqrt(abs(xxx(0)*xxx(1)));
     } else {
     ptep=0.; ptem=0.; Mpt=0.;
     }


// end pt fit
     mmiss2=(My_initstate-My_finalstate)*(My_initstate-My_finalstate);
     ptmiss=My_finalstate.Pt();
   
     TLorentzVector MyNegPion = MyElectron;
     TLorentzVector MyPosPion = MyPositron;
     float mpion=0.13957;
     float ppp=MyElectron.P();
     MyNegPion.SetE(sqrt(mpion*mpion+ppp*ppp));
     ppp=MyPositron.P();
     MyPosPion.SetE(sqrt(mpion*mpion+ppp*ppp));
     missep=(My_initstate-MyProton-MyElectron)*(My_initstate-MyProton-MyElectron);
     missem=(My_initstate-MyProton-MyPositron)*(My_initstate-MyProton-MyPositron);
     missep=(MyProton+MyPosPion).M();
     missem=(MyProton+MyNegPion).M();

     float s=2.*ebeam*mass+mass*mass;
     float gcm=(ebeam+mass)/sqrt(s);
     float bcm=ebeam/(ebeam+mass);
     float pcm=ebeam*mass/sqrt(s);
     jpmass=Minv;
     float ejpcm=(s+jpmass*jpmass-mass*mass)/2./sqrt(s);
     float epcm=(s-jpmass*jpmass+mass*mass)/2./sqrt(s);
     float pjpcm=sqrt(ejpcm*ejpcm-jpmass*jpmass);
     float ppcm=sqrt(epcm*epcm-mass*mass);

     float eepcm=gcm*(pep-bcm*pep*cos(thep));
     float pepcm=gcm*(-bcm*pep+pep*cos(thep));
     float eemcm=gcm*(pem-bcm*pem*cos(them));
     float pemcm=gcm*(-bcm*pem+pem*cos(them));

     float e_m=sqrt(pp_m*pp_m+mass*mass);
     float e=sqrt(pp*pp+mass*mass);
     float ejp_m=pep_m+pem_m;
     float ejp=pep+pem;

     cthcm_m=(e_m/gcm-epcm)/bcm/ppcm;
     cthcm=(e/gcm-epcm)/bcm/ppcm;
     cthcm_jpm=(ejp_m/gcm-ejpcm)/bcm/pjpcm;
     cthcm_jp=(ejp/gcm-ejpcm)/bcm/pjpcm;

     tmin=(pcm-ejpcm)*(pcm-ejpcm)-(pcm-pjpcm)*(pcm-pjpcm);

     TVector3 jp_mom=jpsi_IM_KF.Vect();
     TVector3 ep_mom=MyPositron_KF.Vect();
     TVector3 em_mom=MyElectron_KF.Vect();

     float gjp=ejpcm/jpmass;
     float bjp=pjpcm/ejpcm;
     //float pem_ljp=(jp_mom*ep_mom)/jp_mom.Mag();
     //float pem_tjp=sqrt(pem*pem-pem_ljp*pem_ljp);
     cthcm_em=(eemcm/gjp-ejpcm/2.)/bjp/(ejpcm/2.);
     cthcm_ep=(eepcm/gjp-ejpcm/2.)/bjp/(ejpcm/2.);

     event_no=EventNumber;
     run_no=RunNumber;
     fOut->cd();
     //cok short if(((bcalep/pep>0.5&&pbcalep/pep>(1-bcalep/pep))||fcalep/pep>0.7)&&((bcalem/pem>0.5&&pbcalem/pem>(1-bcalem/pem))||fcalem/pem>0.7))
//cold     if(ebeam>8.2&&(fcalep/pep_m>0.7||bcalep/pep_m>0.7)&&(fcalem/pem_m>0.7||bcalem/pem_m>0.7)){ JP->Fill();}
//     if(ebeam>8.2&&(fcalep/pep_m>0.7||bcalep/pep_m>0.7)&&(fcalem/pem_m>0.7||bcalem/pem_m>0.7)){ JP->Fill();}
     if(ebeam>8.2&&(fcalep/pep>0.7||bcalep/pep>0.7)&&(fcalem/pem>0.7||bcalem/pem>0.7)){ JP->Fill();}
     //cold if(Minv>0.95)JP->Fill();
    //check  if(ebeam>8.2&&(fcalep/pep>0.7||bcalep/pep>0.5)&&(fcalem/pem>0.7||bcalem/pem>0.5)){ JP->Fill();}
     
     //cold } //end ebeam max
   
     } //end loop over combos

   return kTRUE;
}

void jpsi_lp_selector::SlaveTerminate()
{
   // The SlaveTerminate() function is called after all entries or objects
   // have been processed. When running with PROOF SlaveTerminate() is called
   // on each slave server.

}

void jpsi_lp_selector::Terminate()
{
   // The Terminate() function is the last function to be called during
   // a query. It always runs on the client, it can be used to present
   // the results graphically or save the results to file.

   /*
   TH1F *ptmiss_y = (TH1F*) ptmiss->ProjectionY();
   for (int i=1;i<91;i++){
      float bnorm=ptmiss_y->GetBinContent(i);
      for (int j=1;j<101;j++){
         int ibin=ptmiss->GetBin(j,i);
         float cont=ptmiss->GetBinContent(ibin);
         if(bnorm>0){
           cont=cont/bnorm;
         } else {
           cont=0.;
         }
         ptmiss_norm->SetBinContent(ibin,cont);
      }
   } 
   */
   
   fOut->Write();
}

void jpsi_lp_selector::Kine(double Eb,TLorentzVector Lem, TLorentzVector Lep, TLorentzVector Lp1)
{

  TLorentzVector Lbeam;
  TLorentzVector Lp;
    Lp.SetPxPyPzE(0., 0., 0., mass);
    Lbeam.SetPxPyPzE(0, 0, Eb, Eb);  
  TLorentzVector L_mis = Lbeam + Lp - Lem - Lep - Lp1;
  TLorentzVector Lg = Lp1 + Lem + Lep - Lp;
  TLorentzVector Lemep = Lem + Lep;

  double Minv = Lemep.M();
  double Eq_prime = Lemep.E();
  double tM = (Lp - Lp1).M2();
  double Egamma = Lg.E();
  double MM2 = L_mis.M2();
  double mis_mom = L_mis.P();
  double px_mis = L_mis.Px();
  double py_mis = L_mis.Py();
  double pz_mis = L_mis.Pz();

  double Q2 = 2*Eb*(mis_mom - pz_mis);

  TLorentzVector Lcm = Lg + Lp;

  TLorentzVector Lp_cm = Lp;
  Lp_cm.Boost( -Lcm.BoostVector() );
  TLorentzVector Lp1_cm = Lp1;
  Lp1_cm.Boost( -Lcm.BoostVector() );
  TLorentzVector Lem_cm = Lem;
  Lem_cm.Boost( -Lcm.BoostVector() );
  TLorentzVector Lep_cm = Lep;
  Lep_cm.Boost( -Lcm.BoostVector() );

  TLorentzVector Lem_eep_cm = Lem;
  Lem_eep_cm.Boost( -Lemep.BoostVector() );
  TLorentzVector Lep_eep_cm = Lep;
  Lep_eep_cm.Boost( -Lemep.BoostVector() );
  TLorentzVector Lp1_eep_cm = Lp1;
  Lp1_eep_cm.Boost( -Lemep.BoostVector() );

  TVector3 TV3_em = Lem_cm.Vect();
  TVector3 TV3_ep = Lep_cm.Vect();
  TVector3 TV3_p = Lp_cm.Vect();
  TVector3 TV3_p1 = Lp1_cm.Vect();
  TVector3 TV3_emep_crs = TV3_ep.Cross(TV3_em);
  TVector3 TV3_pp1_crs = TV3_p.Cross(TV3_p1);

  if( TV3_em.Dot( TV3_pp1_crs ) > 0 )
    {
      phi_cm = TV3_emep_crs.Angle(TV3_pp1_crs)*radian;
    }
  else
    {
      phi_cm = (2*PI - TV3_emep_crs.Angle(TV3_pp1_crs))*radian;
    }

  theta_cm = (PI - Lem_eep_cm.Angle(Lp1_eep_cm.Vect()))*radian;

  double bb = 2*(Lem - Lep).Dot(Lp - Lp1);
  double L = ((Minv*Minv - tM)*(Minv*Minv - tM) - bb*bb)/4.;
  double L0 = Minv*Minv*Minv*Minv*sin(theta_cm/radian)*sin(theta_cm/radian);

}
