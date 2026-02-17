{
ofstream outF("jpsi_run_ev.txt");
int run_no;
int event_no;
//TFile *_file0 = TFile::Open("jpsi_lp_out2_131.root");
TFile *_file0 = TFile::Open("jpsi_lp_out2.root");
TH2F *in_trd_ep = new TH2F("in_trd_ep","Positrons projected at the TRD; X at TRD Plane [cm]; Y at TRD Plane [cm]",1000,-100.,100.,1000,-100.,100.);
TH2F *in_trd_em = new TH2F("in_trd_em","Electrons projected at the TRD; X at TRD Plane [cm]; Y at TRD Plane [cm]",1000,-100.,100.,1000,-100.,100.);
TH2F *all_ep = new TH2F("all_ep","All positrons projected at the TRD plane; X at TRD Plane [cm]; Y at TRD Plane [cm]",1000,-100.,100.,1000,-100.,100.);
TH2F *all_em = new TH2F("all_em","All electrons projected at the TRD plane; X at TRD Plane [cm]; Y at TRD Plane [cm]",1000,-100.,100.,1000,-100.,100.);

TLegend *l0 = new TLegend(0.13,0.79,0.33,0.9);
l0->SetNColumns(2);
l0->SetTextSize(0.035);
in_trd_ep->SetMarkerStyle(38); //crosshair circle
in_trd_ep->SetMarkerColor(2);
in_trd_ep->SetMarkerSize(2);
l0->AddEntry(in_trd_ep,"e^{+} in TRD","p");
all_ep->SetMarkerStyle(24); //open circle
all_ep->SetMarkerColor(2);
all_ep->SetMarkerSize(2);
l0->AddEntry(all_ep,"e^{+}","p");
in_trd_em->SetMarkerStyle(38); //crosshair circle
in_trd_em->SetMarkerColor(4);
in_trd_em->SetMarkerSize(2);
l0->AddEntry(in_trd_em,"e^{-} in TRD","p");
all_em->SetMarkerStyle(24); //open circle
all_em->SetMarkerColor(4);
all_em->SetMarkerSize(2);
l0->AddEntry(all_em,"e^{-}","p");

TTree *JP = (TTree*)_file0->Get("JP");
//JP->SetBranchAddress("run_no",run_no);
//JP->SetBranchAddress("event_no",event_no);
//JP->SetScanFileName("jpsi_run_ev.txt");
TBox *fbox = new TBox(-83.47,-68.6,-11.47,-32.61); //For TRD acceptance visualization
fbox->SetFillStyle(0);
fbox->SetLineWidth(2);
TCanvas *c0 = new TCanvas("c0","c0",1300,1000);
c0->cd();
gPad->SetGridy();
gPad->SetGridx();
gStyle->SetOptStat(0);

//e+ in all acceptance
JP->Draw("sin(phep)*thep*467.4:cos(phep)*thep*467.4>>all_ep","(Minv>3.0555&&Minv<3.1365&&((pep/bcalep<1.173300&&pep/bcalep>0.875700&&bcalep>0&&pbcalep*sin(thep)>0.030000)||(pep/fcalep<1.233800&&pep/fcalep>0.915800&&fcalep>0))&&((pem/bcalem<1.173300&&pem/bcalem>0.875700&&bcalem>0&&pbcalem*sin(them)>0.030000)||(pem/fcalem<1.233800&&pem/fcalem>0.915800&&fcalem>0))&&chi2>0&&chi2<5000.000000&&ebeam>8.20&&ebeam<11.8&&(them>0.034907&&thep>0.034907))&&abs(trf-tbeam)<2.&&thep*467.4<83.","");
//e- in all acceptance
JP->Draw("sin(phem)*them*467.4:cos(phem)*them*467.4>>all_em","(Minv>3.0555&&Minv<3.1365&&((pep/bcalep<1.173300&&pep/bcalep>0.875700&&bcalep>0&&pbcalep*sin(thep)>0.030000)||(pep/fcalep<1.233800&&pep/fcalep>0.915800&&fcalep>0))&&((pem/bcalem<1.173300&&pem/bcalem>0.875700&&bcalem>0&&pbcalem*sin(them)>0.030000)||(pem/fcalem<1.233800&&pem/fcalem>0.915800&&fcalem>0))&&chi2>0&&chi2<5000.000000&&ebeam>8.20&&ebeam<11.8&&(them>0.034907&&thep>0.034907))&&abs(trf-tbeam)<2.&&them*467.4<83.","same");
//e- in TRD
JP->Draw("sin(phem)*them*467.4:cos(phem)*them*467.4>>in_trd_em","(Minv>3.0555&&Minv<3.1365&&((pep/bcalep<1.173300&&pep/bcalep>0.875700&&bcalep>0&&pbcalep*sin(thep)>0.030000)||(pep/fcalep<1.233800&&pep/fcalep>0.915800&&fcalep>0))&&((pem/bcalem<1.173300&&pem/bcalem>0.875700&&bcalem>0&&pbcalem*sin(them)>0.030000)||(pem/fcalem<1.233800&&pem/fcalem>0.915800&&fcalem>0))&&chi2>0&&chi2<5000.000000&&ebeam>8.20&&ebeam<11.8&&(them>0.034907&&thep>0.034907))&&abs(trf-tbeam)<2.&&them*467.4<83.&&cos(phem)*them*467.4<-11.47&&sin(phem)*them*467.4<-32.61","same");
//e+ in TRD
JP->Draw("sin(phep)*thep*467.4:cos(phep)*thep*467.4>>in_trd_ep","(Minv>3.0555&&Minv<3.1365&&((pep/bcalep<1.173300&&pep/bcalep>0.875700&&bcalep>0&&pbcalep*sin(thep)>0.030000)||(pep/fcalep<1.233800&&pep/fcalep>0.915800&&fcalep>0))&&((pem/bcalem<1.173300&&pem/bcalem>0.875700&&bcalem>0&&pbcalem*sin(them)>0.030000)||(pem/fcalem<1.233800&&pem/fcalem>0.915800&&fcalem>0))&&chi2>0&&chi2<5000.000000&&ebeam>8.20&&ebeam<11.8&&(them>0.034907&&thep>0.034907))&&abs(trf-tbeam)<2.&&thep*467.4<83.&&cos(phep)*thep*467.4<-11.47&&sin(phep)*thep*467.4<-32.61","same");
all_ep->SetTitle("e^{+}e^{-} Pairs in All Acceptance at TRD Z Plane");
fbox->Draw("same");
l0->Draw("same");
c0->SaveAs("JPsiAtTRDPlane.pdf");

//JP->Scan("run_no:event_no:ebeam:Minv:pem_m:pep_m:pp_m:fcalem:pem:cos(phem)*them*467.4:sin(phem)*them*467.4:thp*53.7:thep*57.3:them_m*57.3:them*57.3:phem_m*57.3:phem*57.3","(Minv>3.0555&&Minv<3.1365&&((pep/bcalep<1.173300&&pep/bcalep>0.875700&&bcalep>0&&pbcalep*sin(thep)>0.030000)||(pep/fcalep<1.233800&&pep/fcalep>0.915800&&fcalep>0))&&((pem/bcalem<1.173300&&pem/bcalem>0.875700&&bcalem>0&&pbcalem*sin(them)>0.030000)||(pem/fcalem<1.233800&&pem/fcalem>0.915800&&fcalem>0))&&chi2>0&&chi2<5000.000000&&ebeam>8.20&&ebeam<11.8&&(them>0.034907&&thep>0.034907))&&abs(trf-tbeam)<2.&&them*467.4<83.&&cos(phem)*them*467.4<-8.&&sin(phem)*them*467.4<-30.","same");
//JP->Scan("run_no:event_no:ebeam:Minv:pep_m:pem_m:pp_m:fcalep:pep:cos(phep)*thep*467.4:sin(phep)*thep*467.4:thep_m*57.3:thep*57.3:phep_m*57.3:phep*57.3","(Minv>3.0555&&Minv<3.1365&&((pep/bcalep<1.173300&&pep/bcalep>0.875700&&bcalep>0&&pbcalep*sin(thep)>0.030000)||(pep/fcalep<1.233800&&pep/fcalep>0.915800&&fcalep>0))&&((pem/bcalem<1.173300&&pem/bcalem>0.875700&&bcalem>0&&pbcalem*sin(them)>0.030000)||(pem/fcalem<1.233800&&pem/fcalem>0.915800&&fcalem>0))&&chi2>0&&chi2<5000.000000&&ebeam>8.20&&ebeam<11.8&&(them>0.034907&&thep>0.034907))&&abs(trf-tbeam)<2.&&thep*467.4<83.&&cos(phep)*thep*467.4<-8.&&sin(phep)*thep*467.4<-30.","same");
//JP->Scan("run_no:event_no","(Minv>3.0555&&Minv<3.1365&&((pep/bcalep<1.173300&&pep/bcalep>0.875700&&bcalep>0&&pbcalep*sin(thep)>0.030000)||(pep/fcalep<1.233800&&pep/fcalep>0.915800&&fcalep>0))&&((pem/bcalem<1.173300&&pem/bcalem>0.875700&&bcalem>0&&pbcalem*sin(them)>0.030000)||(pem/fcalem<1.233800&&pem/fcalem>0.915800&&fcalem>0))&&chi2>0&&chi2<5000.000000&&ebeam>8.20&&ebeam<11.8&&(them>0.034907&&thep>0.034907))&&abs(trf-tbeam)<2.&&them*467.4<83","");
//JP->Scan("run_no:event_no","(Minv>3.0555&&Minv<3.1365&&((pep/bcalep<1.173300&&pep/bcalep>0.875700&&bcalep>0&&pbcalep*sin(thep)>0.030000)||(pep/fcalep<1.233800&&pep/fcalep>0.915800&&fcalep>0))&&((pem/bcalem<1.173300&&pem/bcalem>0.875700&&bcalem>0&&pbcalem*sin(them)>0.030000)||(pem/fcalem<1.233800&&pem/fcalem>0.915800&&fcalem>0))&&chi2>0&&chi2<5000.000000&&ebeam>8.20&&ebeam<11.8&&(them>0.034907&&thep>0.034907))&&abs(trf-tbeam)<2.&&thep*467.4<83.","");
}
