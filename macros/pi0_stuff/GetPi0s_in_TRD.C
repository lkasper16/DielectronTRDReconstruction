{
bool saveText=true;
int run_no;
int event_no;
//TFile *_file0 = TFile::Open("jpsi_lp_out_test.root");
TFile *_file0 = TFile::Open("pi0_lp_out_test.root");
TH2F *in_trd_ep = new TH2F("in_trd_ep","Positrons projected at the TRD; X at TRD Plane [cm]; Y at TRD Plane [cm]",1000,-100.,100.,1000,-100.,100.);
TH2F *in_trd_em = new TH2F("in_trd_em","Electrons projected at the TRD; X at TRD Plane [cm]; Y at TRD Plane [cm]",1000,-100.,100.,1000,-100.,100.);
TH2F *all_ep = new TH2F("all_ep","All positrons projected at the TRD plane; X at TRD Plane [cm]; Y at TRD Plane [cm]",1000,-100.,100.,1000,-100.,100.);
TH2F *all_em = new TH2F("all_em","All electrons projected at the TRD plane; X at TRD Plane [cm]; Y at TRD Plane [cm]",1000,-100.,100.,1000,-100.,100.);

//int n_ep=0, n_em=0;

TLegend *l0 = new TLegend(0.7,0.79,0.9,0.9);
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

TTree *Pi0 = (TTree*)_file0->Get("Pi0");
TBox *fbox = new TBox(-83.47,-68.6,-11.47,-32.61); //For TRD acceptance visualization
fbox->SetFillStyle(0);
fbox->SetLineWidth(2);
TCanvas *c0 = new TCanvas("c0","c0",1300,1000);
c0->cd();
gPad->SetGridy();
gPad->SetGridx();
gStyle->SetOptStat(0);



//e+ in all acceptance
Pi0->Draw("sin(phep)*thep*467.4:cos(phep)*thep*467.4>>all_ep","(Cinv>0.0849768&&Cinv<0.1849768&&((pep/bcalep<1.173300&&pep/bcalep>0.875700&&bcalep>0&&pbcalep*sin(thep)>0.030000)||(pep/fcalep<1.233800&&pep/fcalep>0.915800&&fcalep>0))&&((pem/bcalem<1.173300&&pem/bcalem>0.875700&&bcalem>0&&pbcalem*sin(them)>0.030000)||(pem/fcalem<1.233800&&pem/fcalem>0.915800&&fcalem>0))&&chi2>0&&chi2<5000.000000&&ebeam>8.20&&ebeam<11.8&&(them>0.034907&&thep>0.034907))&&abs(trf-tbeam)<2.&&thep*467.4<83.","");

//e- in all acceptance
Pi0->Draw("sin(phem)*them*467.4:cos(phem)*them*467.4>>all_em","(Cinv>0.0849768&&Cinv<0.1849768&&((pep/bcalep<1.173300&&pep/bcalep>0.875700&&bcalep>0&&pbcalep*sin(thep)>0.030000)||(pep/fcalep<1.233800&&pep/fcalep>0.915800&&fcalep>0))&&((pem/bcalem<1.173300&&pem/bcalem>0.875700&&bcalem>0&&pbcalem*sin(them)>0.030000)||(pem/fcalem<1.233800&&pem/fcalem>0.915800&&fcalem>0))&&chi2>0&&chi2<5000.000000&&ebeam>8.20&&ebeam<11.8&&(them>0.034907&&thep>0.034907))&&abs(trf-tbeam)<2.&&them*467.4<83.","same");

//e- in TRD
Pi0->Draw("sin(phem)*them*467.4:cos(phem)*them*467.4>>in_trd_em","(Cinv>0.0849768&&Cinv<0.1849768&&((pep/bcalep<1.173300&&pep/bcalep>0.875700&&bcalep>0&&pbcalep*sin(thep)>0.030000)||(pep/fcalep<1.233800&&pep/fcalep>0.915800&&fcalep>0))&&((pem/bcalem<1.173300&&pem/bcalem>0.875700&&bcalem>0&&pbcalem*sin(them)>0.030000)||(pem/fcalem<1.233800&&pem/fcalem>0.915800&&fcalem>0))&&chi2>0&&chi2<5000.000000&&ebeam>8.20&&ebeam<11.8&&(them>0.034907&&thep>0.034907))&&abs(trf-tbeam)<2.&&them*467.4<83.&&cos(phem)*them*467.4<-11.47&&cos(phem)*them*467.4>-83.47&&sin(phem)*them*467.4<-32.61&&sin(phem)*them*467.4>-68.6","same");
//Pi0->Draw("sin(phem)*them*467.4:cos(phem)*them*467.4>>in_trd_em","cos(phem)*them*467.4<-11.47&&cos(phem)*them*467.4>-83.47&&sin(phem)*them*467.4<-32.61&&sin(phem)*them*467.4>-68.6","same");

//e+ in TRD
Pi0->Draw("sin(phep)*thep*467.4:cos(phep)*thep*467.4>>in_trd_ep","(Cinv>0.0849768&&Cinv<0.1849768&&((pep/bcalep<1.173300&&pep/bcalep>0.875700&&bcalep>0&&pbcalep*sin(thep)>0.030000)||(pep/fcalep<1.233800&&pep/fcalep>0.915800&&fcalep>0))&&((pem/bcalem<1.173300&&pem/bcalem>0.875700&&bcalem>0&&pbcalem*sin(them)>0.030000)||(pem/fcalem<1.233800&&pem/fcalem>0.915800&&fcalem>0))&&chi2>0&&chi2<5000.000000&&ebeam>8.20&&ebeam<11.8&&(them>0.034907&&thep>0.034907))&&abs(trf-tbeam)<2.&&thep*467.4<83.&&cos(phep)*thep*467.4<-11.47&&cos(phep)*thep*467.4>-83.47&&sin(phep)*thep*467.4<-32.61&&sin(phep)*thep*467.4>-68.6","same");
//Pi0->Draw("sin(phem)*them*467.4:cos(phem)*them*467.4>>in_trd_em","cos(phem)*them*467.4<-11.47&&cos(phem)*them*467.4>-83.47&&sin(phem)*them*467.4<-32.61&&sin(phem)*them*467.4>-68.6","same");

all_ep->SetTitle("e^{+}e^{-} Pairs in All Acceptance at TRD Z Plane");
fbox->Draw("same");
l0->Draw("same");
c0->SaveAs("Pi0AtTRDPlane.pdf");

if (saveText) {
	ofstream outF("pi0_run_ev.txt");
	Pi0->SetScanField(0);
	
	//e- in TRD
	Pi0->Scan("run_no:event_no:ebeam:Cinv:pem_m:pep_m:pp_m:fcalem:pem:cos(phem)*them*467.4:sin(phem)*them*467.4:thp*53.7:thep*57.3:them_m*57.3:them*57.3:phem_m*57.3:phem*57.3","(Cinv>0.0849768&&Cinv<0.1849768&&((pep/bcalep<1.173300&&pep/bcalep>0.875700&&bcalep>0&&pbcalep*sin(thep)>0.030000)||(pep/fcalep<1.233800&&pep/fcalep>0.915800&&fcalep>0))&&((pem/bcalem<1.173300&&pem/bcalem>0.875700&&bcalem>0&&pbcalem*sin(them)>0.030000)||(pem/fcalem<1.233800&&pem/fcalem>0.915800&&fcalem>0))&&chi2>0&&chi2<5000.000000&&ebeam>8.20&&ebeam<11.8&&(them>0.034907&&thep>0.034907))&&abs(trf-tbeam)<2.&&them*467.4<83.&&cos(phem)*them*467.4<-11.47&&cos(phem)*them*467.4>-83.47&&sin(phem)*them*467.4<-32.61&&sin(phem)*them*467.4>-68.6","");
	
	//e+ in TRD
	Pi0->Scan("run_no:event_no:ebeam:Cinv:pep_m:pem_m:pp_m:fcalep:pep:cos(phep)*thep*467.4:sin(phep)*thep*467.4:thep_m*57.3:thep*57.3:phep_m*57.3:phep*57.3","(Cinv>0.0849768&&Cinv<0.1849768&&((pep/bcalep<1.173300&&pep/bcalep>0.875700&&bcalep>0&&pbcalep*sin(thep)>0.030000)||(pep/fcalep<1.233800&&pep/fcalep>0.915800&&fcalep>0))&&((pem/bcalem<1.173300&&pem/bcalem>0.875700&&bcalem>0&&pbcalem*sin(them)>0.030000)||(pem/fcalem<1.233800&&pem/fcalem>0.915800&&fcalem>0))&&chi2>0&&chi2<5000.000000&&ebeam>8.20&&ebeam<11.8&&(them>0.034907&&thep>0.034907))&&abs(trf-tbeam)<2.&&thep*467.4<83.&&cos(phep)*thep*467.4<-11.47&&cos(phep)*thep*467.4>-83.47&&sin(phep)*thep*467.4<-32.61&&sin(phep)*thep*467.4>-68.6","");
	
	//e- in all acceptance
	//Pi0->Scan("run_no:event_no","(Cinv>0.0849768&&Cinv<0.1849768&&((pep/bcalep<1.173300&&pep/bcalep>0.875700&&bcalep>0&&pbcalep*sin(thep)>0.030000)||(pep/fcalep<1.233800&&pep/fcalep>0.915800&&fcalep>0))&&((pem/bcalem<1.173300&&pem/bcalem>0.875700&&bcalem>0&&pbcalem*sin(them)>0.030000)||(pem/fcalem<1.233800&&pem/fcalem>0.915800&&fcalem>0))&&chi2>0&&chi2<5000.000000&&ebeam>8.20&&ebeam<11.8&&(them>0.034907&&thep>0.034907))&&abs(trf-tbeam)<2.&&them*467.4<83","");
	
	//e+ in all acceptance
	//Pi0->Scan("run_no:event_no","(Cinv>0.0849768&&Cinv<0.1849768&&((pep/bcalep<1.173300&&pep/bcalep>0.875700&&bcalep>0&&pbcalep*sin(thep)>0.030000)||(pep/fcalep<1.233800&&pep/fcalep>0.915800&&fcalep>0))&&((pem/bcalem<1.173300&&pem/bcalem>0.875700&&bcalem>0&&pbcalem*sin(them)>0.030000)||(pem/fcalem<1.233800&&pem/fcalem>0.915800&&fcalem>0))&&chi2>0&&chi2<5000.000000&&ebeam>8.20&&ebeam<11.8&&(them>0.034907&&thep>0.034907))&&abs(trf-tbeam)<2.&&thep*467.4<83.","");
}
}
