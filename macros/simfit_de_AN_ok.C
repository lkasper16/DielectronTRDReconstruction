{
#define pmass 0.9383  
TFile *fin = TFile::Open("s3de6_roo.root","update");

TH1F *dminva;
TH1F *dminvb;
TH1F *dminv;

#include "RooChi2Var.h"
#include "RooFitResult.h"
using namespace RooFit;
TCanvas *myc;
myc = new TCanvas("myc", "Event", 800, 800);
char htit[128];
//  myc->Divide(4,4);
gROOT->SetStyle("Plain");
//gStyle->SetOptStat(0);

// initialization
RooAbsReal::defaultIntegratorConfig()->method1D().setLabel("RooAdaptiveGaussKronrodIntegrator1D");
RooFitResult *rf;

// E-slices data
TH1F *dejpsi = (TH1F*)fin->Get("dejpsi");
dejpsi->SetName("dejpsi");
TAxis *dx = (TAxis*)dejpsi->GetXaxis();
int nebins=dx->GetNbins();
float e1=dx->GetXmin();
float e2=dx->GetXmax();
float step=(e2-e1)/nebins;
e2-=step;
cout<<" nebins,e1,e2= "<<nebins<<" "<<e1<<" "<<e2<<endl;

//for (int i=nebins-1;i>=0;i--) // start e-loop
//i=6;

RooWorkspace* w = new RooWorkspace("w");
// observable(s)
RooRealVar mass("mass", "M(e^{+}e^{-}), GeV", 2.5, 3.25);
//RooRealVar mcsigma("mcsigma","mcsigma",0.005,0.025);
//RooRealVar mcsigma("mcsigma","mcsigma",ephi->GetBinContent(i+1));
mass.setBins(60);

// classes of data - we use these to tag the different data sets and PDFs
RooCategory dtype("dtype", "dtype");
dtype.defineType("data",1);
dtype.defineType("acc",2);     // accidentals

w->import(RooArgSet(mass,dtype));
//w->import(RooArgSet(mcsigma));

// accidentals
w->factory("Chebychev::acc_shape(mass,a0[0.,-1.e4,1.e4])");
// background and signals
w->factory("Chebychev::bkgd(mass,{Cbkgd[0.,-100.,100.]})");
//cok free sigma w->factory("Gaussian::jpsi(mass,mean[3.091,3.05,3.15],sigma[0.013,0.005,0.017])");
w->factory("Gaussian::jpsi(mass,mean[3.083,3.08,3.10],sigma[0.015,0.005,0.025])");

// build model for data:  accidental shape + background + signals
w->factory("SUM::model(nacc[50,0,1e4]*acc_shape,  Nbkgd[50,0,1e4]*bkgd, njpsi[50,0,1e4]*jpsi)");
w->factory("SUM::my_model(Nbkgd[50,0,1e4]*bkgd, Njpsi[50,0,1e4]*jpsi)");
//w->factory("SUM::model(f1[0.2,0.,1]*acc_shape,  f2[0.4,0.,1]*bkgd, jpsi)");

// simultaneous fit of data and accidental spectrum
w->factory("SIMUL::smodel( dtype, data=model, acc=acc_shape )");

char hnam[128];
char hnama[128];
char hnamb[128];

//  float ebeam=dx->GetBinCenter(i+1);
//  float mlim=sqrt(2.*ebeam*pmass+pmass*pmass)-pmass;
//  cout<<" ibin, minv limit= "<<i<<" "<<mlim<<endl;
float mup=3.25;
float mlow=2.5;
sprintf(hnam,"dminv");
sprintf(hnama,"dminva");
sprintf(hnamb,"dminvb");
dminva = (TH1F*)fin->Get(hnama);
dminvb = (TH1F*)fin->Get(hnamb);
dminv = (TH1F*)fin->Get(hnam);

dminva->Rebin(2);
dminvb->Rebin(2);
dminv->Rebin(2);

dminv->Add(dminva,dminvb,1.,-1.);

RooDataHist *data = new RooDataHist("data","data",RooArgSet(mass),dminva);
RooDataHist *data_acc = new RooDataHist("data_acc","data_acc",RooArgSet(mass),dminvb);
RooDataHist *data_sub = new RooDataHist("data_sub","data_sub",RooArgSet(mass),dminva);

// combine the two sets
RooDataHist combData("combData","combined data",mass,RooFit::Index(dtype),RooFit::Import("data",*data),
  	RooFit::Import("acc",*data_acc)) ;

// do the fit
//const RooFitResult  *rf = w->pdf("smodel")->fitTo(combData,SumW2Error(kFALSE),Hesse(1),Range(mlow,mup));
rf = w->pdf("my_model")->fitTo(*data_sub,SumW2Error(0),Range(mlow,mup),Hesse(0),Save(1));
//w->pdf("my_model")->Print("v");
//const RooFitResult  *rf = w->pdf("my_model")->fitTo(*data_sub,Range(2.9,mup),SumW2Error(0),Hesse(1),Extended(),Verbose(0),PrintLevel(-1000));

float njp = w->var("Njpsi")->getVal();
float enjp = w->var("Njpsi")->getError();
//cedit free sigma float sjp = w->var("sigma")->getVal();
//cedit free sigma float esjp = w->var("sigma")->getError();
//float sjp = w->var("sigma")->getVal();
//float esjp = w->var("sigma")->getError();
float mjp = w->var("mean")->getVal();
float emjp = w->var("mean")->getError();
//float sjp;
//float esjp;

// make some plots
sprintf(htit," J/#psi ");

RooPlot* massframe = mass.frame(RooFit::Title(htit)) ;
//data->plotOn(massframe) ;
data_sub->plotOn(massframe) ;
//w->pdf("model")->plotOn(massframe);
w->pdf("my_model")->paramOn(massframe);
w->pdf("my_model")->plotOn(massframe);

RooChi2Var chi2("chi2","chi2",*(w->pdf("my_model")),*data_sub) ;
cout <<" ........... chi2=........"<<chi2.getVal() << endl ;
//rf->Print("v");

massframe->Draw();
myc->Update();
//cedit free sigma dephi->SetBinContent(i+1,sjp);
//cedit free sigma dephi->SetBinError(i+1,esjp);
//ephi->SetBinContent(i+1,sjp);

//dminv->Write();
dejpsi->Write();
//cedit free sigma dephi->Write();
dephi->Write();

myc->SaveAs("TRD_JPsi.pdf");

}
