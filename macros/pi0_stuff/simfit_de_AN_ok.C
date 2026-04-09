{
#define pmass 0.9383  
TFile *fin = TFile::Open("s3de6_xi_roo.root","update");

TH1F *dminva;
TH1F *dminvb;
TH1F *dminv;
#include "RooChi2Var.h"
#include "RooFitResult.h"
using namespace RooFit;
TCanvas *myc;
myc = new TCanvas("myc", "Event", 800, 800);
char htit[128];
gROOT->SetStyle("Plain");

// initialization
RooAbsReal::defaultIntegratorConfig()->method1D().setLabel("RooAdaptiveGaussKronrodIntegrator1D");
RooFitResult *rf;

// E-slices data
/*
TH1F *dejpsi = (TH1F*)fin->Get("dejpsi");
dejpsi->SetName("dejpsi");
TAxis *dx = (TAxis*)dejpsi->GetXaxis();
int nebins=dx->GetNbins();
float e1=dx->GetXmin();
float e2=dx->GetXmax();
float step=(e2-e1)/nebins;
e2-=step;
cout<<" nebins,e1,e2= "<<nebins<<" "<<e1<<" "<<e2<<endl;
*/

RooWorkspace* w = new RooWorkspace("w");
// observable(s)
RooRealVar mass("mass", "M(e^{+}e^{-}#gamma), GeV", 0.0, 0.2);
//RooRealVar mcsigma("mcsigma","mcsigma",0.005,0.025);
mass.setBins(30); //60

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
//w->factory("Gaussian::jpsi(mass,mean[3.083,3.08,3.10],sigma[0.015,0.005,0.025])"); //?? How are these determined?
w->factory("Gaussian::pi0(mass,mean[0.14,0.12,0.165],sigma[0.015,0.005,0.025])");

// build model for data:  accidental shape + background + signals
w->factory("SUM::model(nacc[50,0,1e4]*acc_shape,  Nbkgd[50,0,1e4]*bkgd, npi0[50,0,1e4]*pi0)");
w->factory("SUM::my_model(Nbkgd[50,0,1e4]*bkgd, Npi0[50,0,1e4]*pi0)");

// simultaneous fit of data and accidental spectrum
w->factory("SIMUL::smodel( dtype, data=model, acc=acc_shape )");

char hnam[128];
char hnama[128];
char hnamb[128];

float mup=0.2;
float mlow=0.0;
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
RooDataHist combData("combData","combined data",mass,RooFit::Index(dtype),RooFit::Import("data",*data),RooFit::Import("acc",*data_acc));

// do the fit
rf = w->pdf("my_model")->fitTo(*data_sub,SumW2Error(0),Range(mlow,mup),Hesse(0),Save(1));

float njp = w->var("Npi0")->getVal();
float enjp = w->var("Npi0")->getError();
float mjp = w->var("mean")->getVal();
float emjp = w->var("mean")->getError();

// make some plots
sprintf(htit," #pi^{0} ");

RooPlot* massframe = mass.frame(RooFit::Title(htit)) ;
//data->plotOn(massframe) ;
data_sub->plotOn(massframe) ;
w->pdf("my_model")->paramOn(massframe);
w->pdf("my_model")->plotOn(massframe);

//RooChi2Var chi2("chi2","chi2",*(w->pdf("my_model")),*data_sub) ;
//cout <<" ........... chi2=........"<<chi2.getVal() << endl ;

massframe->Draw();
myc->Update();

//dminv->Write();
//dejpsi->Write();
//dephi->Write();

myc->SaveAs("TRD_Pi0.pdf");

}
