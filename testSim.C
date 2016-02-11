#include "TDirectory.h"
#include "TPad.h"
#include "TFile.h"
#include "TH1F.h"
#include "TTree.h"
#include "THStack.h"
#include "TH2.h"
#include "TF1.h"
#include "TLine.h"
#include "TCut.h"
#include "TROOT.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TGraphErrors.h"
#include "HiggsAnalysis/CombinedLimit/interface/HZZ2L2QRooPdfs.h"
#include "ZZAnalysis/AnalysisStep/interface/Category.h"
#include "Math/GenVector/LorentzVector.h"
#include "Math/GenVector/PtEtaPhiM4D.h"
#include <stdlib.h>
#include <stdio.h>
#include <cmath>

#include "Math/GenVector/LorentzVector.h"
#include "Math/GenVector/PtEtaPhiM4D.h"

#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooDataHist.h"
#include "RooGaussian.h"
#include "RooCBShape.h"
#include "RooSimultaneous.h"
#include "RooCategory.h"
#include "RooWorkspace.h"
#include "RooHist.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "RooPlot.h"
#include "TAxis.h"
#include "TPaveText.h"
#include "RooAddPdf.h"
#include "RooBreitWigner.h"
#include "RooFitResult.h"
#include "RooFFTConvPdf.h"
#include "RooAddition.h"
#include "RooMinuit.h"
#include "Math/MinimizerOptions.h"
#include <iomanip>
#include "RooAbsCollection.h"
#include "RooWorkspace.h"

using namespace RooFit;


#include <string>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <sstream>
    void testSim()
                            
 {

int massBin[40] ={125,150,200,250,300,350,400,450,500,550,600,650,700,750,800,850,900,950,1000,1050,1100};
int maxMassBin =21;
int channels = 1;
int inputfiles[40]={115,120,124,125,126,130,135,140,145,150,155,160,165,170,175,180,190,210,230,250,270,300,350,400,450,500,550,600,700,800,900,1000};
int Nfiles=32;

 // ------ root settings ---------
  gROOT->Reset();
  gROOT->SetStyle("Plain");
  gStyle->SetPadGridX(kFALSE);
  gStyle->SetPadGridY(kFALSE);
  gStyle->SetOptStat("iourme");
  gStyle->SetOptFit(11);
  gStyle->SetPadLeftMargin(0.14);
  gStyle->SetPadRightMargin(0.06);
  // ------------------------------

  ROOT::Math::MinimizerOptions::SetDefaultTolerance( 1.E-7);

  float m4l,genM;
  Short_t z1flav, z2flav;
  float weight;

  double xMin,xMax;
  xMin = -100;
  xMax = 100;
  cout << "Fit range: [" << xMin << " , " << xMax << "]." << endl;

  RooRealVar x("reso","m_{reco}-m_{true}",0.,xMin,xMax,"GeV");
  x.setBins(70);
  RooRealVar w("myW","myW",1.0,-1000.,1000.);
  RooCategory massrc("massrc","massrc");
  char tempmass[100];


  for (int i=0; i<maxMassBin; i++) {
    sprintf(tempmass,"mh%d",massBin[i]);
    massrc.defineType(tempmass,massBin[i]);
  }

  RooArgSet ntupleVarSet(x,w,massrc);
  RooDataSet dataset("resoM","resoM",ntupleVarSet,WeightVar("myW"));

  TString inputDir = "root://lxcms03://data3/Higgs/160121/";

  vector<TString> files;

  char inputfile[100];
  for (int i=0; i<Nfiles; i++) {
    sprintf(inputfile,"ggH%d/ZZ4lAnalysis.root",inputfiles[i]);
    files.push_back(inputfile);
  }

TChain *ggTree = new TChain("ZZTree/candTree");

for (vector<TString>::const_iterator file = files.begin(); file!=files.end(); ++file) {
ggTree->Add(inputDir+(*file));
}

    int  nentries = ggTree->GetEntries();

    //--- ggTree part
    ggTree->SetBranchAddress("ZZMass",&m4l);
    ggTree->SetBranchAddress("GenHMass",&genM);
    ggTree->SetBranchAddress("Z1Flav",&z1flav);
    ggTree->SetBranchAddress("Z2Flav",&z2flav);
    ggTree->SetBranchAddress("genHEPMCweight",&weight);


    for(int k=0; k<nentries; k++){
      ggTree->GetEvent(k);

      if(channels==0 && z1flav*z2flav != 28561) continue;
      if(channels==1 && z1flav*z2flav != 14641) continue;
      if(channels==2 && z1flav*z2flav != 20449) continue;
      if (weight <= 0 ) cout << "Warning! Negative weight events" << endl;

     for (int i=0; i<maxMassBin; i++) {
      ntupleVarSet.setCatIndex("massrc",massBin[i]);
      ntupleVarSet.setRealValue("reso",m4l-genM);
      ntupleVarSet.setRealValue("myW",weight);
      if(x.getVal()>xMin && x.getVal()<xMax&&genM>(massBin[i]*0.99)&&genM<(massBin[i]*1.01)) 
       dataset.add(ntupleVarSet, weight);

      //--------

    }
  }

  cout << "dataset n entries: " << dataset.sumEntries() << endl;

  RooSimultaneous rs("rs","rs",massrc);
  RooFormulaVar* mean[40];
  RooFormulaVar* sigma[40];
  RooFormulaVar* a1[40];
  RooFormulaVar* a2[40];
  RooFormulaVar* n1[40];
  RooFormulaVar* n2[40];

  RooDoubleCB* DCBall[40];

  RooRealVar mean_p0("mean_p0","mean_p0",0.,-1., 1.) ;
  RooRealVar sigma_p0("sigma_p0","sigma_p0",1.5, 0, 30);
  RooRealVar a1_p0("a1_p0","a1_p0", 1.46, 0.5, 5.);
  RooRealVar n1_p0("n1_p0","n_p0", 1.92, 0, 10);
  RooRealVar a2_p0("a2_p0","a2_p0", 1.46, 1, 10);
  RooRealVar n2_p0("n2_p0","n2_p0", 20., 1., 50.);

  RooDoubleCB DCBall_400("DCBall","Double Crystal ball",x,mean_p0,sigma_p0,a1_p0,n1_p0,a2_p0,n2_p0);
  RooDataSet* dataset_400= (RooDataSet*)dataset.reduce("massrc == massrc::mh400");
  RooFitResult* fitres_400 = (RooFitResult*)DCBall_400.fitTo(*dataset_400,SumW2Error(1),Range(xMin,xMax),Strategy(2),NumCPU(8),Save(true));
  RooArgSet * params400 = DCBall_400.getParameters(x);
  params400->writeToFile("prefit400_param.txt") ;

  RooRealVar mean_p1("mean_p1","mean_p1",0, -1.5, 1.5);
  RooRealVar sigma_p1("sigma_p1","sigma_p1",0,-0.5, 0.5);
  RooRealVar a1_p1("a1_p1","a1_p1", 0, -0.5, 0.5);
  RooRealVar n1_p1("n1_p1","n_p1", 0, -0.5, 0.5);
  RooRealVar a2_p1("a2_p1","a2_p1", 0, -0.5, 0.5);
  RooRealVar n2_p1("n2_p1","n2_p1", 0, -0.5, 0.5);


  for (int i=0; i<maxMassBin; i++) {
    char formulamass[200];
    int massdiff = massBin[i] - 400;
    sprintf(formulamass,"@0+(@1*%d)",massdiff);
    mean[i] = new RooFormulaVar("mean_CB","mean_CB",formulamass,RooArgList(mean_p0,mean_p1));
    sigma[i] = new RooFormulaVar("sigma_CB","sigma_CB",formulamass,RooArgList(sigma_p0,sigma_p1));
    a1[i] = new RooFormulaVar("a1_CB","a1_CB",formulamass,RooArgList(a1_p0,a1_p1));
    n1[i] = new RooFormulaVar("n1_CB","n1_CB",formulamass,RooArgList(n1_p0,n1_p1));
    a2[i] = new RooFormulaVar("a2_CB","a2_CB",formulamass,RooArgList(a2_p0,a2_p1));
    n2[i] = new RooFormulaVar("n2_CB","n2_CB",formulamass,RooArgList(n2_p0,n2_p1));
    DCBall[i] = new RooDoubleCB("DCBall","Double Crystal ball",x,*mean[i],*sigma[i],*a1[i],*n1[i],*a2[i],*n2[i]);
    char tempmass[100];
    sprintf(tempmass,"mh%d",massBin[i]);
    char tempmass2[100];
    sprintf(tempmass2,"massrc == massrc::mh%d",massBin[i]);
    rs.addPdf(*DCBall[i], tempmass);
  }


  RooFitResult *fitressim = (RooFitResult*)rs.fitTo(dataset,SumW2Error(1),Range(xMin,xMax),Strategy(2),NumCPU(8),Save(true));

  mean_p1.setConstant(kTRUE);
  sigma_p1.setConstant(kTRUE);
  a1_p1.setConstant(kTRUE);
  n1_p1.setConstant(kTRUE);
  a2_p1.setConstant(kTRUE);
  n2_p1.setConstant(kTRUE);

  RooArgSet * params2 = rs.getParameters(RooArgList(x,massrc));
  params2->writeToFile("simul_param.txt") ;

  for (int i=0; i<maxMassBin; i++) {
    TCanvas *c1 = new TCanvas("c1","c1",725,725);
    RooPlot* xframe = x.frame() ;
    xframe->SetTitle("");
    xframe->SetName("m4lplot");
    xframe->GetYaxis()->SetTitleOffset(1.5);
    char tempmass[100];
    sprintf(tempmass,"massrc == massrc::mh%d",massBin[i]);
    char tempmass2[100];
    sprintf(tempmass2,"mh%d",massBin[i]);
   

    int col;
    if(channels==0) col=kOrange+7;
    if(channels==1) col=kAzure+2;
    if(channels==2) col=kGreen+3;

    dataset.plotOn(xframe,DataError(RooAbsData::SumW2), MarkerStyle(kOpenCircle), MarkerSize(1.1), Cut(tempmass) );
    rs.plotOn(xframe,LineColor(col),Slice(massrc,tempmass2),ProjWData(massrc,dataset));
    c1->cd();
    xframe->Draw();
    char filename[100];
    sprintf(filename,"Resolution_MH%d.png",massBin[i]);
    c1->SaveAs(filename);
  }
}
