#include <TH1.h>
#include <TH2F.h>
#include <TLegend.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <THStack.h>
#include <TPostScript.h>
#include <TPaveText.h>
#include <TGraph.h>
#include <TMultiGraph.h>
#include <TChain.h>
#include <TBranch.h>
#include <TFile.h>

#include "TMath.h"
#include "TF1.h"

using namespace std;

const TString StorageDirPrefix="file:/home/cms/jruizalv/work/CMSSW_5_3_9_patch2/src/Extractors/PatExtractor/bin/";
const int NOS=5;
//const TString Samples[NOS] = {"Diboson.root", "SingleTop.root", "TTJets.root", "Zjets.root", "Wjets.root", "QCD.root", "Signal.root"};
const TString Samples[NOS] = {"Diboson.root", "SingleTop.root", "TTJets.root", "QCD.root", "Signal.root"};

void AllVariablesPlotter()
{
  TFile *CurrentFile[NOS];
  TH1F *TprimeHistos[NOS];
  TH1F *LeadingJetPT[NOS];
  TH1F *THT[NOS];
  TH1F *DRHjets[NOS];
  TH1F *DRWjets[NOS];
  TH1F *Hpt[NOS];
  TH1F *Tpt[NOS];
  TH1F *DRWH[NOS];
  TH1F *DPHjets[NOS];
  TH1F *DPWjets[NOS];
  TH1F *DPTjets[NOS];
  TH1F *HiggsMass[NOS];
  TH1F *RelHT[NOS];
  TH1F *DRTH[NOS];
  TH1F *PtNormalizedMass[NOS];
  
  for (int i=0; i<NOS; i++)
    {
      CurrentFile[i] = new TFile(StorageDirPrefix + Samples[i], "READ");
      if ( CurrentFile[i]->IsOpen() ) printf( Samples[i] + " File opened successfully\n");
      if (i==0) 
	{
	  TprimeHistos[0]->SetDefaultSumw2(); TprimeHistos[0]= (TH1F*)gDirectory->Get("TprimeMass0");
	  LeadingJetPT[0]->SetDefaultSumw2(); LeadingJetPT[0]= (TH1F*)gDirectory->Get("jet1_pt0");
	  THT[0]->SetDefaultSumw2(); THT[0]= (TH1F*)gDirectory->Get("THT0");
	}
      else if (i==1) 
	{
	  TprimeHistos[1]->SetDefaultSumw2(); TprimeHistos[1]= (TH1F*)gDirectory->Get("TprimeMass3");
	  LeadingJetPT[1]->SetDefaultSumw2(); LeadingJetPT[1]= (TH1F*)gDirectory->Get("jet1_pt3");
	  THT[1]->SetDefaultSumw2(); THT[1]= (TH1F*)gDirectory->Get("THT3");
	}
      else if (i==2) 
	{
	  TprimeHistos[2]->SetDefaultSumw2(); TprimeHistos[2]= (TH1F*)gDirectory->Get("TprimeMass9");
	  LeadingJetPT[2]->SetDefaultSumw2(); LeadingJetPT[2]= (TH1F*)gDirectory->Get("jet1_pt9");
	  THT[2]->SetDefaultSumw2(); THT[2]= (TH1F*)gDirectory->Get("THT9");
	}
      else if (i==3) 
	{
	  TprimeHistos[3]->SetDefaultSumw2(); TprimeHistos[3]= (TH1F*)gDirectory->Get("TprimeMass10");
	  LeadingJetPT[3]->SetDefaultSumw2(); LeadingJetPT[3]= (TH1F*)gDirectory->Get("jet1_pt10");
	  THT[3]->SetDefaultSumw2(); THT[3]= (TH1F*)gDirectory->Get("THT10");
	}
      else if (i==4) 
	{
	  TprimeHistos[4]->SetDefaultSumw2(); TprimeHistos[4]= (TH1F*)gDirectory->Get("TprimeMass15");
	  LeadingJetPT[4]->SetDefaultSumw2(); LeadingJetPT[4]= (TH1F*)gDirectory->Get("jet1_pt15");
	  THT[4]->SetDefaultSumw2(); THT[4]= (TH1F*)gDirectory->Get("THT15");
	}
    }
  
  THStack *BKGandSignal = new THStack("BKGandSignal", "BKG and signal for 5 jets; M_{5j} GeV; Events");
  THStack *BKGLJPT = new THStack("BKGandSignalLJPT", "BKG and signal for Leading Jet PT; pT(j_{1}) GeV; Events");
  THStack *BKGTHT = new THStack("BKGandSignalTHT", "BKG and signal for HT; HT GeV; Events");
  for (int i=0; i<NOS; i++)
    {
      BKGandSignal->Add(TprimeHistos[i]);
      if (i!=NOS-1)
	{
	  BKGLJPT->Add(LeadingJetPT[i]);
	  BKGTHT->Add(THT[i]);
	}
    }

  //PostScript Plotting
  char FileName[100];
  sprintf(FileName,"Single_Tprime.ps");
  TPostScript *ps = new TPostScript(FileName,111);
  TCanvas *MyPlot = new TCanvas("MyPlot","Single t prime to top Higgs with backgrounds",600,800);
  TPaveText *p1;
  p1 = new TPaveText(3,0.5,8,3.5);
  char dateandtime[50];
  sprintf(dateandtime,"date: %s, time: %s",__DATE__,__TIME__);
  p1->AddText(dateandtime);
  //First Page
  MyPlot->Clear();
  MyPlot->cd(1);
  ps->NewPage();
  BKGandSignal->Draw("hist");
  TprimeHistos[4]->Draw("histsame");
  gPad->SetLogy();
  gPad->RedrawAxis("");
  gPad->Update();
  MyPlot->Update();
  //Second Page
  MyPlot->Clear();
  MyPlot->cd(1);
  ps->NewPage();
  BKGLJPT->Draw("hist");
  LeadingJetPT[4]->Draw("histsame");
  gPad->SetLogy();
  gPad->RedrawAxis("");
  gPad->Update();
  MyPlot->Update();
  //Third Page
  MyPlot->Clear();
  MyPlot->cd(1);
  ps->NewPage();
  BKGTHT->Draw("hist");
  THT[4]->Draw("histsame");
  gPad->SetLogy();
  gPad->RedrawAxis("");
  gPad->Update();
  MyPlot->Update();

}
