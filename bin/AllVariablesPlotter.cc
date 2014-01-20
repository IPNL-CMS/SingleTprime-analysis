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
  TH1F *RelativeMass[NOS];
  TH1F *MotherPtNormalizedMass[NOS];
  
  for (int i=0; i<NOS; i++)
    {
      CurrentFile[i] = new TFile(StorageDirPrefix + Samples[i], "READ");
      if ( CurrentFile[i]->IsOpen() ) printf( Samples[i] + " File opened successfully\n");
      if (i==0) 
	{
	  TprimeHistos[0]->SetDefaultSumw2(); TprimeHistos[0]= (TH1F*)gDirectory->Get("TprimeMass0");
	  LeadingJetPT[0]->SetDefaultSumw2(); LeadingJetPT[0]= (TH1F*)gDirectory->Get("jet1_pt0");
	  THT[0]->SetDefaultSumw2(); THT[0]= (TH1F*)gDirectory->Get("THT0");
	  DRHjets[0]->SetDefaultSumw2(); DRHjets[0]= (TH1F*)gDirectory->Get("DeltaR_of_Higgs_Jets0");
	  DRWjets[0]->SetDefaultSumw2(); DRWjets[0]= (TH1F*)gDirectory->Get("DeltaR_of_W_Jets0");
	  Hpt[0]->SetDefaultSumw2(); Hpt[0]= (TH1F*)gDirectory->Get("HPt0");
	  Tpt[0]->SetDefaultSumw2(); Tpt[0]= (TH1F*)gDirectory->Get("TPt0");
	  DRWH[0]->SetDefaultSumw2(); DRWH[0]= (TH1F*)gDirectory->Get("DeltaR_of_W_Higgs0");
	  DPHjets[0]->SetDefaultSumw2(); DPHjets[0]= (TH1F*)gDirectory->Get("DeltaPhi_of_Higgs_jets0");
	  DPWjets[0]->SetDefaultSumw2(); DPWjets[0]= (TH1F*)gDirectory->Get("DeltaPhi_of_W_jets0");
	  DPTjets[0]->SetDefaultSumw2(); DPTjets[0]= (TH1F*)gDirectory->Get("DeltaPhi_of_T_jet0");
	  HiggsMass[0]->SetDefaultSumw2(); HiggsMass[0]= (TH1F*)gDirectory->Get("HM0");
	  RelHT[0]->SetDefaultSumw2(); RelHT[0]= (TH1F*)gDirectory->Get("RelHT0");
	  DRTH[0]->SetDefaultSumw2(); DRTH[0]= (TH1F*)gDirectory->Get("DeltaR_of_Top_Higgs0");
	  PtNormalizedMass[0]->SetDefaultSumw2(); PtNormalizedMass[0]= (TH1F*)gDirectory->Get("PT_Normalized_Mass0");
	  RelativeMass[0]->SetDefaultSumw2(); PtNormalizedMass[0]= (TH1F*)gDirectory->Get("Relative_Mass0");
	  MotherPtNormalizedMass[0]->SetDefaultSumw2(); MotherPtNormalizedMass[0]= (TH1F*)gDirectory->Get("Mother_PT_Normalized_Mass0");
	}
      else if (i==1) 
	{
	  TprimeHistos[1]->SetDefaultSumw2(); TprimeHistos[1]= (TH1F*)gDirectory->Get("TprimeMass3");
	  LeadingJetPT[1]->SetDefaultSumw2(); LeadingJetPT[1]= (TH1F*)gDirectory->Get("jet1_pt3");
	  THT[1]->SetDefaultSumw2(); THT[1]= (TH1F*)gDirectory->Get("THT3");
	  DRHjets[1]->SetDefaultSumw2(); DRHjets[1]= (TH1F*)gDirectory->Get("DeltaR_of_Higgs_Jets3");
	  DRWjets[1]->SetDefaultSumw2(); DRWjets[1]= (TH1F*)gDirectory->Get("DeltaR_of_W_Jets3");
	  Hpt[1]->SetDefaultSumw2(); Hpt[1]= (TH1F*)gDirectory->Get("HPt3");
	  Tpt[1]->SetDefaultSumw2(); Tpt[1]= (TH1F*)gDirectory->Get("TPt3");
	  DRWH[1]->SetDefaultSumw2(); DRWH[1]= (TH1F*)gDirectory->Get("DeltaR_of_W_Higgs3");
	  DPHjets[1]->SetDefaultSumw2(); DPHjets[1]= (TH1F*)gDirectory->Get("DeltaPhi_of_Higgs_jets3");
	  DPWjets[1]->SetDefaultSumw2(); DPWjets[1]= (TH1F*)gDirectory->Get("DeltaPhi_of_W_jets3");
	  DPTjets[1]->SetDefaultSumw2(); DPTjets[1]= (TH1F*)gDirectory->Get("DeltaPhi_of_T_jet3");
	  HiggsMass[1]->SetDefaultSumw2(); HiggsMass[1]= (TH1F*)gDirectory->Get("HM3");
	  RelHT[1]->SetDefaultSumw2(); RelHT[1]= (TH1F*)gDirectory->Get("RelHT3");
	  DRTH[1]->SetDefaultSumw2(); DRTH[1]= (TH1F*)gDirectory->Get("DeltaR_of_Top_Higgs3");
	  PtNormalizedMass[1]->SetDefaultSumw2(); PtNormalizedMass[1]= (TH1F*)gDirectory->Get("PT_Normalized_Mass3");
	  RelativeMass[1]->SetDefaultSumw2(); PtNormalizedMass[1]= (TH1F*)gDirectory->Get("Relative_Mass3");
	  MotherPtNormalizedMass[1]->SetDefaultSumw2(); MotherPtNormalizedMass[1]= (TH1F*)gDirectory->Get("Mother_PT_Normalized_Mass3");
	}
      else if (i==2) 
	{
	  TprimeHistos[2]->SetDefaultSumw2(); TprimeHistos[2]= (TH1F*)gDirectory->Get("TprimeMass9");
	  LeadingJetPT[2]->SetDefaultSumw2(); LeadingJetPT[2]= (TH1F*)gDirectory->Get("jet1_pt9");
	  THT[2]->SetDefaultSumw2(); THT[2]= (TH1F*)gDirectory->Get("THT9");
	  DRHjets[2]->SetDefaultSumw2(); DRHjets[2]= (TH1F*)gDirectory->Get("DeltaR_of_Higgs_Jets9");
	  DRWjets[2]->SetDefaultSumw2(); DRWjets[2]= (TH1F*)gDirectory->Get("DeltaR_of_W_Jets9");
	  Hpt[2]->SetDefaultSumw2(); Hpt[2]= (TH1F*)gDirectory->Get("HPt9");
	  Tpt[2]->SetDefaultSumw2(); Tpt[2]= (TH1F*)gDirectory->Get("TPt9");
	  DRWH[2]->SetDefaultSumw2(); DRWH[2]= (TH1F*)gDirectory->Get("DeltaR_of_W_Higgs9");
	  DPHjets[2]->SetDefaultSumw2(); DPHjets[2]= (TH1F*)gDirectory->Get("DeltaPhi_of_Higgs_jets9");
	  DPWjets[2]->SetDefaultSumw2(); DPWjets[2]= (TH1F*)gDirectory->Get("DeltaPhi_of_W_jets9");
	  DPTjets[2]->SetDefaultSumw2(); DPTjets[2]= (TH1F*)gDirectory->Get("DeltaPhi_of_T_jet9");
	  HiggsMass[2]->SetDefaultSumw2(); HiggsMass[2]= (TH1F*)gDirectory->Get("HM9");
	  RelHT[2]->SetDefaultSumw2(); RelHT[2]= (TH1F*)gDirectory->Get("RelHT9");
	  DRTH[2]->SetDefaultSumw2(); DRTH[2]= (TH1F*)gDirectory->Get("DeltaR_of_Top_Higgs9");
	  PtNormalizedMass[2]->SetDefaultSumw2(); PtNormalizedMass[2]= (TH1F*)gDirectory->Get("PT_Normalized_Mass9");
	  RelativeMass[2]->SetDefaultSumw2(); PtNormalizedMass[2]= (TH1F*)gDirectory->Get("Relative_Mass9");
	  MotherPtNormalizedMass[2]->SetDefaultSumw2(); MotherPtNormalizedMass[2]= (TH1F*)gDirectory->Get("Mother_PT_Normalized_Mass9");
	}
      else if (i==3) 
	{
	  TprimeHistos[3]->SetDefaultSumw2(); TprimeHistos[3]= (TH1F*)gDirectory->Get("TprimeMass10");
	  LeadingJetPT[3]->SetDefaultSumw2(); LeadingJetPT[3]= (TH1F*)gDirectory->Get("jet1_pt10");
	  THT[3]->SetDefaultSumw2(); THT[3]= (TH1F*)gDirectory->Get("THT10");
	  DRHjets[3]->SetDefaultSumw2(); DRHjets[3]= (TH1F*)gDirectory->Get("DeltaR_of_Higgs_Jets10");
	  DRWjets[3]->SetDefaultSumw2(); DRWjets[3]= (TH1F*)gDirectory->Get("DeltaR_of_W_Jets10");
	  Hpt[3]->SetDefaultSumw2(); Hpt[3]= (TH1F*)gDirectory->Get("HPt10");
	  Tpt[3]->SetDefaultSumw2(); Tpt[3]= (TH1F*)gDirectory->Get("TPt10");
	  DRWH[3]->SetDefaultSumw2(); DRWH[3]= (TH1F*)gDirectory->Get("DeltaR_of_W_Higgs10");
	  DPHjets[3]->SetDefaultSumw2(); DPHjets[3]= (TH1F*)gDirectory->Get("DeltaPhi_of_Higgs_jets10");
	  DPWjets[3]->SetDefaultSumw2(); DPWjets[3]= (TH1F*)gDirectory->Get("DeltaPhi_of_W_jets10");
	  DPTjets[3]->SetDefaultSumw2(); DPTjets[3]= (TH1F*)gDirectory->Get("DeltaPhi_of_T_jet10");
	  HiggsMass[3]->SetDefaultSumw2(); HiggsMass[3]= (TH1F*)gDirectory->Get("HM10");
	  RelHT[3]->SetDefaultSumw2(); RelHT[3]= (TH1F*)gDirectory->Get("RelHT10");
	  DRTH[3]->SetDefaultSumw2(); DRTH[3]= (TH1F*)gDirectory->Get("DeltaR_of_Top_Higgs10");
	  PtNormalizedMass[3]->SetDefaultSumw2(); PtNormalizedMass[3]= (TH1F*)gDirectory->Get("PT_Normalized_Mass10");
	  RelativeMass[3]->SetDefaultSumw2(); PtNormalizedMass[3]= (TH1F*)gDirectory->Get("Relative_Mass10");
	  MotherPtNormalizedMass[3]->SetDefaultSumw2(); MotherPtNormalizedMass[3]= (TH1F*)gDirectory->Get("Mother_PT_Normalized_Mass10");
	}
      else if (i==4) 
	{
	  TprimeHistos[4]->SetDefaultSumw2(); TprimeHistos[4]= (TH1F*)gDirectory->Get("TprimeMass15");
	  LeadingJetPT[4]->SetDefaultSumw2(); LeadingJetPT[4]= (TH1F*)gDirectory->Get("jet1_pt15");
	  THT[4]->SetDefaultSumw2(); THT[4]= (TH1F*)gDirectory->Get("THT15");
	  DRHjets[4]->SetDefaultSumw2(); DRHjets[4]= (TH1F*)gDirectory->Get("DeltaR_of_Higgs_Jets15");
	  DRWjets[4]->SetDefaultSumw2(); DRWjets[4]= (TH1F*)gDirectory->Get("DeltaR_of_W_Jets15");
	  Hpt[4]->SetDefaultSumw2(); Hpt[4]= (TH1F*)gDirectory->Get("HPt15");
	  Tpt[4]->SetDefaultSumw2(); Tpt[4]= (TH1F*)gDirectory->Get("TPt15");
	  DRWH[4]->SetDefaultSumw2(); DRWH[4]= (TH1F*)gDirectory->Get("DeltaR_of_W_Higgs15");
	  DPHjets[4]->SetDefaultSumw2(); DPHjets[4]= (TH1F*)gDirectory->Get("DeltaPhi_of_Higgs_jets15");
	  DPWjets[4]->SetDefaultSumw2(); DPWjets[4]= (TH1F*)gDirectory->Get("DeltaPhi_of_W_jets15");
	  DPTjets[4]->SetDefaultSumw2(); DPTjets[4]= (TH1F*)gDirectory->Get("DeltaPhi_of_T_jet15");
	  HiggsMass[4]->SetDefaultSumw2(); HiggsMass[4]= (TH1F*)gDirectory->Get("HM15");
	  RelHT[4]->SetDefaultSumw2(); RelHT[4]= (TH1F*)gDirectory->Get("RelHT15");
	  DRTH[4]->SetDefaultSumw2(); DRTH[4]= (TH1F*)gDirectory->Get("DeltaR_of_Top_Higgs15");
	  PtNormalizedMass[4]->SetDefaultSumw2(); PtNormalizedMass[4]= (TH1F*)gDirectory->Get("PT_Normalized_Mass15");
	  RelativeMass[4]->SetDefaultSumw2(); PtNormalizedMass[4]= (TH1F*)gDirectory->Get("Relative_Mass15");
	  MotherPtNormalizedMass[4]->SetDefaultSumw2(); MotherPtNormalizedMass[4]= (TH1F*)gDirectory->Get("Mother_PT_Normalized_Mass15");
	}
    }
  
  THStack *BKGandSignal = new THStack("BKGandSignal", "BKG and signal for 5 jets; M_{5j} GeV; Events");
  THStack *BKGLJPT = new THStack("BKGLJPT", "BKG for Leading Jet PT; pT(j_{1}) GeV; Events");
  THStack *BKGTHT = new THStack("BKGTHT", "BKG for HT; HT GeV; Events");
  THStack *BKGDRHjets= new THStack("BKGDRHjets", "BKG for DRHjets; #Delta R((jj)^{H}); Events");
  THStack *BKGDRWjets= new THStack("BKGDRWjets", "BKG for DRWjets; #Delta R((jj)^{W}); Events");
  THStack *BKGHpt= new THStack("BKGHpt", "BKG for Hpt; p_{T}(H) GeV; Events");
  THStack *BKGTpt= new THStack("BKGTpt", "BKG for Tpt; p_{T}(t) GeV; Events");
  THStack *BKGDRWH= new THStack("BKGDRWH", "BKG for DRWH; #Delta R(WH); Events");
  THStack *BKGDPHjets= new THStack("BKGDPHjets", "BKG for DPHjets; #Delta #phi((jj)^{H}); Events");
  THStack *BKGDPWjets= new THStack("BKGDPWjets", "BKG for DPWjets; #Delta #phi((jj)^{W}); Events");
  THStack *BKGDPTjets= new THStack("BKGDPTjets", "BKG for DPTjets; #Delta #phi(Wj^{t}); Events");
  THStack *BKGHiggsMass= new THStack("BKGHiggsMass", "BKG for HiggsMass; M_{H} GeV; Events");
  THStack *BKGRelHT= new THStack("BKGRelHT", "BKG for RelHT; (p_{T}(H)+p_{T}(t))/H_{T}; Events");
  THStack *BKGDRTH= new THStack("BKGDRTH", "BKG for DRTH; #Delta R(TH); Events");
  THStack *BKGPtNormalizedMass= new THStack("BKGPtNormalizedMass", "BKG for PtNormalizedMass; PTNM; Events");
  THStack *BKGRelativeMass= new THStack("BKGRelativeMass", "BKG for RelativeMass; RM; Events");
  THStack *BKGMotherPtNormalizedMass= new THStack("BKGMotherPtNormalizedMass", "BKG for MotherPtNormalizedMass; PTNM; Events");
  for (int i=0; i<NOS; i++)
    {
      BKGandSignal->Add(TprimeHistos[i]);
      if (i!=NOS-1)
	{
	  BKGLJPT->Add(LeadingJetPT[i]);
	  BKGTHT->Add(THT[i]);
	  BKGDRHjets->Add(DRHjets[i]);
	  BKGDRWjets->Add(DRWjets[i]);
	  BKGHpt->Add(Hpt[i]);
	  BKGTpt->Add(Tpt[i]);
	  BKGDRWH->Add(DRWH[i]);
	  BKGDPHjets->Add(DPHjets[i]);
	  BKGDPWjets->Add(DPWjets[i]);
	  BKGDPTjets->Add(DPTjets[i]);
	  BKGHiggsMass->Add(HiggsMass[i]);
	  BKGRelHT->Add(RelHT[i]);
	  BKGDRTH->Add(DRTH[i]);
	  BKGPtNormalizedMass->Add(PtNormalizedMass[i]);
	  BKGRelativeMass->Add(RelativeMass[i]);
	  BKGMotherPtNormalizedMass->Add(MotherPtNormalizedMass[i]);
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
  //Fourth Page
  MyPlot->Clear();
  MyPlot->cd(1);
  ps->NewPage();
  BKGDRHjets->Draw("hist");
  DRHjets[4]->Draw("histsame");
  gPad->SetLogy();
  gPad->RedrawAxis("");
  gPad->Update();
  MyPlot->Update();
  //5th Page
  MyPlot->Clear();
  MyPlot->cd(1);
  ps->NewPage();
  BKGDRWjets->Draw("hist");
  DRWjets[4]->Draw("histsame");
  gPad->SetLogy();
  gPad->RedrawAxis("");
  gPad->Update();
  MyPlot->Update();
  //6th Page
  MyPlot->Clear();
  MyPlot->cd(1);
  ps->NewPage();
  BKGHpt->Draw("hist");
  Hpt[4]->Draw("histsame");
  gPad->SetLogy();
  gPad->RedrawAxis("");
  gPad->Update();
  MyPlot->Update();
  //7th Page
  MyPlot->Clear();
  MyPlot->cd(1);
  ps->NewPage();
  BKGTpt->Draw("hist");
  Tpt[4]->Draw("histsame");
  gPad->SetLogy();
  gPad->RedrawAxis("");
  gPad->Update();
  MyPlot->Update();
  //8th Page
  MyPlot->Clear();
  MyPlot->cd(1);
  ps->NewPage();
  BKGDRWH->Draw("hist");
  DRWH[4]->Draw("histsame");
  gPad->SetLogy();
  gPad->RedrawAxis("");
  gPad->Update();
  MyPlot->Update();
  //9th Page
  MyPlot->Clear();
  MyPlot->cd(1);
  ps->NewPage();
  BKGDPHjets->Draw("hist");
  DPHjets[4]->Draw("histsame");
  gPad->SetLogy();
  gPad->RedrawAxis("");
  gPad->Update();
  MyPlot->Update();
  //10th Page
  MyPlot->Clear();
  MyPlot->cd(1);
  ps->NewPage();
  BKGDPWjets->Draw("hist");
  DPWjets[4]->Draw("histsame");
  gPad->SetLogy();
  gPad->RedrawAxis("");
  gPad->Update();
  MyPlot->Update();
  //11th Page
  MyPlot->Clear();
  MyPlot->cd(1);
  ps->NewPage();
  BKGDPTjets->Draw("hist");
  DPTjets[4]->Draw("histsame");
  gPad->SetLogy();
  gPad->RedrawAxis("");
  gPad->Update();
  MyPlot->Update();
  //12th Page
  MyPlot->Clear();
  MyPlot->cd(1);
  ps->NewPage();
  BKGHiggsMass->Draw("hist");
  HiggsMass[4]->Draw("histsame");
  gPad->SetLogy();
  gPad->RedrawAxis("");
  gPad->Update();
  MyPlot->Update();
  //13th Page
  MyPlot->Clear();
  MyPlot->cd(1);
  ps->NewPage();
  BKGRelHT->Draw("hist");
  RelHT[4]->Draw("histsame");
  gPad->SetLogy();
  gPad->RedrawAxis("");
  gPad->Update();
  MyPlot->Update();
  //14th Page
  MyPlot->Clear();
  MyPlot->cd(1);
  ps->NewPage();
  BKGDRTH->Draw("hist");
  DRTH[4]->Draw("histsame");
  gPad->SetLogy();
  gPad->RedrawAxis("");
  gPad->Update();
  MyPlot->Update();
  //15th Page
  MyPlot->Clear();
  MyPlot->cd(1);
  ps->NewPage();
  BKGPtNormalizedMass->Draw("hist");
  PtNormalizedMass[4]->Draw("histsame");
  gPad->SetLogy();
  gPad->RedrawAxis("");
  gPad->Update();
  MyPlot->Update();
  //16th Page
  MyPlot->Clear();
  MyPlot->cd(1);
  ps->NewPage();
  BKGRelativeMass->Draw("hist");
  RelativeMass[4]->Draw("histsame");
  gPad->SetLogy();
  gPad->RedrawAxis("");
  gPad->Update();
  MyPlot->Update();
  //17th Page
  MyPlot->Clear();
  MyPlot->cd(1);
  ps->NewPage();
  BKGMotherPtNormalizedMass->Draw("hist");
  MotherPtNormalizedMass[4]->Draw("histsame");
  gPad->SetLogy();
  gPad->RedrawAxis("");
  gPad->Update();
  MyPlot->Update();
}
