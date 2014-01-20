#include <TFile.h>
#include <TTree.h>
#include <TNtuple.h>
#include <TString.h>
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

#include "TMath.h"
#include "TF1.h"
#include "FilesNamesNoCuts.C"

using namespace std;

// Global Parameters
// Location of root files

const TString MainFolder = "file:/gridgroup/cms/jruizalv/Extracted_MC/NoCuts-Partial/";

const TString StorageDirPrefix = "file:/home/cms/jruizalv/work/CMSSW_5_3_9_patch2/src/Extractors/PatExtractor/test/";

const int NumberOfHistos=20;
const TString Histos[NumberOfHistos] = {"Cut_0", "Cut_1", "Cut_2", "Cut_3", "Cut_4", "Cut_5", "Cut_6", "Cut_7", "Cut_8", "Cut_9", "Cut_10", "Cut_11", "Cut_12", "Cut_13", "Cut_14", "Cut_15", "Cut_16", "Cut_17", "Cut_18", "Cut_19"};

//Computing weights (Everything in pb)
const float XS = 0.15; 
/*float T_s_xs=2.82;
float T_t_xs=47;
float T_tW_xs=10.7;
float Tbar_s_xs=1.57;
float Tbar_t_xs=25;
float Tbar_tW_xs=10.7;
float QCD_120_xs=156293.3;
float QCD_170_xs=34138.15;
float QCD_300_xs=1759.549;
float QCD_470_xs=113.8791;
float QCD_600_xs=26.9921;
float QCD_800_xs=3.550036;
float WW_xs=54.838;
float WZ_xs=33.72;
float ZZ_xs=8.258;
float TTbar_xs=234.0;
float Signal_xs=0.15;
float DYToBB_xs=3840.86;
float DYToCC_xs=3060.099;
float Wjets_VBF_xs=1205.0;*/

float Lumi=20000;

void Signal()
{
  TH1F *FiveJetsMass;
  TH1F *LeadingJetPT;
  TH1F *THT;
  TH1F *DRHjets;
  TH1F *DRWjets;
  TH1F *Hpt;
  TH1F *Tpt;
  TH1F *DRWH;
  TH1F *DPHjets;
  TH1F *DPWjets;
  TH1F *DPTjets;
  TH1F *HiggsMass;
  TH1F *RelHT;
  TH1F *DRTH;
  TH1F *PtNormalizedMass;
  TH1F *RelativeMass;
  TH1F *MotherPtNormalizedMass;
  int EntriePerSample;
  bool SurvivalMarker;

  TChain CutsChain("cuts");
  TChain AnalysisChain("stp");
  TH1F *ALLCuts[NumberOfHistos];
  CutsChain.Add(StorageDirPrefix + SamplesSignal[0]);
  AnalysisChain.Add(StorageDirPrefix + SamplesSignal[0]);
  
  int EntriesPerCut[NumberOfHistos];
  float PassedPerCut[NumberOfHistos];
  
  for (int j=0; j<NumberOfHistos; j++)
    {
      CutsChain.Draw(Histos[j] + " >> " + Histos[j]);
      ALLCuts[j] =  (TH1F*)gDirectory->Get(Histos[j]);
      EntriesPerCut[j]=ALLCuts[j]->GetEntries();
      PassedPerCut[j]=ALLCuts[j]->GetBinContent(ALLCuts[j]->GetXaxis()->FindBin(1));
      if (j==0) cout << "Cut " << j << ": Total " << EntriesPerCut[j] << " Passed " << PassedPerCut[j] << " Efficiency " << PassedPerCut[j]/EntriesPerCut[j] << endl;
      else cout << "Cut " << j << ": Total " << EntriesPerCut[j] << " Passed " << PassedPerCut[j] << " Efficiency " << PassedPerCut[j]/PassedPerCut[j-1] << endl;
    }
  
  EntriePerSample=EntriesPerCut[0];
  cout << PassedPerCut[NumberOfHistos-1] << endl;
  gPad->Close();
  if (PassedPerCut[NumberOfHistos-1]!=0)
    {
      string A1 = Form("Reconstructed_Tprime.M() >> TprimeMass%i(60,400,1600)",15);
      string A2 = Form("TprimeMass%i",15);
      AnalysisChain.Draw(A1.c_str());
      FiveJetsMass = (TH1F*)gDirectory->Get(A2.c_str());
      gPad->Close();
      SurvivalMarker=true;
      /////////////////////
      //Saving Histograms//
      /////////////////////
      string LJPT1 = Form("jet1_pt >> jet1_pt%i(100,10,500)",15);
      string LJPT2 = Form("jet1_pt%i",15);
      AnalysisChain.Draw(LJPT1.c_str());
      LeadingJetPT = (TH1F*)gDirectory->Get(LJPT2.c_str());
      gPad->Close();
      //
      string THT1 = Form("THT >> THT%i(100,600,1400)",15);
      string THT2 = Form("THT%i",15);
      AnalysisChain.Draw(THT1.c_str());
      THT = (TH1F*)gDirectory->Get(THT2.c_str());
      gPad->Close(); 
      //
      string DRHJ1 = Form("DeltaR_of_Higgs_Jets >> DeltaR_of_Higgs_Jets%i(40,0.4,7)",15);
      string DRHJ2 = Form("DeltaR_of_Higgs_Jets%i",15);
      AnalysisChain.Draw(DRHJ1.c_str());
      DRHjets = (TH1F*)gDirectory->Get(DRHJ2.c_str());
      gPad->Close();
      //
      string DRWJ1 = Form("DeltaR_of_W_Jets >> DeltaR_of_W_Jets%i(40,0.4,7)",15);
      string DRWJ2 = Form("DeltaR_of_W_Jets%i",15);
      AnalysisChain.Draw(DRWJ1.c_str(), "DeltaR_of_W_Jets<3");
      DRWjets = (TH1F*)gDirectory->Get(DRWJ2.c_str());
      gPad->Close();
      //
      string HPT1 = Form("Reconstructed_Higgs.Pt() >> HPt%i(200,10,800)",15);
      string HPT2 = Form("HPt%i",15);
      AnalysisChain.Draw(HPT1.c_str());
      Hpt = (TH1F*)gDirectory->Get(HPT2.c_str());
      gPad->Close();
      //
      string TPT1 = Form("Reconstructed_Top.Pt() >> TPt%i(200,10,800)",15);
      string TPT2 = Form("TPt%i",15);
      AnalysisChain.Draw(TPT1.c_str());
      Tpt = (TH1F*)gDirectory->Get(TPT2.c_str());
      gPad->Close();
      //
      string DRWH1 = Form("DeltaR_of_W_Higgs >> DeltaR_of_W_Higgs%i(40,0.4,7)",15);
      string DRWH2 = Form("DeltaR_of_W_Higgs%i",15);
      AnalysisChain.Draw(DRWH1.c_str());
      DRWH = (TH1F*)gDirectory->Get(DRWH2.c_str());
      gPad->Close();
      //
      string DPHJ1 = Form("TMath::Abs(First_Higgs_Jet.Phi()-Second_Higgs_Jet.Phi()) >> DeltaPhi_of_Higgs_jets%i(40,0.0,3.5)",15);
      string DPHJ2 = Form("DeltaPhi_of_Higgs_jets%i",15);
      AnalysisChain.Draw(DPHJ1.c_str());
      DPHjets = (TH1F*)gDirectory->Get(DPHJ2.c_str());
      gPad->Close();
      //
      string DPWJ1 = Form("TMath::Abs(First_W_Jet.Phi()-Second_W_Jet.Phi()) >> DeltaPhi_of_W_jets%i(40,0.0,3.5)",15);
      string DPWJ2 = Form("DeltaPhi_of_W_jets%i",15);
      AnalysisChain.Draw(DPWJ1.c_str());
      DPWjets = (TH1F*)gDirectory->Get(DPWJ2.c_str());
      gPad->Close();
      //
      string DPTJ1 = Form("TMath::Abs(Top_Jet.Phi()-Reconstructed_W.Phi()) >> DeltaPhi_of_T_jet%i(40,0.0,3.5)",15);
      string DPTJ2 = Form("DeltaPhi_of_T_jet%i",15);
      AnalysisChain.Draw(DPTJ1.c_str());
      DPTjets = (TH1F*)gDirectory->Get(DPTJ2.c_str());
      gPad->Close();
      //
      string HM1 = Form("Reconstructed_Higgs.M() >> HM%i(30,60,180)",15);
      string HM2 = Form("HM%i",15);
      AnalysisChain.Draw(HM1.c_str());
      HiggsMass = (TH1F*)gDirectory->Get(HM2.c_str());
      gPad->Close();
      //
      string RHT1 = Form("Relative_THT >> RelHT%i(20,0,1)",15);
      string RHT2 = Form("RelHT%i",15);
      AnalysisChain.Draw(RHT1.c_str());
      RelHT = (TH1F*)gDirectory->Get(RHT2.c_str());
      gPad->Close();
      //
      string DRTH1 = Form("DeltaR_of_Top_Higgs >> DeltaR_of_Top_Higgs%i(40,0.4,7)",15);
      string DRTH2 = Form("DeltaR_of_Top_Higgs%i",15);
      AnalysisChain.Draw(DRTH1.c_str());
      DRTH = (TH1F*)gDirectory->Get(DRTH2.c_str());
      gPad->Close();
      //
      string PTNM1 = Form("PT_Normalized_Mass >> PT_Normalized_Mass%i(40,0.4,7)",15);
      string PTNM2 = Form("PT_Normalized_Mass%i",15);
      AnalysisChain.Draw(PTNM1.c_str());
      PtNormalizedMass = (TH1F*)gDirectory->Get(PTNM2.c_str());
      gPad->Close();
      //
      string RM1 = Form("Relative_Mass >> Relative_Mass%i(20,0.0,1)",15);
      string RM2 = Form("Relative_Mass%i",15);
      AnalysisChain.Draw(RM1.c_str());
      RelativeMass = (TH1F*)gDirectory->Get(RM2.c_str());
      gPad->Close();
      //
      string MPTNM1 = Form("Mother_PT_Normalized_Mass >> Mother_PT_Normalized_Mass%i(40,0.4,7)",15);
      string MPTNM2 = Form("Mother_PT_Normalized_Mass%i",15);
      AnalysisChain.Draw(MPTNM1.c_str());
      MotherPtNormalizedMass = (TH1F*)gDirectory->Get(MPTNM2.c_str());
      gPad->Close();
    }
  
  //Making Stack
  
  THStack *BKGandSignal = new THStack("BKGandSignal", "BKG and signal for 5 jets; M_{5j} GeV; Events");
  TLegend* BKGandSignallegend = new TLegend(0.75,0.55,0.90,0.9);
  
  FiveJetsMass->Scale(Lumi*XS/EntriePerSample);
  LeadingJetPT->Scale(Lumi*XS/EntriePerSample);
  THT->Scale(Lumi*XS/EntriePerSample);
  DRHjets->Scale(Lumi*XS/EntriePerSample);
  DRWjets->Scale(Lumi*XS/EntriePerSample);
  Hpt->Scale(Lumi*XS/EntriePerSample);
  Tpt->Scale(Lumi*XS/EntriePerSample);
  DRWH->Scale(Lumi*XS/EntriePerSample);
  DPHjets->Scale(Lumi*XS/EntriePerSample);
  DPWjets->Scale(Lumi*XS/EntriePerSample);
  DPTjets->Scale(Lumi*XS/EntriePerSample);
  HiggsMass->Scale(Lumi*XS/EntriePerSample);
  RelHT->Scale(Lumi*XS/EntriePerSample);
  DRTH->Scale(Lumi*XS/EntriePerSample);
  PtNormalizedMass->Scale(Lumi*XS/EntriePerSample);
  RelativeMass->Scale(Lumi*XS/EntriePerSample);
  MotherPtNormalizedMass->Scale(Lumi*XS/EntriePerSample);
  
  //Settings for signal
  if (SurvivalMarker)
    {
      //TprimeMass
      FiveJetsMass->SetFillColor(kSpring);
      FiveJetsMass->SetFillStyle(3444);
      FiveJetsMass->SetLineWidth(3);
      TFile f("Signal.root", "RECREATE");
      FiveJetsMass->Write();
      //LJPT
      LeadingJetPT->SetFillColor(kSpring);
      LeadingJetPT->SetFillStyle(3444);
      LeadingJetPT->SetLineWidth(3);
      LeadingJetPT->Write();
      //THT
      THT->SetFillColor(kSpring);
      THT->SetFillStyle(3444);
      THT->SetLineWidth(3);
      THT->Write();
      //
      DRHjets->SetFillColor(kSpring);
      DRHjets->SetFillStyle(3444);
      DRHjets->SetLineWidth(3);
      DRHjets->Write();
      //
      DRWjets->SetFillColor(kSpring);
      DRWjets->SetFillStyle(3444);
      DRWjets->SetLineWidth(3);
      DRWjets->Write();
      //
      Hpt->SetFillColor(kSpring);
      Hpt->SetFillStyle(3444);
      Hpt->SetLineWidth(3);
      Hpt->Write();
      //
      Tpt->SetFillColor(kSpring);
      Tpt->SetFillStyle(3444);
      Tpt->SetLineWidth(3);
      Tpt->Write();
      //
      DRWH->SetFillColor(kSpring);
      DRWH->SetFillStyle(3444);
      DRWH->SetLineWidth(3);
      DRWH->Write();
      //
      DPHjets->SetFillColor(kSpring);
      DPHjets->SetFillStyle(3444);
      DPHjets->SetLineWidth(3);
      DPHjets->Write();
      //
      DPWjets->SetFillColor(kSpring);
      DPWjets->SetFillStyle(3444);
      DPWjets->SetLineWidth(3);
      DPWjets->Write();
      //
      DPTjets->SetFillColor(kSpring);
      DPTjets->SetFillStyle(3444);
      DPTjets->SetLineWidth(3);
      DPTjets->Write();
      //
      HiggsMass->SetFillColor(kSpring);
      HiggsMass->SetFillStyle(3444);
      HiggsMass->SetLineWidth(3);
      HiggsMass->Write();
      //
      RelHT->SetFillColor(kSpring);
      RelHT->SetFillStyle(3444);
      RelHT->SetLineWidth(3);
      RelHT->Write();
      //
      DRTH->SetFillColor(kSpring);
      DRTH->SetFillStyle(3444);
      DRTH->SetLineWidth(3);
      DRTH->Write();
      //
      PtNormalizedMass->SetFillColor(kSpring);
      PtNormalizedMass->SetFillStyle(3444);
      PtNormalizedMass->SetLineWidth(3);
      PtNormalizedMass->Write();
      //
      RelativeMass->SetFillColor(kSpring);
      RelativeMass->SetFillStyle(3444);
      RelativeMass->SetLineWidth(3);
      RelativeMass->Write();
      //
      MotherPtNormalizedMass->SetFillColor(kSpring);
      MotherPtNormalizedMass->SetFillStyle(3444);
      MotherPtNormalizedMass->SetLineWidth(3);
      MotherPtNormalizedMass->Write();
    }
  BKGandSignal->Add(FiveJetsMass);
  BKGandSignallegend->AddEntry(FiveJetsMass);
}

