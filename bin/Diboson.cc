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

const TString StorageDirPrefix[NumberOfProcesses]={MainFolder + "ZZ/", MainFolder + "WZ/", MainFolder + "WW/", MainFolder + "T-s/", MainFolder + "T-t/",MainFolder + "T-tw/", MainFolder + "Tbar-s/", MainFolder + "Tbar-t/", MainFolder + "Tbar-tw/", MainFolder + "TTJets/", MainFolder + "QCD_PT_170_300/", MainFolder + "QCD_PT_300_470/", MainFolder + "QCD_PT_470_600/", MainFolder + "QCD_PT_600_800/", MainFolder + "QCD_PT_800_1000/", "file:/home/cms/jruizalv/work/CMSSW_5_3_9_patch2/src/Extractors/PatExtractor/test/"};

const int NumberOfHistos=20;
const TString Histos[NumberOfHistos] = {"Cut_0", "Cut_1", "Cut_2", "Cut_3", "Cut_4", "Cut_5", "Cut_6", "Cut_7", "Cut_8", "Cut_9", "Cut_10", "Cut_11", "Cut_12", "Cut_13", "Cut_14", "Cut_15", "Cut_16", "Cut_17", "Cut_18", "Cut_19"};

//Computing weights (Everything in pb)
const float XS[NumberOfProcesses]={8.258, 33.72, 54.838, 2.82, 47, 10.7, 1.57, 25, 10.7, 234.0, 34138.15, 1759.549, 113.8791, 26.9921, 3.550036, 0.15}; 
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

void Diboson()
{
  int ParcialTestMax=3; //Don't change
  int ParcialTestMin=1; //Don't change (When ZZ is available put 0 instead of 1)
  TH1F *FiveJetsMass[NumberOfProcesses];
  TH1F *LeadingJetPT[NumberOfProcesses];
  TH1F *THT[NumberOfProcesses];
  TH1F *DRHjets[NumberOfProcesses];
  TH1F *DRWjets[NumberOfProcesses];
  TH1F *Hpt[NumberOfProcesses];
  TH1F *Tpt[NumberOfProcesses];
  TH1F *DRWH[NumberOfProcesses];
  TH1F *DPHjets[NumberOfProcesses];
  TH1F *DPWjets[NumberOfProcesses];
  TH1F *DPTjets[NumberOfProcesses];
  TH1F *HiggsMass[NumberOfProcesses];
  TH1F *RelHT[NumberOfProcesses];
  TH1F *DRTH[NumberOfProcesses];
  TH1F *PtNormalizedMass[NumberOfProcesses];
  TH1F *RelativeMass[NumberOfProcesses];
  TH1F *MotherPtNormalizedMass[NumberOfProcesses];
  int EntriePerSample[NumberOfProcesses];
  bool SurvivalMarker[NumberOfProcesses];

  for (int k=ParcialTestMin; k<ParcialTestMax; k++)
    {
      TChain CutsChain("cuts");
      TChain AnalysisChain("stp");
      TH1F *ALLCuts[NumberOfHistos];
      for (int i=0; i<NumberOfSamples[k]; i++)
	{
	  if (k==0)
	    {
	      CutsChain.Add(StorageDirPrefix[k] + SamplesZZ[i]);
	      AnalysisChain.Add(StorageDirPrefix[k] + SamplesZZ[i]);
	    }
	  else if (k==1)
	    {
	      CutsChain.Add(StorageDirPrefix[k] + SamplesWZ[i]);
	      AnalysisChain.Add(StorageDirPrefix[k] + SamplesWZ[i]);
	    }
	  else if (k==2)
	    {
	      CutsChain.Add(StorageDirPrefix[k] + SamplesWW[i]);
	      AnalysisChain.Add(StorageDirPrefix[k] + SamplesWW[i]);
	    }
	}
      
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
      
      EntriePerSample[k]=EntriesPerCut[0];
      cout << PassedPerCut[NumberOfHistos-1] << endl;
      gPad->Close();
      if (PassedPerCut[NumberOfHistos-1]!=0)
	{
	  string A1 = Form("Reconstructed_Tprime.M() >> TprimeMass%i(60,400,1600)",k);
	  string A2 = Form("TprimeMass%i",k);
	  AnalysisChain.Draw(A1.c_str());
	  FiveJetsMass[k] = (TH1F*)gDirectory->Get(A2.c_str());
	  gPad->Close();
	  SurvivalMarker[k]=true;
	  /////////////////////
	  //Saving Histograms//
	  /////////////////////
	  string LJPT1 = Form("jet1_pt >> jet1_pt%i(100,10,500)",k);
	  string LJPT2 = Form("jet1_pt%i",k);
	  AnalysisChain.Draw(LJPT1.c_str());
	  LeadingJetPT[k] = (TH1F*)gDirectory->Get(LJPT2.c_str());
	  gPad->Close();
	  //
	  string THT1 = Form("THT >> THT%i(100,600,1400)",k);
	  string THT2 = Form("THT%i",k);
	  AnalysisChain.Draw(THT1.c_str());
	  THT[k] = (TH1F*)gDirectory->Get(THT2.c_str());
	  gPad->Close(); 
	  //
	  string DRHJ1 = Form("DeltaR_of_Higgs_Jets >> DeltaR_of_Higgs_Jets%i(40,0.4,7)",k);
	  string DRHJ2 = Form("DeltaR_of_Higgs_Jets%i",k);
	  AnalysisChain.Draw(DRHJ1.c_str());
	  DRHjets[k] = (TH1F*)gDirectory->Get(DRHJ2.c_str());
	  gPad->Close();
	  //
	  string DRWJ1 = Form("DeltaR_of_W_Jets >> DeltaR_of_W_Jets%i(40,0.4,7)",k);
	  string DRWJ2 = Form("DeltaR_of_W_Jets%i",k);
	  AnalysisChain.Draw(DRWJ1.c_str(), "DeltaR_of_W_Jets<3");
	  DRWjets[k] = (TH1F*)gDirectory->Get(DRWJ2.c_str());
	  gPad->Close();
	  //
	  string HPT1 = Form("Reconstructed_Higgs.Pt() >> HPt%i(200,10,800)",k);
	  string HPT2 = Form("HPt%i",k);
	  AnalysisChain.Draw(HPT1.c_str());
	  Hpt[k] = (TH1F*)gDirectory->Get(HPT2.c_str());
	  gPad->Close();
	  //
	  string TPT1 = Form("Reconstructed_Top.Pt() >> TPt%i(200,10,800)",k);
	  string TPT2 = Form("TPt%i",k);
	  AnalysisChain.Draw(TPT1.c_str());
	  Tpt[k] = (TH1F*)gDirectory->Get(TPT2.c_str());
	  gPad->Close();
	  //
	  string DRWH1 = Form("DeltaR_of_W_Higgs >> DeltaR_of_W_Higgs%i(40,0.4,7)",k);
	  string DRWH2 = Form("DeltaR_of_W_Higgs%i",k);
	  AnalysisChain.Draw(DRWH1.c_str());
	  DRWH[k] = (TH1F*)gDirectory->Get(DRWH2.c_str());
	  gPad->Close();
	  //
	  string DPHJ1 = Form("TMath::Abs(First_Higgs_Jet.Phi()-Second_Higgs_Jet.Phi()) >> DeltaPhi_of_Higgs_jets%i(40,0.0,3.5)",k);
	  string DPHJ2 = Form("DeltaPhi_of_Higgs_jets%i",k);
	  AnalysisChain.Draw(DPHJ1.c_str());
	  DPHjets[k] = (TH1F*)gDirectory->Get(DPHJ2.c_str());
	  gPad->Close();
	  //
	  string DPWJ1 = Form("TMath::Abs(First_W_Jet.Phi()-Second_W_Jet.Phi()) >> DeltaPhi_of_W_jets%i(40,0.0,3.5)",k);
	  string DPWJ2 = Form("DeltaPhi_of_W_jets%i",k);
	  AnalysisChain.Draw(DPWJ1.c_str());
	  DPWjets[k] = (TH1F*)gDirectory->Get(DPWJ2.c_str());
	  gPad->Close();
	  //
	  string DPTJ1 = Form("TMath::Abs(Top_Jet.Phi()-Reconstructed_W.Phi()) >> DeltaPhi_of_T_jet%i(40,0.0,3.5)",k);
	  string DPTJ2 = Form("DeltaPhi_of_T_jet%i",k);
	  AnalysisChain.Draw(DPTJ1.c_str());
	  DPTjets[k] = (TH1F*)gDirectory->Get(DPTJ2.c_str());
	  gPad->Close();
	  //
	  string HM1 = Form("Reconstructed_Higgs.M() >> HM%i(30,60,180)",k);
	  string HM2 = Form("HM%i",k);
	  AnalysisChain.Draw(HM1.c_str());
	  HiggsMass[k] = (TH1F*)gDirectory->Get(HM2.c_str());
	  gPad->Close();
	  //
	  string RHT1 = Form("Relative_THT >> RelHT%i(20,0,1)",k);
	  string RHT2 = Form("RelHT%i",k);
	  AnalysisChain.Draw(RHT1.c_str());
	  RelHT[k] = (TH1F*)gDirectory->Get(RHT2.c_str());
	  gPad->Close();
	  //
	  string DRTH1 = Form("DeltaR_of_Top_Higgs >> DeltaR_of_Top_Higgs%i(40,0.4,7)",k);
	  string DRTH2 = Form("DeltaR_of_Top_Higgs%i",k);
	  AnalysisChain.Draw(DRTH1.c_str());
	  DRTH[k] = (TH1F*)gDirectory->Get(DRTH2.c_str());
	  gPad->Close();
	  //
	  string PTNM1 = Form("PT_Normalized_Mass >> PT_Normalized_Mass%i(40,0.4,7)",k);
	  string PTNM2 = Form("PT_Normalized_Mass%i",k);
	  AnalysisChain.Draw(PTNM1.c_str());
	  PtNormalizedMass[k] = (TH1F*)gDirectory->Get(PTNM2.c_str());
	  gPad->Close();
	  //
	  string RM1 = Form("Relative_Mass >> Relative_Mass%i(20,0.0,1)",k);
	  string RM2 = Form("Relative_Mass%i",k);
	  AnalysisChain.Draw(RM1.c_str());
	  RelativeMass[k] = (TH1F*)gDirectory->Get(RM2.c_str());
	  gPad->Close();
	  //
	  string MPTNM1 = Form("Mother_PT_Normalized_Mass >> Mother_PT_Normalized_Mass%i(40,0.4,7)",k);
	  string MPTNM2 = Form("Mother_PT_Normalized_Mass%i",k);
	  AnalysisChain.Draw(MPTNM1.c_str());
	  MotherPtNormalizedMass[k] = (TH1F*)gDirectory->Get(MPTNM2.c_str());
	  gPad->Close();
	}

    }

  //TprimeMass histo recollection
  TH1F *DibosonTprimeMass=(TH1F*)gDirectory->Get("TprimeMass1");
  //LeadingJetPT histo recollection
  TH1F *DibosonTprimeMassLJPT=(TH1F*)gDirectory->Get("jet1_pt1");
  //THT histo recollection
  TH1F *DibosonTprimeMassTHT=(TH1F*)gDirectory->Get("THT1");
  //DR Higgs jets histo recollection
  TH1F *DibosonTprimeMassDRHJ=(TH1F*)gDirectory->Get("DeltaR_of_Higgs_Jets1");
  //DR W jets histo recollection
  TH1F *DibosonTprimeMassDRWJ=(TH1F*)gDirectory->Get("DeltaR_of_W_Jets1");
  //Higgs PT histo recollection
  TH1F *DibosonTprimeMassHPT=(TH1F*)gDirectory->Get("HPt1");
  //Top PT histo recollection
  TH1F *DibosonTprimeMassTPT=(TH1F*)gDirectory->Get("TPt1");
  //DR W Higgs histo recollection
  TH1F *DibosonTprimeMassDRWH=(TH1F*)gDirectory->Get("DeltaR_of_W_Higgs1");
  //DP Higgs jets histo recollection
  TH1F *DibosonTprimeMassDPHJ=(TH1F*)gDirectory->Get("DeltaPhi_of_Higgs_jets1");
  //DP W jets histo recollection
  TH1F *DibosonTprimeMassDPWJ=(TH1F*)gDirectory->Get("DeltaPhi_of_W_jets1");
  //DP Top jet histo recollection
  TH1F *DibosonTprimeMassDPTJ=(TH1F*)gDirectory->Get("DeltaPhi_of_T_jet1");
  //Higgs Mass histo recollection
  TH1F *DibosonTprimeMassHM=(TH1F*)gDirectory->Get("HM1");
  //Relative HT histo recollection
  TH1F *DibosonTprimeMassRelHT=(TH1F*)gDirectory->Get("RelHT1");
  //DR Top Higgs histo recollection
  TH1F *DibosonTprimeMassDRTH=(TH1F*)gDirectory->Get("DeltaR_of_Top_Higgs1");
  //PT Normalized Mass histo recollection
  TH1F *DibosonTprimeMassPTNM=(TH1F*)gDirectory->Get("PT_Normalized_Mass1");
  //Relative Mass histo recollection
  TH1F *DibosonTprimeMassRM=(TH1F*)gDirectory->Get("Relative_Mass1");
  //Mother PT Normalized Mass histo recollection
  TH1F *DibosonTprimeMassMPTNM=(TH1F*)gDirectory->Get("Mother_PT_Normalized_Mass1");
  for (int k=ParcialTestMin; k<ParcialTestMax; k++)
    {
      if (!SurvivalMarker[k]) continue;
      FiveJetsMass[k]->Scale(Lumi*XS[k]/EntriePerSample[k]);      
      LeadingJetPT[k]->Scale(Lumi*XS[k]/EntriePerSample[k]);
      THT[k]->Scale(Lumi*XS[k]/EntriePerSample[k]);
      DRHjets[k]->Scale(Lumi*XS[k]/EntriePerSample[k]);
      DRWjets[k]->Scale(Lumi*XS[k]/EntriePerSample[k]);
      Hpt[k]->Scale(Lumi*XS[k]/EntriePerSample[k]);
      Tpt[k]->Scale(Lumi*XS[k]/EntriePerSample[k]);
      DRWH[k]->Scale(Lumi*XS[k]/EntriePerSample[k]);
      DPHjets[k]->Scale(Lumi*XS[k]/EntriePerSample[k]);
      DPWjets[k]->Scale(Lumi*XS[k]/EntriePerSample[k]);
      DPTjets[k]->Scale(Lumi*XS[k]/EntriePerSample[k]);
      HiggsMass[k]->Scale(Lumi*XS[k]/EntriePerSample[k]);
      RelHT[k]->Scale(Lumi*XS[k]/EntriePerSample[k]);
      DRTH[k]->Scale(Lumi*XS[k]/EntriePerSample[k]);
      PtNormalizedMass[k]->Scale(Lumi*XS[k]/EntriePerSample[k]);
      RelativeMass[k]->Scale(Lumi*XS[k]/EntriePerSample[k]);
      MotherPtNormalizedMass[k]->Scale(Lumi*XS[k]/EntriePerSample[k]);
      //Settings for DiBoson
      if (k<=2 && k>=0) 
	{
	  FiveJetsMass[k]->SetFillColor(kWhite); 
	  LeadingJetPT[k]->SetFillColor(kWhite); 
	  THT[k]->SetFillColor(kWhite); 
	  DRHjets[k]->SetFillColor(kWhite);
	  DRWjets[k]->SetFillColor(kWhite);
	  Hpt[k]->SetFillColor(kWhite);
	  Tpt[k]->SetFillColor(kWhite);
	  DRWH[k]->SetFillColor(kWhite);
	  DPHjets[k]->SetFillColor(kWhite);
	  DPWjets[k]->SetFillColor(kWhite);
	  DPTjets[k]->SetFillColor(kWhite);
	  HiggsMass[k]->SetFillColor(kWhite);
	  RelHT[k]->SetFillColor(kWhite);
	  DRTH[k]->SetFillColor(kWhite);
	  PtNormalizedMass[k]->SetFillColor(kWhite);
	  RelativeMass[k]->SetFillColor(kWhite);
	  MotherPtNormalizedMass[k]->SetFillColor(kWhite);
	  if (k!=0) 
	    {
	      DibosonTprimeMass->Add(FiveJetsMass[k]);
	      DibosonTprimeMassLJPT->Add(LeadingJetPT[k]);
	      DibosonTprimeMassTHT->Add(THT[k]);
	      DibosonTprimeMassDRHJ->Add(DRHjets[k]);
	      DibosonTprimeMassDRWJ->Add(DRWjets[k]);
	      DibosonTprimeMassHPT->Add(Hpt[k]);
	      DibosonTprimeMassTPT->Add(Tpt[k]);
	      DibosonTprimeMassDRWH->Add(DRWH[k]);
	      DibosonTprimeMassDPHJ->Add(DPHjets[k]);
	      DibosonTprimeMassDPWJ->Add(DPWjets[k]);
	      DibosonTprimeMassDPTJ->Add(DPTjets[k]);
	      DibosonTprimeMassHM->Add(HiggsMass[k]);
	      DibosonTprimeMassRelHT->Add(RelHT[k]);
	      DibosonTprimeMassDRTH->Add(DRTH[k]);
	      DibosonTprimeMassPTNM->Add(PtNormalizedMass[k]);
	      DibosonTprimeMassRM->Add(RelativeMass[k]);
	      DibosonTprimeMassMPTNM->Add(MotherPtNormalizedMass[k]);
	    }
	  if (k==2)
	    {
	      TFile f("Diboson.root", "RECREATE"); 
	      DibosonTprimeMass->Write();
	      DibosonTprimeMassLJPT->Write();
	      DibosonTprimeMassTHT->Write();
	      DibosonTprimeMassDRHJ->Write();
	      DibosonTprimeMassDRWJ->Write();
	      DibosonTprimeMassHPT->Write();
	      DibosonTprimeMassTPT->Write();
	      DibosonTprimeMassDRWH->Write();
	      DibosonTprimeMassDPHJ->Write();
	      DibosonTprimeMassDPWJ->Write();
	      DibosonTprimeMassDPTJ->Write();
	      DibosonTprimeMassHM->Write();
	      DibosonTprimeMassRelHT->Write();
	      DibosonTprimeMassDRTH->Write();
	      DibosonTprimeMassPTNM->Write();
	      DibosonTprimeMassRM->Write();
	      DibosonTprimeMassMPTNM->Write();
	    }
	}

    }
  
}
