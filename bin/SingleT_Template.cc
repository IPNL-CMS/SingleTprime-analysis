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

using namespace std;

// Global Parameters
// Location of root files

const int NumberOfProcesses=16;

const TString MainFolder = "file:/afs/cern.ch/work/j/jruizalv/private/Analysis/CMSSW_5_3_9_patch3/src/Extractors/PatExtractor/bin/WithTrigger/SUFFIX/";

const int NumberOfHistos=23;
const TString Histos[NumberOfHistos] = {"Cut_0", "Cut_1", "Cut_2", "Cut_3", "Cut_chi2", "Cut_4", "Cut_5", "Cut_6", "Cut_7", "Cut_8", "Cut_9", "Cut_10", "Cut_11", "Cut_12", "Cut_13", 
"Cut_14", "Cut_15", "Cut_16", "Cut_17", "Cut_18", "Cut_19", "Cut_20", "Cut_21"};

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

void SingleT()
{
  int ParcialTestMax=9;
  int ParcialTestMin=4;
  TH1F *FiveJetsMass[NumberOfProcesses];
  TH1F *LeadingJetPT[NumberOfProcesses];
  TH1F *Leading2JetPT[NumberOfProcesses];
  TH1F *Leading3JetPT[NumberOfProcesses];
  TH1F *Leading4JetPT[NumberOfProcesses];
  TH1F *Leading5JetPT[NumberOfProcesses];
  TH1F *Leading6JetPT[NumberOfProcesses];
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
  TH1F *NumberOfTops[NumberOfProcesses];
  TH1F *HiggsMassOverTopMass[NumberOfProcesses];
  TH1F *HiggsTopAsymmetry[NumberOfProcesses];  
  TH1F *ThirdLooseBtag[NumberOfProcesses];
  TH1F *TopMass[NumberOfProcesses];
  TH1F *Chi2[NumberOfProcesses];
  TH1F *UQuarkContent[NumberOfProcesses];
  TH1F *DQuarkContent[NumberOfProcesses];
  TH1F *SQuarkContent[NumberOfProcesses];
  TH1F *CQuarkContent[NumberOfProcesses];
  TH1F *BQuarkContent[NumberOfProcesses];
  int EntriePerSample[NumberOfProcesses];
  bool SurvivalMarker[NumberOfProcesses];

  for (int k=ParcialTestMin; k<ParcialTestMax; k++)
    {
      TChain CutsChain("cuts");
      TChain AnalysisChain("stp");
      TH1F *ALLCuts[NumberOfHistos];
      if (k==3)
	{
	  CutsChain.Add(MainFolder + "T-s_Full_analyzed.root");
	  AnalysisChain.Add(MainFolder + "T-s_Full_analyzed.root");
	}
      else if (k==4)
	{
	  CutsChain.Add(MainFolder + "T-t_Full_analyzed.root");
	  AnalysisChain.Add(MainFolder + "T-t_Full_analyzed.root");
	}
      else if (k==5)
	{
	  CutsChain.Add(MainFolder + "T-tw_Full_analyzed.root");
	  AnalysisChain.Add(MainFolder + "T-tw_Full_analyzed.root");
	}
      else if (k==6)
	{
	  CutsChain.Add(MainFolder + "Tbar-s_Full_analyzed.root");
	  AnalysisChain.Add(MainFolder + "Tbar-s_Full_analyzed.root");
	}
      else if (k==7)
	{
	  CutsChain.Add(MainFolder + "Tbar-t_Full_analyzed.root");
	  AnalysisChain.Add(MainFolder + "Tbar-t_Full_analyzed.root");
	}
      else if (k==8)
	{
	  CutsChain.Add(MainFolder + "Tbar-tw_Full_analyzed.root");
	  AnalysisChain.Add(MainFolder + "Tbar-tw_Full_analyzed.root");
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
      SurvivalMarker[k]=false;
      if (PassedPerCut[NumberOfHistos-1]!=0)
	{
	  /////////////////////
	  //Saving Histograms//
	  /////////////////////
	  string A1 = Form("Reconstructed_Tprime.M() >> TprimeMass%i(60,400,1600)",k);
	  string A2 = Form("TprimeMass%i",k);
	  AnalysisChain.Draw(A1.c_str());
	  FiveJetsMass[k] = (TH1F*)gDirectory->Get(A2.c_str());
	  gPad->Close();
	  SurvivalMarker[k]=true;
	  //
	  string LJPT1 = Form("jet1_pt >> jet1_pt%i(50,20,1000)",k);
	  string LJPT2 = Form("jet1_pt%i",k);
	  AnalysisChain.Draw(LJPT1.c_str());
	  LeadingJetPT[k] = (TH1F*)gDirectory->Get(LJPT2.c_str());
	  gPad->Close();
	  //
	  string SLJPT1 = Form("jet2_pt >> jet2_pt%i(50,20,1000)",k);
	  string SLJPT2 = Form("jet2_pt%i",k);
	  AnalysisChain.Draw(SLJPT1.c_str());
	  Leading2JetPT[k] = (TH1F*)gDirectory->Get(SLJPT2.c_str());
	  gPad->Close();
	  //
	  string SSLJPT1 = Form("jet3_pt >> jet3_pt%i(35,20,700)",k);
	  string SSLJPT2 = Form("jet3_pt%i",k);
	  AnalysisChain.Draw(SSLJPT1.c_str());
	  Leading3JetPT[k] = (TH1F*)gDirectory->Get(SSLJPT2.c_str());
	  gPad->Close();
	  //
	  string SSSLJPT1 = Form("jet4_pt >> jet4_pt%i(20,20,400)",k);
	  string SSSLJPT2 = Form("jet4_pt%i",k);
	  AnalysisChain.Draw(SSSLJPT1.c_str());
	  Leading4JetPT[k] = (TH1F*)gDirectory->Get(SSSLJPT2.c_str());
	  gPad->Close();
	  //
	  string SSSSLJPT1 = Form("jet5_pt >> jet5_pt%i(15,20,300)",k);
	  string SSSSLJPT2 = Form("jet5_pt%i",k);
	  AnalysisChain.Draw(SSSSLJPT1.c_str());
	  Leading5JetPT[k] = (TH1F*)gDirectory->Get(SSSSLJPT2.c_str());
	  gPad->Close();
	  //
	  string SSSSSLJPT1 = Form("jet6_pt >> jet6_pt%i(10,20,200)",k);
	  string SSSSSLJPT2 = Form("jet6_pt%i",k);
	  AnalysisChain.Draw(SSSSSLJPT1.c_str());
	  Leading6JetPT[k] = (TH1F*)gDirectory->Get(SSSSSLJPT2.c_str());
	  gPad->Close();
	  //
	  string THT1 = Form("THT >> THT%i(65,300,1600)",k);
	  string THT2 = Form("THT%i",k);
	  AnalysisChain.Draw(THT1.c_str());
	  THT[k] = (TH1F*)gDirectory->Get(THT2.c_str());
	  gPad->Close(); 
	  //
	  string DRHJ1 = Form("DeltaR_of_Higgs_Jets >> DeltaR_of_Higgs_Jets%i(65,0.5,7)",k);
	  string DRHJ2 = Form("DeltaR_of_Higgs_Jets%i",k);
	  AnalysisChain.Draw(DRHJ1.c_str());
	  DRHjets[k] = (TH1F*)gDirectory->Get(DRHJ2.c_str());
	  gPad->Close();
	  //
	  string DRWJ1 = Form("DeltaR_of_W_Jets >> DeltaR_of_W_Jets%i(65,0.5,7)",k);
	  string DRWJ2 = Form("DeltaR_of_W_Jets%i",k);
	  AnalysisChain.Draw(DRWJ1.c_str());
	  DRWjets[k] = (TH1F*)gDirectory->Get(DRWJ2.c_str());
	  gPad->Close();
	  //
	  string HPT1 = Form("Reconstructed_Higgs.Pt() >> HPt%i(40,10,800)",k);
	  string HPT2 = Form("HPt%i",k);
	  AnalysisChain.Draw(HPT1.c_str());
	  Hpt[k] = (TH1F*)gDirectory->Get(HPT2.c_str());
	  gPad->Close();
	  //
	  string TPT1 = Form("Reconstructed_Top.Pt() >> TPt%i(40,10,800)",k);
	  string TPT2 = Form("TPt%i",k);
	  AnalysisChain.Draw(TPT1.c_str());
	  Tpt[k] = (TH1F*)gDirectory->Get(TPT2.c_str());
	  gPad->Close();
	  //
	  string DRWH1 = Form("DeltaR_of_W_Higgs >> DeltaR_of_W_Higgs%i(65,0.5,7)",k);
	  string DRWH2 = Form("DeltaR_of_W_Higgs%i",k);
	  AnalysisChain.Draw(DRWH1.c_str());
	  DRWH[k] = (TH1F*)gDirectory->Get(DRWH2.c_str());
	  gPad->Close();
	  //
	  string DPHJ1 = Form("TMath::Abs(First_Higgs_Jet.Phi()-Second_Higgs_Jet.Phi()) >> DeltaPhi_of_Higgs_jets%i(30,0.0,3.0)",k);
	  string DPHJ2 = Form("DeltaPhi_of_Higgs_jets%i",k);
	  AnalysisChain.Draw(DPHJ1.c_str());
	  DPHjets[k] = (TH1F*)gDirectory->Get(DPHJ2.c_str());
	  gPad->Close();
	  //
	  string DPWJ1 = Form("TMath::Abs(First_W_Jet.Phi()-Second_W_Jet.Phi()) >> DeltaPhi_of_W_jets%i(30,0.0,3.0)",k);
	  string DPWJ2 = Form("DeltaPhi_of_W_jets%i",k);
	  AnalysisChain.Draw(DPWJ1.c_str());
	  DPWjets[k] = (TH1F*)gDirectory->Get(DPWJ2.c_str());
	  gPad->Close();
	  //
	  string DPTJ1 = Form("TMath::Abs(Top_Jet.Phi()-Reconstructed_W.Phi()) >> DeltaPhi_of_T_jet%i(30,0.0,3.0)",k);
	  string DPTJ2 = Form("DeltaPhi_of_T_jet%i",k);
	  AnalysisChain.Draw(DPTJ1.c_str());
	  DPTjets[k] = (TH1F*)gDirectory->Get(DPTJ2.c_str());
	  gPad->Close();
	  //
	  string HM1 = Form("Reconstructed_Higgs.M() >> HM%i(36,60,180)",k);
	  string HM2 = Form("HM%i",k);
	  AnalysisChain.Draw(HM1.c_str());
	  HiggsMass[k] = (TH1F*)gDirectory->Get(HM2.c_str());
	  gPad->Close();
	  //
	  string RHT1 = Form("Relative_THT >> RelHT%i(30,0,1)",k);
	  string RHT2 = Form("RelHT%i",k);
	  AnalysisChain.Draw(RHT1.c_str());
	  RelHT[k] = (TH1F*)gDirectory->Get(RHT2.c_str());
	  gPad->Close();
	  //
	  string DRTH1 = Form("DeltaR_of_Top_Higgs >> DeltaR_of_Top_Higgs%i(65,0.5,7)",k);
	  string DRTH2 = Form("DeltaR_of_Top_Higgs%i",k);
	  AnalysisChain.Draw(DRTH1.c_str());
	  DRTH[k] = (TH1F*)gDirectory->Get(DRTH2.c_str());
	  gPad->Close();
	  //
	  string PTNM1 = Form("PT_Normalized_Mass >> PT_Normalized_Mass%i(70,0.4,5)",k);
	  string PTNM2 = Form("PT_Normalized_Mass%i",k);
	  AnalysisChain.Draw(PTNM1.c_str());
	  PtNormalizedMass[k] = (TH1F*)gDirectory->Get(PTNM2.c_str());
	  gPad->Close();
	  //
	  string RM1 = Form("Relative_Mass >> Relative_Mass%i(30,0.0,1)",k);
	  string RM2 = Form("Relative_Mass%i",k);
	  AnalysisChain.Draw(RM1.c_str());
	  RelativeMass[k] = (TH1F*)gDirectory->Get(RM2.c_str());
	  gPad->Close();
	  //
	  string MPTNM1 = Form("Mother_PT_Normalized_Mass >> Mother_PT_Normalized_Mass%i(25,0.0,50)",k);
	  string MPTNM2 = Form("Mother_PT_Normalized_Mass%i",k);
	  AnalysisChain.Draw(MPTNM1.c_str());
	  MotherPtNormalizedMass[k] = (TH1F*)gDirectory->Get(MPTNM2.c_str());
	  gPad->Close();
	  //
	  string NTops1 = Form("Number_of_Tops >> Number_of_Tops%i(8,0.0,8)",k);
	  string NTops2 = Form("Number_of_Tops%i",k);
	  AnalysisChain.Draw(NTops1.c_str());
	  NumberOfTops[k] = (TH1F*)gDirectory->Get(NTops2.c_str());
	  gPad->Close();
	  ///////////NEW VARIABLES//////////////////
	  string HMTM1 = Form("Reconstructed_Higgs.M()/Reconstructed_Top.M() >> HMoverTM%i(30,0,1.)",k);
	  string HMTM2 = Form("HMoverTM%i",k);
	  AnalysisChain.Draw(HMTM1.c_str());
	  HiggsMassOverTopMass[k] = (TH1F*)gDirectory->Get(HMTM2.c_str());
	  gPad->Close();
	  //
	  string HTA1 = Form("((Reconstructed_Higgs.Pt()/Reconstructed_Higgs.M())-(Reconstructed_Top.Pt()/Reconstructed_Top.M()))/((Reconstructed_Higgs.Pt()/Reconstructed_Higgs.M())+(Reconstructed_Top.Pt()/Reconstructed_Top.M())) >> HTAsym%i(12,0,1.)",k);
	  string HTA2 = Form("HTAsym%i",k);
	  AnalysisChain.Draw(HTA1.c_str());
	  HiggsTopAsymmetry[k] = (TH1F*)gDirectory->Get(HTA2.c_str());
	  gPad->Close();
	  //
	  string TLBT1 = Form("Number_of_Loose_and_non_med_b_tags >> TLBTag%i(10,0,10)",k);
	  string TLBT2 = Form("TLBTag%i",k);
	  AnalysisChain.Draw(TLBT1.c_str());
	  ThirdLooseBtag[k] = (TH1F*)gDirectory->Get(TLBT2.c_str());
	  gPad->Close();
	  //
	  string TM1 = Form("Reconstructed_Top.M() >> TMass%i(40,100,300)",k);
	  string TM2 = Form("TMass%i",k);
	  AnalysisChain.Draw(TM1.c_str());
	  TopMass[k] = (TH1F*)gDirectory->Get(TM2.c_str());
	  gPad->Close();
	  //
	  string C21 = Form("ChiSquaredSorting >> ChiSq%i(100,0,1000)",k);
	  string C22 = Form("ChiSq%i",k);
	  AnalysisChain.Draw(C21.c_str());
	  Chi2[k] = (TH1F*)gDirectory->Get(C22.c_str());
	  gPad->Close();
	  //
	  string UQ1 = Form("U Quark Content >> UQC%i(10,0,10)",k);
	  string UQ2 = Form("UQC%i",k);
	  AnalysisChain.Draw(UQ1.c_str());
	  UQuarkContent[k] = (TH1F*)gDirectory->Get(UQ2.c_str());
	  gPad->Close();
	  //
	  string DQ1 = Form("D Quark Content >> DQC%i(10,0,10)",k);
	  string DQ2 = Form("DQC%i",k);
	  AnalysisChain.Draw(DQ1.c_str());
	  DQuarkContent[k] = (TH1F*)gDirectory->Get(DQ2.c_str());
	  gPad->Close();
	  //
	  string SQ1 = Form("S Quark Content >> SQC%i(10,0,10)",k);
	  string SQ2 = Form("SQC%i",k);
	  AnalysisChain.Draw(SQ1.c_str());
	  SQuarkContent[k] = (TH1F*)gDirectory->Get(SQ2.c_str());
	  gPad->Close();
	  //
	  string CQ1 = Form("C Quark Content >> CQC%i(10,0,10)",k);
	  string CQ2 = Form("CQC%i",k);
	  AnalysisChain.Draw(CQ1.c_str());
	  CQuarkContent[k] = (TH1F*)gDirectory->Get(CQ2.c_str());
	  gPad->Close();
	  //
	  string BQ1 = Form("B Quark Content >> BQC%i(10,0,10)",k);
	  string BQ2 = Form("BQC%i",k);
	  AnalysisChain.Draw(BQ1.c_str());
	  BQuarkContent[k] = (TH1F*)gDirectory->Get(BQ2.c_str());
	  gPad->Close();
	}

    }

  //TprimeMass histo recollection
  TH1F *SingleTopTprimeMass=(TH1F*)gDirectory->Get("TprimeMass4");
  //LeadingJetPT histo recollection
  TH1F *SingleTopTprimeMassLJPT=(TH1F*)gDirectory->Get("jet1_pt4");
  TH1F *SingleTopTprimeMassL2JPT=(TH1F*)gDirectory->Get("jet2_pt4");
  TH1F *SingleTopTprimeMassL3JPT=(TH1F*)gDirectory->Get("jet3_pt4");
  TH1F *SingleTopTprimeMassL4JPT=(TH1F*)gDirectory->Get("jet4_pt4");
  TH1F *SingleTopTprimeMassL5JPT=(TH1F*)gDirectory->Get("jet5_pt4");
  TH1F *SingleTopTprimeMassL6JPT=(TH1F*)gDirectory->Get("jet6_pt4");
  //THT histo recollection
  TH1F *SingleTopTprimeMassTHT=(TH1F*)gDirectory->Get("THT4");
  //DR Higgs jets histo recollection
  TH1F *SingleTopTprimeMassDRHJ=(TH1F*)gDirectory->Get("DeltaR_of_Higgs_Jets4");
  //DR W jets histo recollection
  TH1F *SingleTopTprimeMassDRWJ=(TH1F*)gDirectory->Get("DeltaR_of_W_Jets4");
  //Higgs PT histo recollection
  TH1F *SingleTopTprimeMassHPT=(TH1F*)gDirectory->Get("HPt4");
  //Top PT histo recollection
  TH1F *SingleTopTprimeMassTPT=(TH1F*)gDirectory->Get("TPt4");
  //DR W Higgs histo recollection
  TH1F *SingleTopTprimeMassDRWH=(TH1F*)gDirectory->Get("DeltaR_of_W_Higgs4");
  //DP Higgs jets histo recollection
  TH1F *SingleTopTprimeMassDPHJ=(TH1F*)gDirectory->Get("DeltaPhi_of_Higgs_jets4");
  //DP W jets histo recollection
  TH1F *SingleTopTprimeMassDPWJ=(TH1F*)gDirectory->Get("DeltaPhi_of_W_jets4");
  //DP Top jet histo recollection
  TH1F *SingleTopTprimeMassDPTJ=(TH1F*)gDirectory->Get("DeltaPhi_of_T_jet4");
  //Higgs Mass histo recollection
  TH1F *SingleTopTprimeMassHM=(TH1F*)gDirectory->Get("HM4");
  //Relative HT histo recollection
  TH1F *SingleTopTprimeMassRelHT=(TH1F*)gDirectory->Get("RelHT4");
  //DR Top Higgs histo recollection
  TH1F *SingleTopTprimeMassDRTH=(TH1F*)gDirectory->Get("DeltaR_of_Top_Higgs4");
  //PT Normalized Mass histo recollection
  TH1F *SingleTopTprimeMassPTNM=(TH1F*)gDirectory->Get("PT_Normalized_Mass4");
  //Relative Mass histo recollection
  TH1F *SingleTopTprimeMassRM=(TH1F*)gDirectory->Get("Relative_Mass4");
  //Mother PT Normalized Mass histo recollection
  TH1F *SingleTopTprimeMassMPTNM=(TH1F*)gDirectory->Get("Mother_PT_Normalized_Mass4");
  //Number Of Tops Mass histo recollection
  TH1F *SingleTopTprimeMassNT=(TH1F*)gDirectory->Get("Number_of_Tops4");
  //Higgs Mass over Top Mass histo recollection
  TH1F *SingleTopTprimeMassHMTM=(TH1F*)gDirectory->Get("HMoverTM4");
  //Higgs Top Asym histo recollection
  TH1F *SingleTopTprimeMassHTAsym=(TH1F*)gDirectory->Get("HTAsym4");
  //Number of loose and non medium b-tag histo recollection
  TH1F *SingleTopTprimeMassTLBT=(TH1F*)gDirectory->Get("TLBTag4");
  //Top mass histo recollection
  TH1F *SingleTopTprimeMassTM=(TH1F*)gDirectory->Get("TMass4");
  //Chi2 histo recollection
  TH1F *SingleTopTprimeMassChi2=(TH1F*)gDirectory->Get("ChiSq4");
  //U quark content histo recollection
  TH1F *SingleTopTprimeMassUQC=(TH1F*)gDirectory->Get("UQC4");
  //D quark content histo recollection
  TH1F *SingleTopTprimeMassDQC=(TH1F*)gDirectory->Get("DQC4");
  //S quark content histo recollection
  TH1F *SingleTopTprimeMassSQC=(TH1F*)gDirectory->Get("SQC4");
  //C quark content histo recollection
  TH1F *SingleTopTprimeMassCQC=(TH1F*)gDirectory->Get("CQC4");
  //B quark content histo recollection
  TH1F *SingleTopTprimeMassBQC=(TH1F*)gDirectory->Get("BQC4");
  for (int k=ParcialTestMin; k<ParcialTestMax; k++)  
    {
      if (!SurvivalMarker[k]) continue;
      cout << k << endl;
      FiveJetsMass[k]->Scale(Lumi*XS[k]/EntriePerSample[k]);
      LeadingJetPT[k]->Scale(Lumi*XS[k]/EntriePerSample[k]);
      Leading2JetPT[k]->Scale(Lumi*XS[k]/EntriePerSample[k]);
      Leading3JetPT[k]->Scale(Lumi*XS[k]/EntriePerSample[k]);
      Leading4JetPT[k]->Scale(Lumi*XS[k]/EntriePerSample[k]);
      Leading5JetPT[k]->Scale(Lumi*XS[k]/EntriePerSample[k]);
      Leading6JetPT[k]->Scale(Lumi*XS[k]/EntriePerSample[k]);
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
      NumberOfTops[k]->Scale(Lumi*XS[k]/EntriePerSample[k]); 
      HiggsMassOverTopMass[k]->Scale(Lumi*XS[k]/EntriePerSample[k]);
      HiggsTopAsymmetry[k]->Scale(Lumi*XS[k]/EntriePerSample[k]);
      ThirdLooseBtag[k]->Scale(Lumi*XS[k]/EntriePerSample[k]);
      TopMass[k]->Scale(Lumi*XS[k]/EntriePerSample[k]);
      Chi2[k]->Scale(Lumi*XS[k]/EntriePerSample[k]);
      UQuarkContent[k]->Scale(Lumi*XS[k]/EntriePerSample[k]);
      DQuarkContent[k]->Scale(Lumi*XS[k]/EntriePerSample[k]);
      SQuarkContent[k]->Scale(Lumi*XS[k]/EntriePerSample[k]);
      CQuarkContent[k]->Scale(Lumi*XS[k]/EntriePerSample[k]);
      BQuarkContent[k]->Scale(Lumi*XS[k]/EntriePerSample[k]);      
      //Settings for Single Top
      if (k<=8 && k>=3) 
	{
	  FiveJetsMass[k]->SetFillColor(kBlack);
	  FiveJetsMass[k]->SetFillStyle(3305);
	  LeadingJetPT[k]->SetFillColor(kBlack);
	  LeadingJetPT[k]->SetFillStyle(3305);
	  Leading2JetPT[k]->SetFillColor(kBlack);
	  Leading2JetPT[k]->SetFillStyle(3305);
	  Leading3JetPT[k]->SetFillColor(kBlack);
	  Leading3JetPT[k]->SetFillStyle(3305);
	  Leading4JetPT[k]->SetFillColor(kBlack);
	  Leading4JetPT[k]->SetFillStyle(3305);
	  Leading5JetPT[k]->SetFillColor(kBlack);
	  Leading5JetPT[k]->SetFillStyle(3305);
	  Leading6JetPT[k]->SetFillColor(kBlack);
	  Leading6JetPT[k]->SetFillStyle(3305);
	  THT[k]->SetFillColor(kBlack);
	  THT[k]->SetFillStyle(3305);
	  DRHjets[k]->SetFillColor(kBlack);
	  DRHjets[k]->SetFillStyle(3305);
	  DRWjets[k]->SetFillColor(kBlack);
	  DRWjets[k]->SetFillStyle(3305);
	  Hpt[k]->SetFillColor(kBlack);
	  Hpt[k]->SetFillStyle(3305);
	  Tpt[k]->SetFillColor(kBlack);
	  Tpt[k]->SetFillStyle(3305);
	  DRWH[k]->SetFillColor(kBlack);
	  DRWH[k]->SetFillStyle(3305);
	  DPHjets[k]->SetFillColor(kBlack);
	  DPHjets[k]->SetFillStyle(3305);
	  DPWjets[k]->SetFillColor(kBlack);
	  DPWjets[k]->SetFillStyle(3305);
	  DPTjets[k]->SetFillColor(kBlack);
	  DPTjets[k]->SetFillStyle(3305);
	  HiggsMass[k]->SetFillColor(kBlack);
	  HiggsMass[k]->SetFillStyle(3305);
	  RelHT[k]->SetFillColor(kBlack);
	  RelHT[k]->SetFillStyle(3305);
	  DRTH[k]->SetFillColor(kBlack);
	  DRTH[k]->SetFillStyle(3305);
	  PtNormalizedMass[k]->SetFillColor(kBlack);
	  PtNormalizedMass[k]->SetFillStyle(3305);
	  RelativeMass[k]->SetFillColor(kBlack);
	  RelativeMass[k]->SetFillStyle(3305);
	  MotherPtNormalizedMass[k]->SetFillColor(kBlack);
	  MotherPtNormalizedMass[k]->SetFillStyle(3305);
	  NumberOfTops[k]->SetFillColor(kBlack);
	  NumberOfTops[k]->SetFillStyle(3305);
	  HiggsMassOverTopMass[k]->SetFillColor(kBlack);
	  HiggsMassOverTopMass[k]->SetFillStyle(3305);
	  HiggsTopAsymmetry[k]->SetFillColor(kBlack);
	  HiggsTopAsymmetry[k]->SetFillStyle(3305);
	  ThirdLooseBtag[k]->SetFillColor(kBlack);
	  ThirdLooseBtag[k]->SetFillStyle(3305);
	  TopMass[k]->SetFillColor(kBlack);
	  TopMass[k]->SetFillStyle(3305);
	  Chi2[k]->SetFillColor(kBlack);
	  Chi2[k]->SetFillStyle(3305);
	  UQuarkContent[k]->SetFillColor(kBlack);
	  UQuarkContent[k]->SetFillStyle(3305);
	  DQuarkContent[k]->SetFillColor(kBlack);
	  DQuarkContent[k]->SetFillStyle(3305);
	  SQuarkContent[k]->SetFillColor(kBlack);
	  SQuarkContent[k]->SetFillStyle(3305);
	  CQuarkContent[k]->SetFillColor(kBlack);
	  CQuarkContent[k]->SetFillStyle(3305);
	  BQuarkContent[k]->SetFillColor(kBlack);
	  BQuarkContent[k]->SetFillStyle(3305);	  
	  if (k!=ParcialTestMin) 
	    {
	      SingleTopTprimeMass->Add(FiveJetsMass[k]); 
	      SingleTopTprimeMassLJPT->Add(LeadingJetPT[k]);
	      SingleTopTprimeMassL2JPT->Add(Leading2JetPT[k]);
	      SingleTopTprimeMassL3JPT->Add(Leading3JetPT[k]);
	      SingleTopTprimeMassL4JPT->Add(Leading4JetPT[k]);
	      SingleTopTprimeMassL5JPT->Add(Leading5JetPT[k]);
	      SingleTopTprimeMassL6JPT->Add(Leading6JetPT[k]);
	      SingleTopTprimeMassTHT->Add(THT[k]);
	      SingleTopTprimeMassDRHJ->Add(DRHjets[k]);
	      SingleTopTprimeMassDRWJ->Add(DRWjets[k]);
	      SingleTopTprimeMassHPT->Add(Hpt[k]);
	      SingleTopTprimeMassTPT->Add(Tpt[k]);
	      SingleTopTprimeMassDRWH->Add(DRWH[k]);
	      SingleTopTprimeMassDPHJ->Add(DPHjets[k]);
	      SingleTopTprimeMassDPWJ->Add(DPWjets[k]);
	      SingleTopTprimeMassDPTJ->Add(DPTjets[k]);
	      SingleTopTprimeMassHM->Add(HiggsMass[k]);
	      SingleTopTprimeMassRelHT->Add(RelHT[k]);
	      SingleTopTprimeMassDRTH->Add(DRTH[k]);
	      SingleTopTprimeMassPTNM->Add(PtNormalizedMass[k]);
	      SingleTopTprimeMassRM->Add(RelativeMass[k]);
	      SingleTopTprimeMassMPTNM->Add(MotherPtNormalizedMass[k]);
	      SingleTopTprimeMassNT->Add(NumberOfTops[k]);
	      SingleTopTprimeMassHMTM->Add(HiggsMassOverTopMass[k]);
	      SingleTopTprimeMassHTAsym->Add(HiggsTopAsymmetry[k]);
	      SingleTopTprimeMassTLBT->Add(ThirdLooseBtag[k]);
	      SingleTopTprimeMassTM->Add(TopMass[k]);
	      SingleTopTprimeMassChi2->Add(Chi2[k]);
	      SingleTopTprimeMassUQC->Add(UQuarkContent[k]);
	      SingleTopTprimeMassDQC->Add(DQuarkContent[k]);
	      SingleTopTprimeMassSQC->Add(SQuarkContent[k]);
	      SingleTopTprimeMassCQC->Add(CQuarkContent[k]);
	      SingleTopTprimeMassBQC->Add(BQuarkContent[k]);
	    }
	  if (k==8)
	    {
	      TFile f("SingleTop.root", "RECREATE");
	      SingleTopTprimeMass->Write();
	      SingleTopTprimeMassLJPT->Write();
	      SingleTopTprimeMassL2JPT->Write();
	      SingleTopTprimeMassL3JPT->Write();
	      SingleTopTprimeMassL4JPT->Write();
	      SingleTopTprimeMassL5JPT->Write();
	      SingleTopTprimeMassL6JPT->Write();
	      SingleTopTprimeMassTHT->Write();
	      SingleTopTprimeMassDRHJ->Write();
	      SingleTopTprimeMassDRWJ->Write();
	      SingleTopTprimeMassHPT->Write();
	      SingleTopTprimeMassTPT->Write();
	      SingleTopTprimeMassDRWH->Write();
	      SingleTopTprimeMassDPHJ->Write();
	      SingleTopTprimeMassDPWJ->Write();
	      SingleTopTprimeMassDPTJ->Write();
	      SingleTopTprimeMassHM->Write();
	      SingleTopTprimeMassRelHT->Write();
	      SingleTopTprimeMassDRTH->Write();
	      SingleTopTprimeMassPTNM->Write();
	      SingleTopTprimeMassRM->Write();
	      SingleTopTprimeMassMPTNM->Write();
	      SingleTopTprimeMassNT->Write();
	      SingleTopTprimeMassHMTM->Write();
	      SingleTopTprimeMassHTAsym->Write();
	      SingleTopTprimeMassTLBT->Write();
	      SingleTopTprimeMassTM->Write();
	      SingleTopTprimeMassChi2->Write();
	      SingleTopTprimeMassUQC->Write();
	      SingleTopTprimeMassDQC->Write();
	      SingleTopTprimeMassSQC->Write();
	      SingleTopTprimeMassCQC->Write();
	      SingleTopTprimeMassBQC->Write();
	    }
	}
    }

  exit(0);
  
}
