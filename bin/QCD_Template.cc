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

void QCD()
{
  int ParcialTestMax=NumberOfProcesses-1; //Don't change!!!
  int ParcialTestMin=10; //Don't change!!!
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
  int EntriePerSample[NumberOfProcesses];
  bool SurvivalMarker[NumberOfProcesses];

  for (int k=ParcialTestMin; k<ParcialTestMax; k++)
    {
      TChain CutsChain("cuts");
      TChain AnalysisChain("stp");
      TH1F *ALLCuts[NumberOfHistos];

      if (k==10)
	{
	  CutsChain.Add(MainFolder + "QCD_PT_170_300_Full_analyzed.root");
	  AnalysisChain.Add(MainFolder + "QCD_PT_170_300_Full_analyzed.root");
	}
      else if (k==11)
	{
	  CutsChain.Add(MainFolder + "QCD_PT_300_470_Full_analyzed.root");
	  AnalysisChain.Add(MainFolder + "QCD_PT_300_470_Full_analyzed.root");
	}
      else if (k==12)
	{
	  CutsChain.Add(MainFolder + "QCD_PT_470_600_Full_analyzed.root");
	  AnalysisChain.Add(MainFolder + "QCD_PT_470_600_Full_analyzed.root");
	}
      else if (k==13)
	{
	  cout << "1 Marker" << endl; 
	  CutsChain.Add(MainFolder + "QCD_PT_600_800_Full_analyzed.root");
	  AnalysisChain.Add(MainFolder + "QCD_PT_600_800_Full_analyzed.root");
	  cout << "2 Marker" << endl;
	}
      else if (k==14)
	{
	  CutsChain.Add(MainFolder + "QCD_PT_800_1000_Full_analyzed.root");
	  AnalysisChain.Add(MainFolder + "QCD_PT_800_1000_Full_analyzed.root");
	}
      
      int EntriesPerCut[NumberOfHistos];
      float PassedPerCut[NumberOfHistos];
      cout << "3 Marker" << endl;
      
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
	  string TLBT1 = Form("Number_of_Loose_and_non_med_b_tags >> TLBTag%i(10,0,10.)",k);
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
	}

    }

  //TprimeMass histo recollection
  TH1F *QCDTprimeMass=(TH1F*)gDirectory->Get("TprimeMass10");  
  //LeadingJetPT histo recollection
  TH1F *QCDTprimeMassLJPT=(TH1F*)gDirectory->Get("jet1_pt10");  
  //LeadingJetPT histo recollection
  TH1F *QCDTprimeMassL2JPT=(TH1F*)gDirectory->Get("jet2_pt10");
  //LeadingJetPT histo recollection
  TH1F *QCDTprimeMassL3JPT=(TH1F*)gDirectory->Get("jet3_pt10");
  //LeadingJetPT histo recollection
  TH1F *QCDTprimeMassL4JPT=(TH1F*)gDirectory->Get("jet4_pt10");
  //LeadingJetPT histo recollection
  TH1F *QCDTprimeMassL5JPT=(TH1F*)gDirectory->Get("jet5_pt10");
  //LeadingJetPT histo recollection
  TH1F *QCDTprimeMassL6JPT=(TH1F*)gDirectory->Get("jet6_pt10");
  //THT histo recollection
  TH1F *QCDTprimeMassTHT=(TH1F*)gDirectory->Get("THT10");  
  //DR Higgs jets histo recollection
  TH1F *QCDTprimeMassDRHJ=(TH1F*)gDirectory->Get("DeltaR_of_Higgs_Jets10");
  //DR W jets histo recollection
  TH1F *QCDTprimeMassDRWJ=(TH1F*)gDirectory->Get("DeltaR_of_W_Jets10");
  //Higgs PT histo recollection
  TH1F *QCDTprimeMassHPT=(TH1F*)gDirectory->Get("HPt10");
  //Top PT histo recollection
  TH1F *QCDTprimeMassTPT=(TH1F*)gDirectory->Get("TPt10");
  //DR W Higgs histo recollection
  TH1F *QCDTprimeMassDRWH=(TH1F*)gDirectory->Get("DeltaR_of_W_Higgs10");
  //DP Higgs jets histo recollection
  TH1F *QCDTprimeMassDPHJ=(TH1F*)gDirectory->Get("DeltaPhi_of_Higgs_jets10");
  //DP W jets histo recollection
  TH1F *QCDTprimeMassDPWJ=(TH1F*)gDirectory->Get("DeltaPhi_of_W_jets10");
  //DP Top jet histo recollection
  TH1F *QCDTprimeMassDPTJ=(TH1F*)gDirectory->Get("DeltaPhi_of_T_jet10");
  //Higgs Mass histo recollection
  TH1F *QCDTprimeMassHM=(TH1F*)gDirectory->Get("HM10");
  //Relative HT histo recollection
  TH1F *QCDTprimeMassRelHT=(TH1F*)gDirectory->Get("RelHT10");
  //DR Top Higgs histo recollection
  TH1F *QCDTprimeMassDRTH=(TH1F*)gDirectory->Get("DeltaR_of_Top_Higgs10");
  //PT Normalized Mass histo recollection
  TH1F *QCDTprimeMassPTNM=(TH1F*)gDirectory->Get("PT_Normalized_Mass10");
  //Relative Mass histo recollection
  TH1F *QCDTprimeMassRM=(TH1F*)gDirectory->Get("Relative_Mass10");
  //Mother PT Normalized Mass histo recollection
  TH1F *QCDTprimeMassMPTNM=(TH1F*)gDirectory->Get("Mother_PT_Normalized_Mass10");
  //Number of Tops histo recollection
  TH1F *QCDTprimeMassNT=(TH1F*)gDirectory->Get("Number_of_Tops10");
  //Higgs Mass over Top Mass histo recollection
  TH1F *QCDTprimeMassHMTM=(TH1F*)gDirectory->Get("HMoverTM10");
  //Higgs Top Asym histo recollection
  TH1F *QCDTprimeMassHTAsym=(TH1F*)gDirectory->Get("HTAsym10");
  //Number of loose and non medium b-tag histo recollection
  TH1F *QCDTprimeMassTLBT=(TH1F*)gDirectory->Get("TLBTag10");
  //Top Mass histo recollection
  TH1F *QCDTprimeMassTM=(TH1F*)gDirectory->Get("TMass10");
  //Chi2 histo recollection
  TH1F *QCDTprimeMassChi2=(TH1F*)gDirectory->Get("ChiSq10");
  for (int k=ParcialTestMin; k<ParcialTestMax; k++)
    {
      if (!SurvivalMarker[k]) continue;
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
      //Settings for QCD
      if (k>=10 && k<=14)
	{
	  FiveJetsMass[k]->SetFillColor(kViolet); 
	  LeadingJetPT[k]->SetFillColor(kViolet); 
	  Leading2JetPT[k]->SetFillColor(kViolet);
	  Leading3JetPT[k]->SetFillColor(kViolet);
	  Leading4JetPT[k]->SetFillColor(kViolet);
	  Leading5JetPT[k]->SetFillColor(kViolet);
	  Leading6JetPT[k]->SetFillColor(kViolet);
	  THT[k]->SetFillColor(kViolet); 
	  DRHjets[k]->SetFillColor(kViolet);
	  DRWjets[k]->SetFillColor(kViolet);
	  Hpt[k]->SetFillColor(kViolet);
	  Tpt[k]->SetFillColor(kViolet);
	  DRWH[k]->SetFillColor(kViolet);
	  DPHjets[k]->SetFillColor(kViolet);
	  DPWjets[k]->SetFillColor(kViolet);
	  DPTjets[k]->SetFillColor(kViolet);
	  HiggsMass[k]->SetFillColor(kViolet);
	  RelHT[k]->SetFillColor(kViolet);
	  DRTH[k]->SetFillColor(kViolet);
	  PtNormalizedMass[k]->SetFillColor(kViolet);
	  RelativeMass[k]->SetFillColor(kViolet);
	  MotherPtNormalizedMass[k]->SetFillColor(kViolet);
	  NumberOfTops[k]->SetFillColor(kViolet);
	  HiggsMassOverTopMass[k]->SetFillColor(kViolet);
	  HiggsTopAsymmetry[k]->SetFillColor(kViolet);
	  ThirdLooseBtag[k]->SetFillColor(kViolet);
	  TopMass[k]->SetFillColor(kViolet);
	  Chi2[k]->SetFillColor(kViolet);
	  if (k!=10) 
	    {
	      QCDTprimeMass->Add(FiveJetsMass[k]);
	      QCDTprimeMassLJPT->Add(LeadingJetPT[k]);
	      QCDTprimeMassL2JPT->Add(Leading2JetPT[k]);
	      QCDTprimeMassL3JPT->Add(Leading3JetPT[k]);
	      QCDTprimeMassL4JPT->Add(Leading4JetPT[k]);
	      QCDTprimeMassL5JPT->Add(Leading5JetPT[k]);
	      QCDTprimeMassL6JPT->Add(Leading6JetPT[k]);
	      QCDTprimeMassTHT->Add(THT[k]);
	      QCDTprimeMassDRHJ->Add(DRHjets[k]);
	      QCDTprimeMassDRWJ->Add(DRWjets[k]);
	      QCDTprimeMassHPT->Add(Hpt[k]);
	      QCDTprimeMassTPT->Add(Tpt[k]);
	      QCDTprimeMassDRWH->Add(DRWH[k]);
	      QCDTprimeMassDPHJ->Add(DPHjets[k]);
	      QCDTprimeMassDPWJ->Add(DPWjets[k]);
	      QCDTprimeMassDPTJ->Add(DPTjets[k]);
	      QCDTprimeMassHM->Add(HiggsMass[k]);
	      QCDTprimeMassRelHT->Add(RelHT[k]);
	      QCDTprimeMassDRTH->Add(DRTH[k]);
	      QCDTprimeMassPTNM->Add(PtNormalizedMass[k]);
	      QCDTprimeMassRM->Add(RelativeMass[k]);
	      QCDTprimeMassMPTNM->Add(MotherPtNormalizedMass[k]);
	      QCDTprimeMassNT->Add(NumberOfTops[k]);
	      QCDTprimeMassHMTM->Add(HiggsMassOverTopMass[k]);
	      QCDTprimeMassHTAsym->Add(HiggsTopAsymmetry[k]);
	      QCDTprimeMassTLBT->Add(ThirdLooseBtag[k]);
	      QCDTprimeMassTM->Add(TopMass[k]);
	      QCDTprimeMassChi2->Add(Chi2[k]);
	    }
	  if (k==14)
	    {
	      TFile f("QCD.root", "RECREATE"); 
	      QCDTprimeMass->Write();
	      QCDTprimeMassLJPT->Write();
	      QCDTprimeMassL2JPT->Write();
	      QCDTprimeMassL3JPT->Write();
	      QCDTprimeMassL4JPT->Write();
	      QCDTprimeMassL5JPT->Write();
	      QCDTprimeMassL6JPT->Write();
	      QCDTprimeMassTHT->Write();
	      QCDTprimeMassDRHJ->Write();
	      QCDTprimeMassDRWJ->Write();
	      QCDTprimeMassHPT->Write();
	      QCDTprimeMassTPT->Write();
	      QCDTprimeMassDRWH->Write();
	      QCDTprimeMassDPHJ->Write();
	      QCDTprimeMassDPWJ->Write();
	      QCDTprimeMassDPTJ->Write();
	      QCDTprimeMassHM->Write();
	      QCDTprimeMassRelHT->Write();
	      QCDTprimeMassDRTH->Write();
	      QCDTprimeMassPTNM->Write();
	      QCDTprimeMassRM->Write();
	      QCDTprimeMassMPTNM->Write();
	      QCDTprimeMassNT->Write();
	      QCDTprimeMassHMTM->Write();
	      QCDTprimeMassHTAsym->Write();
	      QCDTprimeMassTLBT->Write();
	      QCDTprimeMassTM->Write();
	      QCDTprimeMassChi2->Write();
	    }
	}
    }

  exit(0);
  
}
