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
#include "TF2.h"
#include "PU_Reweighting.h"

using namespace std;

// Global Parameters
// Location of root files

//const int NumberOfProcesses=19;

const TString MainFolder = "file:/sps/cms/ruizalva/CMSSW_5_3_9_patch3/src/Extractors/PatExtractor/bin/WithTrigger/SUFFIX/";

const int NumberOfHistos=23;
const TString Histos[NumberOfHistos] = {"Cut_0", "Cut_1", "Cut_2", "Cut_3", "Cut_chi2", "Cut_4", "Cut_5", "Cut_6", "Cut_7", "Cut_8", "Cut_9", "Cut_10", "Cut_11", "Cut_12", "Cut_13", 
"Cut_14", "Cut_15", "Cut_16", "Cut_17", "Cut_18", "Cut_19", "Cut_20", "Cut_21"};

//Computing weights (Everything in pb)
//const float XS[NumberOfProcesses]={8.258, 33.72, 54.838, 2.82, 47, 10.7, 1.57, 25, 10.7, 234.0, 156293.3, 34138.15, 1759.549, 113.8791, 26.9921, 3.550036, 0.15, 3840.86, 3060.099}; 
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

//float Lumi=20000;

float PUR_function(int TI) //function with input the Number of True Interactions
{

  if (TI>=PUBins) return 1;
  else return PU_weight[TI];

}

void DY()
{
  int ParcialTestMax=21;
  int ParcialTestMin=19;
  TH1F *FiveJetsMass[NumberOfProcesses];
  TH1F *FiveJetsMass3L[NumberOfProcesses];
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
  TH1F *CSVLB[NumberOfProcesses];
  TH1F *CSVMB[NumberOfProcesses];
  TH1F *CSVTB[NumberOfProcesses];
  TH1F *Vtcs[NumberOfProcesses];   
  TH2F *HptTpt[6][NumberOfProcesses]; 
  TH2F *HTCSVM[NumberOfProcesses];  
  TH1F *HiggsMassCuts[6][NumberOfProcesses];
  TH1F *HiggsMassReversedHptTpt[6][NumberOfProcesses];
  TH1F *TprimeMassReversedHptTpt[6][NumberOfProcesses];    
  TH1F *WMassFromHiggs[6][NumberOfProcesses];  
  TH1F *TopMassFromHiggs[NumberOfProcesses];

  //ADDITIONAL INFO ON JETS  
  TH1F *JET1ETA[NumberOfProcesses];   
  TH1F *JET2ETA[NumberOfProcesses]; 
  TH1F *JET3ETA[NumberOfProcesses]; 
  TH1F *JET4ETA[NumberOfProcesses];
  TH1F *JET5ETA[NumberOfProcesses];
  TH1F *JET6ETA[NumberOfProcesses];
  TH1F *JET1PHI[NumberOfProcesses];   
  TH1F *JET2PHI[NumberOfProcesses]; 
  TH1F *JET3PHI[NumberOfProcesses]; 
  TH1F *JET4PHI[NumberOfProcesses];
  TH1F *JET5PHI[NumberOfProcesses];
  TH1F *JET6PHI[NumberOfProcesses];
  TH1F *JETMULTI[NumberOfProcesses];

  int EntriePerSample[NumberOfProcesses];
  bool SurvivalMarker[NumberOfProcesses];
  int FirstGoodMarker=ParcialTestMin;
  int LastGoodMarker=ParcialTestMax;
  int MarkerCounter=0;

  for (int k=ParcialTestMin; k<ParcialTestMax; k++)
    {
      TChain CutsChain("cuts");
      TChain AnalysisChain("stp");
      TChain AnalysisChain3L("stp3L");
      TH1F *ALLCuts[NumberOfHistos];
      if (k==19)
	{
	  CutsChain.Add(MainFolder + "DYToBB_Full_analyzed.root");
	  AnalysisChain.Add(MainFolder + "DYToBB_Full_analyzed.root");
	  AnalysisChain3L.Add(MainFolder + "DYToBB_Full_analyzed.root");
	}
      else if (k==20)
	{
	  CutsChain.Add(MainFolder + "DYToCC_Full_analyzed.root");
	  AnalysisChain.Add(MainFolder + "DYToCC_Full_analyzed.root");
	  AnalysisChain3L.Add(MainFolder + "DYToCC_Full_analyzed.root");
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
	  AnalysisChain.Draw(A1.c_str(),finalcut);
	  FiveJetsMass[k] = (TH1F*)gDirectory->Get(A2.c_str());
	  gPad->Close();
	  SurvivalMarker[k]=true;
	  if (MarkerCounter==0) FirstGoodMarker=k;
	  MarkerCounter++;
	  LastGoodMarker=k;
	  //
	  string LJPT1 = Form("jet1_pt >> jet1_pt%i(50,20,1000)",k);
	  string LJPT2 = Form("jet1_pt%i",k);
	  AnalysisChain.Draw(LJPT1.c_str(),"(PUR_function(Number_True_Interactions))*weight");
	  LeadingJetPT[k] = (TH1F*)gDirectory->Get(LJPT2.c_str());
	  gPad->Close();
	  //
	  string SLJPT1 = Form("jet2_pt >> jet2_pt%i(50,20,1000)",k);
	  string SLJPT2 = Form("jet2_pt%i",k);
	  AnalysisChain.Draw(SLJPT1.c_str(),"(PUR_function(Number_True_Interactions))*weight");
	  Leading2JetPT[k] = (TH1F*)gDirectory->Get(SLJPT2.c_str());
	  gPad->Close();
	  //
	  string SSLJPT1 = Form("jet3_pt >> jet3_pt%i(35,20,700)",k);
	  string SSLJPT2 = Form("jet3_pt%i",k);
	  AnalysisChain.Draw(SSLJPT1.c_str(),"(PUR_function(Number_True_Interactions))*weight");
	  Leading3JetPT[k] = (TH1F*)gDirectory->Get(SSLJPT2.c_str());
	  gPad->Close();
	  //
	  string SSSLJPT1 = Form("jet4_pt >> jet4_pt%i(20,20,400)",k);
	  string SSSLJPT2 = Form("jet4_pt%i",k);
	  AnalysisChain.Draw(SSSLJPT1.c_str(),"(PUR_function(Number_True_Interactions))*weight");
	  Leading4JetPT[k] = (TH1F*)gDirectory->Get(SSSLJPT2.c_str());
	  gPad->Close();
	  //
	  string SSSSLJPT1 = Form("jet5_pt >> jet5_pt%i(15,20,300)",k);
	  string SSSSLJPT2 = Form("jet5_pt%i",k);
	  AnalysisChain.Draw(SSSSLJPT1.c_str(),"(PUR_function(Number_True_Interactions))*weight");
	  Leading5JetPT[k] = (TH1F*)gDirectory->Get(SSSSLJPT2.c_str());
	  gPad->Close();
	  //
	  string SSSSSLJPT1 = Form("jet6_pt >> jet6_pt%i(10,20,200)",k);
	  string SSSSSLJPT2 = Form("jet6_pt%i",k);
	  AnalysisChain.Draw(SSSSSLJPT1.c_str(),"(PUR_function(Number_True_Interactions))*weight");
	  Leading6JetPT[k] = (TH1F*)gDirectory->Get(SSSSSLJPT2.c_str());
	  gPad->Close();
	  //
	  string THT1 = Form("THT >> THT%i(65,300,1600)",k);
	  string THT2 = Form("THT%i",k);
	  AnalysisChain.Draw(THT1.c_str(),"PUR_function(Number_True_Interactions)*weight*(THT>630)"); //PUR_function(Number_True_Interactions)*(THT>630)");
	  THT[k] = (TH1F*)gDirectory->Get(THT2.c_str());
	  gPad->Close(); 
	  //
	  string DRHJ1 = Form("DeltaR_of_Higgs_Jets >> DeltaR_of_Higgs_Jets%i(65,0.5,7)",k);
	  string DRHJ2 = Form("DeltaR_of_Higgs_Jets%i",k);
	  AnalysisChain.Draw(DRHJ1.c_str(),"(PUR_function(Number_True_Interactions))*weight");
	  DRHjets[k] = (TH1F*)gDirectory->Get(DRHJ2.c_str());
	  gPad->Close();
	  //
	  string DRWJ1 = Form("DeltaR_of_W_Jets >> DeltaR_of_W_Jets%i(65,0.5,7)",k);
	  string DRWJ2 = Form("DeltaR_of_W_Jets%i",k);
	  AnalysisChain.Draw(DRWJ1.c_str(),"(PUR_function(Number_True_Interactions))*weight");
	  DRWjets[k] = (TH1F*)gDirectory->Get(DRWJ2.c_str());
	  gPad->Close();
	  //
	  string HPT1 = Form("Reconstructed_Higgs.Pt() >> HPt%i(40,10,800)",k);
	  string HPT2 = Form("HPt%i",k);
	  AnalysisChain.Draw(HPT1.c_str(),"(PUR_function(Number_True_Interactions))*weight");
	  Hpt[k] = (TH1F*)gDirectory->Get(HPT2.c_str());
	  gPad->Close();
	  //
	  string TPT1 = Form("Reconstructed_Top.Pt() >> TPt%i(40,10,800)",k);
	  string TPT2 = Form("TPt%i",k);
	  AnalysisChain.Draw(TPT1.c_str(),"(PUR_function(Number_True_Interactions))*weight");
	  Tpt[k] = (TH1F*)gDirectory->Get(TPT2.c_str());
	  gPad->Close();
	  //
	  string DRWH1 = Form("DeltaR_of_W_Higgs >> DeltaR_of_W_Higgs%i(65,0.5,7)",k);
	  string DRWH2 = Form("DeltaR_of_W_Higgs%i",k);
	  AnalysisChain.Draw(DRWH1.c_str(),"(PUR_function(Number_True_Interactions))*weight");
	  DRWH[k] = (TH1F*)gDirectory->Get(DRWH2.c_str());
	  gPad->Close();
	  //
	  string DPHJ1 = Form("TMath::Abs(First_Higgs_Jet.Phi()-Second_Higgs_Jet.Phi()) >> DeltaPhi_of_Higgs_jets%i(30,0.0,3.0)",k);
	  string DPHJ2 = Form("DeltaPhi_of_Higgs_jets%i",k);
	  AnalysisChain.Draw(DPHJ1.c_str(),"(PUR_function(Number_True_Interactions))*weight");
	  DPHjets[k] = (TH1F*)gDirectory->Get(DPHJ2.c_str());
	  gPad->Close();
	  //
	  string DPWJ1 = Form("TMath::Abs(First_W_Jet.Phi()-Second_W_Jet.Phi()) >> DeltaPhi_of_W_jets%i(30,0.0,3.0)",k);
	  string DPWJ2 = Form("DeltaPhi_of_W_jets%i",k);
	  AnalysisChain.Draw(DPWJ1.c_str(),"(PUR_function(Number_True_Interactions))*weight");
	  DPWjets[k] = (TH1F*)gDirectory->Get(DPWJ2.c_str());
	  gPad->Close();
	  //
	  string DPTJ1 = Form("TMath::Abs(Top_Jet.Phi()-Reconstructed_W.Phi()) >> DeltaPhi_of_T_jet%i(30,0.0,3.0)",k);
	  string DPTJ2 = Form("DeltaPhi_of_T_jet%i",k);
	  AnalysisChain.Draw(DPTJ1.c_str(),"(PUR_function(Number_True_Interactions))*weight");
	  DPTjets[k] = (TH1F*)gDirectory->Get(DPTJ2.c_str());
	  gPad->Close();
	  //
	  string HM1 = Form("Reconstructed_Higgs.M() >> HM%i(36,60,180)",k);
	  string HM2 = Form("HM%i",k);
	  AnalysisChain.Draw(HM1.c_str(),"(PUR_function(Number_True_Interactions))*weight");
	  HiggsMass[k] = (TH1F*)gDirectory->Get(HM2.c_str());
	  gPad->Close();
	  //
	  string RHT1 = Form("Relative_THT >> RelHT%i(30,0,1)",k);
	  string RHT2 = Form("RelHT%i",k);
	  AnalysisChain.Draw(RHT1.c_str(),"(PUR_function(Number_True_Interactions))*weight");
	  RelHT[k] = (TH1F*)gDirectory->Get(RHT2.c_str());
	  gPad->Close();
	  //
	  string DRTH1 = Form("DeltaR_of_Top_Higgs >> DeltaR_of_Top_Higgs%i(65,0.5,7)",k);
	  string DRTH2 = Form("DeltaR_of_Top_Higgs%i",k);
	  AnalysisChain.Draw(DRTH1.c_str(),"(PUR_function(Number_True_Interactions))*weight");
	  DRTH[k] = (TH1F*)gDirectory->Get(DRTH2.c_str());
	  gPad->Close();
	  //
	  string PTNM1 = Form("PT_Normalized_Mass >> PT_Normalized_Mass%i(70,0.4,5)",k);
	  string PTNM2 = Form("PT_Normalized_Mass%i",k);
	  AnalysisChain.Draw(PTNM1.c_str(),"(PUR_function(Number_True_Interactions))*weight");
	  PtNormalizedMass[k] = (TH1F*)gDirectory->Get(PTNM2.c_str());
	  gPad->Close();
	  //
	  string RM1 = Form("Relative_Mass >> Relative_Mass%i(30,0.0,1)",k);
	  string RM2 = Form("Relative_Mass%i",k);
	  AnalysisChain.Draw(RM1.c_str(),"(PUR_function(Number_True_Interactions))*weight");
	  RelativeMass[k] = (TH1F*)gDirectory->Get(RM2.c_str());
	  gPad->Close();
	  //
	  string MPTNM1 = Form("Mother_PT_Normalized_Mass >> Mother_PT_Normalized_Mass%i(25,0.0,50)",k);
	  string MPTNM2 = Form("Mother_PT_Normalized_Mass%i",k);
	  AnalysisChain.Draw(MPTNM1.c_str(),"(PUR_function(Number_True_Interactions))*weight");
	  MotherPtNormalizedMass[k] = (TH1F*)gDirectory->Get(MPTNM2.c_str());
	  gPad->Close();
	  //
	  string NTops1 = Form("Number_of_Tops >> Number_of_Tops%i(8,0.0,8)",k);
	  string NTops2 = Form("Number_of_Tops%i",k);
	  AnalysisChain.Draw(NTops1.c_str(),"(PUR_function(Number_True_Interactions))*weight");
	  NumberOfTops[k] = (TH1F*)gDirectory->Get(NTops2.c_str());
	  gPad->Close();
	  ///////////NEW VARIABLES//////////////////
	  string HMTM1 = Form("Reconstructed_Higgs.M()/Reconstructed_Top.M() >> HMoverTM%i(30,0,1.)",k);
	  string HMTM2 = Form("HMoverTM%i",k);
	  AnalysisChain.Draw(HMTM1.c_str(),"(PUR_function(Number_True_Interactions))*weight");
	  HiggsMassOverTopMass[k] = (TH1F*)gDirectory->Get(HMTM2.c_str());
	  gPad->Close();
	  //
	  string HTA1 = Form("((Reconstructed_Higgs.Pt()/Reconstructed_Higgs.M())-(Reconstructed_Top.Pt()/Reconstructed_Top.M()))/((Reconstructed_Higgs.Pt()/Reconstructed_Higgs.M())+(Reconstructed_Top.Pt()/Reconstructed_Top.M())) >> HTAsym%i(12,0,1.)",k);
	  string HTA2 = Form("HTAsym%i",k);
	  AnalysisChain.Draw(HTA1.c_str(),"(PUR_function(Number_True_Interactions))*weight");
	  HiggsTopAsymmetry[k] = (TH1F*)gDirectory->Get(HTA2.c_str());
	  gPad->Close();
	  //
	  string TLBT1 = Form("Number_of_Loose_and_non_med_b_tags >> TLBTag%i(10,0,10)",k);
	  string TLBT2 = Form("TLBTag%i",k);
	  AnalysisChain.Draw(TLBT1.c_str(),"(PUR_function(Number_True_Interactions))*weight");
	  ThirdLooseBtag[k] = (TH1F*)gDirectory->Get(TLBT2.c_str());
	  gPad->Close();
	  //
	  string TM1 = Form("Reconstructed_Top.M() >> TMass%i(60,100,700)",k);
	  string TM2 = Form("TMass%i",k);
	  AnalysisChain.Draw(TM1.c_str(),"(PUR_function(Number_True_Interactions))*weight");
	  TopMass[k] = (TH1F*)gDirectory->Get(TM2.c_str());
	  gPad->Close();
	  //
	  string C21 = Form("ChiSquaredSorting >> ChiSq%i(100,0,1000)",k);
	  string C22 = Form("ChiSq%i",k);
	  AnalysisChain.Draw(C21.c_str(),"(PUR_function(Number_True_Interactions))*weight");
	  Chi2[k] = (TH1F*)gDirectory->Get(C22.c_str());
	  gPad->Close();
	  //
	  string UQ1 = Form("U Quark Content >> UQC%i(10,0,10)",k);
	  string UQ2 = Form("UQC%i",k);
	  AnalysisChain.Draw(UQ1.c_str(),"(PUR_function(Number_True_Interactions))*weight");
	  UQuarkContent[k] = (TH1F*)gDirectory->Get(UQ2.c_str());
	  gPad->Close();
	  //
	  string DQ1 = Form("D Quark Content >> DQC%i(10,0,10)",k);
	  string DQ2 = Form("DQC%i",k);
	  AnalysisChain.Draw(DQ1.c_str(),"(PUR_function(Number_True_Interactions))*weight");
	  DQuarkContent[k] = (TH1F*)gDirectory->Get(DQ2.c_str());
	  gPad->Close();
	  //
	  string SQ1 = Form("S Quark Content >> SQC%i(10,0,10)",k);
	  string SQ2 = Form("SQC%i",k);
	  AnalysisChain.Draw(SQ1.c_str(),"(PUR_function(Number_True_Interactions))*weight");
	  SQuarkContent[k] = (TH1F*)gDirectory->Get(SQ2.c_str());
	  gPad->Close();
	  //
	  string CQ1 = Form("C Quark Content >> CQC%i(10,0,10)",k);
	  string CQ2 = Form("CQC%i",k);
	  AnalysisChain.Draw(CQ1.c_str(),"(PUR_function(Number_True_Interactions))*weight");
	  CQuarkContent[k] = (TH1F*)gDirectory->Get(CQ2.c_str());
	  gPad->Close();
	  //
	  string BQ1 = Form("B Quark Content >> BQC%i(10,0,10)",k);
	  string BQ2 = Form("BQC%i",k);
	  AnalysisChain.Draw(BQ1.c_str(),"(PUR_function(Number_True_Interactions))*weight");
	  BQuarkContent[k] = (TH1F*)gDirectory->Get(BQ2.c_str());
	  gPad->Close();
	  //B-tagging Working point	  
	  string BTL1 = Form("Number_CSVLbtagged_jets >> CSVL%i(10,0,10)",k);
	  string BTL2 = Form("CSVL%i",k);
	  AnalysisChain.Draw(BTL1.c_str(),"(PUR_function(Number_True_Interactions))*weight");
	  CSVLB[k] = (TH1F*)gDirectory->Get(BTL2.c_str());
	  gPad->Close();
	  //	  
	  string BTM1 = Form("Number_CSVMbtagged_jets >> CSVM%i(10,0,10)",k);
	  string BTM2 = Form("CSVM%i",k);
	  AnalysisChain.Draw(BTM1.c_str(),"(PUR_function(Number_True_Interactions))*weight");
	  CSVMB[k] = (TH1F*)gDirectory->Get(BTM2.c_str());
	  gPad->Close();
	  //	  
	  string BTT1 = Form("Number_CSVTbtagged_jets >> CSVT%i(10,0,10)",k);
	  string BTT2 = Form("CSVT%i",k);
	  AnalysisChain.Draw(BTT1.c_str(),"(PUR_function(Number_True_Interactions))*weight");
	  CSVTB[k] = (TH1F*)gDirectory->Get(BTT2.c_str());
	  gPad->Close();
	  //	  
	  string VT1 = Form("Vertices >> VTX%i(40,1,41)",k);
	  string VT2 = Form("VTX%i",k);
	  AnalysisChain.Draw(VT1.c_str(),"PUR_function(Number_True_Interactions)*weight*(THT>630)"); //PUR_function(Number_True_Interactions)*(THT>630)");
	  Vtcs[k] = (TH1F*)gDirectory->Get(VT2.c_str());
	  gPad->Close();
	  //
	  string HTCSVM1 = Form("Number_CSVMbtagged_jets : THT >> HT_CSVM%i(65,300,1600,10)",k);
	  //cout << HTCSVM1 << endl;
	  //string HTCSVM1 = Form("Number_CSVMbtagged_jets : Number_CSVLbtagged_jets >> HT_CSVM%i(10,0,10,10,0,10)",k);
	  string HTCSVM2 = Form("HT_CSVM%i",k);
	  //AnalysisChain.Draw("THT:Number_CSVMbtagged_jets");
	  AnalysisChain.Draw(HTCSVM1.c_str(),"(PUR_function(Number_True_Interactions))*weight");
	  HTCSVM[k] = (TH2F*)gDirectory->Get(HTCSVM2.c_str());
	  cout << "Correlation Factor: " << HTCSVM[k]->GetCorrelationFactor() << endl;
	  gPad->Close();
	  //////////////////////////////////////
	  /////////ABCD BKG Estimation//////////
	  //////////////////////////////////////
	  HM1 = Form("Reconstructed_Higgs.M() >> HMBE0%i(36,60,180)",k);
	  HM2 = Form("HMBE0%i",k);
	  AnalysisChain.Draw(HM1.c_str(),"Reconstructed_Higgs.Pt()<200 || Reconstructed_Top.Pt()<300");
	  HiggsMassReversedHptTpt[0][k] = (TH1F*)gDirectory->Get(HM2.c_str());
	  gPad->Close();
	  //	  
	  A1 = Form("Reconstructed_Tprime.M() >> TprimeMassBE0%i(60,400,1600)",k);
	  A2 = Form("TprimeMassBE0%i",k);
	  AnalysisChain.Draw(A1.c_str(),"Reconstructed_Higgs.Pt()<200 || Reconstructed_Top.Pt()<300");
	  TprimeMassReversedHptTpt[0][k] = (TH1F*)gDirectory->Get(A2.c_str());
	  gPad->Close();
	  //	  
	  HM1 = Form("Reconstructed_Higgs.M() >> HMBE1%i(36,60,180)",k);
	  HM2 = Form("HMBE1%i",k);
	  AnalysisChain.Draw(HM1.c_str(),"(Reconstructed_Higgs.Pt()<200 || Reconstructed_Top.Pt()<300) && (DeltaR_of_W_Higgs>=2.2 && DeltaR_of_W_Higgs<=3.5)");
	  HiggsMassReversedHptTpt[1][k] = (TH1F*)gDirectory->Get(HM2.c_str());
	  gPad->Close();
	  //	  
	  A1 = Form("Reconstructed_Tprime.M() >> TprimeMassBE1%i(60,400,1600)",k);
	  A2 = Form("TprimeMassBE1%i",k);
	  AnalysisChain.Draw(A1.c_str(),"(Reconstructed_Higgs.Pt()<200 || Reconstructed_Top.Pt()<300) && (DeltaR_of_W_Higgs>=2.2 && DeltaR_of_W_Higgs<=3.5)");
	  TprimeMassReversedHptTpt[1][k] = (TH1F*)gDirectory->Get(A2.c_str());
	  gPad->Close();
	  //	  
	  HM1 = Form("Reconstructed_Higgs.M() >> HMBE2%i(36,60,180)",k);
	  HM2 = Form("HMBE2%i",k);
	  AnalysisChain.Draw(HM1.c_str(),"(Reconstructed_Higgs.Pt()<200 || Reconstructed_Top.Pt()<300) && (DeltaR_of_W_Higgs>=2.2 && DeltaR_of_W_Higgs<=3.5) && (TMath::Abs(First_Higgs_Jet.Phi()-Second_Higgs_Jet.Phi())<=1.2 && TMath::Abs(Top_Jet.Phi()-Reconstructed_W.Phi())<=1.2)");
	  HiggsMassReversedHptTpt[2][k] = (TH1F*)gDirectory->Get(HM2.c_str());
	  gPad->Close();
	  //	  
	  A1 = Form("Reconstructed_Tprime.M() >> TprimeMassBE2%i(60,400,1600)",k);
	  A2 = Form("TprimeMassBE2%i",k);
	  AnalysisChain.Draw(A1.c_str(),"(Reconstructed_Higgs.Pt()<200 || Reconstructed_Top.Pt()<300) && (DeltaR_of_W_Higgs>=2.2 && DeltaR_of_W_Higgs<=3.5) && (TMath::Abs(First_Higgs_Jet.Phi()-Second_Higgs_Jet.Phi())<=1.2 && TMath::Abs(Top_Jet.Phi()-Reconstructed_W.Phi())<=1.2)");
	  TprimeMassReversedHptTpt[2][k] = (TH1F*)gDirectory->Get(A2.c_str());
	  gPad->Close();
	  //	  
	  HM1 = Form("Reconstructed_Higgs.M() >> HMBE3%i(36,60,180)",k);
	  HM2 = Form("HMBE3%i",k);
	  AnalysisChain.Draw(HM1.c_str(),"(Reconstructed_Higgs.Pt()<200 || Reconstructed_Top.Pt()<300) && (DeltaR_of_W_Higgs>=2.2 && DeltaR_of_W_Higgs<=3.5) && (TMath::Abs(First_Higgs_Jet.Phi()-Second_Higgs_Jet.Phi())<=1.2 && TMath::Abs(Top_Jet.Phi()-Reconstructed_W.Phi())<=1.2) && TMath::Abs(First_W_Jet.Phi()-Second_W_Jet.Phi())<=2.0");
	  HiggsMassReversedHptTpt[3][k] = (TH1F*)gDirectory->Get(HM2.c_str());
	  gPad->Close();
	  //	  
	  A1 = Form("Reconstructed_Tprime.M() >> TprimeMassBE3%i(60,400,1600)",k);
	  A2 = Form("TprimeMassBE3%i",k);
	  AnalysisChain.Draw(A1.c_str(),"(Reconstructed_Higgs.Pt()<200 || Reconstructed_Top.Pt()<300) && (DeltaR_of_W_Higgs>=2.2 && DeltaR_of_W_Higgs<=3.5) && (TMath::Abs(First_Higgs_Jet.Phi()-Second_Higgs_Jet.Phi())<=1.2 && TMath::Abs(Top_Jet.Phi()-Reconstructed_W.Phi())<=1.2) && TMath::Abs(First_W_Jet.Phi()-Second_W_Jet.Phi())<=2.0");
	  TprimeMassReversedHptTpt[3][k] = (TH1F*)gDirectory->Get(A2.c_str());
	  gPad->Close();
	  //	  
	  HM1 = Form("Reconstructed_Higgs.M() >> HMBE4%i(36,60,180)",k);
	  HM2 = Form("HMBE4%i",k);
	  AnalysisChain.Draw(HM1.c_str(),"(Reconstructed_Higgs.Pt()<200 || Reconstructed_Top.Pt()<300) && (DeltaR_of_W_Higgs>=2.2 && DeltaR_of_W_Higgs<=3.5) && (TMath::Abs(First_Higgs_Jet.Phi()-Second_Higgs_Jet.Phi())<=1.2 && TMath::Abs(Top_Jet.Phi()-Reconstructed_W.Phi())<=1.2) && TMath::Abs(First_W_Jet.Phi()-Second_W_Jet.Phi())<=2.0 && Relative_THT>=0.65");
	  HiggsMassReversedHptTpt[4][k] = (TH1F*)gDirectory->Get(HM2.c_str());
	  gPad->Close();
	  //	  
	  A1 = Form("Reconstructed_Tprime.M() >> TprimeMassBE4%i(60,400,1600)",k);
	  A2 = Form("TprimeMassBE4%i",k);
	  AnalysisChain.Draw(A1.c_str(),"(Reconstructed_Higgs.Pt()<200 || Reconstructed_Top.Pt()<300) && (DeltaR_of_W_Higgs>=2.2 && DeltaR_of_W_Higgs<=3.5) && (TMath::Abs(First_Higgs_Jet.Phi()-Second_Higgs_Jet.Phi())<=1.2 && TMath::Abs(Top_Jet.Phi()-Reconstructed_W.Phi())<=1.2) && TMath::Abs(First_W_Jet.Phi()-Second_W_Jet.Phi())<=2.0 && Relative_THT>=0.65");
	  TprimeMassReversedHptTpt[4][k] = (TH1F*)gDirectory->Get(A2.c_str());
	  //	  
	  HM1 = Form("Reconstructed_Higgs.M() >> HMBE5%i(36,60,180)",k);
	  HM2 = Form("HMBE5%i",k);
	  AnalysisChain.Draw(HM1.c_str(),"(Reconstructed_Higgs.Pt()<200 || Reconstructed_Top.Pt()<300) && (DeltaR_of_W_Higgs>=2.2 && DeltaR_of_W_Higgs<=3.5) && (TMath::Abs(First_Higgs_Jet.Phi()-Second_Higgs_Jet.Phi())<=1.2 && TMath::Abs(Top_Jet.Phi()-Reconstructed_W.Phi())<=1.2) && TMath::Abs(First_W_Jet.Phi()-Second_W_Jet.Phi())<=2.0 && Relative_THT>=0.65 && (Relative_Mass>=0.3 && Relative_Mass<=0.5)");
	  HiggsMassReversedHptTpt[5][k] = (TH1F*)gDirectory->Get(HM2.c_str());
	  gPad->Close();
	  //	  
	  A1 = Form("Reconstructed_Tprime.M() >> TprimeMassBE5%i(60,400,1600)",k);
	  A2 = Form("TprimeMassBE5%i",k);
	  AnalysisChain.Draw(A1.c_str(),"(Reconstructed_Higgs.Pt()<200 || Reconstructed_Top.Pt()<300) && (DeltaR_of_W_Higgs>=2.2 && DeltaR_of_W_Higgs<=3.5) && (TMath::Abs(First_Higgs_Jet.Phi()-Second_Higgs_Jet.Phi())<=1.2 && TMath::Abs(Top_Jet.Phi()-Reconstructed_W.Phi())<=1.2) && TMath::Abs(First_W_Jet.Phi()-Second_W_Jet.Phi())<=2.0 && Relative_THT>=0.65 && (Relative_Mass>=0.3 && Relative_Mass<=0.5)");
	  TprimeMassReversedHptTpt[5][k] = (TH1F*)gDirectory->Get(A2.c_str());
	  gPad->Close();
	  //
	  string HPTTPT1 = Form("Reconstructed_Higgs.Pt():Reconstructed_Top.Pt() >> HPtTPt0%i(40,10,800,40,10,800)",k);
	  string HPTTPT2 = Form("HPtTPt0%i",k);
	  AnalysisChain.Draw(HPTTPT1.c_str());
	  HptTpt[0][k] = (TH2F*)gDirectory->Get(HPTTPT2.c_str());
	  cout << "Correlation Factor: " << HptTpt[0][k]->GetCorrelationFactor() << endl;
	  gPad->Close();
	  //
	  HPTTPT1 = Form("Reconstructed_Higgs.Pt():Reconstructed_Top.Pt() >> HPtTPt1%i(40,10,800,40,10,800)",k);
	  HPTTPT2 = Form("HPtTPt1%i",k);
	  AnalysisChain.Draw(HPTTPT1.c_str(),"(DeltaR_of_W_Higgs>=2.2 && DeltaR_of_W_Higgs<=3.5)");
	  HptTpt[1][k] = (TH2F*)gDirectory->Get(HPTTPT2.c_str());
	  cout << "Correlation Factor: " << HptTpt[1][k]->GetCorrelationFactor() << endl;
	  gPad->Close();
	  //
	  HPTTPT1 = Form("Reconstructed_Higgs.Pt():Reconstructed_Top.Pt() >> HPtTPt2%i(40,10,800,40,10,800)",k);
	  HPTTPT2 = Form("HPtTPt2%i",k);
	  AnalysisChain.Draw(HPTTPT1.c_str(),"(DeltaR_of_W_Higgs>=2.2 && DeltaR_of_W_Higgs<=3.5) && (TMath::Abs(First_Higgs_Jet.Phi()-Second_Higgs_Jet.Phi())<=1.2 && TMath::Abs(Top_Jet.Phi()-Reconstructed_W.Phi())<=1.2)");
	  HptTpt[2][k] = (TH2F*)gDirectory->Get(HPTTPT2.c_str());
	  cout << "Correlation Factor: " << HptTpt[2][k]->GetCorrelationFactor() << endl;
	  gPad->Close();
	  //
	  HPTTPT1 = Form("Reconstructed_Higgs.Pt():Reconstructed_Top.Pt() >> HPtTPt3%i(40,10,800,40,10,800)",k);
	  HPTTPT2 = Form("HPtTPt3%i",k);
	  AnalysisChain.Draw(HPTTPT1.c_str(),"(DeltaR_of_W_Higgs>=2.2 && DeltaR_of_W_Higgs<=3.5) && (TMath::Abs(First_Higgs_Jet.Phi()-Second_Higgs_Jet.Phi())<=1.2 && TMath::Abs(Top_Jet.Phi()-Reconstructed_W.Phi())<=1.2) && TMath::Abs(First_W_Jet.Phi()-Second_W_Jet.Phi())<=2.0");
	  HptTpt[3][k] = (TH2F*)gDirectory->Get(HPTTPT2.c_str());
	  cout << "Correlation Factor: " << HptTpt[3][k]->GetCorrelationFactor() << endl;
	  gPad->Close();
	  //
	  HPTTPT1 = Form("Reconstructed_Higgs.Pt():Reconstructed_Top.Pt() >> HPtTPt4%i(40,10,800,40,10,800)",k);
	  HPTTPT2 = Form("HPtTPt4%i",k);
	  AnalysisChain.Draw(HPTTPT1.c_str(),"(DeltaR_of_W_Higgs>=2.2 && DeltaR_of_W_Higgs<=3.5) && (TMath::Abs(First_Higgs_Jet.Phi()-Second_Higgs_Jet.Phi())<=1.2 && TMath::Abs(Top_Jet.Phi()-Reconstructed_W.Phi())<=1.2) && TMath::Abs(First_W_Jet.Phi()-Second_W_Jet.Phi())<=2.0 && Relative_THT>=0.65");
	  HptTpt[4][k] = (TH2F*)gDirectory->Get(HPTTPT2.c_str());
	  cout << "Correlation Factor: " << HptTpt[4][k]->GetCorrelationFactor() << endl;
	  gPad->Close();
	  //
	  HPTTPT1 = Form("Reconstructed_Higgs.Pt():Reconstructed_Top.Pt() >> HPtTPt5%i(40,10,800,40,10,800)",k);
	  HPTTPT2 = Form("HPtTPt5%i",k);
	  AnalysisChain.Draw(HPTTPT1.c_str(),"(DeltaR_of_W_Higgs>=2.2 && DeltaR_of_W_Higgs<=3.5) && (TMath::Abs(First_Higgs_Jet.Phi()-Second_Higgs_Jet.Phi())<=1.2 && TMath::Abs(Top_Jet.Phi()-Reconstructed_W.Phi())<=1.2) && TMath::Abs(First_W_Jet.Phi()-Second_W_Jet.Phi())<=2.0 && Relative_THT>=0.65 && (Relative_Mass>=0.3 && Relative_Mass<=0.5)");
	  HptTpt[5][k] = (TH2F*)gDirectory->Get(HPTTPT2.c_str());
	  cout << "Correlation Factor: " << HptTpt[5][k]->GetCorrelationFactor() << endl;
	  gPad->Close();
	  HM1 = Form("Reconstructed_Higgs.M() >> HM0%i(36,60,180)",k);
	  HM2 = Form("HM0%i",k);
	  AnalysisChain.Draw(HM1.c_str(),"Reconstructed_Higgs.Pt()>200 && Reconstructed_Top.Pt()>200");
	  HiggsMassCuts[0][k] = (TH1F*)gDirectory->Get(HM2.c_str());
	  gPad->Close();
	  //	  
	  HM1 = Form("Reconstructed_Higgs.M() >> HM1%i(36,60,180)",k);
	  HM2 = Form("HM1%i",k);
	  AnalysisChain.Draw(HM1.c_str(),"(Reconstructed_Higgs.Pt()>200 && Reconstructed_Top.Pt()>200) && (DeltaR_of_W_Higgs>=2.2 && DeltaR_of_W_Higgs<=3.5)");
	  HiggsMassCuts[1][k] = (TH1F*)gDirectory->Get(HM2.c_str());
	  gPad->Close();
	  //	  
	  HM1 = Form("Reconstructed_Higgs.M() >> HM2%i(36,60,180)",k);
	  HM2 = Form("HM2%i",k);
	  AnalysisChain.Draw(HM1.c_str(),"(Reconstructed_Higgs.Pt()>200 && Reconstructed_Top.Pt()>200) && (DeltaR_of_W_Higgs>=2.2 && DeltaR_of_W_Higgs<=3.5) && (TMath::Abs(First_Higgs_Jet.Phi()-Second_Higgs_Jet.Phi())<=1.2 && TMath::Abs(Top_Jet.Phi()-Reconstructed_W.Phi())<=1.2)");
	  HiggsMassCuts[2][k] = (TH1F*)gDirectory->Get(HM2.c_str());
	  gPad->Close();
	  //	  
	  HM1 = Form("Reconstructed_Higgs.M() >> HM3%i(36,60,180)",k);
	  HM2 = Form("HM3%i",k);
	  AnalysisChain.Draw(HM1.c_str(),"(Reconstructed_Higgs.Pt()>200 && Reconstructed_Top.Pt()>200) && (DeltaR_of_W_Higgs>=2.2 && DeltaR_of_W_Higgs<=3.5) && (TMath::Abs(First_Higgs_Jet.Phi()-Second_Higgs_Jet.Phi())<=1.2 && TMath::Abs(Top_Jet.Phi()-Reconstructed_W.Phi())<=1.2) && TMath::Abs(First_W_Jet.Phi()-Second_W_Jet.Phi())<=2.0");
	  HiggsMassCuts[3][k] = (TH1F*)gDirectory->Get(HM2.c_str());
	  gPad->Close();
	  //
	  HM1 = Form("Reconstructed_Higgs.M() >> HM4%i(36,60,180)",k);
	  HM2 = Form("HM4%i",k);
	  AnalysisChain.Draw(HM1.c_str(),"(Reconstructed_Higgs.Pt()>200 && Reconstructed_Top.Pt()>200) && (DeltaR_of_W_Higgs>=2.2 && DeltaR_of_W_Higgs<=3.5) && (TMath::Abs(First_Higgs_Jet.Phi()-Second_Higgs_Jet.Phi())<=1.2 && TMath::Abs(Top_Jet.Phi()-Reconstructed_W.Phi())<=1.2) && TMath::Abs(First_W_Jet.Phi()-Second_W_Jet.Phi())<=2.0 && Relative_THT>=0.65");
	  HiggsMassCuts[4][k] = (TH1F*)gDirectory->Get(HM2.c_str());
	  gPad->Close();
	  //
	  HM1 = Form("Reconstructed_Higgs.M() >> HM5%i(36,60,180)",k);
	  HM2 = Form("HM5%i",k);
	  AnalysisChain.Draw(HM1.c_str(),"(Reconstructed_Higgs.Pt()>200 && Reconstructed_Top.Pt()>200) && (DeltaR_of_W_Higgs>=2.2 && DeltaR_of_W_Higgs<=3.5) && (TMath::Abs(First_Higgs_Jet.Phi()-Second_Higgs_Jet.Phi())<=1.2 && TMath::Abs(Top_Jet.Phi()-Reconstructed_W.Phi())<=1.2) && TMath::Abs(First_W_Jet.Phi()-Second_W_Jet.Phi())<=2.0 && Relative_THT>=0.65 && (Relative_Mass>=0.3 && Relative_Mass<=0.5)");
	  HiggsMassCuts[5][k] = (TH1F*)gDirectory->Get(HM2.c_str());
	  gPad->Close();
	  /*///////////////////////
	  //Sideband estimation//
	  ///////////////////////
	  string WFH1 = Form("W_From_Matching.M() >> WMFH0%i(44,60.0,500)",k);
	  string WFH2 = Form("WMFH0%i",k);
	  AnalysisChain.Draw(WFH1.c_str(),"Reconstructed_Higgs.Pt()>200 && Reconstructed_Top.Pt()>200");
	  WMassFromHiggs[0][k] = (TH1F*)gDirectory->Get(WFH2.c_str());
	  gPad->Close();
	  //
	  WFH1 = Form("W_From_Matching.M() >> WMFH1%i(44,60.0,500)",k);
	  WFH2 = Form("WMFH1%i",k);
	  AnalysisChain.Draw(WFH1.c_str(),"(Reconstructed_Higgs.Pt()>200 && Reconstructed_Top.Pt()>200) && (DeltaR_of_W_Higgs>=2.2 && DeltaR_of_W_Higgs<=3.5)");
	  WMassFromHiggs[1][k] = (TH1F*)gDirectory->Get(WFH2.c_str());
	  gPad->Close();
	  //
	  WFH1 = Form("W_From_Matching.M() >> WMFH2%i(44,60.0,500)",k);
	  WFH2 = Form("WMFH2%i",k);
	  AnalysisChain.Draw(WFH1.c_str(),"(Reconstructed_Higgs.Pt()>200 && Reconstructed_Top.Pt()>200) && (DeltaR_of_W_Higgs>=2.2 && DeltaR_of_W_Higgs<=3.5) && (TMath::Abs(First_Higgs_Jet.Phi()-Second_Higgs_Jet.Phi())<=1.2 && TMath::Abs(Top_Jet.Phi()-Reconstructed_W.Phi())<=1.2)");
	  WMassFromHiggs[2][k] = (TH1F*)gDirectory->Get(WFH2.c_str());
	  gPad->Close();
	  //
	  WFH1 = Form("W_From_Matching.M() >> WMFH3%i(44,60.0,500)",k);
	  WFH2 = Form("WMFH3%i",k);
	  AnalysisChain.Draw(WFH1.c_str(),"(Reconstructed_Higgs.Pt()>200 && Reconstructed_Top.Pt()>200) && (DeltaR_of_W_Higgs>=2.2 && DeltaR_of_W_Higgs<=3.5) && (TMath::Abs(First_Higgs_Jet.Phi()-Second_Higgs_Jet.Phi())<=1.2 && TMath::Abs(Top_Jet.Phi()-Reconstructed_W.Phi())<=1.2) && TMath::Abs(First_W_Jet.Phi()-Second_W_Jet.Phi())<=2.0");
	  WMassFromHiggs[3][k] = (TH1F*)gDirectory->Get(WFH2.c_str());
	  gPad->Close();
	  //
	  WFH1 = Form("W_From_Matching.M() >> WMFH4%i(44,60.0,500)",k);
	  WFH2 = Form("WMFH4%i",k);
	  AnalysisChain.Draw(WFH1.c_str(),"(Reconstructed_Higgs.Pt()>200 && Reconstructed_Top.Pt()>200) && (DeltaR_of_W_Higgs>=2.2 && DeltaR_of_W_Higgs<=3.5) && (TMath::Abs(First_Higgs_Jet.Phi()-Second_Higgs_Jet.Phi())<=1.2 && TMath::Abs(Top_Jet.Phi()-Reconstructed_W.Phi())<=1.2) && TMath::Abs(First_W_Jet.Phi()-Second_W_Jet.Phi())<=2.0 && Relative_THT>=0.65");
	  WMassFromHiggs[4][k] = (TH1F*)gDirectory->Get(WFH2.c_str());
	  gPad->Close();
	  //
	  WFH1 = Form("W_From_Matching.M() >> WMFH5%i(44,60.0,500)",k);
	  WFH2 = Form("WMFH5%i",k);
	  AnalysisChain.Draw(WFH1.c_str(),"(Reconstructed_Higgs.Pt()>200 && Reconstructed_Top.Pt()>200) && (DeltaR_of_W_Higgs>=2.2 && DeltaR_of_W_Higgs<=3.5) && (TMath::Abs(First_Higgs_Jet.Phi()-Second_Higgs_Jet.Phi())<=1.2 && TMath::Abs(Top_Jet.Phi()-Reconstructed_W.Phi())<=1.2) && TMath::Abs(First_W_Jet.Phi()-Second_W_Jet.Phi())<=2.0 && Relative_THT>=0.65 && (Relative_Mass>=0.3 && Relative_Mass<=0.5)");
	  WMassFromHiggs[5][k] = (TH1F*)gDirectory->Get(WFH2.c_str());
	  gPad->Close();*/
          ////////
          //////////
          A1 = Form("Top_From_Higgs_Chi2.M() >> TopFromHiggsChi2Mass%i(95,50,1000)",k);
          A2 = Form("TopFromHiggsChi2Mass%i",k);
          AnalysisChain.Draw(A1.c_str(),"(PUR_function(Number_True_Interactions))*weight");
          TopMassFromHiggs[k] = (TH1F*)gDirectory->Get(A2.c_str());
          gPad->Close();
	  //////////////////////////////////////////////////////////////////////////////
	  //////////////////////////////////////////////////////////////////////////////
	  //////////////////////////////////////////////////////////////////////////////
	  A1 = Form("Reconstructed_Tprime3L.M() >> TprimeMass3L%i(60,400,1600)",k);
	  A2 = Form("TprimeMass3L%i",k);
	  AnalysisChain3L.Draw(A1.c_str()); //,"(Reconstructed_Top3L.M()>=230)");
	  FiveJetsMass3L[k] = (TH1F*)gDirectory->Get(A2.c_str());
	  gPad->Close();
	  //
	  /*A1 = Form("Reconstructed_Tprime3L.M() >> TprimeMass3L_L%i(60,400,1600)",k);
	  A2 = Form("TprimeMass3L_L%i",k);
	  AnalysisChain3L.Draw(A1.c_str(),"(Reconstructed_Top3L.M()<230)");
	  FiveJetsMass3L_L[k] = (TH1F*)gDirectory->Get(A2.c_str());
	  gPad->Close();
	  //	
	  A1 = Form("Reconstructed_Tprime3T.M() >> TprimeMass3T%i(60,400,1600)",k);
	  A2 = Form("TprimeMass3T%i",k);
	  AnalysisChain3T.Draw(A1.c_str()); //,"(Reconstructed_Top3T.M()>=230)");
	  FiveJetsMass3T[k] = (TH1F*)gDirectory->Get(A2.c_str());
	  gPad->Close();  
	  //	
	  A1 = Form("Reconstructed_Tprime3T.M() >> TprimeMass3T_L%i(60,400,1600)",k);
	  A2 = Form("TprimeMass3T_L%i",k);
	  AnalysisChain3T.Draw(A1.c_str()); //,"(Reconstructed_Top3T.M()<230)");
	  FiveJetsMass3T_L[k] = (TH1F*)gDirectory->Get(A2.c_str());
	  gPad->Close();  */
	  /////////////////////
	  ///////////////////// ADDTIONAL INFO ON JETS
	  /////////////////////
	  A1 = Form("JET1.Eta() >> jet1_eta%i(100,-5,5)",k);
	  A2 = Form("jet1_eta%i",k);
	  AnalysisChain.Draw(A1.c_str(),finalcut);
	  JET1ETA[k] = (TH1F*)gDirectory->Get(A2.c_str());
	  gPad->Close();
	  //
	  A1 = Form("JET2.Eta() >> jet2_eta%i(100,-5,5)",k);
	  A2 = Form("jet2_eta%i",k);
	  AnalysisChain.Draw(A1.c_str(),finalcut);
	  JET2ETA[k] = (TH1F*)gDirectory->Get(A2.c_str());
	  gPad->Close();
	  //
	  A1 = Form("JET3.Eta() >> jet3_eta%i(100,-5,5)",k);
	  A2 = Form("jet3_eta%i",k);
	  AnalysisChain.Draw(A1.c_str(),finalcut);
	  JET3ETA[k] = (TH1F*)gDirectory->Get(A2.c_str());
	  gPad->Close();
	  //
	  A1 = Form("JET4.Eta() >> jet4_eta%i(100,-5,5)",k);
	  A2 = Form("jet4_eta%i",k);
	  AnalysisChain.Draw(A1.c_str(),finalcut);
	  JET4ETA[k] = (TH1F*)gDirectory->Get(A2.c_str());
	  gPad->Close();
	  //
	  A1 = Form("JET5.Eta() >> jet5_eta%i(100,-5,5)",k);
	  A2 = Form("jet5_eta%i",k);
	  AnalysisChain.Draw(A1.c_str(),finalcut);
	  JET5ETA[k] = (TH1F*)gDirectory->Get(A2.c_str());
	  gPad->Close();
	  //
	  A1 = Form("JET6.Eta() >> jet6_eta%i(100,-5,5)",k);
	  A2 = Form("jet6_eta%i",k);
	  AnalysisChain.Draw(A1.c_str(),finalcut);
	  JET6ETA[k] = (TH1F*)gDirectory->Get(A2.c_str());
	  gPad->Close();
	  //
	  A1 = Form("JET1.Phi() >> jet1_phi%i(60,-3,3)",k);
	  A2 = Form("jet1_phi%i",k);
	  AnalysisChain.Draw(A1.c_str(),finalcut);
	  JET1PHI[k] = (TH1F*)gDirectory->Get(A2.c_str());
	  gPad->Close();
	  //
	  A1 = Form("JET2.Phi() >> jet2_phi%i(60,-3,3)",k);
	  A2 = Form("jet2_phi%i",k);
	  AnalysisChain.Draw(A1.c_str(),finalcut);
	  JET2PHI[k] = (TH1F*)gDirectory->Get(A2.c_str());
	  gPad->Close();
	  //
	  A1 = Form("JET3.Phi() >> jet3_phi%i(60,-3,3)",k);
	  A2 = Form("jet3_phi%i",k);
	  AnalysisChain.Draw(A1.c_str(),finalcut);
	  JET3PHI[k] = (TH1F*)gDirectory->Get(A2.c_str());
	  gPad->Close();
	  //
	  A1 = Form("JET4.Phi() >> jet4_phi%i(60,-3,3)",k);
	  A2 = Form("jet4_phi%i",k);
	  AnalysisChain.Draw(A1.c_str(),finalcut);
	  JET4PHI[k] = (TH1F*)gDirectory->Get(A2.c_str());
	  gPad->Close();
	  //
	  A1 = Form("JET5.Phi() >> jet5_phi%i(60,-3,3)",k);
	  A2 = Form("jet5_phi%i",k);
	  AnalysisChain.Draw(A1.c_str(),finalcut);
	  JET5PHI[k] = (TH1F*)gDirectory->Get(A2.c_str());
	  gPad->Close();
	  //
	  A1 = Form("JET6.Phi() >> jet6_phi%i(60,-3,3)",k);
	  A2 = Form("jet6_phi%i",k);
	  AnalysisChain.Draw(A1.c_str(),finalcut);
	  JET6PHI[k] = (TH1F*)gDirectory->Get(A2.c_str());
	  gPad->Close();
	  //
	  A1 = Form("Number_Jets >> Num_jets%i(20,0,20)",k);
	  A2 = Form("Num_jets%i",k);
	  AnalysisChain.Draw(A1.c_str(),finalcut);
	  JETMULTI[k] = (TH1F*)gDirectory->Get(A2.c_str());
	  gPad->Close();
	}

    }

  //TprimeMass histo recollection
  string A=Form("TprimeMass%i",FirstGoodMarker);
  TH1F *SingleTopTprimeMass=(TH1F*)gDirectory->Get(A.c_str());  
  //LeadingJetPT histo recollection
  A=Form("jet1_pt%i",FirstGoodMarker);
  TH1F *SingleTopTprimeMassLJPT=(TH1F*)gDirectory->Get(A.c_str());  
  //LeadingJetPT histo recollection
  A=Form("jet2_pt%i",FirstGoodMarker);
  TH1F *SingleTopTprimeMassL2JPT=(TH1F*)gDirectory->Get(A.c_str());
  //LeadingJetPT histo recollection
  A=Form("jet3_pt%i",FirstGoodMarker);
  TH1F *SingleTopTprimeMassL3JPT=(TH1F*)gDirectory->Get(A.c_str());
  //LeadingJetPT histo recollection
  A=Form("jet4_pt%i",FirstGoodMarker);
  TH1F *SingleTopTprimeMassL4JPT=(TH1F*)gDirectory->Get(A.c_str());
  //LeadingJetPT histo recollection
  A=Form("jet5_pt%i",FirstGoodMarker);
  TH1F *SingleTopTprimeMassL5JPT=(TH1F*)gDirectory->Get(A.c_str());
  //LeadingJetPT histo recollection
  A=Form("jet6_pt%i",FirstGoodMarker);
  TH1F *SingleTopTprimeMassL6JPT=(TH1F*)gDirectory->Get(A.c_str());
  //THT histo recollection
  A=Form("THT%i",FirstGoodMarker);
  TH1F *SingleTopTprimeMassTHT=(TH1F*)gDirectory->Get(A.c_str());  
  //DR Higgs jets histo recollection
  A=Form("DeltaR_of_Higgs_Jets%i",FirstGoodMarker);
  TH1F *SingleTopTprimeMassDRHJ=(TH1F*)gDirectory->Get(A.c_str());
  //DR W jets histo recollection
  A=Form("DeltaR_of_W_Jets%i",FirstGoodMarker);
  TH1F *SingleTopTprimeMassDRWJ=(TH1F*)gDirectory->Get(A.c_str());
  //Higgs PT histo recollection
  A=Form("HPt%i",FirstGoodMarker);
  TH1F *SingleTopTprimeMassHPT=(TH1F*)gDirectory->Get(A.c_str());
  //Top PT histo recollection
  A=Form("TPt%i",FirstGoodMarker);
  TH1F *SingleTopTprimeMassTPT=(TH1F*)gDirectory->Get(A.c_str());
  //DR W Higgs histo recollection
  A=Form("DeltaR_of_W_Higgs%i",FirstGoodMarker);
  TH1F *SingleTopTprimeMassDRWH=(TH1F*)gDirectory->Get(A.c_str());
  //DP Higgs jets histo recollection
  A=Form("DeltaPhi_of_Higgs_jets%i",FirstGoodMarker);
  TH1F *SingleTopTprimeMassDPHJ=(TH1F*)gDirectory->Get(A.c_str());
  //DP W jets histo recollection
  A=Form("DeltaPhi_of_W_jets%i",FirstGoodMarker);
  TH1F *SingleTopTprimeMassDPWJ=(TH1F*)gDirectory->Get(A.c_str());
  //DP Top jet histo recollection
  A=Form("DeltaPhi_of_T_jet%i",FirstGoodMarker);
  TH1F *SingleTopTprimeMassDPTJ=(TH1F*)gDirectory->Get(A.c_str());
  //Higgs Mass histo recollection
  A=Form("HM%i",FirstGoodMarker);
  TH1F *SingleTopTprimeMassHM=(TH1F*)gDirectory->Get(A.c_str());
  //Relative HT histo recollection
  A=Form("RelHT%i",FirstGoodMarker);
  TH1F *SingleTopTprimeMassRelHT=(TH1F*)gDirectory->Get(A.c_str());
  //DR Top Higgs histo recollection
  A=Form("DeltaR_of_Top_Higgs%i",FirstGoodMarker);
  TH1F *SingleTopTprimeMassDRTH=(TH1F*)gDirectory->Get(A.c_str());
  //PT Normalized Mass histo recollection
  A=Form("PT_Normalized_Mass%i",FirstGoodMarker);
  TH1F *SingleTopTprimeMassPTNM=(TH1F*)gDirectory->Get(A.c_str());
  //Relative Mass histo recollection
  A=Form("Relative_Mass%i",FirstGoodMarker);
  TH1F *SingleTopTprimeMassRM=(TH1F*)gDirectory->Get(A.c_str());
  //Mother PT Normalized Mass histo recollection
  A=Form("Mother_PT_Normalized_Mass%i",FirstGoodMarker);
  TH1F *SingleTopTprimeMassMPTNM=(TH1F*)gDirectory->Get(A.c_str());
  //Number of Tops histo recollection
  A=Form("Number_of_Tops%i",FirstGoodMarker);
  TH1F *SingleTopTprimeMassNT=(TH1F*)gDirectory->Get(A.c_str());
  //Higgs Mass over Top Mass histo recollection
  A=Form("HMoverTM%i",FirstGoodMarker);
  TH1F *SingleTopTprimeMassHMTM=(TH1F*)gDirectory->Get(A.c_str());
  //Higgs Top Asym histo recollection
  A=Form("HTAsym%i",FirstGoodMarker);
  TH1F *SingleTopTprimeMassHTAsym=(TH1F*)gDirectory->Get(A.c_str());
  //Number of loose and non medium b-tag histo recollection
  A=Form("TLBTag%i",FirstGoodMarker);
  TH1F *SingleTopTprimeMassTLBT=(TH1F*)gDirectory->Get(A.c_str());
  //Top Mass histo recollection
  A=Form("TMass%i",FirstGoodMarker);
  TH1F *SingleTopTprimeMassTM=(TH1F*)gDirectory->Get(A.c_str());
  //Chi2 histo recollection
  A=Form("ChiSq%i",FirstGoodMarker);
  TH1F *SingleTopTprimeMassChi2=(TH1F*)gDirectory->Get(A.c_str());
  //U quark content histo recollection
  A=Form("UQC%i",FirstGoodMarker);
  TH1F *SingleTopTprimeMassUQC=(TH1F*)gDirectory->Get(A.c_str());
  //D quark content histo recollection
  A=Form("DQC%i",FirstGoodMarker);
  TH1F *SingleTopTprimeMassDQC=(TH1F*)gDirectory->Get(A.c_str());
  //S quark content histo recollection
  A=Form("SQC%i",FirstGoodMarker);
  TH1F *SingleTopTprimeMassSQC=(TH1F*)gDirectory->Get(A.c_str());
  //C quark content histo recollection
  A=Form("CQC%i",FirstGoodMarker);
  TH1F *SingleTopTprimeMassCQC=(TH1F*)gDirectory->Get(A.c_str());
  //B quark content histo recollection
  A=Form("BQC%i",FirstGoodMarker);
  TH1F *SingleTopTprimeMassBQC=(TH1F*)gDirectory->Get(A.c_str());
  //CSVL tagging
  A=Form("CSVL%i",FirstGoodMarker);
  TH1F *SingleTopTprimeMassCSVL=(TH1F*)gDirectory->Get(A.c_str());
  //CSVM tagging
  A=Form("CSVM%i",FirstGoodMarker);
  TH1F *SingleTopTprimeMassCSVM=(TH1F*)gDirectory->Get(A.c_str());
  //CSVT tagging
  A=Form("CSVT%i",FirstGoodMarker);
  TH1F *SingleTopTprimeMassCSVT=(TH1F*)gDirectory->Get(A.c_str());
  //Vertices
  A=Form("VTX%i",FirstGoodMarker);
  TH1F *SingleTopTprimeMassVTX=(TH1F*)gDirectory->Get(A.c_str());
  ///////////BKG Estimation///////////////
  TH1F *QCDTprimeMassHMBE[6];
  TH1F *QCDTprimeMassBE[6];
  //Higgs Mass histo recollection
  A=Form("HMBE0%i",FirstGoodMarker);
  QCDTprimeMassHMBE[0]=(TH1F*)gDirectory->Get(A.c_str());
  //TprimeMass histo recollection
  A=Form("TprimeMassBE0%i",FirstGoodMarker);
  QCDTprimeMassBE[0]=(TH1F*)gDirectory->Get(A.c_str());
  //
  A=Form("HMBE1%i",FirstGoodMarker);
  QCDTprimeMassHMBE[1]=(TH1F*)gDirectory->Get(A.c_str());
  //
  A=Form("TprimeMassBE1%i",FirstGoodMarker);
  QCDTprimeMassBE[1]=(TH1F*)gDirectory->Get(A.c_str());
  //
  A=Form("HMBE2%i",FirstGoodMarker);
  QCDTprimeMassHMBE[2]=(TH1F*)gDirectory->Get(A.c_str());
  //
  A=Form("TprimeMassBE2%i",FirstGoodMarker);
  QCDTprimeMassBE[2]=(TH1F*)gDirectory->Get(A.c_str());
  //
  A=Form("HMBE3%i",FirstGoodMarker);
  QCDTprimeMassHMBE[3]=(TH1F*)gDirectory->Get(A.c_str());
  //
  A=Form("TprimeMassBE3%i",FirstGoodMarker);
  QCDTprimeMassBE[3]=(TH1F*)gDirectory->Get(A.c_str());
  //
  A=Form("HMBE4%i",FirstGoodMarker);
  QCDTprimeMassHMBE[4]=(TH1F*)gDirectory->Get(A.c_str());
  //
  A=Form("TprimeMassBE4%i",FirstGoodMarker);
  QCDTprimeMassBE[4]=(TH1F*)gDirectory->Get(A.c_str());
  //
  A=Form("HMBE5%i",FirstGoodMarker);
  QCDTprimeMassHMBE[5]=(TH1F*)gDirectory->Get(A.c_str());
  //
  A=Form("TprimeMassBE5%i",FirstGoodMarker);
  QCDTprimeMassBE[5]=(TH1F*)gDirectory->Get(A.c_str());
  //HTCSVM
  A=Form("HT_CSVM%i",FirstGoodMarker);
  TH2F *SingleTop_HTCSVM=(TH2F*)gDirectory->Get(A.c_str());
  //HptToppt
  TH2F *SingleTop_HptTpt[6];
  A=Form("HPtTPt0%i",FirstGoodMarker);
  SingleTop_HptTpt[0]=(TH2F*)gDirectory->Get(A.c_str());
  //HptToppt
  A=Form("HPtTPt1%i",FirstGoodMarker);
  SingleTop_HptTpt[1]=(TH2F*)gDirectory->Get(A.c_str());
  //HptToppt
  A=Form("HPtTPt2%i",FirstGoodMarker);
  SingleTop_HptTpt[2]=(TH2F*)gDirectory->Get(A.c_str());
  //HptToppt
  A=Form("HPtTPt3%i",FirstGoodMarker);
  SingleTop_HptTpt[3]=(TH2F*)gDirectory->Get(A.c_str());
  //HptToppt
  A=Form("HPtTPt4%i",FirstGoodMarker);
  SingleTop_HptTpt[4]=(TH2F*)gDirectory->Get(A.c_str());
  //HptToppt
  A=Form("HPtTPt5%i",FirstGoodMarker);
  SingleTop_HptTpt[5]=(TH2F*)gDirectory->Get(A.c_str());
  //Higgs Mass histo recollection
  TH1F *QCDTprimeMassHMCuts[6];
  A=Form("HM0%i",FirstGoodMarker);
  QCDTprimeMassHMCuts[0]=(TH1F*)gDirectory->Get(A.c_str());
  //
  A=Form("HM1%i",FirstGoodMarker);
  QCDTprimeMassHMCuts[1]=(TH1F*)gDirectory->Get(A.c_str());
  //
  A=Form("HM2%i",FirstGoodMarker);
  QCDTprimeMassHMCuts[2]=(TH1F*)gDirectory->Get(A.c_str());
  //
  A=Form("HM3%i",FirstGoodMarker);
  QCDTprimeMassHMCuts[3]=(TH1F*)gDirectory->Get(A.c_str());
  //
  A=Form("HM4%i",FirstGoodMarker);
  QCDTprimeMassHMCuts[4]=(TH1F*)gDirectory->Get(A.c_str());
  //
  A=Form("HM5%i",FirstGoodMarker);
  QCDTprimeMassHMCuts[5]=(TH1F*)gDirectory->Get(A.c_str()); 
  /*//W Mass from Higgs histo recollection
  TH1F *QCDTprimeMassWMFHCuts[6];
  A=Form("WMFH0%i",FirstGoodMarker);
  QCDTprimeMassWMFHCuts[0]=(TH1F*)gDirectory->Get(A.c_str()); 
  //
  A=Form("WMFH1%i",FirstGoodMarker);
  QCDTprimeMassWMFHCuts[1]=(TH1F*)gDirectory->Get(A.c_str()); 
  //
  A=Form("WMFH2%i",FirstGoodMarker);
  QCDTprimeMassWMFHCuts[2]=(TH1F*)gDirectory->Get(A.c_str()); 
  //
  A=Form("WMFH3%i",FirstGoodMarker);
  QCDTprimeMassWMFHCuts[3]=(TH1F*)gDirectory->Get(A.c_str()); 
  //
  A=Form("WMFH4%i",FirstGoodMarker);
  QCDTprimeMassWMFHCuts[4]=(TH1F*)gDirectory->Get(A.c_str()); 
  //
  A=Form("WMFH5%i",FirstGoodMarker);
  QCDTprimeMassWMFHCuts[5]=(TH1F*)gDirectory->Get(A.c_str()); */
  //   
  //TopMassFromHiggs histo recollection
  A=Form("TopFromHiggsChi2Mass%i",FirstGoodMarker);
  TH1F *QCDTopMassFromHiggs=(TH1F*)gDirectory->Get(A.c_str()); 
  //
  A=Form("jet1_eta%i",FirstGoodMarker);
  TH1F *QCDJET1ETA=(TH1F*)gDirectory->Get(A.c_str());
  QCDJET1ETA->Scale(Lumi*XS[FirstGoodMarker]*0.5/ProcessedEvents[FirstGoodMarker]);
  A=Form("jet2_eta%i",FirstGoodMarker);
  TH1F *QCDJET2ETA=(TH1F*)gDirectory->Get(A.c_str());
  QCDJET2ETA->Scale(Lumi*XS[FirstGoodMarker]*0.5/ProcessedEvents[FirstGoodMarker]);
  A=Form("jet3_eta%i",FirstGoodMarker);
  TH1F *QCDJET3ETA=(TH1F*)gDirectory->Get(A.c_str());
  QCDJET3ETA->Scale(Lumi*XS[FirstGoodMarker]*0.5/ProcessedEvents[FirstGoodMarker]);
  A=Form("jet4_eta%i",FirstGoodMarker);
  TH1F *QCDJET4ETA=(TH1F*)gDirectory->Get(A.c_str());
  QCDJET4ETA->Scale(Lumi*XS[FirstGoodMarker]*0.5/ProcessedEvents[FirstGoodMarker]);
  A=Form("jet5_eta%i",FirstGoodMarker);
  TH1F *QCDJET5ETA=(TH1F*)gDirectory->Get(A.c_str());
  QCDJET5ETA->Scale(Lumi*XS[FirstGoodMarker]*0.5/ProcessedEvents[FirstGoodMarker]);
  A=Form("jet6_eta%i",FirstGoodMarker);
  TH1F *QCDJET6ETA=(TH1F*)gDirectory->Get(A.c_str());
  QCDJET6ETA->Scale(Lumi*XS[FirstGoodMarker]*0.5/ProcessedEvents[FirstGoodMarker]);
  A=Form("jet1_phi%i",FirstGoodMarker);
  TH1F *QCDJET1PHI=(TH1F*)gDirectory->Get(A.c_str());
  QCDJET1PHI->Scale(Lumi*XS[FirstGoodMarker]*0.5/ProcessedEvents[FirstGoodMarker]);
  A=Form("jet2_phi%i",FirstGoodMarker);
  TH1F *QCDJET2PHI=(TH1F*)gDirectory->Get(A.c_str());
  QCDJET2PHI->Scale(Lumi*XS[FirstGoodMarker]*0.5/ProcessedEvents[FirstGoodMarker]);
  A=Form("jet3_phi%i",FirstGoodMarker);
  TH1F *QCDJET3PHI=(TH1F*)gDirectory->Get(A.c_str());
  QCDJET3PHI->Scale(Lumi*XS[FirstGoodMarker]*0.5/ProcessedEvents[FirstGoodMarker]);
  A=Form("jet4_phi%i",FirstGoodMarker);
  TH1F *QCDJET4PHI=(TH1F*)gDirectory->Get(A.c_str());
  QCDJET4PHI->Scale(Lumi*XS[FirstGoodMarker]*0.5/ProcessedEvents[FirstGoodMarker]);
  A=Form("jet5_phi%i",FirstGoodMarker);
  TH1F *QCDJET5PHI=(TH1F*)gDirectory->Get(A.c_str());
  QCDJET5PHI->Scale(Lumi*XS[FirstGoodMarker]*0.5/ProcessedEvents[FirstGoodMarker]);
  A=Form("jet6_phi%i",FirstGoodMarker);
  TH1F *QCDJET6PHI=(TH1F*)gDirectory->Get(A.c_str());
  QCDJET6PHI->Scale(Lumi*XS[FirstGoodMarker]*0.5/ProcessedEvents[FirstGoodMarker]);
  A=Form("Num_jets%i",FirstGoodMarker);
  TH1F *QCDJETMULTI=(TH1F*)gDirectory->Get(A.c_str());
  QCDJETMULTI->Scale(Lumi*XS[FirstGoodMarker]*0.5/ProcessedEvents[FirstGoodMarker]);
  for (int k=ParcialTestMin; k<ParcialTestMax; k++)  
    {
      if (!SurvivalMarker[k]) continue;
      cout << k << endl;
      FiveJetsMass[k]->Scale(Lumi*XS[k]/ProcessedEvents[k]);
      LeadingJetPT[k]->Scale(Lumi*XS[k]/ProcessedEvents[k]);
      Leading2JetPT[k]->Scale(Lumi*XS[k]/ProcessedEvents[k]);
      Leading3JetPT[k]->Scale(Lumi*XS[k]/ProcessedEvents[k]);
      Leading4JetPT[k]->Scale(Lumi*XS[k]/ProcessedEvents[k]);
      Leading5JetPT[k]->Scale(Lumi*XS[k]/ProcessedEvents[k]);
      Leading6JetPT[k]->Scale(Lumi*XS[k]/ProcessedEvents[k]);
      THT[k]->Scale(Lumi*XS[k]/ProcessedEvents[k]);
      DRHjets[k]->Scale(Lumi*XS[k]/ProcessedEvents[k]);
      DRWjets[k]->Scale(Lumi*XS[k]/ProcessedEvents[k]);
      Hpt[k]->Scale(Lumi*XS[k]/ProcessedEvents[k]);
      Tpt[k]->Scale(Lumi*XS[k]/ProcessedEvents[k]);
      DRWH[k]->Scale(Lumi*XS[k]/ProcessedEvents[k]);
      DPHjets[k]->Scale(Lumi*XS[k]/ProcessedEvents[k]);
      DPWjets[k]->Scale(Lumi*XS[k]/ProcessedEvents[k]);
      DPTjets[k]->Scale(Lumi*XS[k]/ProcessedEvents[k]);
      HiggsMass[k]->Scale(Lumi*XS[k]/ProcessedEvents[k]);
      RelHT[k]->Scale(Lumi*XS[k]/ProcessedEvents[k]);
      DRTH[k]->Scale(Lumi*XS[k]/ProcessedEvents[k]);
      PtNormalizedMass[k]->Scale(Lumi*XS[k]/ProcessedEvents[k]);
      RelativeMass[k]->Scale(Lumi*XS[k]/ProcessedEvents[k]);
      MotherPtNormalizedMass[k]->Scale(Lumi*XS[k]/ProcessedEvents[k]);   
      NumberOfTops[k]->Scale(Lumi*XS[k]/ProcessedEvents[k]); 
      HiggsMassOverTopMass[k]->Scale(Lumi*XS[k]/ProcessedEvents[k]);
      HiggsTopAsymmetry[k]->Scale(Lumi*XS[k]/ProcessedEvents[k]);
      ThirdLooseBtag[k]->Scale(Lumi*XS[k]/ProcessedEvents[k]);
      TopMass[k]->Scale(Lumi*XS[k]/ProcessedEvents[k]);
      Chi2[k]->Scale(Lumi*XS[k]/ProcessedEvents[k]);
      UQuarkContent[k]->Scale(Lumi*XS[k]/ProcessedEvents[k]);
      DQuarkContent[k]->Scale(Lumi*XS[k]/ProcessedEvents[k]);
      SQuarkContent[k]->Scale(Lumi*XS[k]/ProcessedEvents[k]);
      CQuarkContent[k]->Scale(Lumi*XS[k]/ProcessedEvents[k]);
      BQuarkContent[k]->Scale(Lumi*XS[k]/ProcessedEvents[k]);  
      CSVLB[k]->Scale(Lumi*XS[k]/ProcessedEvents[k]);
      CSVMB[k]->Scale(Lumi*XS[k]/ProcessedEvents[k]);
      CSVTB[k]->Scale(Lumi*XS[k]/ProcessedEvents[k]);  
      Vtcs[k]->Scale(Lumi*XS[k]/ProcessedEvents[k]);  
      for (int i=0; i<6; i++) HiggsMassReversedHptTpt[i][k]->Scale(Lumi*XS[k]/ProcessedEvents[k]);
      for (int i=0; i<6; i++) TprimeMassReversedHptTpt[i][k]->Scale(Lumi*XS[k]/ProcessedEvents[k]);
      for (int i=0; i<6; i++) HiggsMassCuts[i][k]->Scale(Lumi*XS[k]/ProcessedEvents[k]);
      /*for (int i=0; i<6; i++) WMassFromHiggs[i][k]->Scale(Lumi*XS[k]/ProcessedEvents[k]);*/
      TopMassFromHiggs[k]->Scale(Lumi*XS[k]/ProcessedEvents[k]);
      JET1ETA[k]->Scale(Lumi*XS[k]/ProcessedEvents[k]);
      JET2ETA[k]->Scale(Lumi*XS[k]/ProcessedEvents[k]);
      JET3ETA[k]->Scale(Lumi*XS[k]/ProcessedEvents[k]);
      JET4ETA[k]->Scale(Lumi*XS[k]/ProcessedEvents[k]);
      JET5ETA[k]->Scale(Lumi*XS[k]/ProcessedEvents[k]);
      JET6ETA[k]->Scale(Lumi*XS[k]/ProcessedEvents[k]);
      JET1PHI[k]->Scale(Lumi*XS[k]/ProcessedEvents[k]);
      JET2PHI[k]->Scale(Lumi*XS[k]/ProcessedEvents[k]);
      JET3PHI[k]->Scale(Lumi*XS[k]/ProcessedEvents[k]);
      JET4PHI[k]->Scale(Lumi*XS[k]/ProcessedEvents[k]);
      JET5PHI[k]->Scale(Lumi*XS[k]/ProcessedEvents[k]);
      JET6PHI[k]->Scale(Lumi*XS[k]/ProcessedEvents[k]);
      JETMULTI[k]->Scale(Lumi*XS[k]/ProcessedEvents[k]);
      //Settings for Single Top
      if (k<=ParcialTestMax && k>=ParcialTestMin) 
	{
	  FiveJetsMass[k]->SetFillColor(kBlue);
	  LeadingJetPT[k]->SetFillColor(kBlue);
	  Leading2JetPT[k]->SetFillColor(kBlue);
	  Leading3JetPT[k]->SetFillColor(kBlue);
	  Leading4JetPT[k]->SetFillColor(kBlue);
	  Leading5JetPT[k]->SetFillColor(kBlue);
	  Leading6JetPT[k]->SetFillColor(kBlue);
	  THT[k]->SetFillColor(kBlue);
	  DRHjets[k]->SetFillColor(kBlue);
	  DRWjets[k]->SetFillColor(kBlue);
	  Hpt[k]->SetFillColor(kBlue);
	  Tpt[k]->SetFillColor(kBlue);
	  DRWH[k]->SetFillColor(kBlue);
	  DPHjets[k]->SetFillColor(kBlue);
	  DPWjets[k]->SetFillColor(kBlue);
	  DPTjets[k]->SetFillColor(kBlue);
	  HiggsMass[k]->SetFillColor(kBlue);
	  RelHT[k]->SetFillColor(kBlue);
	  DRTH[k]->SetFillColor(kBlue);
	  PtNormalizedMass[k]->SetFillColor(kBlue);
	  RelativeMass[k]->SetFillColor(kBlue);
	  MotherPtNormalizedMass[k]->SetFillColor(kBlue);
	  NumberOfTops[k]->SetFillColor(kBlue);
	  HiggsMassOverTopMass[k]->SetFillColor(kBlue);
	  HiggsTopAsymmetry[k]->SetFillColor(kBlue);
	  ThirdLooseBtag[k]->SetFillColor(kBlue);
	  TopMass[k]->SetFillColor(kBlue);
	  Chi2[k]->SetFillColor(kBlue);
	  UQuarkContent[k]->SetFillColor(kBlue);
	  DQuarkContent[k]->SetFillColor(kBlue);
	  SQuarkContent[k]->SetFillColor(kBlue);
	  CQuarkContent[k]->SetFillColor(kBlue);
	  BQuarkContent[k]->SetFillColor(kBlue);
	  CSVLB[k]->SetFillColor(kBlue);
	  CSVMB[k]->SetFillColor(kBlue);
	  CSVTB[k]->SetFillColor(kBlue);
	  Vtcs[k]->SetFillColor(kBlue);
	  for (int i=0; i<6; i++) HiggsMassReversedHptTpt[i][k]->SetFillColor(kBlue);
	  for (int i=0; i<6; i++) TprimeMassReversedHptTpt[i][k]->SetFillColor(kBlue);
	  for (int i=0; i<6; i++) HiggsMassCuts[i][k]->SetFillColor(kBlue);
	  /*for (int i=0; i<6; i++) WMassFromHiggs[i][k]->SetFillColor(kBlue);*/
          
	  TopMassFromHiggs[k]->SetFillColor(kBlue); 
	  JET1ETA[k]->SetFillColor(kBlue);
	  JET2ETA[k]->SetFillColor(kBlue);
	  JET3ETA[k]->SetFillColor(kBlue);
	  JET4ETA[k]->SetFillColor(kBlue);
	  JET5ETA[k]->SetFillColor(kBlue);
	  JET6ETA[k]->SetFillColor(kBlue);
	  JET1PHI[k]->SetFillColor(kBlue);
	  JET2PHI[k]->SetFillColor(kBlue);
	  JET3PHI[k]->SetFillColor(kBlue);
	  JET4PHI[k]->SetFillColor(kBlue);
	  JET5PHI[k]->SetFillColor(kBlue);
	  JET6PHI[k]->SetFillColor(kBlue);
	  JETMULTI[k]->SetFillColor(kBlue);	
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
	      SingleTopTprimeMassCSVL->Add(CSVLB[k]);
	      SingleTopTprimeMassCSVM->Add(CSVMB[k]);
	      SingleTopTprimeMassCSVT->Add(CSVTB[k]);
	      SingleTopTprimeMassVTX->Add(Vtcs[k]);
	      for (int i=0; i<6; i++) SingleTop_HptTpt[i]->Add(HptTpt[i][k]);
	      SingleTop_HTCSVM->Add(HTCSVM[k]);
	      for (int i=0; i<6; i++) QCDTprimeMassHMBE[i]->Add(HiggsMassReversedHptTpt[i][k]);
	      for (int i=0; i<6; i++) QCDTprimeMassBE[i]->Add(TprimeMassReversedHptTpt[i][k]);
	      for (int i=0; i<6; i++) QCDTprimeMassHMCuts[i]->Add(HiggsMassCuts[i][k]);
	      /*for (int i=0; i<6; i++) QCDTprimeMassWMFHCuts[i]->Add(WMassFromHiggs[i][k]);*/
              QCDTopMassFromHiggs->Add(TopMassFromHiggs[k]);
	      QCDJET1ETA->Add(JET1ETA[k]);
	      QCDJET2ETA->Add(JET2ETA[k]);
	      QCDJET3ETA->Add(JET3ETA[k]);
	      QCDJET4ETA->Add(JET4ETA[k]);
	      QCDJET5ETA->Add(JET5ETA[k]);
	      QCDJET6ETA->Add(JET6ETA[k]);
	      QCDJET1PHI->Add(JET1PHI[k]);
	      QCDJET2PHI->Add(JET2PHI[k]);
	      QCDJET3PHI->Add(JET3PHI[k]);
	      QCDJET4PHI->Add(JET4PHI[k]);
	      QCDJET5PHI->Add(JET5PHI[k]);
	      QCDJET6PHI->Add(JET6PHI[k]);
	      QCDJETMULTI->Add(JETMULTI[k]);
	    }
	  if (k==LastGoodMarker)
	    {
	      TFile f("DY.root", "RECREATE");
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
	      SingleTopTprimeMassCSVL->Write();
	      SingleTopTprimeMassCSVM->Write();
	      SingleTopTprimeMassCSVT->Write();
	      SingleTopTprimeMassVTX->Write();
	      for (int i=0; i<6; i++) SingleTop_HptTpt[i]->Write();
	      SingleTop_HTCSVM->Write();
	      for (int i=0; i<6; i++) QCDTprimeMassHMBE[i]->Write();
	      for (int i=0; i<6; i++) QCDTprimeMassBE[i]->Write();
	      for (int i=0; i<6; i++) QCDTprimeMassHMCuts[i]->Write();
	      /*for (int i=0; i<6; i++) QCDTprimeMassWMFHCuts[i]->Write();*/              
              QCDTopMassFromHiggs->Write();
	      QCDJET1ETA->Write();
	      QCDJET2ETA->Write();
	      QCDJET3ETA->Write();
	      QCDJET4ETA->Write();
	      QCDJET5ETA->Write();
	      QCDJET6ETA->Write();
	      QCDJET1PHI->Write();
	      QCDJET2PHI->Write();
	      QCDJET3PHI->Write();
	      QCDJET4PHI->Write();
	      QCDJET5PHI->Write();
	      QCDJET6PHI->Write();
	      QCDJETMULTI->Write();
	    }
	}
    }

  double RegionA=0; double RegionB=0; double RegionC=0; double RegionD=0;
  RegionA=SingleTop_HptTpt[0]->Integral(SingleTop_HptTpt[0]->GetXaxis()->FindBin(200.0),SingleTop_HptTpt[0]->GetNbinsX(),SingleTop_HptTpt[0]->GetYaxis()->FindBin(200.0),SingleTop_HptTpt[0]->GetNbinsY());
  RegionB=SingleTop_HptTpt[0]->Integral(SingleTop_HptTpt[0]->GetXaxis()->FindBin(200.0),SingleTop_HptTpt[0]->GetNbinsX(),SingleTop_HptTpt[0]->GetYaxis()->FindBin(0.0),SingleTop_HptTpt[0]->GetYaxis()->FindBin(200.0));
  RegionC=SingleTop_HptTpt[0]->Integral(SingleTop_HptTpt[0]->GetXaxis()->FindBin(0.0),SingleTop_HptTpt[0]->GetXaxis()->FindBin(200.0),SingleTop_HptTpt[0]->GetYaxis()->FindBin(200.0),SingleTop_HptTpt[0]->GetNbinsY());
  RegionD=SingleTop_HptTpt[0]->Integral(SingleTop_HptTpt[0]->GetXaxis()->FindBin(0.0),SingleTop_HptTpt[0]->GetXaxis()->FindBin(200.0),SingleTop_HptTpt[0]->GetYaxis()->FindBin(0.0),SingleTop_HptTpt[0]->GetYaxis()->FindBin(200.0));

  cout << "ABCD Method Info for Hpt Top pt:" << endl;
  cout << "Number of events on region A (signal enriched) " <<  RegionA << endl;
  cout << "Number of events on region B " <<  RegionB << endl;
  cout << "Number of events on region C " <<  RegionC << endl;
  cout << "Number of events on region D " <<  RegionD << endl;

  double Region2A=0; double Region2B=0; double Region2C=0; double Region2D=0;
  Region2A=SingleTop_HTCSVM->Integral(SingleTop_HTCSVM->GetXaxis()->FindBin(630.0),SingleTop_HTCSVM->GetNbinsX(),SingleTop_HTCSVM->GetYaxis()->FindBin(3.0),SingleTop_HTCSVM->GetNbinsY());
  Region2B=SingleTop_HTCSVM->Integral(SingleTop_HTCSVM->GetXaxis()->FindBin(630.0),SingleTop_HTCSVM->GetNbinsX(),SingleTop_HTCSVM->GetYaxis()->FindBin(0.0),SingleTop_HTCSVM->GetYaxis()->FindBin(3.0));
  Region2C=SingleTop_HTCSVM->Integral(SingleTop_HTCSVM->GetXaxis()->FindBin(0.0),SingleTop_HTCSVM->GetXaxis()->FindBin(630.0),SingleTop_HTCSVM->GetYaxis()->FindBin(3.0),SingleTop_HTCSVM->GetNbinsY());
  Region2D=SingleTop_HTCSVM->Integral(SingleTop_HTCSVM->GetXaxis()->FindBin(0.0),SingleTop_HTCSVM->GetXaxis()->FindBin(630.0),SingleTop_HTCSVM->GetYaxis()->FindBin(0.0),SingleTop_HTCSVM->GetYaxis()->FindBin(3.0));

  cout << "ABCD Method Info for HT and CSVM b-tags:" << endl;
  cout << "Number of events on region A (signal enriched) " <<  Region2A << endl;
  cout << "Number of events on region B " <<  Region2B << endl;
  cout << "Number of events on region C " <<  Region2C << endl;
  cout << "Number of events on region D " <<  Region2D << endl;

  exit(0);
  
}
