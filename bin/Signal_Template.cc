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

const TString MainFolder = "file:/afs/cern.ch/work/j/jruizalv/private/Analysis/CMSSW_5_3_9_patch3/src/Extractors/PatExtractor/test/";

const int NumberOfHistos=23;
const TString Histos[NumberOfHistos] = {"Cut_0", "Cut_1", "Cut_2", "Cut_3", "Cut_chi2", "Cut_4", "Cut_5", "Cut_6", "Cut_7", "Cut_8", "Cut_9", "Cut_10", "Cut_11", "Cut_12", "Cut_13", 
"Cut_14", "Cut_15", "Cut_16", "Cut_17", "Cut_18", "Cut_19", "Cut_20", "Cut_21"};

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
  TH1F *Leading2JetPT;
  TH1F *Leading3JetPT;
  TH1F *Leading4JetPT;
  TH1F *Leading5JetPT;
  TH1F *Leading6JetPT;
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
  TH1F *NumberOfTops;
  TH1F *HiggsMassOverTopMass;
  TH1F *HiggsTopAsymmetry;
  TH1F *ThirdLooseBtag;  
  TH2F *HptTpt;
  TH1F *TopMass;
  TH1F *Chi2;
  
  int EntriePerSample;
  bool SurvivalMarker;

  TChain CutsChain("cuts");
  TChain AnalysisChain("stp");
  TH1F *ALLCuts[NumberOfHistos];
  CutsChain.Add(MainFolder + "Signal_SUFFIX_analyzed.root");
  AnalysisChain.Add(MainFolder + "Signal_SUFFIX_analyzed.root");
  
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

  //Cuts String
  //"jet1_pt>150 && THT>630 && (Reconstructed_Higgs.Pt()>200 && Reconstructed_Top.Pt()>300) && (DeltaR_of_W_Higgs>2.2 && DeltaR_of_W_Higgs<3.5) && (TMath::Abs(First_Higgs_Jet.Phi()-Second_Higgs_Jet.Phi())<2.0 && TMath::Abs(Top_Jet.Phi()-Reconstructed_W.Phi())<3.3) && TMath::Abs(First_W_Jet.Phi()-Second_W_Jet.Phi())<2.3 && (Reconstructed_Higgs.M()>100 && Reconstructed_Higgs.M()<135) && Relative_THT>0.65 && (DeltaR_of_Top_Higgs>2.8 && DeltaR_of_Top_Higgs<3.3) && PT_Normalized_Mass<0.7 && Mother_PT_Normalized_Mass<10 && Number_of_Tops<2"
  string CurrentCut = "jet1_pt>0";
  //////////////////////
  //Activate Only One!//
  //////////////////////
  if (false) CurrentCut = "jet1_pt>=150"; //Cut1 Only!

  if (false) CurrentCut = "jet1_pt>=150 && THT>=630 && (Reconstructed_Higgs.Pt()>=200 && Reconstructed_Top.Pt()>=300)"; //Up to Cut2 Only!

  if (false) CurrentCut = "jet1_pt>=150 && THT>=630 && (Reconstructed_Higgs.Pt()>=200 && Reconstructed_Top.Pt()>=300) && (DeltaR_of_W_Higgs>=2.2 && DeltaR_of_W_Higgs<=3.5)"; //Up to Cut6 Only!

  if (false) CurrentCut = "jet1_pt>=150 && THT>=630 && (Reconstructed_Higgs.Pt()>=200 && Reconstructed_Top.Pt()>=300) && (DeltaR_of_W_Higgs>=2.2 && DeltaR_of_W_Higgs<=3.5) && (TMath::Abs(First_Higgs_Jet.Phi()-Second_Higgs_Jet.Phi())<=2.0 && TMath::Abs(Top_Jet.Phi()-Reconstructed_W.Phi())<=3.3)"; //Up to Cut7 Only!

  if (false) CurrentCut = "jet1_pt>=150 && THT>=630 && (Reconstructed_Higgs.Pt()>=200 && Reconstructed_Top.Pt()>=300) && (DeltaR_of_W_Higgs>=2.2 && DeltaR_of_W_Higgs<=3.5) && (TMath::Abs(First_Higgs_Jet.Phi()-Second_Higgs_Jet.Phi())<=2.0 && TMath::Abs(Top_Jet.Phi()-Reconstructed_W.Phi())<=3.3) && TMath::Abs(First_W_Jet.Phi()-Second_W_Jet.Phi())<=2.3"; //Up to Cut8 Only!

  if (false) CurrentCut = "jet1_pt>=150 && THT>=630 && (Reconstructed_Higgs.Pt()>=200 && Reconstructed_Top.Pt()>=300) && (DeltaR_of_W_Higgs>=2.2 && DeltaR_of_W_Higgs<=3.5) && (TMath::Abs(First_Higgs_Jet.Phi()-Second_Higgs_Jet.Phi())<=2.0 && TMath::Abs(Top_Jet.Phi()-Reconstructed_W.Phi())<=3.3) && TMath::Abs(First_W_Jet.Phi()-Second_W_Jet.Phi())<=2.3 && (Reconstructed_Higgs.M()>=100 && Reconstructed_Higgs.M()<=135)"; //Up to Cut10 Only!

  if (false) CurrentCut = "jet1_pt>=150 && THT>=630 && (Reconstructed_Higgs.Pt()>=200 && Reconstructed_Top.Pt()>=300) && (DeltaR_of_W_Higgs>=2.2 && DeltaR_of_W_Higgs<=3.5) && (TMath::Abs(First_Higgs_Jet.Phi()-Second_Higgs_Jet.Phi())<=2.0 && TMath::Abs(Top_Jet.Phi()-Reconstructed_W.Phi())<=3.3) && TMath::Abs(First_W_Jet.Phi()-Second_W_Jet.Phi())<=2.3 && (Reconstructed_Higgs.M()>=100 && Reconstructed_Higgs.M()<=135) && Relative_THT>=0.6)"; //Up to Cut11 Only!

  if (false) CurrentCut = "jet1_pt>=150 && THT>=630 && (Reconstructed_Higgs.Pt()>=200 && Reconstructed_Top.Pt()>=300) && (DeltaR_of_W_Higgs>=2.2 && DeltaR_of_W_Higgs<=3.5) && (TMath::Abs(First_Higgs_Jet.Phi()-Second_Higgs_Jet.Phi())<=2.0 && TMath::Abs(Top_Jet.Phi()-Reconstructed_W.Phi())<=3.3) && TMath::Abs(First_W_Jet.Phi()-Second_W_Jet.Phi())<=2.3 && (Reconstructed_Higgs.M()>=100 && Reconstructed_Higgs.M()<=135) && Relative_THT>=0.65 && (DeltaR_of_Top_Higgs>=2.8 && DeltaR_of_Top_Higgs<=3.3)"; //Up to Cut12 Only!

  if (false) CurrentCut = "jet1_pt>=150 && THT>=630 && (Reconstructed_Higgs.Pt()>=200 && Reconstructed_Top.Pt()>=300) && (DeltaR_of_W_Higgs>=2.2 && DeltaR_of_W_Higgs<=3.5) && (TMath::Abs(First_Higgs_Jet.Phi()-Second_Higgs_Jet.Phi())<=2.0 && TMath::Abs(Top_Jet.Phi()-Reconstructed_W.Phi())<=3.3) && TMath::Abs(First_W_Jet.Phi()-Second_W_Jet.Phi())<=2.3 && (Reconstructed_Higgs.M()>=100 && Reconstructed_Higgs.M()<=135) && Relative_THT>=0.65 && (DeltaR_of_Top_Higgs>=2.8 && DeltaR_of_Top_Higgs<=3.3) && Mother_PT_Normalized_Mass<=10"; //To try 1

  if (false) CurrentCut = "jet1_pt>=150 && THT>=630 && (Reconstructed_Higgs.Pt()>=200 && Reconstructed_Top.Pt()>=300) && (DeltaR_of_W_Higgs>=2.2 && DeltaR_of_W_Higgs<=3.5) && (TMath::Abs(First_Higgs_Jet.Phi()-Second_Higgs_Jet.Phi())<=2.0 && TMath::Abs(Top_Jet.Phi()-Reconstructed_W.Phi())<=3.3) && TMath::Abs(First_W_Jet.Phi()-Second_W_Jet.Phi())<=2.3 && (Reconstructed_Higgs.M()>=100 && Reconstructed_Higgs.M()<=135) && Relative_THT>=0.65 && (DeltaR_of_Top_Higgs>=2.8 && DeltaR_of_Top_Higgs<=3.3) && Number_of_Tops<2"; //To trye 2
  
  if (PassedPerCut[NumberOfHistos-1]!=0)
    {	  
      /////////////////////
      //Saving Histograms//
      /////////////////////
      string A1 = Form("Reconstructed_Tprime.M() >> TprimeMass%i(60,400,1600)",15);
      string A2 = Form("TprimeMass%i",15);
      AnalysisChain.Draw(A1.c_str(),CurrentCut.c_str());
      FiveJetsMass = (TH1F*)gDirectory->Get(A2.c_str());
      gPad->Close();
      SurvivalMarker=true;
      //
      string LJPT1 = Form("jet1_pt >> jet1_pt%i(50,20,1000)",15);
      string LJPT2 = Form("jet1_pt%i",15);
      AnalysisChain.Draw(LJPT1.c_str(),CurrentCut.c_str());
      LeadingJetPT = (TH1F*)gDirectory->Get(LJPT2.c_str());
      gPad->Close();
      //
      string SLJPT1 = Form("jet2_pt >> jet2_pt%i(50,20,1000)",15);
      string SLJPT2 = Form("jet2_pt%i",15);
      AnalysisChain.Draw(SLJPT1.c_str(),CurrentCut.c_str());
      Leading2JetPT = (TH1F*)gDirectory->Get(SLJPT2.c_str());
      gPad->Close();
      //
      string SSLJPT1 = Form("jet3_pt >> jet3_pt%i(35,20,700)",15);
      string SSLJPT2 = Form("jet3_pt%i",15);
      AnalysisChain.Draw(SSLJPT1.c_str(),CurrentCut.c_str());
      Leading3JetPT = (TH1F*)gDirectory->Get(SSLJPT2.c_str());
      gPad->Close();
      //
      string SSSLJPT1 = Form("jet4_pt >> jet4_pt%i(20,20,400)",15);
      string SSSLJPT2 = Form("jet4_pt%i",15);
      AnalysisChain.Draw(SSSLJPT1.c_str(),CurrentCut.c_str());
      Leading4JetPT = (TH1F*)gDirectory->Get(SSSLJPT2.c_str());
      gPad->Close();
      //
      string SSSSLJPT1 = Form("jet5_pt >> jet5_pt%i(15,20,300)",15);
      string SSSSLJPT2 = Form("jet5_pt%i",15);
      AnalysisChain.Draw(SSSSLJPT1.c_str(),CurrentCut.c_str());
      Leading5JetPT = (TH1F*)gDirectory->Get(SSSSLJPT2.c_str());
      gPad->Close();
      //
      string SSSSSLJPT1 = Form("jet6_pt >> jet6_pt%i(10,20,200)",15);
      string SSSSSLJPT2 = Form("jet6_pt%i",15);
      AnalysisChain.Draw(SSSSSLJPT1.c_str(),CurrentCut.c_str());
      Leading6JetPT = (TH1F*)gDirectory->Get(SSSSSLJPT2.c_str());
      gPad->Close();
      //
      string THT1 = Form("THT >> THT%i(65,300,1600)",15);
      string THT2 = Form("THT%i",15);
      AnalysisChain.Draw(THT1.c_str(),CurrentCut.c_str());
      THT = (TH1F*)gDirectory->Get(THT2.c_str());
      gPad->Close(); 
      //
      string DRHJ1 = Form("DeltaR_of_Higgs_Jets >> DeltaR_of_Higgs_Jets%i(65,0.5,7)",15);
      string DRHJ2 = Form("DeltaR_of_Higgs_Jets%i",15);
      AnalysisChain.Draw(DRHJ1.c_str(),CurrentCut.c_str());
      DRHjets = (TH1F*)gDirectory->Get(DRHJ2.c_str());
      gPad->Close();
      //
      string DRWJ1 = Form("DeltaR_of_W_Jets >> DeltaR_of_W_Jets%i(65,0.5,7)",15);
      string DRWJ2 = Form("DeltaR_of_W_Jets%i",15);
      AnalysisChain.Draw(DRWJ1.c_str(),CurrentCut.c_str());
      DRWjets = (TH1F*)gDirectory->Get(DRWJ2.c_str());
      gPad->Close();
      //
      string HPT1 = Form("Reconstructed_Higgs.Pt() >> HPt%i(40,10,800)",15);
      string HPT2 = Form("HPt%i",15);
      AnalysisChain.Draw(HPT1.c_str(),CurrentCut.c_str());
      Hpt = (TH1F*)gDirectory->Get(HPT2.c_str());
      gPad->Close();
      //
      string TPT1 = Form("Reconstructed_Top.Pt() >> TPt%i(40,10,800)",15);
      string TPT2 = Form("TPt%i",15);
      AnalysisChain.Draw(TPT1.c_str(),CurrentCut.c_str());
      Tpt = (TH1F*)gDirectory->Get(TPT2.c_str());
      gPad->Close();
      //
      string DRWH1 = Form("DeltaR_of_W_Higgs >> DeltaR_of_W_Higgs%i(65,0.5,7)",15);
      string DRWH2 = Form("DeltaR_of_W_Higgs%i",15);
      AnalysisChain.Draw(DRWH1.c_str(),CurrentCut.c_str());
      DRWH = (TH1F*)gDirectory->Get(DRWH2.c_str());
      gPad->Close();
      //
      string DPHJ1 = Form("TMath::Abs(First_Higgs_Jet.Phi()-Second_Higgs_Jet.Phi()) >> DeltaPhi_of_Higgs_jets%i(30,0.0,3.0)",15);
      string DPHJ2 = Form("DeltaPhi_of_Higgs_jets%i",15);
      AnalysisChain.Draw(DPHJ1.c_str(),CurrentCut.c_str());
      DPHjets = (TH1F*)gDirectory->Get(DPHJ2.c_str());
      gPad->Close();
      //
      string DPWJ1 = Form("TMath::Abs(First_W_Jet.Phi()-Second_W_Jet.Phi()) >> DeltaPhi_of_W_jets%i(30,0.0,3.0)",15);
      string DPWJ2 = Form("DeltaPhi_of_W_jets%i",15);
      AnalysisChain.Draw(DPWJ1.c_str(),CurrentCut.c_str());
      DPWjets = (TH1F*)gDirectory->Get(DPWJ2.c_str());
      gPad->Close();
      //
      string DPTJ1 = Form("TMath::Abs(Top_Jet.Phi()-Reconstructed_W.Phi()) >> DeltaPhi_of_T_jet%i(30,0.0,3.0)",15);
      string DPTJ2 = Form("DeltaPhi_of_T_jet%i",15);
      AnalysisChain.Draw(DPTJ1.c_str(),CurrentCut.c_str());
      DPTjets = (TH1F*)gDirectory->Get(DPTJ2.c_str());
      gPad->Close();
      //
      string HM1 = Form("Reconstructed_Higgs.M() >> HM%i(36,60,180)",15);
      string HM2 = Form("HM%i",15);
      AnalysisChain.Draw(HM1.c_str(),CurrentCut.c_str());
      HiggsMass = (TH1F*)gDirectory->Get(HM2.c_str());
      gPad->Close();
      //
      string RHT1 = Form("Relative_THT >> RelHT%i(30,0,1)",15);
      string RHT2 = Form("RelHT%i",15);
      AnalysisChain.Draw(RHT1.c_str(),CurrentCut.c_str());
      RelHT = (TH1F*)gDirectory->Get(RHT2.c_str());
      gPad->Close();
      //
      string DRTH1 = Form("DeltaR_of_Top_Higgs >> DeltaR_of_Top_Higgs%i(65,0.5,7)",15);
      string DRTH2 = Form("DeltaR_of_Top_Higgs%i",15);
      AnalysisChain.Draw(DRTH1.c_str(),CurrentCut.c_str());
      DRTH = (TH1F*)gDirectory->Get(DRTH2.c_str());
      gPad->Close();
      //
      string PTNM1 = Form("PT_Normalized_Mass >> PT_Normalized_Mass%i(70,0.4,5)",15);
      string PTNM2 = Form("PT_Normalized_Mass%i",15);
      AnalysisChain.Draw(PTNM1.c_str(),CurrentCut.c_str());
      PtNormalizedMass = (TH1F*)gDirectory->Get(PTNM2.c_str());
      gPad->Close();
      //
      string RM1 = Form("Relative_Mass >> Relative_Mass%i(30,0.0,1)",15);
      string RM2 = Form("Relative_Mass%i",15);
      AnalysisChain.Draw(RM1.c_str(),CurrentCut.c_str());
      RelativeMass = (TH1F*)gDirectory->Get(RM2.c_str());
      gPad->Close();
      //
      string MPTNM1 = Form("Mother_PT_Normalized_Mass >> Mother_PT_Normalized_Mass%i(25,0.0,50)",15);
      string MPTNM2 = Form("Mother_PT_Normalized_Mass%i",15);
      AnalysisChain.Draw(MPTNM1.c_str(),CurrentCut.c_str());
      MotherPtNormalizedMass = (TH1F*)gDirectory->Get(MPTNM2.c_str());
      gPad->Close();
      //
      string NTops1 = Form("Number_of_Tops >> Number_of_Tops%i(8,0.0,8)",15);
      string NTops2 = Form("Number_of_Tops%i",15);
      AnalysisChain.Draw(NTops1.c_str(),CurrentCut.c_str());
      NumberOfTops = (TH1F*)gDirectory->Get(NTops2.c_str());
      gPad->Close();
      //
      string HPTTPT1 = Form("Reconstructed_Higgs.Pt():Reconstructed_Top.Pt() >> HPtTPt%i(40,10,800,40,10,800)",15);
      string HPTTPT2 = Form("HPtTPt%i",15);
      AnalysisChain.Draw(HPTTPT1.c_str(),CurrentCut.c_str());
      HptTpt = (TH2F*)gDirectory->Get(HPTTPT2.c_str());
      gPad->Close();
      ///////////NEW VARIABLES//////////////////
      string HMTM1 = Form("Reconstructed_Higgs.M()/Reconstructed_Top.M() >> HMoverTM%i(30,0,1.)",15);
      string HMTM2 = Form("HMoverTM%i",15);
      AnalysisChain.Draw(HMTM1.c_str());
      HiggsMassOverTopMass = (TH1F*)gDirectory->Get(HMTM2.c_str());
      gPad->Close();
      //
      string HTA1 = Form("((Reconstructed_Higgs.Pt()/Reconstructed_Higgs.M())-(Reconstructed_Top.Pt()/Reconstructed_Top.M()))/((Reconstructed_Higgs.Pt()/Reconstructed_Higgs.M())+(Reconstructed_Top.Pt()/Reconstructed_Top.M())) >> HTAsym%i(12,0,1.)",15);
      string HTA2 = Form("HTAsym%i",15);
      AnalysisChain.Draw(HTA1.c_str());
      HiggsTopAsymmetry = (TH1F*)gDirectory->Get(HTA2.c_str());
      gPad->Close();
      //
      string TLBT1 = Form("Number_of_Loose_and_non_med_b_tags >> TLBTag%i(10,0,10)",15);
      string TLBT2 = Form("TLBTag%i",15);
      AnalysisChain.Draw(TLBT1.c_str());
      ThirdLooseBtag = (TH1F*)gDirectory->Get(TLBT2.c_str());
      gPad->Close();
      //
      string TM1 = Form("Reconstructed_Top.M() >> TMass%i(40,100,300)",15);
      string TM2 = Form("TMass%i",15);
      AnalysisChain.Draw(TM1.c_str());
      TopMass = (TH1F*)gDirectory->Get(TM2.c_str());
      gPad->Close();
      //
      string C21 = Form("ChiSquaredSorting >> ChiSq%i(100,0,1000)",15);
      string C22 = Form("ChiSq%i",15);
      AnalysisChain.Draw(C21.c_str());
      Chi2 = (TH1F*)gDirectory->Get(C22.c_str());
      gPad->Close();
    }
    
  FiveJetsMass->Scale(Lumi*XS/EntriePerSample);
  LeadingJetPT->Scale(Lumi*XS/EntriePerSample);
  Leading2JetPT->Scale(Lumi*XS/EntriePerSample);
  Leading3JetPT->Scale(Lumi*XS/EntriePerSample);
  Leading4JetPT->Scale(Lumi*XS/EntriePerSample);
  Leading5JetPT->Scale(Lumi*XS/EntriePerSample);
  Leading6JetPT->Scale(Lumi*XS/EntriePerSample);
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
  NumberOfTops->Scale(Lumi*XS/EntriePerSample);
  HptTpt->Scale(Lumi*XS/EntriePerSample);
  HiggsMassOverTopMass->Scale(Lumi*XS/EntriePerSample);
  HiggsTopAsymmetry->Scale(Lumi*XS/EntriePerSample);
  ThirdLooseBtag->Scale(Lumi*XS/EntriePerSample);
  TopMass->Scale(Lumi*XS/EntriePerSample);
  Chi2->Scale(Lumi*XS/EntriePerSample);
  
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
      Leading2JetPT->SetFillColor(kSpring);
      Leading2JetPT->SetFillStyle(3444);
      Leading2JetPT->SetLineWidth(3);
      Leading2JetPT->Write();
      Leading3JetPT->SetFillColor(kSpring);
      Leading3JetPT->SetFillStyle(3444);
      Leading3JetPT->SetLineWidth(3);
      Leading3JetPT->Write();
      Leading4JetPT->SetFillColor(kSpring);
      Leading4JetPT->SetFillStyle(3444);
      Leading4JetPT->SetLineWidth(3);
      Leading4JetPT->Write();
      Leading5JetPT->SetFillColor(kSpring);
      Leading5JetPT->SetFillStyle(3444);
      Leading5JetPT->SetLineWidth(3);
      Leading5JetPT->Write();
      Leading6JetPT->SetFillColor(kSpring);
      Leading6JetPT->SetFillStyle(3444);
      Leading6JetPT->SetLineWidth(3);
      Leading6JetPT->Write();
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
      //
      NumberOfTops->SetFillColor(kSpring);
      NumberOfTops->SetFillStyle(3444);
      NumberOfTops->SetLineWidth(3);
      NumberOfTops->Write();
      //
      HptTpt->SetFillColor(kSpring);
      HptTpt->SetFillStyle(3444);
      HptTpt->SetLineWidth(3);
      HptTpt->Write();
      //
      HiggsMassOverTopMass->SetFillColor(kSpring);
      HiggsMassOverTopMass->SetFillStyle(3444);
      HiggsMassOverTopMass->SetLineWidth(3);
      HiggsMassOverTopMass->Write();
      //
      HiggsTopAsymmetry->SetFillColor(kSpring);
      HiggsTopAsymmetry->SetFillStyle(3444);
      HiggsTopAsymmetry->SetLineWidth(3);
      HiggsTopAsymmetry->Write();
      //
      ThirdLooseBtag->SetFillColor(kSpring);
      ThirdLooseBtag->SetFillStyle(3444);
      ThirdLooseBtag->SetLineWidth(3);
      ThirdLooseBtag->Write();
      //
      TopMass->SetFillColor(kSpring);
      TopMass->SetFillStyle(3444);
      TopMass->SetLineWidth(3);
      TopMass->Write();
      //
      Chi2->SetFillColor(kSpring);
      Chi2->SetFillStyle(3444);
      Chi2->SetLineWidth(3);
      Chi2->Write();
    }

  exit(0);

}

