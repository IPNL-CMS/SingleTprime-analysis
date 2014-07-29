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

//const int NumberOfProcesses=17;

const TString MainFolder = "file:/sps/cms/ruizalva/CMSSW_5_3_9_patch3/src/Extractors/PatExtractor/test/";

const int NumberOfHistos=23;
const TString Histos[NumberOfHistos] = {"Cut_0", "Cut_1", "Cut_2", "Cut_3", "Cut_chi2", "Cut_4", "Cut_5", "Cut_6", "Cut_7", "Cut_8", "Cut_9", "Cut_10", "Cut_11", "Cut_12", "Cut_13", 
"Cut_14", "Cut_15", "Cut_16", "Cut_17", "Cut_18", "Cut_19", "Cut_20", "Cut_21"};

//Computing weights (Everything in pb)
const float Xs = XS[18]; 
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

//float Lumi=20000.;

float PUR_function(int TI) //function with input the Number of True Interactions
{

  if (TI>=PUBins) return 1;
  else return PU_weight[TI];

}

double ABs(double X, double Y)
{
  return TMath::Abs(X-Y);
}

void Signal()
{
  TH1F *FiveJetsMass;
  TH1F *FiveJetsMass3L;
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
  TH1F *TopMass;
  TH1F *Chi2;
  TH1F *UQuarkContent;
  TH1F *DQuarkContent;
  TH1F *SQuarkContent;
  TH1F *CQuarkContent;
  TH1F *BQuarkContent;
  TH1F *CSVLB;
  TH1F *CSVMB;
  TH1F *CSVTB;
  TH2F *HptTpt[6];
  TH1F *Vtcs; 
  TH2F *HTCSVM; 
  TH1F *HiggsMassCuts[6]; 
  TH1F *HiggsMassReversedHptTpt[6];
  TH1F *TprimeMassReversedHptTpt[6];   
  TH1F *WMassFromHiggs[6];   
  TH1F *WMassFromHiggsChi2[6];
  TH1F *TprimeNotWFromHiggs;
  TH1F *FiveJetsMassBE; //5 jets mass inside the last cut
  TH1F *FiveJetsMassLC; //5 jets mass outside last cut
  TH1F *TopMassFromHiggs;
  TH1F *FiveJetsMassLCoverBE; //5 jets mass ratio plot in vs out the cut
  TH1F *HMBE;
  TH1F *HMLC;
  TH1F *HMLCoverBE;
  TH1F *FJLC_QCDBE;
  TH1F *FJLC_QCDLC;
  TH1F *FJLC_QCDLCoverBE;

  //First top study for estimation of overall QCD
  TH1F *TM_LC_HMW; //Top mass after ttbar estimation cut inside the Higgs mass window

  //HT estimation
  TH2F *DPWJ_HT;

  //ADDITIONAL INFO ON JETS  
  TH1F *JET1ETA;   
  TH1F *JET2ETA; 
  TH1F *JET3ETA; 
  TH1F *JET4ETA;
  TH1F *JET5ETA;
  TH1F *JET6ETA;
  TH1F *JET1PHI;   
  TH1F *JET2PHI; 
  TH1F *JET3PHI; 
  TH1F *JET4PHI;
  TH1F *JET5PHI;
  TH1F *JET6PHI;
  TH1F *JETMULTI;
  
  int EntriePerSample;
  bool SurvivalMarker;

  TChain CutsChain("cuts");
  TChain AnalysisChain("stp");
  TChain AnalysisChain3L("stp3L");
  TH1F *ALLCuts[NumberOfHistos];
  CutsChain.Add(MainFolder + "Signal_SUFFIX_analyzed.root");
  AnalysisChain.Add(MainFolder + "Signal_SUFFIX_analyzed.root");
  AnalysisChain3L.Add(MainFolder + "Signal_SUFFIX_analyzed.root");
  
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
      string A1 = Form("Reconstructed_Tprime.M() >> TprimeMass%i(60,400,1600)",18);
      string A2 = Form("TprimeMass%i",18);
      AnalysisChain.Draw(A1.c_str(),finalcut);
      FiveJetsMass = (TH1F*)gDirectory->Get(A2.c_str());
      gPad->Close();
      SurvivalMarker=true;
      //
      string LJPT1 = Form("jet1_pt >> jet1_pt%i(50,20,1000)",18);
      string LJPT2 = Form("jet1_pt%i",18);
      AnalysisChain.Draw(LJPT1.c_str(),"(PUR_function(Number_True_Interactions))*weight");
      LeadingJetPT = (TH1F*)gDirectory->Get(LJPT2.c_str());
      gPad->Close();
      //
      A1 = Form("Reconstructed_Tprime.M() >> TprimeMassBkgE%i(60,400,1600)",18);
      A2 = Form("TprimeMassBkgE%i",18);
      AnalysisChain.Draw(A1.c_str(),"Top_From_Higgs_Chi2.M()>=140 && Top_From_Higgs_Chi2.M()<=230");
      FiveJetsMassBE = (TH1F*)gDirectory->Get(A2.c_str());
      gPad->Close();
      //
      A1 = Form("Reconstructed_Tprime.M() >> TprimeMassLC%i(60,400,1600)",18);
      A2 = Form("TprimeMassLC%i",18);
      AnalysisChain.Draw(A1.c_str(),"Top_From_Higgs_Chi2.M()<140 || Top_From_Higgs_Chi2.M()>230");
      FiveJetsMassLC = (TH1F*)gDirectory->Get(A2.c_str());
      gPad->Close();
      //
      A1 = Form("Reconstructed_Tprime.M() >> TprimeMassLCoverBE%i(60,400,1600)",18);
      A2 = Form("TprimeMassLCoverBE%i",18);
      AnalysisChain.Draw(A1.c_str(),"Top_From_Higgs_Chi2.M()<140 || Top_From_Higgs_Chi2.M()>230");
      FiveJetsMassLCoverBE = (TH1F*)gDirectory->Get(A2.c_str());
      gPad->Close();
      //
      string HMBEstim1 = Form("Reconstructed_Higgs.M() >> HiggsMassBE%i(36,60,180)",18);
      string HMBEstim2 = Form("HiggsMassBE%i",18);
      AnalysisChain.Draw(HMBEstim1.c_str(),"Top_From_Higgs_Chi2.M()>=140 && Top_From_Higgs_Chi2.M()<=230");
      HMBE = (TH1F*)gDirectory->Get(HMBEstim2.c_str());
      gPad->Close();
      //
      HMBEstim1 = Form("Reconstructed_Higgs.M() >> HiggsMassLC%i(36,60,180)",18);
      HMBEstim2 = Form("HiggsMassLC%i",18);
      AnalysisChain.Draw(HMBEstim1.c_str(),"Top_From_Higgs_Chi2.M()<140 || Top_From_Higgs_Chi2.M()>230");
      HMLC = (TH1F*)gDirectory->Get(HMBEstim2.c_str());
      gPad->Close();
      //
      HMBEstim1 = Form("Reconstructed_Higgs.M() >> HiggsMassLCoverBE%i(36,60,180)",18);
      HMBEstim2 = Form("HiggsMassLCoverBE%i",18);
      AnalysisChain.Draw(HMBEstim1.c_str(),"Top_From_Higgs_Chi2.M()<140 || Top_From_Higgs_Chi2.M()>230");
      HMLCoverBE = (TH1F*)gDirectory->Get(HMBEstim2.c_str());
      gPad->Close();
      //
      HMBEstim1 = Form("Reconstructed_Tprime.M() >> FJLCQCDBE%i(60,400,1600)",18);
      HMBEstim2 = Form("FJLCQCDBE%i",18);
      AnalysisChain.Draw(HMBEstim1.c_str(),"(Top_From_Higgs_Chi2.M()<140 || Top_From_Higgs_Chi2.M()>230) && (Reconstructed_Higgs.M()<110 || Reconstructed_Higgs.M()>140)");
      FJLC_QCDBE = (TH1F*)gDirectory->Get(HMBEstim2.c_str());
      gPad->Close();
      //
      HMBEstim1 = Form("Reconstructed_Tprime.M() >> FJLCQCDLC%i(60,400,1600)",18);
      HMBEstim2 = Form("FJLCQCDLC%i",18);
      AnalysisChain.Draw(HMBEstim1.c_str(),"(Top_From_Higgs_Chi2.M()<140 || Top_From_Higgs_Chi2.M()>230) && (Reconstructed_Higgs.M()>=110 && Reconstructed_Higgs.M()<=140)");
      FJLC_QCDLC = (TH1F*)gDirectory->Get(HMBEstim2.c_str());
      gPad->Close();
      //
      HMBEstim1 = Form("Reconstructed_Tprime.M() >> FJLCQCDLCoverBE%i(60,400,1600)",18);
      HMBEstim2 = Form("FJLCQCDLCoverBE%i",18);
      AnalysisChain.Draw(HMBEstim1.c_str(),"(Top_From_Higgs_Chi2.M()<140 || Top_From_Higgs_Chi2.M()>230) && (Reconstructed_Higgs.M()>=110 && Reconstructed_Higgs.M()<=140)");
      FJLC_QCDLCoverBE = (TH1F*)gDirectory->Get(HMBEstim2.c_str());
      gPad->Close();
      //////////
      //////////
      A1 = Form("Top_From_Higgs_Chi2.M() >> TopFromHiggsChi2Mass%i(95,50,1000)",18);
      A2 = Form("TopFromHiggsChi2Mass%i",18);
      AnalysisChain.Draw(A1.c_str(),"(PUR_function(Number_True_Interactions))*weight");
      TopMassFromHiggs = (TH1F*)gDirectory->Get(A2.c_str());
      gPad->Close();
      //////////
      //////////
      string TM1_LC = Form("Reconstructed_Top.M() >> TMass_LC_HMW%i(60,100,700)",18);
      string TM2_LC = Form("TMass_LC_HMW%i",18);
      AnalysisChain.Draw(TM1_LC.c_str(),"(Top_From_Higgs_Chi2.M()<140 || Top_From_Higgs_Chi2.M()>230) && (Reconstructed_Higgs.M()>=110 && Reconstructed_Higgs.M()<=140)");
      TM_LC_HMW = (TH1F*)gDirectory->Get(TM2_LC.c_str());
      gPad->Close();
      //////////
      //////////
      string SLJPT1 = Form("jet2_pt >> jet2_pt%i(50,20,1000)",18);
      string SLJPT2 = Form("jet2_pt%i",18);
      AnalysisChain.Draw(SLJPT1.c_str(),"(PUR_function(Number_True_Interactions))*weight");
      Leading2JetPT = (TH1F*)gDirectory->Get(SLJPT2.c_str());
      gPad->Close();
      //
      string SSLJPT1 = Form("jet3_pt >> jet3_pt%i(35,20,700)",18);
      string SSLJPT2 = Form("jet3_pt%i",18);
      AnalysisChain.Draw(SSLJPT1.c_str(),"(PUR_function(Number_True_Interactions))*weight");
      Leading3JetPT = (TH1F*)gDirectory->Get(SSLJPT2.c_str());
      gPad->Close();
      //
      string SSSLJPT1 = Form("jet4_pt >> jet4_pt%i(20,20,400)",18);
      string SSSLJPT2 = Form("jet4_pt%i",18);
      AnalysisChain.Draw(SSSLJPT1.c_str(),"(PUR_function(Number_True_Interactions))*weight");
      Leading4JetPT = (TH1F*)gDirectory->Get(SSSLJPT2.c_str());
      gPad->Close();
      //
      string SSSSLJPT1 = Form("jet5_pt >> jet5_pt%i(15,20,300)",18);
      string SSSSLJPT2 = Form("jet5_pt%i",18);
      AnalysisChain.Draw(SSSSLJPT1.c_str(),"(PUR_function(Number_True_Interactions))*weight");
      Leading5JetPT = (TH1F*)gDirectory->Get(SSSSLJPT2.c_str());
      gPad->Close();
      //
      string SSSSSLJPT1 = Form("jet6_pt >> jet6_pt%i(10,20,200)",18);
      string SSSSSLJPT2 = Form("jet6_pt%i",18);
      AnalysisChain.Draw(SSSSSLJPT1.c_str(),"(PUR_function(Number_True_Interactions))*weight");
      Leading6JetPT = (TH1F*)gDirectory->Get(SSSSSLJPT2.c_str());
      gPad->Close();
      //
      string THT1 = Form("THT >> THT%i(65,300,1600)",18);
      string THT2 = Form("THT%i",18);
      AnalysisChain.Draw(THT1.c_str(),"(PUR_function(Number_True_Interactions))*weight");
      THT = (TH1F*)gDirectory->Get(THT2.c_str());
      gPad->Close(); 
      //
      string DRHJ1 = Form("DeltaR_of_Higgs_Jets >> DeltaR_of_Higgs_Jets%i(65,0.5,7)",18);
      string DRHJ2 = Form("DeltaR_of_Higgs_Jets%i",18);
      AnalysisChain.Draw(DRHJ1.c_str(),"(PUR_function(Number_True_Interactions))*weight");
      DRHjets = (TH1F*)gDirectory->Get(DRHJ2.c_str());
      gPad->Close();
      //
      string DRWJ1 = Form("DeltaR_of_W_Jets >> DeltaR_of_W_Jets%i(65,0.5,7)",18);
      string DRWJ2 = Form("DeltaR_of_W_Jets%i",18);
      AnalysisChain.Draw(DRWJ1.c_str(),"(PUR_function(Number_True_Interactions))*weight");
      DRWjets = (TH1F*)gDirectory->Get(DRWJ2.c_str());
      gPad->Close();
      //
      string HPT1 = Form("Reconstructed_Higgs.Pt() >> HPt%i(40,10,800)",18);
      string HPT2 = Form("HPt%i",18);
      AnalysisChain.Draw(HPT1.c_str(),"(PUR_function(Number_True_Interactions))*weight");
      Hpt = (TH1F*)gDirectory->Get(HPT2.c_str());
      gPad->Close();
      //
      string TPT1 = Form("Reconstructed_Top.Pt() >> TPt%i(40,10,800)",18);
      string TPT2 = Form("TPt%i",18);
      AnalysisChain.Draw(TPT1.c_str(),"(PUR_function(Number_True_Interactions))*weight");
      Tpt = (TH1F*)gDirectory->Get(TPT2.c_str());
      gPad->Close();
      //
      string DRWH1 = Form("DeltaR_of_W_Higgs >> DeltaR_of_W_Higgs%i(65,0.5,7)",18);
      string DRWH2 = Form("DeltaR_of_W_Higgs%i",18);
      AnalysisChain.Draw(DRWH1.c_str(),"(PUR_function(Number_True_Interactions))*weight");
      DRWH = (TH1F*)gDirectory->Get(DRWH2.c_str());
      gPad->Close();
      //
      string DPHJ1 = Form("TMath::Abs(First_Higgs_Jet.Phi()-Second_Higgs_Jet.Phi()) >> DeltaPhi_of_Higgs_jets%i(30,0.0,3.0)",18);
      string DPHJ2 = Form("DeltaPhi_of_Higgs_jets%i",18);
      AnalysisChain.Draw(DPHJ1.c_str(),"(PUR_function(Number_True_Interactions))*weight");
      DPHjets = (TH1F*)gDirectory->Get(DPHJ2.c_str());
      gPad->Close();
      //
      string DPWJ1 = Form("TMath::Abs(First_W_Jet.Phi()-Second_W_Jet.Phi()) >> DeltaPhi_of_W_jets%i(30,0.0,3.0)",18);
      string DPWJ2 = Form("DeltaPhi_of_W_jets%i",18);
      AnalysisChain.Draw(DPWJ1.c_str(),"(PUR_function(Number_True_Interactions))*weight");
      DPWjets = (TH1F*)gDirectory->Get(DPWJ2.c_str());
      gPad->Close();
      //
      string DPTJ1 = Form("TMath::Abs(Top_Jet.Phi()-Reconstructed_W.Phi()) >> DeltaPhi_of_T_jet%i(30,0.0,3.0)",18);
      string DPTJ2 = Form("DeltaPhi_of_T_jet%i",18);
      AnalysisChain.Draw(DPTJ1.c_str(),"(PUR_function(Number_True_Interactions))*weight");
      DPTjets = (TH1F*)gDirectory->Get(DPTJ2.c_str());
      gPad->Close();
      //
      string HM1 = Form("Reconstructed_Higgs.M() >> HM%i(36,60,180)",18);
      string HM2 = Form("HM%i",18);
      AnalysisChain.Draw(HM1.c_str(),"(PUR_function(Number_True_Interactions))*weight");
      HiggsMass = (TH1F*)gDirectory->Get(HM2.c_str());
      gPad->Close();
      //
      string RHT1 = Form("Relative_THT >> RelHT%i(30,0,1)",18);
      string RHT2 = Form("RelHT%i",18);
      AnalysisChain.Draw(RHT1.c_str(),"(PUR_function(Number_True_Interactions))*weight");
      RelHT = (TH1F*)gDirectory->Get(RHT2.c_str());
      gPad->Close();
      //
      string DRTH1 = Form("DeltaR_of_Top_Higgs >> DeltaR_of_Top_Higgs%i(65,0.5,7)",18);
      string DRTH2 = Form("DeltaR_of_Top_Higgs%i",18);
      AnalysisChain.Draw(DRTH1.c_str(),"(PUR_function(Number_True_Interactions))*weight");
      DRTH = (TH1F*)gDirectory->Get(DRTH2.c_str());
      gPad->Close();
      //
      string PTNM1 = Form("PT_Normalized_Mass >> PT_Normalized_Mass%i(70,0.4,5)",18);
      string PTNM2 = Form("PT_Normalized_Mass%i",18);
      AnalysisChain.Draw(PTNM1.c_str(),"(PUR_function(Number_True_Interactions))*weight");
      PtNormalizedMass = (TH1F*)gDirectory->Get(PTNM2.c_str());
      gPad->Close();
      //
      string RM1 = Form("Relative_Mass >> Relative_Mass%i(30,0.0,1)",18);
      string RM2 = Form("Relative_Mass%i",18);
      AnalysisChain.Draw(RM1.c_str(),"(PUR_function(Number_True_Interactions))*weight");
      RelativeMass = (TH1F*)gDirectory->Get(RM2.c_str());
      gPad->Close();
      //
      string MPTNM1 = Form("Mother_PT_Normalized_Mass >> Mother_PT_Normalized_Mass%i(25,0.0,50)",18);
      string MPTNM2 = Form("Mother_PT_Normalized_Mass%i",18);
      AnalysisChain.Draw(MPTNM1.c_str(),"(PUR_function(Number_True_Interactions))*weight");
      MotherPtNormalizedMass = (TH1F*)gDirectory->Get(MPTNM2.c_str());
      gPad->Close();
      //
      string NTops1 = Form("Number_of_Tops >> Number_of_Tops%i(8,0.0,8)",18);
      string NTops2 = Form("Number_of_Tops%i",18);
      AnalysisChain.Draw(NTops1.c_str(),"(PUR_function(Number_True_Interactions))*weight");
      NumberOfTops = (TH1F*)gDirectory->Get(NTops2.c_str());
      gPad->Close();
      ///////////NEW VARIABLES//////////////////
      string HMTM1 = Form("Reconstructed_Higgs.M()/Reconstructed_Top.M() >> HMoverTM%i(30,0,1.)",18);
      string HMTM2 = Form("HMoverTM%i",18);
      AnalysisChain.Draw(HMTM1.c_str(),"(PUR_function(Number_True_Interactions))*weight");
      HiggsMassOverTopMass = (TH1F*)gDirectory->Get(HMTM2.c_str());
      gPad->Close();
      //
      string HTA1 = Form("((Reconstructed_Higgs.Pt()/Reconstructed_Higgs.M())-(Reconstructed_Top.Pt()/Reconstructed_Top.M()))/((Reconstructed_Higgs.Pt()/Reconstructed_Higgs.M())+(Reconstructed_Top.Pt()/Reconstructed_Top.M())) >> HTAsym%i(12,0,1.)",18);
      string HTA2 = Form("HTAsym%i",18);
      AnalysisChain.Draw(HTA1.c_str(),"(PUR_function(Number_True_Interactions))*weight");
      HiggsTopAsymmetry = (TH1F*)gDirectory->Get(HTA2.c_str());
      gPad->Close();
      //
      string TLBT1 = Form("Number_of_Loose_and_non_med_b_tags >> TLBTag%i(10,0,10)",18);
      string TLBT2 = Form("TLBTag%i",18);
      AnalysisChain.Draw(TLBT1.c_str(),"(PUR_function(Number_True_Interactions))*weight");
      ThirdLooseBtag = (TH1F*)gDirectory->Get(TLBT2.c_str());
      gPad->Close();
      //
      string TM1 = Form("Reconstructed_Top.M() >> TMass%i(60,100,700)",18);
      string TM2 = Form("TMass%i",18);
      AnalysisChain.Draw(TM1.c_str(),"(PUR_function(Number_True_Interactions))*weight");
      TopMass = (TH1F*)gDirectory->Get(TM2.c_str());
      gPad->Close();
      //
      string C21 = Form("ChiSquaredSorting >> ChiSq%i(100,0,1000)",18);
      string C22 = Form("ChiSq%i",18);
      AnalysisChain.Draw(C21.c_str(),"(PUR_function(Number_True_Interactions))*weight");
      Chi2 = (TH1F*)gDirectory->Get(C22.c_str());
      gPad->Close();
      //
      string UQ1 = Form("U Quark Content >> UQC%i(10,0,10)",18);
      string UQ2 = Form("UQC%i",18);
      AnalysisChain.Draw(UQ1.c_str(),"(PUR_function(Number_True_Interactions))*weight");
      UQuarkContent = (TH1F*)gDirectory->Get(UQ2.c_str());
      gPad->Close();
      //
      string DQ1 = Form("D Quark Content >> DQC%i(10,0,10)",18);
      string DQ2 = Form("DQC%i",18);
      AnalysisChain.Draw(DQ1.c_str(),"(PUR_function(Number_True_Interactions))*weight");
      DQuarkContent = (TH1F*)gDirectory->Get(DQ2.c_str());
      gPad->Close();
      //
      string SQ1 = Form("S Quark Content >> SQC%i(10,0,10)",18);
      string SQ2 = Form("SQC%i",18);
      AnalysisChain.Draw(SQ1.c_str(),"(PUR_function(Number_True_Interactions))*weight");
      SQuarkContent = (TH1F*)gDirectory->Get(SQ2.c_str());
      gPad->Close();
      //
      string CQ1 = Form("C Quark Content >> CQC%i(10,0,10)",18);
      string CQ2 = Form("CQC%i",18);
      AnalysisChain.Draw(CQ1.c_str(),"(PUR_function(Number_True_Interactions))*weight");
      CQuarkContent = (TH1F*)gDirectory->Get(CQ2.c_str());
      gPad->Close();
      //
      string BQ1 = Form("B Quark Content >> BQC%i(10,0,10)",18);
      string BQ2 = Form("BQC%i",18);
      AnalysisChain.Draw(BQ1.c_str(),"(PUR_function(Number_True_Interactions))*weight");
      BQuarkContent = (TH1F*)gDirectory->Get(BQ2.c_str());
      gPad->Close();
      //B-tagging Working point	  
      string BTL1 = Form("Number_CSVLbtagged_jets >> CSVL%i(10,0,10)",18);
      string BTL2 = Form("CSVL%i",18);
      AnalysisChain.Draw(BTL1.c_str(),"(PUR_function(Number_True_Interactions))*weight");
      CSVLB = (TH1F*)gDirectory->Get(BTL2.c_str());
      gPad->Close();
      //	  
      string BTM1 = Form("Number_CSVMbtagged_jets >> CSVM%i(10,0,10)",18);
      string BTM2 = Form("CSVM%i",18);
      AnalysisChain.Draw(BTM1.c_str(),"(PUR_function(Number_True_Interactions))*weight");
      CSVMB = (TH1F*)gDirectory->Get(BTM2.c_str());
      gPad->Close();
      //	  
      string BTT1 = Form("Number_CSVTbtagged_jets >> CSVT%i(10,0,10)",18);
      string BTT2 = Form("CSVT%i",18);
      AnalysisChain.Draw(BTT1.c_str(),"(PUR_function(Number_True_Interactions))*weight");
      CSVTB = (TH1F*)gDirectory->Get(BTT2.c_str());
      gPad->Close();
      //	  
      string VT1 = Form("Vertices >> VTX%i(40,0,40)",18);
      string VT2 = Form("VTX%i",18);
      AnalysisChain.Draw(VT1.c_str(),"(PUR_function(Number_True_Interactions))*weight");
      Vtcs = (TH1F*)gDirectory->Get(VT2.c_str());
      gPad->Close();
      //
      string HTCSVM1 = Form("THT:Number_CSVMbtagged_jets >> HT_CSVM%i(65,300,1600,10,0,10)",18);
      string HTCSVM2 = Form("HT_CSVM%i",18);
      AnalysisChain.Draw(HTCSVM1.c_str(),"(PUR_function(Number_True_Interactions))*weight");
      HTCSVM = (TH2F*)gDirectory->Get(HTCSVM2.c_str());
      gPad->Close();
      //////////////////////////////////////
      /////////ABCD BKG Estimation//////////
      //////////////////////////////////////
      HM1 = Form("Reconstructed_Higgs.M() >> HMBE0%i(36,60,180)",18);
      HM2 = Form("HMBE0%i",18);
      AnalysisChain.Draw(HM1.c_str(),"Reconstructed_Higgs.Pt()<200 || Reconstructed_Top.Pt()<200");
      HiggsMassReversedHptTpt[0] = (TH1F*)gDirectory->Get(HM2.c_str());
      gPad->Close();
      //	  
      A1 = Form("Reconstructed_Tprime.M() >> TprimeMassBE0%i(60,400,1600)",18);
      A2 = Form("TprimeMassBE0%i",18);
      AnalysisChain.Draw(A1.c_str(),"Reconstructed_Higgs.Pt()<200 || Reconstructed_Top.Pt()<200");
      TprimeMassReversedHptTpt[0] = (TH1F*)gDirectory->Get(A2.c_str());
      gPad->Close();
      //	  
      HM1 = Form("Reconstructed_Higgs.M() >> HMBE1%i(36,60,180)",18);
      HM2 = Form("HMBE1%i",18);
      AnalysisChain.Draw(HM1.c_str(),"(Reconstructed_Higgs.Pt()<200 || Reconstructed_Top.Pt()<200) && (DeltaR_of_W_Higgs>=2.2 && DeltaR_of_W_Higgs<=3.5)");
      HiggsMassReversedHptTpt[1] = (TH1F*)gDirectory->Get(HM2.c_str());
      gPad->Close();
      //	  
      A1 = Form("Reconstructed_Tprime.M() >> TprimeMassBE1%i(60,400,1600)",18);
      A2 = Form("TprimeMassBE1%i",18);
      AnalysisChain.Draw(A1.c_str(),"(Reconstructed_Higgs.Pt()<200 || Reconstructed_Top.Pt()<200) && (DeltaR_of_W_Higgs>=2.2 && DeltaR_of_W_Higgs<=3.5)");
      TprimeMassReversedHptTpt[1] = (TH1F*)gDirectory->Get(A2.c_str());
      gPad->Close();
      //	  
      HM1 = Form("Reconstructed_Higgs.M() >> HMBE2%i(36,60,180)",18);
      HM2 = Form("HMBE2%i",18);
      AnalysisChain.Draw(HM1.c_str(),"(Reconstructed_Higgs.Pt()<200 || Reconstructed_Top.Pt()<200) && (DeltaR_of_W_Higgs>=2.2 && DeltaR_of_W_Higgs<=3.5) && (TMath::Abs(First_Higgs_Jet.Phi()-Second_Higgs_Jet.Phi())<=1.2 && TMath::Abs(Top_Jet.Phi()-Reconstructed_W.Phi())<=1.2)");
      HiggsMassReversedHptTpt[2] = (TH1F*)gDirectory->Get(HM2.c_str());
      gPad->Close();
      //	  
      A1 = Form("Reconstructed_Tprime.M() >> TprimeMassBE2%i(60,400,1600)",18);
      A2 = Form("TprimeMassBE2%i",18);
      AnalysisChain.Draw(A1.c_str(),"(Reconstructed_Higgs.Pt()<200 || Reconstructed_Top.Pt()<200) && (DeltaR_of_W_Higgs>=2.2 && DeltaR_of_W_Higgs<=3.5) && (TMath::Abs(First_Higgs_Jet.Phi()-Second_Higgs_Jet.Phi())<=1.2 && TMath::Abs(Top_Jet.Phi()-Reconstructed_W.Phi())<=1.2)");
      TprimeMassReversedHptTpt[2] = (TH1F*)gDirectory->Get(A2.c_str());
      gPad->Close();
      //	  
      HM1 = Form("Reconstructed_Higgs.M() >> HMBE3%i(36,60,180)",18);
      HM2 = Form("HMBE3%i",18);
      AnalysisChain.Draw(HM1.c_str(),"(Reconstructed_Higgs.Pt()<200 || Reconstructed_Top.Pt()<200) && (DeltaR_of_W_Higgs>=2.2 && DeltaR_of_W_Higgs<=3.5) && (TMath::Abs(First_Higgs_Jet.Phi()-Second_Higgs_Jet.Phi())<=1.2 && TMath::Abs(Top_Jet.Phi()-Reconstructed_W.Phi())<=1.2) && TMath::Abs(First_W_Jet.Phi()-Second_W_Jet.Phi())<=2.0");
      HiggsMassReversedHptTpt[3] = (TH1F*)gDirectory->Get(HM2.c_str());
      gPad->Close();
      //	  
      A1 = Form("Reconstructed_Tprime.M() >> TprimeMassBE3%i(60,400,1600)",18);
      A2 = Form("TprimeMassBE3%i",18);
      AnalysisChain.Draw(A1.c_str(),"(Reconstructed_Higgs.Pt()<200 || Reconstructed_Top.Pt()<200) && (DeltaR_of_W_Higgs>=2.2 && DeltaR_of_W_Higgs<=3.5) && (TMath::Abs(First_Higgs_Jet.Phi()-Second_Higgs_Jet.Phi())<=1.2 && TMath::Abs(Top_Jet.Phi()-Reconstructed_W.Phi())<=1.2) && TMath::Abs(First_W_Jet.Phi()-Second_W_Jet.Phi())<=2.0");
      TprimeMassReversedHptTpt[3] = (TH1F*)gDirectory->Get(A2.c_str());
      gPad->Close();
      //	  
      HM1 = Form("Reconstructed_Higgs.M() >> HMBE4%i(36,60,180)",18);
      HM2 = Form("HMBE4%i",18);
      AnalysisChain.Draw(HM1.c_str(),"(Reconstructed_Higgs.Pt()<200 || Reconstructed_Top.Pt()<200) && (DeltaR_of_W_Higgs>=2.2 && DeltaR_of_W_Higgs<=3.5) && (TMath::Abs(First_Higgs_Jet.Phi()-Second_Higgs_Jet.Phi())<=1.2 && TMath::Abs(Top_Jet.Phi()-Reconstructed_W.Phi())<=1.2) && TMath::Abs(First_W_Jet.Phi()-Second_W_Jet.Phi())<=2.0 && Relative_THT>=0.65");
      HiggsMassReversedHptTpt[4] = (TH1F*)gDirectory->Get(HM2.c_str());
      gPad->Close();
      //	  
      A1 = Form("Reconstructed_Tprime.M() >> TprimeMassBE4%i(60,400,1600)",18);
      A2 = Form("TprimeMassBE4%i",18);
      AnalysisChain.Draw(A1.c_str(),"(Reconstructed_Higgs.Pt()<200 || Reconstructed_Top.Pt()<200) && (DeltaR_of_W_Higgs>=2.2 && DeltaR_of_W_Higgs<=3.5) && (TMath::Abs(First_Higgs_Jet.Phi()-Second_Higgs_Jet.Phi())<=1.2 && TMath::Abs(Top_Jet.Phi()-Reconstructed_W.Phi())<=1.2) && TMath::Abs(First_W_Jet.Phi()-Second_W_Jet.Phi())<=2.0 && Relative_THT>=0.65");
      TprimeMassReversedHptTpt[4] = (TH1F*)gDirectory->Get(A2.c_str());
      //	  
      HM1 = Form("Reconstructed_Higgs.M() >> HMBE5%i(36,60,180)",18);
      HM2 = Form("HMBE5%i",18);
      AnalysisChain.Draw(HM1.c_str(),"(Reconstructed_Higgs.Pt()<200 || Reconstructed_Top.Pt()<200) && (DeltaR_of_W_Higgs>=2.2 && DeltaR_of_W_Higgs<=3.5) && (TMath::Abs(First_Higgs_Jet.Phi()-Second_Higgs_Jet.Phi())<=1.2 && TMath::Abs(Top_Jet.Phi()-Reconstructed_W.Phi())<=1.2) && TMath::Abs(First_W_Jet.Phi()-Second_W_Jet.Phi())<=2.0 && Relative_THT>=0.65 && (Relative_Mass>=0.3 && Relative_Mass<=0.5)");
      HiggsMassReversedHptTpt[5] = (TH1F*)gDirectory->Get(HM2.c_str());
      gPad->Close();
      //	  
      A1 = Form("Reconstructed_Tprime.M() >> TprimeMassBE5%i(60,400,1600)",18);
      A2 = Form("TprimeMassBE5%i",18);
      AnalysisChain.Draw(A1.c_str(),"(Reconstructed_Higgs.Pt()<200 || Reconstructed_Top.Pt()<200) && (DeltaR_of_W_Higgs>=2.2 && DeltaR_of_W_Higgs<=3.5) && (TMath::Abs(First_Higgs_Jet.Phi()-Second_Higgs_Jet.Phi())<=1.2 && TMath::Abs(Top_Jet.Phi()-Reconstructed_W.Phi())<=1.2) && TMath::Abs(First_W_Jet.Phi()-Second_W_Jet.Phi())<=2.0 && Relative_THT>=0.65 && (Relative_Mass>=0.3 && Relative_Mass<=0.5)");
      TprimeMassReversedHptTpt[5] = (TH1F*)gDirectory->Get(A2.c_str());
      gPad->Close();
      //
      string HPTTPT1 = Form("Reconstructed_Higgs.Pt():Reconstructed_Top.Pt() >> HPtTPt0%i(40,10,800,40,10,800)",18);
      string HPTTPT2 = Form("HPtTPt0%i",18);
      AnalysisChain.Draw(HPTTPT1.c_str());
      HptTpt[0] = (TH2F*)gDirectory->Get(HPTTPT2.c_str());
      cout << "Correlation Factor: " << HptTpt[0]->GetCorrelationFactor() << endl;
      gPad->Close();
      //
      HPTTPT1 = Form("Reconstructed_Higgs.Pt():Reconstructed_Top.Pt() >> HPtTPt1%i(40,10,800,40,10,800)",18);
      HPTTPT2 = Form("HPtTPt1%i",18);
      AnalysisChain.Draw(HPTTPT1.c_str(),"(DeltaR_of_W_Higgs>=2.2 && DeltaR_of_W_Higgs<=3.5)");
      HptTpt[1] = (TH2F*)gDirectory->Get(HPTTPT2.c_str());
      cout << "Correlation Factor: " << HptTpt[1]->GetCorrelationFactor() << endl;
      gPad->Close();
      //
      HPTTPT1 = Form("Reconstructed_Higgs.Pt():Reconstructed_Top.Pt() >> HPtTPt2%i(40,10,800,40,10,800)",18);
      HPTTPT2 = Form("HPtTPt2%i",18);
      AnalysisChain.Draw(HPTTPT1.c_str(),"(DeltaR_of_W_Higgs>=2.2 && DeltaR_of_W_Higgs<=3.5) && (TMath::Abs(First_Higgs_Jet.Phi()-Second_Higgs_Jet.Phi())<=1.2 && TMath::Abs(Top_Jet.Phi()-Reconstructed_W.Phi())<=1.2)");
      HptTpt[2] = (TH2F*)gDirectory->Get(HPTTPT2.c_str());
      cout << "Correlation Factor: " << HptTpt[2]->GetCorrelationFactor() << endl;
      gPad->Close();
      //
      HPTTPT1 = Form("Reconstructed_Higgs.Pt():Reconstructed_Top.Pt() >> HPtTPt3%i(40,10,800,40,10,800)",18);
      HPTTPT2 = Form("HPtTPt3%i",18);
      AnalysisChain.Draw(HPTTPT1.c_str(),"(DeltaR_of_W_Higgs>=2.2 && DeltaR_of_W_Higgs<=3.5) && (TMath::Abs(First_Higgs_Jet.Phi()-Second_Higgs_Jet.Phi())<=1.2 && TMath::Abs(Top_Jet.Phi()-Reconstructed_W.Phi())<=1.2) && TMath::Abs(First_W_Jet.Phi()-Second_W_Jet.Phi())<=2.0");
      HptTpt[3] = (TH2F*)gDirectory->Get(HPTTPT2.c_str());
      cout << "Correlation Factor: " << HptTpt[3]->GetCorrelationFactor() << endl;
      gPad->Close();
      //
      HPTTPT1 = Form("Reconstructed_Higgs.Pt():Reconstructed_Top.Pt() >> HPtTPt4%i(40,10,800,40,10,800)",18);
      HPTTPT2 = Form("HPtTPt4%i",18);
      AnalysisChain.Draw(HPTTPT1.c_str(),"(DeltaR_of_W_Higgs>=2.2 && DeltaR_of_W_Higgs<=3.5) && (TMath::Abs(First_Higgs_Jet.Phi()-Second_Higgs_Jet.Phi())<=1.2 && TMath::Abs(Top_Jet.Phi()-Reconstructed_W.Phi())<=1.2) && TMath::Abs(First_W_Jet.Phi()-Second_W_Jet.Phi())<=2.0 && Relative_THT>=0.65");
      HptTpt[4] = (TH2F*)gDirectory->Get(HPTTPT2.c_str());
      cout << "Correlation Factor: " << HptTpt[4]->GetCorrelationFactor() << endl;
      gPad->Close();
      //
      HPTTPT1 = Form("Reconstructed_Higgs.Pt():Reconstructed_Top.Pt() >> HPtTPt5%i(40,10,800,40,10,800)",18);
      HPTTPT2 = Form("HPtTPt5%i",18);
      AnalysisChain.Draw(HPTTPT1.c_str(),"(DeltaR_of_W_Higgs>=2.2 && DeltaR_of_W_Higgs<=3.5) && (TMath::Abs(First_Higgs_Jet.Phi()-Second_Higgs_Jet.Phi())<=1.2 && TMath::Abs(Top_Jet.Phi()-Reconstructed_W.Phi())<=1.2) && TMath::Abs(First_W_Jet.Phi()-Second_W_Jet.Phi())<=2.0 && Relative_THT>=0.65 && (Relative_Mass>=0.3 && Relative_Mass<=0.5)");
      HptTpt[5] = (TH2F*)gDirectory->Get(HPTTPT2.c_str());
      cout << "Correlation Factor: " << HptTpt[5]->GetCorrelationFactor() << endl;
      gPad->Close();
      //
      HM1 = Form("Reconstructed_Higgs.M() >> HM0%i(36,60,180)",18);
      HM2 = Form("HM0%i",18);
      AnalysisChain.Draw(HM1.c_str(),"Reconstructed_Higgs.Pt()>200 && Reconstructed_Top.Pt()>200");
      HiggsMassCuts[0] = (TH1F*)gDirectory->Get(HM2.c_str());
      gPad->Close();
      //	  
      HM1 = Form("Reconstructed_Higgs.M() >> HM1%i(36,60,180)",18);
      HM2 = Form("HM1%i",18);
      AnalysisChain.Draw(HM1.c_str(),"(Reconstructed_Higgs.Pt()>200 && Reconstructed_Top.Pt()>200) && (DeltaR_of_W_Higgs>=2.2 && DeltaR_of_W_Higgs<=3.5)");
      HiggsMassCuts[1] = (TH1F*)gDirectory->Get(HM2.c_str());
      gPad->Close();
      //	  
      HM1 = Form("Reconstructed_Higgs.M() >> HM2%i(36,60,180)",18);
      HM2 = Form("HM2%i",18);
      AnalysisChain.Draw(HM1.c_str(),"(Reconstructed_Higgs.Pt()>200 && Reconstructed_Top.Pt()>200) && (DeltaR_of_W_Higgs>=2.2 && DeltaR_of_W_Higgs<=3.5) && (TMath::Abs(First_Higgs_Jet.Phi()-Second_Higgs_Jet.Phi())<=1.2 && TMath::Abs(Top_Jet.Phi()-Reconstructed_W.Phi())<=1.2)");
      HiggsMassCuts[2] = (TH1F*)gDirectory->Get(HM2.c_str());
      gPad->Close();
      //	  
      HM1 = Form("Reconstructed_Higgs.M() >> HM3%i(36,60,180)",18);
      HM2 = Form("HM3%i",18);
      AnalysisChain.Draw(HM1.c_str(),"(Reconstructed_Higgs.Pt()>200 && Reconstructed_Top.Pt()>200) && (DeltaR_of_W_Higgs>=2.2 && DeltaR_of_W_Higgs<=3.5) && (TMath::Abs(First_Higgs_Jet.Phi()-Second_Higgs_Jet.Phi())<=1.2 && TMath::Abs(Top_Jet.Phi()-Reconstructed_W.Phi())<=1.2) && TMath::Abs(First_W_Jet.Phi()-Second_W_Jet.Phi())<=2.0");
      HiggsMassCuts[3] = (TH1F*)gDirectory->Get(HM2.c_str());
      gPad->Close();
      //
      HM1 = Form("Reconstructed_Higgs.M() >> HM4%i(36,60,180)",18);
      HM2 = Form("HM4%i",18);
      AnalysisChain.Draw(HM1.c_str(),"(Reconstructed_Higgs.Pt()>200 && Reconstructed_Top.Pt()>200) && (DeltaR_of_W_Higgs>=2.2 && DeltaR_of_W_Higgs<=3.5) && (TMath::Abs(First_Higgs_Jet.Phi()-Second_Higgs_Jet.Phi())<=1.2 && TMath::Abs(Top_Jet.Phi()-Reconstructed_W.Phi())<=1.2) && TMath::Abs(First_W_Jet.Phi()-Second_W_Jet.Phi())<=2.0 && Relative_THT>=0.65");
      HiggsMassCuts[4] = (TH1F*)gDirectory->Get(HM2.c_str());
      gPad->Close();
      //
      HM1 = Form("Reconstructed_Higgs.M() >> HM5%i(36,60,180)",18);
      HM2 = Form("HM5%i",18);
      AnalysisChain.Draw(HM1.c_str(),"(Reconstructed_Higgs.Pt()>200 && Reconstructed_Top.Pt()>200) && (DeltaR_of_W_Higgs>=2.2 && DeltaR_of_W_Higgs<=3.5) && (TMath::Abs(First_Higgs_Jet.Phi()-Second_Higgs_Jet.Phi())<=1.2 && TMath::Abs(Top_Jet.Phi()-Reconstructed_W.Phi())<=1.2) && TMath::Abs(First_W_Jet.Phi()-Second_W_Jet.Phi())<=2.0 && Relative_THT>=0.65 && (Relative_Mass>=0.3 && Relative_Mass<=0.5)");
      HiggsMassCuts[5] = (TH1F*)gDirectory->Get(HM2.c_str());
      gPad->Close();
      ///////////////////////
      //Sideband estimation//
      ///////////////////////
      string WFH1 = Form("W_From_Higgs.M() >> WMFH0%i(44,60.0,500)",18);
      string WFH2 = Form("WMFH0%i",18);
      AnalysisChain.Draw(WFH1.c_str(),"Reconstructed_Higgs.Pt()>200 && Reconstructed_Top.Pt()>200");
      WMassFromHiggs[0] = (TH1F*)gDirectory->Get(WFH2.c_str());
      gPad->Close();
      //
      WFH1 = Form("W_From_Higgs.M() >> WMFH1%i(44,60.0,500)",18);
      WFH2 = Form("WMFH1%i",18);
      AnalysisChain.Draw(WFH1.c_str(),"(Reconstructed_Higgs.Pt()>200 && Reconstructed_Top.Pt()>200) && (DeltaR_of_W_Higgs>=2.2 && DeltaR_of_W_Higgs<=3.5)");
      WMassFromHiggs[1] = (TH1F*)gDirectory->Get(WFH2.c_str());
      gPad->Close();
      //
      WFH1 = Form("W_From_Higgs.M() >> WMFH2%i(44,60.0,500)",18);
      WFH2 = Form("WMFH2%i",18);
      AnalysisChain.Draw(WFH1.c_str(),"(Reconstructed_Higgs.Pt()>200 && Reconstructed_Top.Pt()>200) && (DeltaR_of_W_Higgs>=2.2 && DeltaR_of_W_Higgs<=3.5) && (TMath::Abs(First_Higgs_Jet.Phi()-Second_Higgs_Jet.Phi())<=1.2 && TMath::Abs(Top_Jet.Phi()-Reconstructed_W.Phi())<=1.2)");
      WMassFromHiggs[2] = (TH1F*)gDirectory->Get(WFH2.c_str());
      gPad->Close();
      //
      WFH1 = Form("W_From_Higgs.M() >> WMFH3%i(44,60.0,500)",18);
      WFH2 = Form("WMFH3%i",18);
      AnalysisChain.Draw(WFH1.c_str(),"(Reconstructed_Higgs.Pt()>200 && Reconstructed_Top.Pt()>200) && (DeltaR_of_W_Higgs>=2.2 && DeltaR_of_W_Higgs<=3.5) && (TMath::Abs(First_Higgs_Jet.Phi()-Second_Higgs_Jet.Phi())<=1.2 && TMath::Abs(Top_Jet.Phi()-Reconstructed_W.Phi())<=1.2) && TMath::Abs(First_W_Jet.Phi()-Second_W_Jet.Phi())<=2.0");
      WMassFromHiggs[3] = (TH1F*)gDirectory->Get(WFH2.c_str());
      gPad->Close();
      //
      WFH1 = Form("W_From_Higgs.M() >> WMFH4%i(44,60.0,500)",18);
      WFH2 = Form("WMFH4%i",18);
      AnalysisChain.Draw(WFH1.c_str(),"(Reconstructed_Higgs.Pt()>200 && Reconstructed_Top.Pt()>200) && (DeltaR_of_W_Higgs>=2.2 && DeltaR_of_W_Higgs<=3.5) && (TMath::Abs(First_Higgs_Jet.Phi()-Second_Higgs_Jet.Phi())<=1.2 && TMath::Abs(Top_Jet.Phi()-Reconstructed_W.Phi())<=1.2) && TMath::Abs(First_W_Jet.Phi()-Second_W_Jet.Phi())<=2.0 && Relative_THT>=0.65");
      WMassFromHiggs[4] = (TH1F*)gDirectory->Get(WFH2.c_str());
      gPad->Close();
      //
      WFH1 = Form("W_From_Higgs.M() >> WMFH5%i(44,60.0,500)",18);
      WFH2 = Form("WMFH5%i",18);
      AnalysisChain.Draw(WFH1.c_str(),"(Reconstructed_Higgs.Pt()>200 && Reconstructed_Top.Pt()>200) && (DeltaR_of_W_Higgs>=2.2 && DeltaR_of_W_Higgs<=3.5) && (TMath::Abs(First_Higgs_Jet.Phi()-Second_Higgs_Jet.Phi())<=1.2 && TMath::Abs(Top_Jet.Phi()-Reconstructed_W.Phi())<=1.2) && TMath::Abs(First_W_Jet.Phi()-Second_W_Jet.Phi())<=2.0 && Relative_THT>=0.65 && (Relative_Mass>=0.3 && Relative_Mass<=0.5)");
      WMassFromHiggs[5] = (TH1F*)gDirectory->Get(WFH2.c_str());
      gPad->Close();
      //Chi2////////////////
      WFH1 = Form("W_From_Higgs_Chi2.M() >> WMFHC0%i(44,60.0,500)",18);
      WFH2 = Form("WMFHC0%i",18);
      AnalysisChain.Draw(WFH1.c_str(),"Reconstructed_Higgs.Pt()>200 && Reconstructed_Top.Pt()>200");
      WMassFromHiggsChi2[0] = (TH1F*)gDirectory->Get(WFH2.c_str());
      gPad->Close();
      //
      WFH1 = Form("W_From_Higgs_Chi2.M() >> WMFHC1%i(44,60.0,500)",18);
      WFH2 = Form("WMFHC1%i",18);
      AnalysisChain.Draw(WFH1.c_str(),"(Reconstructed_Higgs.Pt()>200 && Reconstructed_Top.Pt()>200) && (DeltaR_of_W_Higgs>=2.2 && DeltaR_of_W_Higgs<=3.5)");
      WMassFromHiggsChi2[1] = (TH1F*)gDirectory->Get(WFH2.c_str());
      gPad->Close();
      //
      WFH1 = Form("W_From_Higgs_Chi2.M() >> WMFHC2%i(44,60.0,500)",18);
      WFH2 = Form("WMFHC2%i",18);
      AnalysisChain.Draw(WFH1.c_str(),"(Reconstructed_Higgs.Pt()>200 && Reconstructed_Top.Pt()>200) && (DeltaR_of_W_Higgs>=2.2 && DeltaR_of_W_Higgs<=3.5) && (TMath::Abs(First_Higgs_Jet.Phi()-Second_Higgs_Jet.Phi())<=1.2 && TMath::Abs(Top_Jet.Phi()-Reconstructed_W.Phi())<=1.2)");
      WMassFromHiggsChi2[2] = (TH1F*)gDirectory->Get(WFH2.c_str());
      gPad->Close();
      //
      WFH1 = Form("W_From_Higgs_Chi2.M() >> WMFHC3%i(44,60.0,500)",18);
      WFH2 = Form("WMFHC3%i",18);
      AnalysisChain.Draw(WFH1.c_str(),"(Reconstructed_Higgs.Pt()>200 && Reconstructed_Top.Pt()>200) && (DeltaR_of_W_Higgs>=2.2 && DeltaR_of_W_Higgs<=3.5) && (TMath::Abs(First_Higgs_Jet.Phi()-Second_Higgs_Jet.Phi())<=1.2 && TMath::Abs(Top_Jet.Phi()-Reconstructed_W.Phi())<=1.2) && TMath::Abs(First_W_Jet.Phi()-Second_W_Jet.Phi())<=2.0");
      WMassFromHiggsChi2[3] = (TH1F*)gDirectory->Get(WFH2.c_str());
      gPad->Close();
      //
      WFH1 = Form("W_From_Higgs_Chi2.M() >> WMFHC4%i(44,60.0,500)",18);
      WFH2 = Form("WMFHC4%i",18);
      AnalysisChain.Draw(WFH1.c_str(),"(Reconstructed_Higgs.Pt()>200 && Reconstructed_Top.Pt()>200) && (DeltaR_of_W_Higgs>=2.2 && DeltaR_of_W_Higgs<=3.5) && (TMath::Abs(First_Higgs_Jet.Phi()-Second_Higgs_Jet.Phi())<=1.2 && TMath::Abs(Top_Jet.Phi()-Reconstructed_W.Phi())<=1.2) && TMath::Abs(First_W_Jet.Phi()-Second_W_Jet.Phi())<=2.0 && Relative_THT>=0.65");
      WMassFromHiggsChi2[4] = (TH1F*)gDirectory->Get(WFH2.c_str());
      gPad->Close();
      //
      WFH1 = Form("W_From_Higgs_Chi2.M() >> WMFHC5%i(44,60.0,500)",18);
      WFH2 = Form("WMFHC5%i",18);
      AnalysisChain.Draw(WFH1.c_str(),"(Reconstructed_Higgs.Pt()>200 && Reconstructed_Top.Pt()>200) && (DeltaR_of_W_Higgs>=2.2 && DeltaR_of_W_Higgs<=3.5) && (TMath::Abs(First_Higgs_Jet.Phi()-Second_Higgs_Jet.Phi())<=1.2 && TMath::Abs(Top_Jet.Phi()-Reconstructed_W.Phi())<=1.2) && TMath::Abs(First_W_Jet.Phi()-Second_W_Jet.Phi())<=2.0 && Relative_THT>=0.65 && (Relative_Mass>=0.3 && Relative_Mass<=0.5)");
      WMassFromHiggsChi2[5] = (TH1F*)gDirectory->Get(WFH2.c_str());
      gPad->Close();
      /////////////////////////////////////////////////////////////
      string TPNWH1 = Form("Reconstructed_Tprime.M() >> TPNWFH%i(60,400,1600)",18);
      string TPNWH2 = Form("TPNWFH%i",18);
      AnalysisChain.Draw(TPNWH1.c_str(),"(Reconstructed_Higgs.Pt()>200 && Reconstructed_Top.Pt()>200) && (DeltaR_of_W_Higgs>=2.2 && DeltaR_of_W_Higgs<=3.5) && (TMath::Abs(First_Higgs_Jet.Phi()-Second_Higgs_Jet.Phi())<=1.2 && TMath::Abs(Top_Jet.Phi()-Reconstructed_W.Phi())<=1.2) && TMath::Abs(First_W_Jet.Phi()-Second_W_Jet.Phi())<=2.0 && Relative_THT>=0.65 && (Relative_Mass>=0.3 && Relative_Mass<=0.5) && (W_From_Higgs.M()<70 || W_From_Higgs.M()>110)");
      TprimeNotWFromHiggs = (TH1F*)gDirectory->Get(TPNWH2.c_str());
      gPad->Close();
      //
      string D1 = Form("ABs(First_W_Jet.Phi(),Second_W_Jet.Phi()):THT >> DPWJ_HT_%i(65,300,1600,60,0.5,3.5)",18);
      string D2 = Form("DPWJ_HT_%i",18);
      AnalysisChain.Draw(D1.c_str()); //,"(Reconstructed_Higgs.Pt()<200 || Reconstructed_Top.Pt()<200) && (DeltaR_of_W_Higgs>=2.2 && DeltaR_of_W_Higgs<=3.5) && (TMath::Abs(First_Higgs_Jet.Phi()-Second_Higgs_Jet.Phi())<=1.2 && TMath::Abs(Top_Jet.Phi()-Reconstructed_W.Phi())<=1.2) && Relative_THT>=0.65 && (Relative_Mass>=0.3 && Relative_Mass<=0.5) && (Top_From_Higgs_Chi2.M()<140 || Top_From_Higgs_Chi2.M()>230) && (Reconstructed_Higgs.M()>=110 && Reconstructed_Higgs.M()<=140)");
      DPWJ_HT = (TH2F*)gDirectory->Get(D2.c_str());
      cout << "Correlation Factor DPWJ and HT: " << DPWJ_HT->GetCorrelationFactor() << " Integral->" << DPWJ_HT->Integral() << endl;
      gPad->Close();
      //////////////////////////////////////////////////////////////////////////////
      //////////////////////////////////////////////////////////////////////////////
      //////////////////////////////////////////////////////////////////////////////
      A1 = Form("Reconstructed_Tprime3L.M() >> TprimeMass3L%i(60,400,1600)",18);
      A2 = Form("TprimeMass3L%i",18);
      AnalysisChain3L.Draw(A1.c_str());
      FiveJetsMass3L = (TH1F*)gDirectory->Get(A2.c_str());
      gPad->Close();
      /////////////////////
      ///////////////////// ADDTIONAL INFO ON JETS
      /////////////////////
      A1 = Form("JET1.Eta() >> jet1_eta%i(100,-5,5)",18);
      A2 = Form("jet1_eta%i",18);
      AnalysisChain.Draw(A1.c_str(),finalcut);
      JET1ETA = (TH1F*)gDirectory->Get(A2.c_str());
      gPad->Close();
      //
      A1 = Form("JET2.Eta() >> jet2_eta%i(100,-5,5)",18);
      A2 = Form("jet2_eta%i",18);
      AnalysisChain.Draw(A1.c_str(),finalcut);
      JET2ETA = (TH1F*)gDirectory->Get(A2.c_str());
      gPad->Close();
      //
      A1 = Form("JET3.Eta() >> jet3_eta%i(100,-5,5)",18);
      A2 = Form("jet3_eta%i",18);
      AnalysisChain.Draw(A1.c_str(),finalcut);
      JET3ETA = (TH1F*)gDirectory->Get(A2.c_str());
      gPad->Close();
      //
      A1 = Form("JET4.Eta() >> jet4_eta%i(100,-5,5)",18);
      A2 = Form("jet4_eta%i",18);
      AnalysisChain.Draw(A1.c_str(),finalcut);
      JET4ETA = (TH1F*)gDirectory->Get(A2.c_str());
      gPad->Close();
      //
      A1 = Form("JET5.Eta() >> jet5_eta%i(100,-5,5)",18);
      A2 = Form("jet5_eta%i",18);
      AnalysisChain.Draw(A1.c_str(),finalcut);
      JET5ETA = (TH1F*)gDirectory->Get(A2.c_str());
      gPad->Close();
      //
      A1 = Form("JET6.Eta() >> jet6_eta%i(100,-5,5)",18);
      A2 = Form("jet6_eta%i",18);
      AnalysisChain.Draw(A1.c_str(),finalcut);
      JET6ETA = (TH1F*)gDirectory->Get(A2.c_str());
      gPad->Close();
      //
      A1 = Form("JET1.Phi() >> jet1_phi%i(60,-3,3)",18);
      A2 = Form("jet1_phi%i",18);
      AnalysisChain.Draw(A1.c_str(),finalcut);
      JET1PHI = (TH1F*)gDirectory->Get(A2.c_str());
      gPad->Close();
      //
      A1 = Form("JET2.Phi() >> jet2_phi%i(60,-3,3)",18);
      A2 = Form("jet2_phi%i",18);
      AnalysisChain.Draw(A1.c_str(),finalcut);
      JET2PHI = (TH1F*)gDirectory->Get(A2.c_str());
      gPad->Close();
      //
      A1 = Form("JET3.Phi() >> jet3_phi%i(60,-3,3)",18);
      A2 = Form("jet3_phi%i",18);
      AnalysisChain.Draw(A1.c_str(),finalcut);
      JET3PHI = (TH1F*)gDirectory->Get(A2.c_str());
      gPad->Close();
      //
      A1 = Form("JET4.Phi() >> jet4_phi%i(60,-3,3)",18);
      A2 = Form("jet4_phi%i",18);
      AnalysisChain.Draw(A1.c_str(),finalcut);
      JET4PHI = (TH1F*)gDirectory->Get(A2.c_str());
      gPad->Close();
      //
      A1 = Form("JET5.Phi() >> jet5_phi%i(60,-3,3)",18);
      A2 = Form("jet5_phi%i",18);
      AnalysisChain.Draw(A1.c_str(),finalcut);
      JET5PHI = (TH1F*)gDirectory->Get(A2.c_str());
      gPad->Close();
      //
      A1 = Form("JET6.Phi() >> jet6_phi%i(60,-3,3)",18);
      A2 = Form("jet6_phi%i",18);
      AnalysisChain.Draw(A1.c_str(),finalcut);
      JET6PHI = (TH1F*)gDirectory->Get(A2.c_str());
      gPad->Close();
      //
      A1 = Form("Number_Jets >> Num_jets%i(20,0,20)",18);
      A2 = Form("Num_jets%i",18);
      AnalysisChain.Draw(A1.c_str(),finalcut);
      JETMULTI = (TH1F*)gDirectory->Get(A2.c_str());
      gPad->Close();
    }
  
  FiveJetsMass->Scale(Lumi*Xs/ProcessedEvents[18]);
  FiveJetsMass3L->Scale(Lumi*Xs/ProcessedEvents[18]);
  LeadingJetPT->Scale(Lumi*Xs/ProcessedEvents[18]);
  Leading2JetPT->Scale(Lumi*Xs/ProcessedEvents[18]);
  Leading3JetPT->Scale(Lumi*Xs/ProcessedEvents[18]);
  Leading4JetPT->Scale(Lumi*Xs/ProcessedEvents[18]);
  Leading5JetPT->Scale(Lumi*Xs/ProcessedEvents[18]);
  Leading6JetPT->Scale(Lumi*Xs/ProcessedEvents[18]);
  THT->Scale(Lumi*Xs/ProcessedEvents[18]);
  DRHjets->Scale(Lumi*Xs/ProcessedEvents[18]);
  DRWjets->Scale(Lumi*Xs/ProcessedEvents[18]);
  Hpt->Scale(Lumi*Xs/ProcessedEvents[18]);
  Tpt->Scale(Lumi*Xs/ProcessedEvents[18]);
  DRWH->Scale(Lumi*Xs/ProcessedEvents[18]);
  DPHjets->Scale(Lumi*Xs/ProcessedEvents[18]);
  DPWjets->Scale(Lumi*Xs/ProcessedEvents[18]);
  DPTjets->Scale(Lumi*Xs/ProcessedEvents[18]);
  HiggsMass->Scale(Lumi*Xs/ProcessedEvents[18]);
  RelHT->Scale(Lumi*Xs/ProcessedEvents[18]);
  DRTH->Scale(Lumi*Xs/ProcessedEvents[18]);
  PtNormalizedMass->Scale(Lumi*Xs/ProcessedEvents[18]);
  RelativeMass->Scale(Lumi*Xs/ProcessedEvents[18]);
  MotherPtNormalizedMass->Scale(Lumi*Xs/ProcessedEvents[18]);
  NumberOfTops->Scale(Lumi*Xs/ProcessedEvents[18]);
  HiggsMassOverTopMass->Scale(Lumi*Xs/ProcessedEvents[18]);
  HiggsTopAsymmetry->Scale(Lumi*Xs/ProcessedEvents[18]);
  ThirdLooseBtag->Scale(Lumi*Xs/ProcessedEvents[18]);
  TopMass->Scale(Lumi*Xs/ProcessedEvents[18]);
  Chi2->Scale(Lumi*Xs/ProcessedEvents[18]);
  UQuarkContent->Scale(Lumi*Xs/ProcessedEvents[18]);
  DQuarkContent->Scale(Lumi*Xs/ProcessedEvents[18]);
  SQuarkContent->Scale(Lumi*Xs/ProcessedEvents[18]);
  CQuarkContent->Scale(Lumi*Xs/ProcessedEvents[18]);
  BQuarkContent->Scale(Lumi*Xs/ProcessedEvents[18]);
  CSVLB->Scale(Lumi*Xs/ProcessedEvents[18]);
  CSVMB->Scale(Lumi*Xs/ProcessedEvents[18]);
  CSVTB->Scale(Lumi*Xs/ProcessedEvents[18]);
  Vtcs->Scale(Lumi*Xs/ProcessedEvents[18]);
  HTCSVM->Scale(Lumi*Xs/ProcessedEvents[18]);
  for (int i=0; i<6; i++) HptTpt[i]->Scale(Lumi*Xs/ProcessedEvents[18]); 
  for (int i=0; i<6; i++) HiggsMassReversedHptTpt[i]->Scale(Lumi*Xs/ProcessedEvents[18]);
  for (int i=0; i<6; i++) TprimeMassReversedHptTpt[i]->Scale(Lumi*Xs/ProcessedEvents[18]);
  for (int i=0; i<6; i++) HiggsMassCuts[i]->Scale(Lumi*Xs/ProcessedEvents[18]);
  for (int i=0; i<6; i++) WMassFromHiggs[i]->Scale(Lumi*Xs/ProcessedEvents[18]);
  for (int i=0; i<6; i++) WMassFromHiggsChi2[i]->Scale(Lumi*Xs/ProcessedEvents[18]);
  TprimeNotWFromHiggs->Scale(Lumi*Xs/ProcessedEvents[18]);
  FiveJetsMassBE->Scale(Lumi*Xs/ProcessedEvents[18]); FiveJetsMassBE->Scale((TopMassFromHiggs->Integral(TopMassFromHiggs->GetXaxis()->FindBin(0.0),TopMassFromHiggs->GetXaxis()->FindBin(140.0))+TopMassFromHiggs->Integral(TopMassFromHiggs->GetXaxis()->FindBin(230.0),TopMassFromHiggs->GetXaxis()->FindBin(10000.0)))/TopMassFromHiggs->Integral(TopMassFromHiggs->GetXaxis()->FindBin(140.0),TopMassFromHiggs->GetXaxis()->FindBin(230.0)));
  FiveJetsMassLC->Scale(Lumi*Xs/ProcessedEvents[18]);
  TopMassFromHiggs->Scale(Lumi*Xs/ProcessedEvents[18]);
  FiveJetsMassLCoverBE->Scale(Lumi*Xs/ProcessedEvents[18]); FiveJetsMassLCoverBE->Divide(FiveJetsMassBE);
  HMBE->Scale(Lumi*Xs/ProcessedEvents[18]);
  HMLC->Scale(Lumi*Xs/ProcessedEvents[18]);
  HMLCoverBE->Scale(Lumi*Xs/ProcessedEvents[18]); HMLCoverBE->Scale(HMBE->Integral()/HMLCoverBE->Integral()); HMLCoverBE->Divide(HMBE);
  FJLC_QCDBE->Scale(Lumi*Xs/ProcessedEvents[18]);
  FJLC_QCDLC->Scale(Lumi*Xs/ProcessedEvents[18]);
  FJLC_QCDLCoverBE->Scale(Lumi*Xs/ProcessedEvents[18]); FJLC_QCDLCoverBE->Scale(FJLC_QCDBE->Integral()/FJLC_QCDLCoverBE->Integral()); FJLC_QCDLCoverBE->Divide(FJLC_QCDBE);

  TM_LC_HMW->Scale(Lumi*Xs/ProcessedEvents[18]);
  DPWJ_HT->Scale(Lumi*Xs/ProcessedEvents[18]);

  JET1ETA->Scale(Lumi*Xs/ProcessedEvents[18]);
  JET2ETA->Scale(Lumi*Xs/ProcessedEvents[18]);
  JET3ETA->Scale(Lumi*Xs/ProcessedEvents[18]);
  JET4ETA->Scale(Lumi*Xs/ProcessedEvents[18]);
  JET5ETA->Scale(Lumi*Xs/ProcessedEvents[18]);
  JET6ETA->Scale(Lumi*Xs/ProcessedEvents[18]);
  JET1PHI->Scale(Lumi*Xs/ProcessedEvents[18]);
  JET2PHI->Scale(Lumi*Xs/ProcessedEvents[18]);
  JET3PHI->Scale(Lumi*Xs/ProcessedEvents[18]);
  JET4PHI->Scale(Lumi*Xs/ProcessedEvents[18]);
  JET5PHI->Scale(Lumi*Xs/ProcessedEvents[18]);
  JET6PHI->Scale(Lumi*Xs/ProcessedEvents[18]);
  JETMULTI->Scale(Lumi*Xs/ProcessedEvents[18]);
  
  //Settings for signal
  if (SurvivalMarker)
    {
      //TprimeMass
      FiveJetsMass->SetFillColor(kSpring);
      FiveJetsMass->SetFillStyle(3444);
      FiveJetsMass->SetLineWidth(3);
      FiveJetsMass3L->SetFillColor(kSpring);
      FiveJetsMass3L->SetFillStyle(3444);
      FiveJetsMass3L->SetLineWidth(3);
      TFile f("Signal.root", "RECREATE");
      FiveJetsMass->Write();
      FiveJetsMass3L->Write();
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
      //HptTpt->SetFillColor(kSpring);
      //HptTpt->SetFillStyle(3444);
      //HptTpt->SetLineWidth(3);
      for (int i=0; i<6; i++) cout << "Correlation Factor: " << HptTpt[i]->GetCorrelationFactor() << endl;
      for (int i=0; i<6; i++) HptTpt[i]->Write();
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
      //
      UQuarkContent->SetFillColor(kSpring);
      UQuarkContent->SetFillStyle(3444);
      UQuarkContent->SetLineWidth(3);
      UQuarkContent->Write();
      //
      DQuarkContent->SetFillColor(kSpring);
      DQuarkContent->SetFillStyle(3444);
      DQuarkContent->SetLineWidth(3);
      DQuarkContent->Write();
      //
      SQuarkContent->SetFillColor(kSpring);
      SQuarkContent->SetFillStyle(3444);
      SQuarkContent->SetLineWidth(3);
      SQuarkContent->Write();
      //
      CQuarkContent->SetFillColor(kSpring);
      CQuarkContent->SetFillStyle(3444);
      CQuarkContent->SetLineWidth(3);
      CQuarkContent->Write();
      //
      BQuarkContent->SetFillColor(kSpring);
      BQuarkContent->SetFillStyle(3444);
      BQuarkContent->SetLineWidth(3);
      BQuarkContent->Write();
      //
      CSVLB->SetFillColor(kSpring);
      CSVLB->SetFillStyle(3444);
      CSVLB->SetLineWidth(3);
      CSVLB->Write();
      //
      CSVMB->SetFillColor(kSpring);
      CSVMB->SetFillStyle(3444);
      CSVMB->SetLineWidth(3);
      CSVMB->Write();
      //
      CSVTB->SetFillColor(kSpring);
      CSVTB->SetFillStyle(3444);
      CSVTB->SetLineWidth(3);
      CSVTB->Write();
      //
      Vtcs->SetFillColor(kSpring);
      Vtcs->SetFillStyle(3444);
      Vtcs->SetLineWidth(3);
      Vtcs->Write();
      //
      //HTCSVM->SetFillColor(kSpring);
      //HTCSVM->SetFillStyle(3444);
      //HTCSVM->SetLineWidth(3);
      HTCSVM->Write();
      //
      for (int i=0; i<6; i++) 
	{
	  HiggsMassReversedHptTpt[i]->SetFillColor(kSpring); HiggsMassReversedHptTpt[i]->SetFillStyle(3444); HiggsMassReversedHptTpt[i]->SetLineWidth(3); HiggsMassReversedHptTpt[i]->Write();
	  TprimeMassReversedHptTpt[i]->SetFillColor(kSpring); TprimeMassReversedHptTpt[i]->SetFillStyle(3444); TprimeMassReversedHptTpt[i]->SetLineWidth(3); TprimeMassReversedHptTpt[i]->Write();
	  HiggsMassCuts[i]->SetFillColor(kSpring); HiggsMassCuts[i]->SetFillStyle(3444); HiggsMassCuts[i]->SetLineWidth(3); HiggsMassCuts[i]->Write();
	  WMassFromHiggs[i]->SetFillColor(kSpring); WMassFromHiggs[i]->SetFillStyle(3444); WMassFromHiggs[i]->Write();
	  WMassFromHiggsChi2[i]->SetFillColor(kSpring); WMassFromHiggsChi2[i]->SetFillStyle(3444); WMassFromHiggsChi2[i]->Write();
	}
      //
      TprimeNotWFromHiggs->SetFillColor(kSpring);
      TprimeNotWFromHiggs->SetFillStyle(3444);
      TprimeNotWFromHiggs->SetLineWidth(3);
      TprimeNotWFromHiggs->Write();
      //
      FiveJetsMassLC->SetFillColor(kSpring);
      FiveJetsMassLC->SetFillStyle(3444);
      FiveJetsMassLC->SetLineWidth(3);
      FiveJetsMassLC->Write();
      //
      FiveJetsMassBE->SetFillColor(kSpring);
      FiveJetsMassBE->SetFillStyle(3444);
      FiveJetsMassBE->SetLineWidth(3);
      FiveJetsMassBE->Write();
      //
      TopMassFromHiggs->SetFillColor(kSpring);
      TopMassFromHiggs->SetFillStyle(3444);
      TopMassFromHiggs->SetLineWidth(3);
      TopMassFromHiggs->Write();
      //
      FiveJetsMassLCoverBE->SetFillColor(kSpring);
      FiveJetsMassLCoverBE->SetFillStyle(3444);
      FiveJetsMassLCoverBE->SetLineWidth(3);
      FiveJetsMassLCoverBE->Write();
      //
      HMBE->SetFillColor(kSpring); HMBE->SetFillStyle(3444); HMBE->SetLineWidth(3); HMBE->Write();
      HMLC->SetFillColor(kSpring); HMLC->SetFillStyle(3444); HMLC->SetLineWidth(3); HMLC->Write();
      HMLCoverBE->Write();
      FJLC_QCDBE->SetFillColor(kSpring); FJLC_QCDBE->SetFillStyle(3444); FJLC_QCDBE->SetLineWidth(3); FJLC_QCDBE->Write();
      FJLC_QCDLC->SetFillColor(kSpring); FJLC_QCDLC->SetFillStyle(3444); FJLC_QCDLC->SetLineWidth(3); FJLC_QCDLC->Write();
      FJLC_QCDLCoverBE->Write();
      //
      TM_LC_HMW->SetFillColor(kSpring);
      TM_LC_HMW->SetFillStyle(3444);
      TM_LC_HMW->SetLineWidth(3);
      TM_LC_HMW->Write();
      //
      DPWJ_HT->Write();
      //
      JET1ETA->SetFillColor(kSpring); JET1ETA->SetFillStyle(3444); JET1ETA->SetLineWidth(3); JET1ETA->Write();
      JET2ETA->SetFillColor(kSpring); JET2ETA->SetFillStyle(3444); JET2ETA->SetLineWidth(3); JET2ETA->Write();
      JET3ETA->SetFillColor(kSpring); JET3ETA->SetFillStyle(3444); JET3ETA->SetLineWidth(3); JET3ETA->Write();
      JET4ETA->SetFillColor(kSpring); JET4ETA->SetFillStyle(3444); JET4ETA->SetLineWidth(3); JET4ETA->Write();
      JET5ETA->SetFillColor(kSpring); JET5ETA->SetFillStyle(3444); JET5ETA->SetLineWidth(3); JET5ETA->Write();
      JET6ETA->SetFillColor(kSpring); JET6ETA->SetFillStyle(3444); JET6ETA->SetLineWidth(3); JET6ETA->Write();
      JET1PHI->SetFillColor(kSpring); JET1PHI->SetFillStyle(3444); JET1PHI->SetLineWidth(3); JET1PHI->Write();
      JET2PHI->SetFillColor(kSpring); JET2PHI->SetFillStyle(3444); JET2PHI->SetLineWidth(3); JET2PHI->Write();
      JET3PHI->SetFillColor(kSpring); JET3PHI->SetFillStyle(3444); JET3PHI->SetLineWidth(3); JET3PHI->Write();
      JET4PHI->SetFillColor(kSpring); JET4PHI->SetFillStyle(3444); JET4PHI->SetLineWidth(3); JET4PHI->Write();
      JET5PHI->SetFillColor(kSpring); JET5PHI->SetFillStyle(3444); JET5PHI->SetLineWidth(3); JET5PHI->Write();
      JET6PHI->SetFillColor(kSpring); JET6PHI->SetFillStyle(3444); JET6PHI->SetLineWidth(3); JET6PHI->Write();
      JETMULTI->SetFillColor(kSpring); JETMULTI->SetFillStyle(3444); JETMULTI->SetLineWidth(3); JETMULTI->Write();
    }

  double RegionA=0; double RegionB=0; double RegionC=0; double RegionD=0;
  RegionA=HptTpt[0]->Integral(HptTpt[0]->GetXaxis()->FindBin(200.0),HptTpt[0]->GetNbinsX(),HptTpt[0]->GetYaxis()->FindBin(200.0),HptTpt[0]->GetNbinsY());
  RegionB=HptTpt[0]->Integral(HptTpt[0]->GetXaxis()->FindBin(200.0),HptTpt[0]->GetNbinsX(),HptTpt[0]->GetYaxis()->FindBin(0.0),HptTpt[0]->GetYaxis()->FindBin(200.0));
  RegionC=HptTpt[0]->Integral(HptTpt[0]->GetXaxis()->FindBin(0.0),HptTpt[0]->GetXaxis()->FindBin(200.0),HptTpt[0]->GetYaxis()->FindBin(200.0),HptTpt[0]->GetNbinsY());
  RegionD=HptTpt[0]->Integral(HptTpt[0]->GetXaxis()->FindBin(0.0),HptTpt[0]->GetXaxis()->FindBin(200.0),HptTpt[0]->GetYaxis()->FindBin(0.0),HptTpt[0]->GetYaxis()->FindBin(200.0));

  cout << "ABCD Method Info for Hpt Top pt:" << endl;
  cout << "Number of events on region A (signal enriched)" <<  RegionA << endl;
  cout << "Number of events on region B " <<  RegionB << endl;
  cout << "Number of events on region C " <<  RegionC << endl;
  cout << "Number of events on region D " <<  RegionD << endl;

  exit(0);

}

