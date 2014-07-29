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
#include "PU_Reweighting.h"

using namespace std;

// Global Parameters
// Location of root files

//const int NumberOfProcesses=17;

const TString MainFolder = "file:/sps/cms/ruizalva/CMSSW_5_3_9_patch3/src/Extractors/PatExtractor/bin/WithTrigger/SUFFIXDATA/";

const int NumberOfHistos=23;
const TString Histos[NumberOfHistos] = {"Cut_0", "Cut_1", "Cut_2", "Cut_3", "Cut_chi2", "Cut_4", "Cut_5", "Cut_6", "Cut_7", "Cut_8", "Cut_9", "Cut_10", "Cut_11", "Cut_12", "Cut_13", 
"Cut_14", "Cut_15", "Cut_16", "Cut_17", "Cut_18", "Cut_19", "Cut_20", "Cut_21"};

double ABs(double X, double Y)
{
  return TMath::Abs(X-Y);
}

void Data()
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
  TH1F *UQuarkContent;
  TH1F *DQuarkContent;
  TH1F *SQuarkContent;
  TH1F *CQuarkContent;
  TH1F *BQuarkContent;
  TH1F *CSVLB;
  TH1F *CSVMB;
  TH1F *CSVTB;
  TH1F *Vtcs;
  TH1F *TopMassFromHiggs;
  TH1F *FiveJetsMassBE;
  
  TH1F *FiveJetsMassLC_HHM;
  TH1F *FiveJetsMassLC_LHM;
  TH1F *FiveJetsMassLC_LHMoverHHM;

  TH1F *TM_LC; TH1F *TM_BE;

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
  TH1F *ALLCuts[NumberOfHistos];
  CutsChain.Add(MainFolder + "Full_Data.root");
  AnalysisChain.Add(MainFolder + "Full_Data.root");
  
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
      /////////////////////
      //Saving Histograms//
      /////////////////////
      string A1 = Form("Reconstructed_Tprime.M() >> TprimeMass(60,400,1600)");
      string A2 = Form("TprimeMass");
      AnalysisChain.Draw(A1.c_str());
      FiveJetsMass = (TH1F*)gDirectory->Get(A2.c_str());
      gPad->Close();
      SurvivalMarker=true;
      //
      string LJPT1 = Form("jet1_pt >> jet1_pt(50,20,1000)");
      string LJPT2 = Form("jet1_pt");
      AnalysisChain.Draw(LJPT1.c_str());
      LeadingJetPT = (TH1F*)gDirectory->Get(LJPT2.c_str());
      gPad->Close();
      //
      string SLJPT1 = Form("jet2_pt >> jet2_pt(50,20,1000)");
      string SLJPT2 = Form("jet2_pt");
      AnalysisChain.Draw(SLJPT1.c_str());
      Leading2JetPT = (TH1F*)gDirectory->Get(SLJPT2.c_str());
      gPad->Close();
      //
      string SSLJPT1 = Form("jet3_pt >> jet3_pt(35,20,700)");
      string SSLJPT2 = Form("jet3_pt");
      AnalysisChain.Draw(SSLJPT1.c_str());
      Leading3JetPT = (TH1F*)gDirectory->Get(SSLJPT2.c_str());
      gPad->Close();
      //
      string SSSLJPT1 = Form("jet4_pt >> jet4_pt(20,20,400)");
      string SSSLJPT2 = Form("jet4_pt");
      AnalysisChain.Draw(SSSLJPT1.c_str());
      Leading4JetPT = (TH1F*)gDirectory->Get(SSSLJPT2.c_str());
      gPad->Close();
      //
      string SSSSLJPT1 = Form("jet5_pt >> jet5_pt(15,20,300)");
      string SSSSLJPT2 = Form("jet5_pt");
      AnalysisChain.Draw(SSSSLJPT1.c_str());
      Leading5JetPT = (TH1F*)gDirectory->Get(SSSSLJPT2.c_str());
      gPad->Close();
      //
      string SSSSSLJPT1 = Form("jet6_pt >> jet6_pt(10,20,200)");
      string SSSSSLJPT2 = Form("jet6_pt");
      AnalysisChain.Draw(SSSSSLJPT1.c_str());
      Leading6JetPT = (TH1F*)gDirectory->Get(SSSSSLJPT2.c_str());
      gPad->Close();
      //
      string THT1 = Form("THT >> THT(65,300,1600)");
      string THT2 = Form("THT");
      AnalysisChain.Draw(THT1.c_str());
      THT = (TH1F*)gDirectory->Get(THT2.c_str());
      gPad->Close(); 
      //
      string DRHJ1 = Form("DeltaR_of_Higgs_Jets >> DeltaR_of_Higgs_Jets(65,0.5,7)");
      string DRHJ2 = Form("DeltaR_of_Higgs_Jets");
      AnalysisChain.Draw(DRHJ1.c_str());
      DRHjets = (TH1F*)gDirectory->Get(DRHJ2.c_str());
      gPad->Close();
      //
      string DRWJ1 = Form("DeltaR_of_W_Jets >> DeltaR_of_W_Jets(65,0.5,7)");
      string DRWJ2 = Form("DeltaR_of_W_Jets");
      AnalysisChain.Draw(DRWJ1.c_str());
      DRWjets = (TH1F*)gDirectory->Get(DRWJ2.c_str());
      gPad->Close();
      //
      string HPT1 = Form("Reconstructed_Higgs.Pt() >> HPt(40,10,800)");
      string HPT2 = Form("HPt");
      AnalysisChain.Draw(HPT1.c_str());
      Hpt = (TH1F*)gDirectory->Get(HPT2.c_str());
      gPad->Close();
      //
      string TPT1 = Form("Reconstructed_Top.Pt() >> TPt(40,10,800)");
      string TPT2 = Form("TPt");
      AnalysisChain.Draw(TPT1.c_str());
      Tpt = (TH1F*)gDirectory->Get(TPT2.c_str());
      gPad->Close();
      //
      string DRWH1 = Form("DeltaR_of_W_Higgs >> DeltaR_of_W_Higgs(65,0.5,7)");
      string DRWH2 = Form("DeltaR_of_W_Higgs");
      AnalysisChain.Draw(DRWH1.c_str());
      DRWH = (TH1F*)gDirectory->Get(DRWH2.c_str());
      gPad->Close();
      //
      string DPHJ1 = Form("TMath::Abs(First_Higgs_Jet.Phi()-Second_Higgs_Jet.Phi()) >> DeltaPhi_of_Higgs_jets(30,0.0,3.0)");
      string DPHJ2 = Form("DeltaPhi_of_Higgs_jets");
      AnalysisChain.Draw(DPHJ1.c_str());
      DPHjets = (TH1F*)gDirectory->Get(DPHJ2.c_str());
      gPad->Close();
      //
      string DPWJ1 = Form("TMath::Abs(First_W_Jet.Phi()-Second_W_Jet.Phi()) >> DeltaPhi_of_W_jets(30,0.0,3.0)");
      string DPWJ2 = Form("DeltaPhi_of_W_jets");
      AnalysisChain.Draw(DPWJ1.c_str());
      DPWjets = (TH1F*)gDirectory->Get(DPWJ2.c_str());
      gPad->Close();
      //
      string DPTJ1 = Form("TMath::Abs(Top_Jet.Phi()-Reconstructed_W.Phi()) >> DeltaPhi_of_T_jet(30,0.0,3.0)");
      string DPTJ2 = Form("DeltaPhi_of_T_jet");
      AnalysisChain.Draw(DPTJ1.c_str());
      DPTjets = (TH1F*)gDirectory->Get(DPTJ2.c_str());
      gPad->Close();
      //
      string HM1 = Form("Reconstructed_Higgs.M() >> HM(36,60,180)");
      string HM2 = Form("HM");
      AnalysisChain.Draw(HM1.c_str());
      HiggsMass = (TH1F*)gDirectory->Get(HM2.c_str());
      gPad->Close();
      //
      string RHT1 = Form("Relative_THT >> RelHT(30,0,1)");
      string RHT2 = Form("RelHT");
      AnalysisChain.Draw(RHT1.c_str());
      RelHT = (TH1F*)gDirectory->Get(RHT2.c_str());
      gPad->Close();
      //
      string DRTH1 = Form("DeltaR_of_Top_Higgs >> DeltaR_of_Top_Higgs(65,0.5,7)");
      string DRTH2 = Form("DeltaR_of_Top_Higgs");
      AnalysisChain.Draw(DRTH1.c_str());
      DRTH = (TH1F*)gDirectory->Get(DRTH2.c_str());
      gPad->Close();
      //
      string PTNM1 = Form("PT_Normalized_Mass >> PT_Normalized_Mass(70,0.4,5)");
      string PTNM2 = Form("PT_Normalized_Mass");
      AnalysisChain.Draw(PTNM1.c_str());
      PtNormalizedMass = (TH1F*)gDirectory->Get(PTNM2.c_str());
      gPad->Close();
      //
      string RM1 = Form("Relative_Mass >> Relative_Mass(30,0.0,1)");
      string RM2 = Form("Relative_Mass");
      AnalysisChain.Draw(RM1.c_str());
      RelativeMass = (TH1F*)gDirectory->Get(RM2.c_str());
      gPad->Close();
      //
      string MPTNM1 = Form("Mother_PT_Normalized_Mass >> Mother_PT_Normalized_Mass(25,0.0,50)");
      string MPTNM2 = Form("Mother_PT_Normalized_Mass");
      AnalysisChain.Draw(MPTNM1.c_str());
      MotherPtNormalizedMass = (TH1F*)gDirectory->Get(MPTNM2.c_str());
      gPad->Close();
      //
      string NTops1 = Form("Number_of_Tops >> Number_of_Tops(8,0.0,8)");
      string NTops2 = Form("Number_of_Tops");
      AnalysisChain.Draw(NTops1.c_str());
      NumberOfTops = (TH1F*)gDirectory->Get(NTops2.c_str());
      gPad->Close();
      //
      string HPTTPT1 = Form("Reconstructed_Higgs.Pt():Reconstructed_Top.Pt() >> HPtTPt(40,10,800,40,10,800)");
      string HPTTPT2 = Form("HPtTPt");
      AnalysisChain.Draw(HPTTPT1.c_str());
      HptTpt = (TH2F*)gDirectory->Get(HPTTPT2.c_str());
      gPad->Close();
      ///////////NEW VARIABLES//////////////////
      string HMTM1 = Form("Reconstructed_Higgs.M()/Reconstructed_Top.M() >> HMoverTM(30,0,1.)");
      string HMTM2 = Form("HMoverTM");
      AnalysisChain.Draw(HMTM1.c_str());
      HiggsMassOverTopMass = (TH1F*)gDirectory->Get(HMTM2.c_str());
      gPad->Close();
      //
      string HTA1 = Form("((Reconstructed_Higgs.Pt()/Reconstructed_Higgs.M())-(Reconstructed_Top.Pt()/Reconstructed_Top.M()))/((Reconstructed_Higgs.Pt()/Reconstructed_Higgs.M())+(Reconstructed_Top.Pt()/Reconstructed_Top.M())) >> HTAsym(12,0,1.)");
      string HTA2 = Form("HTAsym");
      AnalysisChain.Draw(HTA1.c_str());
      HiggsTopAsymmetry = (TH1F*)gDirectory->Get(HTA2.c_str());
      gPad->Close();
      //
      string TLBT1 = Form("Number_of_Loose_and_non_med_b_tags >> TLBTag(10,0,10)");
      string TLBT2 = Form("TLBTag");
      AnalysisChain.Draw(TLBT1.c_str());
      ThirdLooseBtag = (TH1F*)gDirectory->Get(TLBT2.c_str());
      gPad->Close();
      //
      string TM1 = Form("Reconstructed_Top.M() >> TMass(60,100,700)");
      string TM2 = Form("TMass");
      AnalysisChain.Draw(TM1.c_str());
      TopMass = (TH1F*)gDirectory->Get(TM2.c_str());
      gPad->Close();
      //
      string C21 = Form("ChiSquaredSorting >> ChiSq(100,0,1000)");
      string C22 = Form("ChiSq");
      AnalysisChain.Draw(C21.c_str());
      Chi2 = (TH1F*)gDirectory->Get(C22.c_str());
      gPad->Close();
      //
      string UQ1 = Form("U Quark Content >> UQC(10,0,10)");
      string UQ2 = Form("UQC");
      AnalysisChain.Draw(UQ1.c_str());
      UQuarkContent = (TH1F*)gDirectory->Get(UQ2.c_str());
      gPad->Close();
      //
      string DQ1 = Form("D Quark Content >> DQC(10,0,10)");
      string DQ2 = Form("DQC");
      AnalysisChain.Draw(DQ1.c_str());
      DQuarkContent = (TH1F*)gDirectory->Get(DQ2.c_str());
      gPad->Close();
      //
      string SQ1 = Form("S Quark Content >> SQC(10,0,10)");
      string SQ2 = Form("SQC");
      AnalysisChain.Draw(SQ1.c_str());
      SQuarkContent = (TH1F*)gDirectory->Get(SQ2.c_str());
      gPad->Close();
      //
      string CQ1 = Form("C Quark Content >> CQC(10,0,10)");
      string CQ2 = Form("CQC");
      AnalysisChain.Draw(CQ1.c_str());
      CQuarkContent = (TH1F*)gDirectory->Get(CQ2.c_str());
      gPad->Close();
      //
      string BQ1 = Form("B Quark Content >> BQC(10,0,10)");
      string BQ2 = Form("BQC");
      AnalysisChain.Draw(BQ1.c_str());
      BQuarkContent = (TH1F*)gDirectory->Get(BQ2.c_str());
      gPad->Close();
      //B-tagging Working point	  
      string BTL1 = Form("Number_CSVLbtagged_jets >> CSVL(10,0,10)");
      string BTL2 = Form("CSVL");
      AnalysisChain.Draw(BTL1.c_str());
      CSVLB = (TH1F*)gDirectory->Get(BTL2.c_str());
      gPad->Close();
      //	  
      string BTM1 = Form("Number_CSVMbtagged_jets >> CSVM(10,0,10)");
      string BTM2 = Form("CSVM");
      AnalysisChain.Draw(BTM1.c_str());
      CSVMB = (TH1F*)gDirectory->Get(BTM2.c_str());
      gPad->Close();
      //	  
      string BTT1 = Form("Number_CSVTbtagged_jets >> CSVT(10,0,10)");
      string BTT2 = Form("CSVT");
      AnalysisChain.Draw(BTT1.c_str());
      CSVTB = (TH1F*)gDirectory->Get(BTT2.c_str());
      gPad->Close();
      //          
      string VT1 = Form("Vertices >> VTX(40,1,41)");
      string VT2 = Form("VTX");
      AnalysisChain.Draw(VT1.c_str());
      Vtcs = (TH1F*)gDirectory->Get(VT2.c_str());
      gPad->Close();
      //
      A1 = Form("Top_From_Higgs_Chi2.M() >> TopFromHiggsChi2Mass(95,50,1000)");
      A2 = Form("TopFromHiggsChi2Mass");
      AnalysisChain.Draw(A1.c_str(),"(Reconstructed_Higgs.M()>=110 && Reconstructed_Higgs.M()<=140)");
      TopMassFromHiggs = (TH1F*)gDirectory->Get(A2.c_str());
      gPad->Close();
      //
      A1 = Form("Reconstructed_Tprime.M() >> TprimeMassBkgE(60,400,1600)");
      A2 = Form("TprimeMassBkgE");
      AnalysisChain.Draw(A1.c_str(),"Top_From_Higgs_Chi2.M()>=140 && Top_From_Higgs_Chi2.M()<=230");
      FiveJetsMassBE = (TH1F*)gDirectory->Get(A2.c_str());
      gPad->Close();
      /////////////////////////////////
      /////////////////////////////////
      /////////////////////////////////      
      A1 = Form("Reconstructed_Tprime.M() >> TprimeMassLC_HHM(60,400,1600)");
      A2 = Form("TprimeMassLC_HHM");
      AnalysisChain.Draw(A1.c_str(),"(Top_From_Higgs_Chi2.M()<140 || Top_From_Higgs_Chi2.M()>230) && Reconstructed_Higgs.M()>140");
      FiveJetsMassLC_HHM = (TH1F*)gDirectory->Get(A2.c_str());
      gPad->Close();
      //
      A1 = Form("Reconstructed_Tprime.M() >> TprimeMassLC_LHM(60,400,1600)");
      A2 = Form("TprimeMassLC_LHM");
      AnalysisChain.Draw(A1.c_str(),"(Top_From_Higgs_Chi2.M()<140 || Top_From_Higgs_Chi2.M()>230) && Reconstructed_Higgs.M()<110");
      FiveJetsMassLC_LHM = (TH1F*)gDirectory->Get(A2.c_str());
      gPad->Close();
      //
      FiveJetsMassLC_LHMoverHHM=(TH1F*)FiveJetsMassLC_LHM->Clone("FiveJetsMass_LC_LHMoverHHM");
      //A1 = Form("Reconstructed_Tprime.M() >> TprimeMassLC_LHMoverHHM%i(60,400,1600)",9);
      //A2 = Form("TprimeMassLCoverLC_LHMoverHHM%i",9);
      //AnalysisChain.Draw(A1.c_str(),"(Top_From_Higgs_Chi2.M()<140 || Top_From_Higgs_Chi2.M()>230) && Reconstructed_Higgs.M()<110");
      //FiveJetsMassLC_LHMoverHHM = (TH1F*)gDirectory->Get(A2.c_str());
      //gPad->Close();
      //
      string D1 = Form("ABs(First_W_Jet.Phi(),Second_W_Jet.Phi()):THT >> DPWJ_HT(65,300,1600,60,0.5,3.5)");
      string D2 = Form("DPWJ_HT");
      AnalysisChain.Draw(D1.c_str()); //,"(Reconstructed_Higgs.Pt()<200 || Reconstructed_Top.Pt()<200) && (DeltaR_of_W_Higgs>=2.2 && DeltaR_of_W_Higgs<=3.5) && (TMath::Abs(First_Higgs_Jet.Phi()-Second_Higgs_Jet.Phi())<=1.2 && TMath::Abs(Top_Jet.Phi()-Reconstructed_W.Phi())<=1.2) && Relative_THT>=0.65 && (Relative_Mass>=0.3 && Relative_Mass<=0.5) && (Top_From_Higgs_Chi2.M()<140 || Top_From_Higgs_Chi2.M()>230) && (Reconstructed_Higgs.M()>=110 && Reconstructed_Higgs.M()<=140)");
      DPWJ_HT = (TH2F*)gDirectory->Get(D2.c_str());
      cout << "Correlation Factor DPWJ and HT: " << DPWJ_HT->GetCorrelationFactor() << " Integral->" << DPWJ_HT->Integral() << endl;
      gPad->Close();
      //
      string TM1_LC = Form("Reconstructed_Top.M() >> TMass_LC(60,100,700)");
      string TM2_LC = Form("TMass_LC");
      AnalysisChain.Draw(TM1_LC.c_str(),finalcut_data);
      TM_LC = (TH1F*)gDirectory->Get(TM2_LC.c_str());
      gPad->Close();
      //
      string TM1_BE = Form("Reconstructed_Top.M() >> TMass_BE(60,100,700)");
      string TM2_BE = Form("TMass_BE");
      AnalysisChain.Draw(TM1_BE.c_str(),"(Top_From_Higgs_Chi2.M()>=140 && Top_From_Higgs_Chi2.M()<=230) && (Reconstructed_Higgs.M()>=110 && Reconstructed_Higgs.M()<=140)");
      TM_BE = (TH1F*)gDirectory->Get(TM2_BE.c_str());
      gPad->Close();
      /////////////////////
      ///////////////////// ADDTIONAL INFO ON JETS
      /////////////////////
      A1 = Form("JET1.Eta() >> jet1_eta(100,-5,5)");
      A2 = Form("jet1_eta");
      AnalysisChain.Draw(A1.c_str(),finalcut_data);
      JET1ETA = (TH1F*)gDirectory->Get(A2.c_str());
      gPad->Close();
      //
      A1 = Form("JET2.Eta() >> jet2_eta(100,-5,5)");
      A2 = Form("jet2_eta");
      AnalysisChain.Draw(A1.c_str(),finalcut_data);
      JET2ETA = (TH1F*)gDirectory->Get(A2.c_str());
      gPad->Close();
      //
      A1 = Form("JET3.Eta() >> jet3_eta(100,-5,5)");
      A2 = Form("jet3_eta");
      AnalysisChain.Draw(A1.c_str(),finalcut_data);
      JET3ETA = (TH1F*)gDirectory->Get(A2.c_str());
      gPad->Close();
      //
      A1 = Form("JET4.Eta() >> jet4_eta(100,-5,5)");
      A2 = Form("jet4_eta");
      AnalysisChain.Draw(A1.c_str(),finalcut_data);
      JET4ETA = (TH1F*)gDirectory->Get(A2.c_str());
      gPad->Close();
      //
      A1 = Form("JET5.Eta() >> jet5_eta(100,-5,5)");
      A2 = Form("jet5_eta");
      AnalysisChain.Draw(A1.c_str(),finalcut_data);
      JET5ETA = (TH1F*)gDirectory->Get(A2.c_str());
      gPad->Close();
      //
      A1 = Form("JET6.Eta() >> jet6_eta(100,-5,5)");
      A2 = Form("jet6_eta");
      AnalysisChain.Draw(A1.c_str(),finalcut_data);
      JET6ETA = (TH1F*)gDirectory->Get(A2.c_str());
      gPad->Close();
      //
      A1 = Form("JET1.Phi() >> jet1_phi(60,-3,3)");
      A2 = Form("jet1_phi");
      AnalysisChain.Draw(A1.c_str(),finalcut_data);
      JET1PHI = (TH1F*)gDirectory->Get(A2.c_str());
      gPad->Close();
      //
      A1 = Form("JET2.Phi() >> jet2_phi(60,-3,3)");
      A2 = Form("jet2_phi");
      AnalysisChain.Draw(A1.c_str(),finalcut_data);
      JET2PHI = (TH1F*)gDirectory->Get(A2.c_str());
      gPad->Close();
      //
      A1 = Form("JET3.Phi() >> jet3_phi(60,-3,3)");
      A2 = Form("jet3_phi");
      AnalysisChain.Draw(A1.c_str(),finalcut_data);
      JET3PHI = (TH1F*)gDirectory->Get(A2.c_str());
      gPad->Close();
      //
      A1 = Form("JET4.Phi() >> jet4_phi(60,-3,3)");
      A2 = Form("jet4_phi");
      AnalysisChain.Draw(A1.c_str(),finalcut_data);
      JET4PHI = (TH1F*)gDirectory->Get(A2.c_str());
      gPad->Close();
      //
      A1 = Form("JET5.Phi() >> jet5_phi(60,-3,3)");
      A2 = Form("jet5_phi");
      AnalysisChain.Draw(A1.c_str(),finalcut_data);
      JET5PHI = (TH1F*)gDirectory->Get(A2.c_str());
      gPad->Close();
      //
      A1 = Form("JET6.Phi() >> jet6_phi(60,-3,3)");
      A2 = Form("jet6_phi");
      AnalysisChain.Draw(A1.c_str(),finalcut_data);
      JET6PHI = (TH1F*)gDirectory->Get(A2.c_str());
      gPad->Close();
      //
      A1 = Form("Number_Jets >> Num_jets(20,0,20)");
      A2 = Form("Num_jets");
      AnalysisChain.Draw(A1.c_str(),finalcut_data);
      JETMULTI = (TH1F*)gDirectory->Get(A2.c_str());
      gPad->Close();
    }
  
  //Settings for signal
  if (SurvivalMarker)
    {
      //TprimeMass
      FiveJetsMass->SetFillColor(kWhite);
      FiveJetsMass->SetFillStyle(0);
      FiveJetsMass->SetLineWidth(3);
      TFile f("Data.root", "RECREATE");
      FiveJetsMass->Write();
      //LJPT
      LeadingJetPT->SetFillColor(kWhite);
      LeadingJetPT->SetFillStyle(0);
      LeadingJetPT->SetLineWidth(3);
      LeadingJetPT->Write();
      Leading2JetPT->SetFillColor(kWhite);
      Leading2JetPT->SetFillStyle(0);
      Leading2JetPT->SetLineWidth(3);
      Leading2JetPT->Write();
      Leading3JetPT->SetFillColor(kWhite);
      Leading3JetPT->SetFillStyle(0);
      Leading3JetPT->SetLineWidth(3);
      Leading3JetPT->Write();
      Leading4JetPT->SetFillColor(kWhite);
      Leading4JetPT->SetFillStyle(0);
      Leading4JetPT->SetLineWidth(3);
      Leading4JetPT->Write();
      Leading5JetPT->SetFillColor(kWhite);
      Leading5JetPT->SetFillStyle(0);
      Leading5JetPT->SetLineWidth(3);
      Leading5JetPT->Write();
      Leading6JetPT->SetFillColor(kWhite);
      Leading6JetPT->SetFillStyle(0);
      Leading6JetPT->SetLineWidth(3);
      Leading6JetPT->Write();
      //THT
      THT->SetFillColor(kWhite);
      THT->SetFillStyle(0);
      THT->SetLineWidth(3);
      THT->Write();
      //
      DRHjets->SetFillColor(kWhite);
      DRHjets->SetFillStyle(0);
      DRHjets->SetLineWidth(3);
      DRHjets->Write();
      //
      DRWjets->SetFillColor(kWhite);
      DRWjets->SetFillStyle(0);
      DRWjets->SetLineWidth(3);
      DRWjets->Write();
      //
      Hpt->SetFillColor(kWhite);
      Hpt->SetFillStyle(0);
      Hpt->SetLineWidth(3);
      Hpt->Write();
      //
      Tpt->SetFillColor(kWhite);
      Tpt->SetFillStyle(0);
      Tpt->SetLineWidth(3);
      Tpt->Write();
      //
      DRWH->SetFillColor(kWhite);
      DRWH->SetFillStyle(0);
      DRWH->SetLineWidth(3);
      DRWH->Write();
      //
      DPHjets->SetFillColor(kWhite);
      DPHjets->SetFillStyle(0);
      DPHjets->SetLineWidth(3);
      DPHjets->Write();
      //
      DPWjets->SetFillColor(kWhite);
      DPWjets->SetFillStyle(0);
      DPWjets->SetLineWidth(3);
      DPWjets->Write();
      //
      DPTjets->SetFillColor(kWhite);
      DPTjets->SetFillStyle(0);
      DPTjets->SetLineWidth(3);
      DPTjets->Write();
      //
      HiggsMass->SetFillColor(kWhite);
      HiggsMass->SetFillStyle(0);
      HiggsMass->SetLineWidth(3);
      HiggsMass->Write();
      //
      RelHT->SetFillColor(kWhite);
      RelHT->SetFillStyle(0);
      RelHT->SetLineWidth(3);
      RelHT->Write();
      //
      DRTH->SetFillColor(kWhite);
      DRTH->SetFillStyle(0);
      DRTH->SetLineWidth(3);
      DRTH->Write();
      //
      PtNormalizedMass->SetFillColor(kWhite);
      PtNormalizedMass->SetFillStyle(0);
      PtNormalizedMass->SetLineWidth(3);
      PtNormalizedMass->Write();
      //
      RelativeMass->SetFillColor(kWhite);
      RelativeMass->SetFillStyle(0);
      RelativeMass->SetLineWidth(3);
      RelativeMass->Write();
      //
      MotherPtNormalizedMass->SetFillColor(kWhite);
      MotherPtNormalizedMass->SetFillStyle(0);
      MotherPtNormalizedMass->SetLineWidth(3);
      MotherPtNormalizedMass->Write();
      //
      NumberOfTops->SetFillColor(kWhite);
      NumberOfTops->SetFillStyle(0);
      NumberOfTops->SetLineWidth(3);
      NumberOfTops->Write();
      //
      HptTpt->SetFillColor(kWhite);
      HptTpt->SetFillStyle(0);
      HptTpt->SetLineWidth(3);
      HptTpt->Write();
      //
      HiggsMassOverTopMass->SetFillColor(kWhite);
      HiggsMassOverTopMass->SetFillStyle(0);
      HiggsMassOverTopMass->SetLineWidth(3);
      HiggsMassOverTopMass->Write();
      //
      HiggsTopAsymmetry->SetFillColor(kWhite);
      HiggsTopAsymmetry->SetFillStyle(0);
      HiggsTopAsymmetry->SetLineWidth(3);
      HiggsTopAsymmetry->Write();
      //
      ThirdLooseBtag->SetFillColor(kWhite);
      ThirdLooseBtag->SetFillStyle(0);
      ThirdLooseBtag->SetLineWidth(3);
      ThirdLooseBtag->Write();
      //
      TopMass->SetFillColor(kWhite);
      TopMass->SetFillStyle(0);
      TopMass->SetLineWidth(3);
      TopMass->Write();
      //
      Chi2->SetFillColor(kWhite);
      Chi2->SetFillStyle(0);
      Chi2->SetLineWidth(3);
      Chi2->Write();
      //
      UQuarkContent->SetFillColor(kWhite);
      UQuarkContent->SetFillStyle(0);
      UQuarkContent->SetLineWidth(3);
      UQuarkContent->Write();
      //
      DQuarkContent->SetFillColor(kWhite);
      DQuarkContent->SetFillStyle(0);
      DQuarkContent->SetLineWidth(3);
      DQuarkContent->Write();
      //
      SQuarkContent->SetFillColor(kWhite);
      SQuarkContent->SetFillStyle(0);
      SQuarkContent->SetLineWidth(3);
      SQuarkContent->Write();
      //
      CQuarkContent->SetFillColor(kWhite);
      CQuarkContent->SetFillStyle(0);
      CQuarkContent->SetLineWidth(3);
      CQuarkContent->Write();
      //
      BQuarkContent->SetFillColor(kWhite);
      BQuarkContent->SetFillStyle(0);
      BQuarkContent->SetLineWidth(3);
      BQuarkContent->Write();
      //
      CSVLB->SetFillColor(kWhite);
      CSVLB->SetFillStyle(0);
      CSVLB->SetLineWidth(3);
      CSVLB->Write();
      //
      CSVMB->SetFillColor(kWhite);
      CSVMB->SetFillStyle(0);
      CSVMB->SetLineWidth(3);
      CSVMB->Write();
      //
      CSVTB->SetFillColor(kWhite);
      CSVTB->SetFillStyle(0);
      CSVTB->SetLineWidth(3);
      CSVTB->Write();
      //
      Vtcs->SetFillColor(kWhite);
      Vtcs->SetFillStyle(0);
      Vtcs->SetLineWidth(3);
      Vtcs->Write();
      //
      TopMassFromHiggs->SetFillColor(kWhite);
      TopMassFromHiggs->SetFillStyle(0);
      TopMassFromHiggs->SetLineWidth(3);
      TopMassFromHiggs->Write();
      //
      FiveJetsMassBE->SetFillColor(kWhite);
      FiveJetsMassBE->SetFillStyle(0);
      FiveJetsMassBE->SetLineWidth(3);
      FiveJetsMassBE->Write();
      ////////////////////////////////////
      ////////////////////////////////////
      ////////////////////////////////////
      FiveJetsMassLC_HHM->SetFillColor(kWhite);
      FiveJetsMassLC_HHM->SetFillStyle(0);
      FiveJetsMassLC_HHM->SetLineWidth(3);
      FiveJetsMassLC_HHM->Write();
      //
      FiveJetsMassLC_LHM->SetFillColor(kWhite);
      FiveJetsMassLC_LHM->SetFillStyle(0);
      FiveJetsMassLC_LHM->SetLineWidth(3);
      FiveJetsMassLC_LHM->Write();
      //
      FiveJetsMassLC_LHMoverHHM->Scale(FiveJetsMassLC_HHM->Integral()/FiveJetsMassLC_LHM->Integral());
      FiveJetsMassLC_LHMoverHHM->Divide(FiveJetsMassLC_HHM);
      FiveJetsMassLC_LHMoverHHM->Write();
      //
      DPWJ_HT->Write();
      //
      TM_LC->SetFillColor(kWhite);
      TM_LC->SetFillStyle(0);
      TM_LC->SetLineWidth(3);
      TM_LC->Write();
      //
      TM_BE->SetFillColor(kWhite);
      TM_BE->SetFillStyle(0);
      TM_BE->SetLineWidth(3);
      TM_BE->Write();  
      //
      TH1F *TM_LC_ttbar_subs=(TH1F*)TM_LC->Clone("TM_LC_ttbar_subs");
      TM_LC_ttbar_subs->Sumw2(); 
      TM_LC_ttbar_subs->Add(TM_BE,-0.63);
      TM_LC_ttbar_subs->SetFillColor(kWhite);
      TM_LC_ttbar_subs->SetFillStyle(0);
      TM_LC_ttbar_subs->SetLineWidth(3);
      TM_LC_ttbar_subs->Write();
      //
      JET1ETA->SetFillColor(kWhite); JET1ETA->SetFillStyle(0); JET1ETA->SetLineWidth(3); JET1ETA->Write();
      JET2ETA->SetFillColor(kWhite); JET2ETA->SetFillStyle(0); JET2ETA->SetLineWidth(3); JET2ETA->Write();
      JET3ETA->SetFillColor(kWhite); JET3ETA->SetFillStyle(0); JET3ETA->SetLineWidth(3); JET3ETA->Write();
      JET4ETA->SetFillColor(kWhite); JET4ETA->SetFillStyle(0); JET4ETA->SetLineWidth(3); JET4ETA->Write();
      JET5ETA->SetFillColor(kWhite); JET5ETA->SetFillStyle(0); JET5ETA->SetLineWidth(3); JET5ETA->Write();
      JET6ETA->SetFillColor(kWhite); JET6ETA->SetFillStyle(0); JET6ETA->SetLineWidth(3); JET6ETA->Write();
      JET1PHI->SetFillColor(kWhite); JET1PHI->SetFillStyle(0); JET1PHI->SetLineWidth(3); JET1PHI->Write();
      JET2PHI->SetFillColor(kWhite); JET2PHI->SetFillStyle(0); JET2PHI->SetLineWidth(3); JET2PHI->Write();
      JET3PHI->SetFillColor(kWhite); JET3PHI->SetFillStyle(0); JET3PHI->SetLineWidth(3); JET3PHI->Write();
      JET4PHI->SetFillColor(kWhite); JET4PHI->SetFillStyle(0); JET4PHI->SetLineWidth(3); JET4PHI->Write();
      JET5PHI->SetFillColor(kWhite); JET5PHI->SetFillStyle(0); JET5PHI->SetLineWidth(3); JET5PHI->Write();
      JET6PHI->SetFillColor(kWhite); JET6PHI->SetFillStyle(0); JET6PHI->SetLineWidth(3); JET6PHI->Write();
      JETMULTI->SetFillColor(kWhite); JETMULTI->SetFillStyle(0); JETMULTI->SetLineWidth(3); JETMULTI->Write();

    }

  exit(0);

}

