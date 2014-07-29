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

//const int NumberOfProcesses=16;

const TString MainFolder = "file:/sps/cms/ruizalva/CMSSW_5_3_9_patch3/src/Extractors/PatExtractor/bin/WithTrigger/SUFFIX/";

const TString StorageDirPrefix = MainFolder + "TTJets/";

const int NumberOfHistos=23;
const TString Histos[NumberOfHistos] = {"Cut_0", "Cut_1", "Cut_2", "Cut_3", "Cut_chi2", "Cut_4", "Cut_5", "Cut_6", "Cut_7", "Cut_8", "Cut_9", "Cut_10", "Cut_11", "Cut_12", "Cut_13", 
"Cut_14", "Cut_15", "Cut_16", "Cut_17", "Cut_18", "Cut_19", "Cut_20", "Cut_21"};

//Computing weights (Everything in pb)
const float Xs= XS[9];
 
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

// Breit-Wigner Gaussian fitting function definition
double GaussBreitW(Double_t *xx, Double_t *par)
  //----------------------------------------------------------------------
{
  double M = xx[0];
  double M0 = par[1], M02=M0*M0;
  double G0 = 2.49*M0/80.4, G02=G0*G0;
  //double G0 = 2.49*M0/par[3], G02=G0*G0;
  double sigma = par[2], sigma2=sigma*sigma;
  
  double sum = 0;
  for(int i=0 ; i<=15 ; i++)
    {
      double m = i*sigma/5;
      double u = (M+m)*(M+m)-M02;
      double v = (M-m)*(M-m)-M02;
      double BW = 1/(u*u+G02*M02) + 1/(v*v+G02*M02);
      if( i==0 ) BW /= 2;
      sum += BW*exp(-m*m/sigma2/2);
    }
  sum *= G02*M02/31/sqrt(2*acos(-1.))/sigma;
  return par[0]*sum;
}

double bkg(Double_t *x, Double_t *par)
{
  return par[0]+par[1]*x[0];
}

double fit_function(Double_t *x, Double_t *par)
{
  return GaussBreitW(x,&par[2]);
  //return bkg(x,par)+GaussBreitW(x,&par[2]);
}

float PUR_function(int TI) //function with input the Number of True Interactions
{

  if (TI>=PUBins) return 1;
  else return PU_weight[TI];

}

double ABs(double X, double Y)
{
  return TMath::Abs(X-Y);
}

void TTJets()
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
  TH1F *Vtcs;  
  TH1F *CLBT;   
  TH1F *CMBT; 
  TH1F *CTBT;  
  TH1F *FLLBT; 
  TH1F *FLMBT; 
  TH1F *FLTBT; 
  TH1F *FCLBT; 
  TH1F *FCMBT; 
  TH1F *FCTBT; 
  TH2F *HptTpt[6];
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
  TH1F *TopMassFromHiggsLC;
  TH1F *FiveJetsMassLCoverBE; //5 jets mass ratio plot in vs out the cut
  TH1F *HMBE;
  TH1F *HMLC;
  TH1F *HMLCoverBE;
  TH1F *FJLC_QCDBE;
  TH1F *FJLC_QCDLC;
  TH1F *FJLC_QCDLCoverBE;

  //Checking uncorrelated sector of higgs mass and top from higgs mass
  TH1F *FiveJetsMassBE_HHM; //5 jets mass inside the last cut with Higgs mass >140 GeV
  TH1F *FiveJetsMassLC_HHM; //5 jets mass outside last cut with Higgs mass >140 GeV
  TH1F *TopMassFromHiggs_HHM; //Top mass from Higgs with Higgs mass >140
  TH1F *TopMassFromHiggs_HHMoverFHM; //Ratio of the Top mass from Higgs with Higgs mass >140 over Top mass from Higgs with Higgs mass in full range
  TH1F *FiveJetsMassLCoverBE_HHM; //5 jets mass ratio plot in vs out the cut with Higgs mass >140 GeV
  //Checking correlation sector of low higgs mass
  TH1F *FiveJetsMassBE_LHM; 
  TH1F *FiveJetsMassLC_LHM;
  TH1F *FiveJetsMassLCoverBE_LHM;
  TH1F *FiveJetsMassLC_LHMoverHHM;

  //First top study for estimation of overall QCD
  TH1F *TM_LC; TH1F *TM_BE; TH1F *TM_LCoverBE;
  TH1F *TM_LC_HMW; //Top mass after ttbar estimation cut inside the Higgs mass window
  TH1F *TM_BE_HMW;
  TH1F *TM_LCoverBE_HMW;

  //Study on Correlation of Higgs mass and second W mass
  TH1F *HM_W_BE; TH1F *HM_W_LC; TH1F *HM_W_LCoverBE;
  TH1F *FJM_W_BE; TH1F *FJM_W_LC; TH1F *FJM_W_LCoverBE;

  //RM or RHT sideband estimations
  TH1F *FiveJetsMassLC_LRHT; TH1F *FiveJetsMassLC_HRHT;
  TH1F *FiveJetsMassLC_LRM; TH1F *FiveJetsMassLC_HRM; TH1F *FiveJetsMassLC_WRM;
  
  //HT estimation
  TH2F *DPWJ_HT;
  TH2F *DPWJ_HT_BE;

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
  CutsChain.Add(MainFolder + "TTJets_Full_analyzed.root");
  AnalysisChain.Add(MainFolder + "TTJets_Full_analyzed.root");
  AnalysisChain3L.Add(MainFolder + "TTJets_Full_analyzed.root");
        
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
      string A1 = Form("Reconstructed_Tprime.M() >> TprimeMass%i(60,400,1600)",9);
      string A2 = Form("TprimeMass%i",9);
      AnalysisChain.Draw(A1.c_str(),finalcut);
      FiveJetsMass = (TH1F*)gDirectory->Get(A2.c_str());
      gPad->Close();
      SurvivalMarker=true;
      //
      A1 = Form("Reconstructed_Tprime.M() >> TprimeMassBkgE%i(60,400,1600)",9);
      A2 = Form("TprimeMassBkgE%i",9);
      AnalysisChain.Draw(A1.c_str(),"(PUR_function(Number_True_Interactions))*weight*((Top_From_Higgs_Chi2.M()>=140 && Top_From_Higgs_Chi2.M()<=230) && (Reconstructed_Higgs.M()>=110 && Reconstructed_Higgs.M()<=140))");
      FiveJetsMassBE = (TH1F*)gDirectory->Get(A2.c_str());
      gPad->Close();
      //
      A1 = Form("Reconstructed_Tprime.M() >> TprimeMassLC%i(60,400,1600)",9);
      A2 = Form("TprimeMassLC%i",9);
      AnalysisChain.Draw(A1.c_str(),"(PUR_function(Number_True_Interactions))*weight*((Top_From_Higgs_Chi2.M()<140 || Top_From_Higgs_Chi2.M()>230) && (Reconstructed_Higgs.M()>=110 && Reconstructed_Higgs.M()<=140))");
      FiveJetsMassLC = (TH1F*)gDirectory->Get(A2.c_str());
      gPad->Close();
      //
      FiveJetsMassLCoverBE=(TH1F*)FiveJetsMassLC->Clone("FiveJetsMass_LCoverBE");
      //A1 = Form("Reconstructed_Tprime.M() >> TprimeMassLCoverBE%i(60,400,1600)",9);
      //A2 = Form("TprimeMassLCoverBE%i",9);
      //AnalysisChain.Draw(A1.c_str(),"Top_From_Higgs_Chi2.M()<140 || Top_From_Higgs_Chi2.M()>230");
      //FiveJetsMassLCoverBE = (TH1F*)gDirectory->Get(A2.c_str());
      //gPad->Close();
      //
      string HMBEstim1 = Form("Reconstructed_Higgs.M() >> HiggsMassBE%i(36,60,180)",9);
      string HMBEstim2 = Form("HiggsMassBE%i",9);
      AnalysisChain.Draw(HMBEstim1.c_str(),"Top_From_Higgs_Chi2.M()>=140 && Top_From_Higgs_Chi2.M()<=230");
      HMBE = (TH1F*)gDirectory->Get(HMBEstim2.c_str());
      gPad->Close();
      //
      HMBEstim1 = Form("Reconstructed_Higgs.M() >> HiggsMassLC%i(36,60,180)",9);
      HMBEstim2 = Form("HiggsMassLC%i",9);
      AnalysisChain.Draw(HMBEstim1.c_str(),"Top_From_Higgs_Chi2.M()<140 || Top_From_Higgs_Chi2.M()>230");
      HMLC = (TH1F*)gDirectory->Get(HMBEstim2.c_str());
      gPad->Close();
      //
      HMLCoverBE=(TH1F*)HMLC->Clone("HM_LCoverBE");
      //HMBEstim1 = Form("Reconstructed_Higgs.M() >> HiggsMassLCoverBE%i(36,60,180)",9);
      //HMBEstim2 = Form("HiggsMassLCoverBE%i",9);
      //AnalysisChain.Draw(HMBEstim1.c_str(),"Top_From_Higgs_Chi2.M()<140 || Top_From_Higgs_Chi2.M()>230");
      //HMLCoverBE = (TH1F*)gDirectory->Get(HMBEstim2.c_str());
      //gPad->Close();
      //
      HMBEstim1 = Form("Reconstructed_Tprime.M() >> FJLCQCDBE%i(60,400,1600)",9);
      HMBEstim2 = Form("FJLCQCDBE%i",9);
      AnalysisChain.Draw(HMBEstim1.c_str(),"(Top_From_Higgs_Chi2.M()<140 || Top_From_Higgs_Chi2.M()>230) && (Reconstructed_Higgs.M()<110 || Reconstructed_Higgs.M()>140)");
      FJLC_QCDBE = (TH1F*)gDirectory->Get(HMBEstim2.c_str());
      gPad->Close();
      //
      HMBEstim1 = Form("Reconstructed_Tprime.M() >> FJLCQCDLC%i(60,400,1600)",9);
      HMBEstim2 = Form("FJLCQCDLC%i",9);
      AnalysisChain.Draw(HMBEstim1.c_str(),"(Top_From_Higgs_Chi2.M()<140 || Top_From_Higgs_Chi2.M()>230) && (Reconstructed_Higgs.M()>=110 && Reconstructed_Higgs.M()<=140)");
      FJLC_QCDLC = (TH1F*)gDirectory->Get(HMBEstim2.c_str());
      gPad->Close();
      //
      FJLC_QCDLCoverBE=(TH1F*)FJLC_QCDLC->Clone("FJLC_QCDLCoverBE");
      //HMBEstim1 = Form("Reconstructed_Tprime.M() >> FJLCQCDLCoverBE%i(60,400,1600)",9);
      //HMBEstim2 = Form("FJLCQCDLCoverBE%i",9);
      //AnalysisChain.Draw(HMBEstim1.c_str(),"(Top_From_Higgs_Chi2.M()<140 || Top_From_Higgs_Chi2.M()>230) && (Reconstructed_Higgs.M()>=110 && Reconstructed_Higgs.M()<=140)");
      //FJLC_QCDLCoverBE = (TH1F*)gDirectory->Get(HMBEstim2.c_str());
      //gPad->Close();
      //////////
      //////////
      A1 = Form("Top_From_Higgs_Chi2.M() >> TopFromHiggsChi2Mass%i(95,50,1000)",9);
      A2 = Form("TopFromHiggsChi2Mass%i",9);
      AnalysisChain.Draw(A1.c_str(),"(PUR_function(Number_True_Interactions))*weight*(Reconstructed_Higgs.M()>=110 && Reconstructed_Higgs.M()<=140)");
      TopMassFromHiggs = (TH1F*)gDirectory->Get(A2.c_str());
      gPad->Close();
      //
      A1 = Form("Top_From_Higgs_Chi2.M() >> TopFromHiggsChi2MassLC%i(75,50,800)",9);
      A2 = Form("TopFromHiggsChi2MassLC%i",9);
      AnalysisChain.Draw(A1.c_str(),"Number_of_Tops<=2");
      TopMassFromHiggsLC = (TH1F*)gDirectory->Get(A2.c_str());
      gPad->Close();      
      //////////
      //////////
      A1 = Form("Reconstructed_Tprime.M() >> TprimeMassBkgE_HHM%i(60,400,1600)",9);
      A2 = Form("TprimeMassBkgE_HHM%i",9);
      AnalysisChain.Draw(A1.c_str(),"(Top_From_Higgs_Chi2.M()>=140 && Top_From_Higgs_Chi2.M()<=230) && Reconstructed_Higgs.M()>140");
      FiveJetsMassBE_HHM = (TH1F*)gDirectory->Get(A2.c_str());
      gPad->Close();
      //
      A1 = Form("Reconstructed_Tprime.M() >> TprimeMassLC_HHM%i(60,400,1600)",9);
      A2 = Form("TprimeMassLC_HHM%i",9);
      AnalysisChain.Draw(A1.c_str(),"(Top_From_Higgs_Chi2.M()<140 || Top_From_Higgs_Chi2.M()>230) && Reconstructed_Higgs.M()>140");
      FiveJetsMassLC_HHM = (TH1F*)gDirectory->Get(A2.c_str());
      gPad->Close();
      //
      FiveJetsMassLCoverBE_HHM=(TH1F*)FiveJetsMassLC_HHM->Clone("FiveJetsMass_LCoverBE_HHM");
      //A1 = Form("Reconstructed_Tprime.M() >> TprimeMassLCoverBE_HHM%i(60,400,1600)",9);
      //A2 = Form("TprimeMassLCoverBE_HHM%i",9);
      //AnalysisChain.Draw(A1.c_str(),"(Top_From_Higgs_Chi2.M()<140 || Top_From_Higgs_Chi2.M()>230) && Reconstructed_Higgs.M()>140");
      //FiveJetsMassLCoverBE_HHM = (TH1F*)gDirectory->Get(A2.c_str());
      //gPad->Close();
      //
      A1 = Form("Top_From_Higgs_Chi2.M() >> TopFromHiggsChi2Mass_HHM%i(75,50,800)",9);
      A2 = Form("TopFromHiggsChi2Mass_HHM%i",9);
      AnalysisChain.Draw(A1.c_str(),"Reconstructed_Higgs.M()>140");
      TopMassFromHiggs_HHM = (TH1F*)gDirectory->Get(A2.c_str());
      gPad->Close();
      //17
      A1 = Form("Top_From_Higgs_Chi2.M() >> TopFromHiggsChi2Mass_HHMoverFHM%i(75,50,800)",9);
      A2 = Form("TopFromHiggsChi2Mass_HHMoverFHM%i",9);
      AnalysisChain.Draw(A1.c_str(),"Reconstructed_Higgs.M()>140");
      TopMassFromHiggs_HHMoverFHM = (TH1F*)gDirectory->Get(A2.c_str());
      gPad->Close();
      //
      A1 = Form("Reconstructed_Tprime.M() >> TprimeMassBE_LHM%i(60,400,1600)",9);
      A2 = Form("TprimeMassBE_LHM%i",9);
      AnalysisChain.Draw(A1.c_str(),"(Top_From_Higgs_Chi2.M()>=140 && Top_From_Higgs_Chi2.M()<=230) && Reconstructed_Higgs.M()<110");
      FiveJetsMassBE_LHM = (TH1F*)gDirectory->Get(A2.c_str());
      gPad->Close();
      //
      A1 = Form("Reconstructed_Tprime.M() >> TprimeMassLC_LHM%i(60,400,1600)",9);
      A2 = Form("TprimeMassLC_LHM%i",9);
      AnalysisChain.Draw(A1.c_str(),"(Top_From_Higgs_Chi2.M()<140 || Top_From_Higgs_Chi2.M()>230) && Reconstructed_Higgs.M()<110");
      FiveJetsMassLC_LHM = (TH1F*)gDirectory->Get(A2.c_str());
      gPad->Close();
      //
      FiveJetsMassLCoverBE_LHM=(TH1F*)FiveJetsMassLC_LHM->Clone("FiveJetsMass_LCoverBE_LHM");
      //A1 = Form("Reconstructed_Tprime.M() >> TprimeMassLCoverBE_LHM%i(60,400,1600)",9);
      //A2 = Form("TprimeMassLCoverBE_LHM%i",9);
      //AnalysisChain.Draw(A1.c_str(),"(Top_From_Higgs_Chi2.M()<140 || Top_From_Higgs_Chi2.M()>230) && Reconstructed_Higgs.M()<110");
      //FiveJetsMassLCoverBE_LHM = (TH1F*)gDirectory->Get(A2.c_str());
      //gPad->Close();
      //
      FiveJetsMassLC_LHMoverHHM=(TH1F*)FiveJetsMassLC_LHM->Clone("FiveJetsMass_LC_LHMoverHHM");
      //A1 = Form("Reconstructed_Tprime.M() >> TprimeMassLC_LHMoverHHM%i(60,400,1600)",9);
      //A2 = Form("TprimeMassLCoverLC_LHMoverHHM%i",9);
      //AnalysisChain.Draw(A1.c_str(),"(Top_From_Higgs_Chi2.M()<140 || Top_From_Higgs_Chi2.M()>230) && Reconstructed_Higgs.M()<110");
      //FiveJetsMassLC_LHMoverHHM = (TH1F*)gDirectory->Get(A2.c_str());
      //gPad->Close();      
      ///////////////////////////////////////////////
      ///////////////////////////////////////////////
      ///////////////////////////////////////////////
      string TM1_LC = Form("Reconstructed_Top.M() >> TMass_LC%i(60,100,700)",9);
      string TM2_LC = Form("TMass_LC%i",9);
      AnalysisChain.Draw(TM1_LC.c_str(),"(Top_From_Higgs_Chi2.M()<140 || Top_From_Higgs_Chi2.M()>230) && (Reconstructed_Higgs.M()>=110 && Reconstructed_Higgs.M()<=140)");
      TM_LC = (TH1F*)gDirectory->Get(TM2_LC.c_str());
      gPad->Close();
      //
      string TM1_BE = Form("Reconstructed_Top.M() >> TMass_BE%i(60,100,700)",9);
      string TM2_BE = Form("TMass_BE%i",9);
      AnalysisChain.Draw(TM1_BE.c_str(),"(Top_From_Higgs_Chi2.M()>=140 && Top_From_Higgs_Chi2.M()<=230) && (Reconstructed_Higgs.M()>=110 && Reconstructed_Higgs.M()<=140)");
      TM_BE = (TH1F*)gDirectory->Get(TM2_BE.c_str());
      gPad->Close();
      //
      TM_LCoverBE=(TH1F*)TM_LC->Clone("TM__LCoverBE");
      //TM1_LC = Form("Reconstructed_Top.M() >> TMass_LCoverBE%i(60,100,700)",9);
      //TM2_LC = Form("TMass_LCoverBE%i",9);
      //AnalysisChain.Draw(TM1_LC.c_str(),"(Top_From_Higgs_Chi2.M()<140 || Top_From_Higgs_Chi2.M()>230)");
      //TM_LCoverBE = (TH1F*)gDirectory->Get(TM2_LC.c_str());
      //gPad->Close();
      //
      TM1_LC = Form("Reconstructed_Top.M() >> TMass_LC_HMW%i(60,100,700)",9);
      TM2_LC = Form("TMass_LC_HMW%i",9);
      AnalysisChain.Draw(TM1_LC.c_str(),"(Top_From_Higgs_Chi2.M()<140 || Top_From_Higgs_Chi2.M()>230) && (Reconstructed_Higgs.M()>=110 && Reconstructed_Higgs.M()<=140)");
      TM_LC_HMW = (TH1F*)gDirectory->Get(TM2_LC.c_str());
      gPad->Close();
      ///////////////////////////////////////////////
      ////////////////2nd W study////////////////////
      ///////////////////////////////////////////////      
      string WBE1 = Form("Reconstructed_Higgs.M() >> HM_W_LC%i(36,60,180)",9);
      string WBE2 = Form("HM_W_LC%i",9);
      AnalysisChain.Draw(WBE1.c_str(),"(W_From_Higgs_Chi2.M()<60 || W_From_Higgs_Chi2.M()>120)");
      HM_W_LC = (TH1F*)gDirectory->Get(WBE2.c_str());
      gPad->Close();
      //
      HM_W_LCoverBE =(TH1F*)HM_W_LC->Clone("HM_WE_LCoverBE");
      //WBE1 = Form("Reconstructed_Higgs.M() >> HM_W_LCoverBE%i(36,60,180)",9);
      //WBE2 = Form("HM_W_BE%i",9);
      //AnalysisChain.Draw(WBE1.c_str(),"(W_From_Higgs_Chi2.M()<60 || W_From_Higgs_Chi2.M()>120)");
      //HM_W_LCoverBE = (TH1F*)gDirectory->Get(WBE2.c_str());
      //gPad->Close();
      //
      WBE1 = Form("Reconstructed_Higgs.M() >> HM_W_BE%i(36,60,180)",9);
      WBE2 = Form("HM_W_BE%i",9);
      AnalysisChain.Draw(WBE1.c_str(),"(W_From_Higgs_Chi2.M()>60 && W_From_Higgs_Chi2.M()<120)");
      HM_W_BE = (TH1F*)gDirectory->Get(WBE2.c_str());
      gPad->Close();
      //
      WBE1 = Form("Reconstructed_Tprime.M() >> FJM_W_BE%i(60,400,1600)",9);
      WBE2 = Form("FJM_W_BE%i",9);
      AnalysisChain.Draw(WBE1.c_str(),"(W_From_Higgs_Chi2.M()>60 && W_From_Higgs_Chi2.M()<120)");
      FJM_W_BE = (TH1F*)gDirectory->Get(WBE2.c_str());
      gPad->Close();
      //
      WBE1 = Form("Reconstructed_Tprime.M() >> FJM_W_LC%i(60,400,1600)",9);
      WBE2 = Form("FJM_W_LC%i",9);
      AnalysisChain.Draw(WBE1.c_str(),"(W_From_Higgs_Chi2.M()<60 || W_From_Higgs_Chi2.M()>120)");
      FJM_W_LC = (TH1F*)gDirectory->Get(WBE2.c_str());
      gPad->Close();
      //
      FJM_W_LCoverBE =(TH1F*)FJM_W_LC->Clone("FJM_WE_LCoverBE");
      //WBE1 = Form("Reconstructed_Tprime.M() >> FJM_W_LCoverBE%i(60,400,1600)",9);
      //WBE2 = Form("FJM_W_LCoverBE%i",9);
      //AnalysisChain.Draw(WBE1.c_str(),"(W_From_Higgs_Chi2.M()<60 || W_From_Higgs_Chi2.M()>120)");
      //FJM_W_LCoverBE = (TH1F*)gDirectory->Get(WBE2.c_str());
      //gPad->Close();      
      ///////////////////////////////////////////////
      ///////////////////////////////////////////////
      ///////////////////////////////////////////////
      A1 = Form("Reconstructed_Tprime.M() >> TprimeMassLC_LRHT%i(60,400,1600)",9);
      A2 = Form("TprimeMassLC_LRHT%i",9);
      AnalysisChain.Draw(A1.c_str(),"(Reconstructed_Higgs.Pt()>200 && Reconstructed_Top.Pt()>200) && (DeltaR_of_W_Higgs>=2.2 && DeltaR_of_W_Higgs<=3.5) && (TMath::Abs(First_Higgs_Jet.Phi()-Second_Higgs_Jet.Phi())<=1.2 && TMath::Abs(Top_Jet.Phi()-Reconstructed_W.Phi())<=1.2) && TMath::Abs(First_W_Jet.Phi()-Second_W_Jet.Phi())<=2.0 && (Relative_Mass>=0.3 && Relative_Mass<=0.5) && (Reconstructed_Higgs.M()>=110 && Reconstructed_Higgs.M()<=140) && (Top_From_Higgs_Chi2.M()<140 || Top_From_Higgs_Chi2.M()>230) && Relative_THT<0.65");
      FiveJetsMassLC_LRHT = (TH1F*)gDirectory->Get(A2.c_str());
      gPad->Close();
      //
      A1 = Form("Reconstructed_Tprime.M() >> TprimeMassLC_HRHT%i(60,400,1600)",9);
      A2 = Form("TprimeMassLC_HRHT%i",9);
      AnalysisChain.Draw(A1.c_str(),"(Reconstructed_Higgs.Pt()>200 && Reconstructed_Top.Pt()>200) && (DeltaR_of_W_Higgs>=2.2 && DeltaR_of_W_Higgs<=3.5) && (TMath::Abs(First_Higgs_Jet.Phi()-Second_Higgs_Jet.Phi())<=1.2 && TMath::Abs(Top_Jet.Phi()-Reconstructed_W.Phi())<=1.2) && TMath::Abs(First_W_Jet.Phi()-Second_W_Jet.Phi())<=2.0 && (Relative_Mass>=0.3 && Relative_Mass<=0.5) && (Reconstructed_Higgs.M()>=110 && Reconstructed_Higgs.M()<=140) && (Top_From_Higgs_Chi2.M()<140 || Top_From_Higgs_Chi2.M()>230) && Relative_THT>=0.65");
      FiveJetsMassLC_HRHT = (TH1F*)gDirectory->Get(A2.c_str());
      gPad->Close();
      //
      A1 = Form("Reconstructed_Tprime.M() >> TprimeMassLC_LRM%i(60,400,1600)",9);
      A2 = Form("TprimeMassLC_LRM%i",9);
      AnalysisChain.Draw(A1.c_str(),"(Reconstructed_Higgs.Pt()>200 && Reconstructed_Top.Pt()>200) && (DeltaR_of_W_Higgs>=2.2 && DeltaR_of_W_Higgs<=3.5) && (TMath::Abs(First_Higgs_Jet.Phi()-Second_Higgs_Jet.Phi())<=1.2 && TMath::Abs(Top_Jet.Phi()-Reconstructed_W.Phi())<=1.2) && TMath::Abs(First_W_Jet.Phi()-Second_W_Jet.Phi())<=2.0 && Relative_THT>=0.65 && (Reconstructed_Higgs.M()>=110 && Reconstructed_Higgs.M()<=140) && (Top_From_Higgs_Chi2.M()<140 || Top_From_Higgs_Chi2.M()>230) && Relative_Mass<0.3");
      FiveJetsMassLC_LRM = (TH1F*)gDirectory->Get(A2.c_str());
      gPad->Close();
      //
      A1 = Form("Reconstructed_Tprime.M() >> TprimeMassLC_HRM%i(60,400,1600)",9);
      A2 = Form("TprimeMassLC_HRM%i",9);
      AnalysisChain.Draw(A1.c_str(),"(Reconstructed_Higgs.Pt()>200 && Reconstructed_Top.Pt()>200) && (DeltaR_of_W_Higgs>=2.2 && DeltaR_of_W_Higgs<=3.5) && (TMath::Abs(First_Higgs_Jet.Phi()-Second_Higgs_Jet.Phi())<=1.2 && TMath::Abs(Top_Jet.Phi()-Reconstructed_W.Phi())<=1.2) && TMath::Abs(First_W_Jet.Phi()-Second_W_Jet.Phi())<=2.0 && Relative_THT>=0.65 && (Reconstructed_Higgs.M()>=110 && Reconstructed_Higgs.M()<=140) && (Top_From_Higgs_Chi2.M()<140 || Top_From_Higgs_Chi2.M()>230) && Relative_Mass>0.5");
      FiveJetsMassLC_HRM = (TH1F*)gDirectory->Get(A2.c_str());
      gPad->Close();
      //
      A1 = Form("Reconstructed_Tprime.M() >> TprimeMassLC_WRM%i(60,400,1600)",9);
      A2 = Form("TprimeMassLC_WRM%i",9);
      AnalysisChain.Draw(A1.c_str(),"(Reconstructed_Higgs.Pt()>200 && Reconstructed_Top.Pt()>200) && (DeltaR_of_W_Higgs>=2.2 && DeltaR_of_W_Higgs<=3.5) && (TMath::Abs(First_Higgs_Jet.Phi()-Second_Higgs_Jet.Phi())<=1.2 && TMath::Abs(Top_Jet.Phi()-Reconstructed_W.Phi())<=1.2) && TMath::Abs(First_W_Jet.Phi()-Second_W_Jet.Phi())<=2.0 && Relative_THT>=0.65 && (Reconstructed_Higgs.M()>=110 && Reconstructed_Higgs.M()<=140) && (Top_From_Higgs_Chi2.M()<140 || Top_From_Higgs_Chi2.M()>230) && (Relative_Mass>=0.3 && Relative_Mass<=0.5)");
      FiveJetsMassLC_WRM = (TH1F*)gDirectory->Get(A2.c_str());
      gPad->Close();      
      ///////////////////////////////////////////////
      ///////////////////////////////////////////////
      ///////////////////////////////////////////////
      string LJPT1 = Form("jet1_pt >> jet1_pt%i(50,20,1000)",9);
      string LJPT2 = Form("jet1_pt%i",9);
      AnalysisChain.Draw(LJPT1.c_str(),"(PUR_function(Number_True_Interactions))*weight");
      LeadingJetPT = (TH1F*)gDirectory->Get(LJPT2.c_str());
      gPad->Close();
      //
      string SLJPT1 = Form("jet2_pt >> jet2_pt%i(50,20,1000)",9);
      string SLJPT2 = Form("jet2_pt%i",9);
      AnalysisChain.Draw(SLJPT1.c_str(),"(PUR_function(Number_True_Interactions))*weight");
      Leading2JetPT = (TH1F*)gDirectory->Get(SLJPT2.c_str());
      gPad->Close();
      //
      string SSLJPT1 = Form("jet3_pt >> jet3_pt%i(35,20,700)",9);
      string SSLJPT2 = Form("jet3_pt%i",9);
      AnalysisChain.Draw(SSLJPT1.c_str(),"(PUR_function(Number_True_Interactions))*weight");
      Leading3JetPT = (TH1F*)gDirectory->Get(SSLJPT2.c_str());
      gPad->Close();
      //
      string SSSLJPT1 = Form("jet4_pt >> jet4_pt%i(20,20,400)",9);
      string SSSLJPT2 = Form("jet4_pt%i",9);
      AnalysisChain.Draw(SSSLJPT1.c_str(),"(PUR_function(Number_True_Interactions))*weight");
      Leading4JetPT = (TH1F*)gDirectory->Get(SSSLJPT2.c_str());
      gPad->Close();
      //
      string SSSSLJPT1 = Form("jet5_pt >> jet5_pt%i(15,20,300)",9);
      string SSSSLJPT2 = Form("jet5_pt%i",9);
      AnalysisChain.Draw(SSSSLJPT1.c_str(),"(PUR_function(Number_True_Interactions))*weight");
      Leading5JetPT = (TH1F*)gDirectory->Get(SSSSLJPT2.c_str());
      gPad->Close();
      //
      string SSSSSLJPT1 = Form("jet6_pt >> jet6_pt%i(10,20,200)",9);
      string SSSSSLJPT2 = Form("jet6_pt%i",9);
      AnalysisChain.Draw(SSSSSLJPT1.c_str(),"(PUR_function(Number_True_Interactions))*weight");
      Leading6JetPT = (TH1F*)gDirectory->Get(SSSSSLJPT2.c_str());
      gPad->Close();
      //
      string THT1 = Form("THT >> THT%i(65,300,1600)",9);
      string THT2 = Form("THT%i",9);
      AnalysisChain.Draw(THT1.c_str(),"(PUR_function(Number_True_Interactions))*weight","(PUR_function(Number_True_Interactions))*weight*(THT>630)"); //PUR_function(Number_True_Interactions)*(THT>630)");
      THT = (TH1F*)gDirectory->Get(THT2.c_str());
      gPad->Close(); 
      //
      string DRHJ1 = Form("DeltaR_of_Higgs_Jets >> DeltaR_of_Higgs_Jets%i(65,0.5,7)",9);
      string DRHJ2 = Form("DeltaR_of_Higgs_Jets%i",9);
      AnalysisChain.Draw(DRHJ1.c_str(),"(PUR_function(Number_True_Interactions))*weight");
      DRHjets = (TH1F*)gDirectory->Get(DRHJ2.c_str());
      gPad->Close();
      //
      string DRWJ1 = Form("DeltaR_of_W_Jets >> DeltaR_of_W_Jets%i(65,0.5,7)",9);
      string DRWJ2 = Form("DeltaR_of_W_Jets%i",9);
      AnalysisChain.Draw(DRWJ1.c_str(),"(PUR_function(Number_True_Interactions))*weight");
      DRWjets = (TH1F*)gDirectory->Get(DRWJ2.c_str());
      gPad->Close();
      //
      string HPT1 = Form("Reconstructed_Higgs.Pt() >> HPt%i(40,10,800)",9);
      string HPT2 = Form("HPt%i",9);
      AnalysisChain.Draw(HPT1.c_str(),"(PUR_function(Number_True_Interactions))*weight");
      Hpt = (TH1F*)gDirectory->Get(HPT2.c_str());
      gPad->Close();
      //
      string TPT1 = Form("Reconstructed_Top.Pt() >> TPt%i(40,10,800)",9);
      string TPT2 = Form("TPt%i",9);
      AnalysisChain.Draw(TPT1.c_str(),"(PUR_function(Number_True_Interactions))*weight");
      Tpt = (TH1F*)gDirectory->Get(TPT2.c_str());
      gPad->Close();
      //
      string DRWH1 = Form("DeltaR_of_W_Higgs >> DeltaR_of_W_Higgs%i(65,0.5,7)",9);
      string DRWH2 = Form("DeltaR_of_W_Higgs%i",9);
      AnalysisChain.Draw(DRWH1.c_str(),"(PUR_function(Number_True_Interactions))*weight");
      DRWH = (TH1F*)gDirectory->Get(DRWH2.c_str());
      gPad->Close();
      //
      string DPHJ1 = Form("TMath::Abs(First_Higgs_Jet.Phi()-Second_Higgs_Jet.Phi()) >> DeltaPhi_of_Higgs_jets%i(30,0.0,3.0)",9);
      string DPHJ2 = Form("DeltaPhi_of_Higgs_jets%i",9);
      AnalysisChain.Draw(DPHJ1.c_str(),"(PUR_function(Number_True_Interactions))*weight");
      DPHjets = (TH1F*)gDirectory->Get(DPHJ2.c_str());
      gPad->Close();
      //
      string DPWJ1 = Form("TMath::Abs(First_W_Jet.Phi()-Second_W_Jet.Phi()) >> DeltaPhi_of_W_jets%i(30,0.0,3.0)",9);
      string DPWJ2 = Form("DeltaPhi_of_W_jets%i",9);
      AnalysisChain.Draw(DPWJ1.c_str(),"(PUR_function(Number_True_Interactions))*weight");
      DPWjets = (TH1F*)gDirectory->Get(DPWJ2.c_str());
      gPad->Close();
      //
      string DPTJ1 = Form("TMath::Abs(Top_Jet.Phi()-Reconstructed_W.Phi()) >> DeltaPhi_of_T_jet%i(30,0.0,3.0)",9);
      string DPTJ2 = Form("DeltaPhi_of_T_jet%i",9);
      AnalysisChain.Draw(DPTJ1.c_str(),"(PUR_function(Number_True_Interactions))*weight");
      DPTjets = (TH1F*)gDirectory->Get(DPTJ2.c_str());
      gPad->Close();
      //
      string HM1 = Form("Reconstructed_Higgs.M() >> HM%i(36,60,180)",9);
      string HM2 = Form("HM%i",9);
      AnalysisChain.Draw(HM1.c_str(),"(PUR_function(Number_True_Interactions))*weight");
      HiggsMass = (TH1F*)gDirectory->Get(HM2.c_str());
      gPad->Close();
      //
      string RHT1 = Form("Relative_THT >> RelHT%i(30,0,1)",9);
      string RHT2 = Form("RelHT%i",9);
      AnalysisChain.Draw(RHT1.c_str(),"(PUR_function(Number_True_Interactions))*weight");
      RelHT = (TH1F*)gDirectory->Get(RHT2.c_str());
      gPad->Close();
      //
      string DRTH1 = Form("DeltaR_of_Top_Higgs >> DeltaR_of_Top_Higgs%i(65,0.5,7)",9);
      string DRTH2 = Form("DeltaR_of_Top_Higgs%i",9);
      AnalysisChain.Draw(DRTH1.c_str(),"(PUR_function(Number_True_Interactions))*weight");
      DRTH = (TH1F*)gDirectory->Get(DRTH2.c_str());
      gPad->Close();
      //
      string PTNM1 = Form("PT_Normalized_Mass >> PT_Normalized_Mass%i(70,0.4,5)",9);
      string PTNM2 = Form("PT_Normalized_Mass%i",9);
      AnalysisChain.Draw(PTNM1.c_str(),"(PUR_function(Number_True_Interactions))*weight");
      PtNormalizedMass = (TH1F*)gDirectory->Get(PTNM2.c_str());
      gPad->Close();
      //
      string RM1 = Form("Relative_Mass >> Relative_Mass%i(30,0.0,1)",9);
      string RM2 = Form("Relative_Mass%i",9);
      AnalysisChain.Draw(RM1.c_str(),"(PUR_function(Number_True_Interactions))*weight");
      RelativeMass = (TH1F*)gDirectory->Get(RM2.c_str());
      gPad->Close();
      //
      string MPTNM1 = Form("Mother_PT_Normalized_Mass >> Mother_PT_Normalized_Mass%i(25,0.0,50)",9);
      string MPTNM2 = Form("Mother_PT_Normalized_Mass%i",9);
      AnalysisChain.Draw(MPTNM1.c_str(),"(PUR_function(Number_True_Interactions))*weight");
      MotherPtNormalizedMass = (TH1F*)gDirectory->Get(MPTNM2.c_str());
      gPad->Close();
      //
      string NTops1 = Form("Number_of_Tops >> Number_of_Tops%i(8,0.0,8)",9);
      string NTops2 = Form("Number_of_Tops%i",9);
      AnalysisChain.Draw(NTops1.c_str(),"(PUR_function(Number_True_Interactions))*weight");
      NumberOfTops = (TH1F*)gDirectory->Get(NTops2.c_str());
      gPad->Close();
      ///////////NEW VARIABLES//////////////////
      string HMTM1 = Form("Reconstructed_Higgs.M()/Reconstructed_Top.M() >> HMoverTM%i(30,0,1.)",9);
      string HMTM2 = Form("HMoverTM%i",9);
      AnalysisChain.Draw(HMTM1.c_str(),"(PUR_function(Number_True_Interactions))*weight");
      HiggsMassOverTopMass = (TH1F*)gDirectory->Get(HMTM2.c_str());
      gPad->Close();
      //
      string HTA1 = Form("((Reconstructed_Higgs.Pt()/Reconstructed_Higgs.M())-(Reconstructed_Top.Pt()/Reconstructed_Top.M()))/((Reconstructed_Higgs.Pt()/Reconstructed_Higgs.M())+(Reconstructed_Top.Pt()/Reconstructed_Top.M())) >> HTAsym%i(12,0,1.)",9);
      string HTA2 = Form("HTAsym%i",9);
      AnalysisChain.Draw(HTA1.c_str(),"(PUR_function(Number_True_Interactions))*weight");
      HiggsTopAsymmetry = (TH1F*)gDirectory->Get(HTA2.c_str());
      gPad->Close();
      //
      string TLBT1 = Form("Number_of_Loose_and_non_med_b_tags >> TLBTag%i(10,0,10)",9);
      string TLBT2 = Form("TLBTag%i",9);
      AnalysisChain.Draw(TLBT1.c_str(),"(PUR_function(Number_True_Interactions))*weight");
      ThirdLooseBtag = (TH1F*)gDirectory->Get(TLBT2.c_str());
      gPad->Close();
      //
      string TM1 = Form("Reconstructed_Top.M() >> TMass%i(60,100,700)",9);
      string TM2 = Form("TMass%i",9);
      AnalysisChain.Draw(TM1.c_str(),"(PUR_function(Number_True_Interactions))*weight");
      TopMass = (TH1F*)gDirectory->Get(TM2.c_str());
      gPad->Close();
      //
      string C21 = Form("ChiSquaredSorting >> ChiSq%i(100,0,1000)",9);
      string C22 = Form("ChiSq%i",9);
      AnalysisChain.Draw(C21.c_str(),"(PUR_function(Number_True_Interactions))*weight");
      Chi2 = (TH1F*)gDirectory->Get(C22.c_str());
      gPad->Close();
      //
      string UQ1 = Form("U Quark Content >> UQC%i(10,0,10)",9);
      string UQ2 = Form("UQC%i",9);
      AnalysisChain.Draw(UQ1.c_str(),"(PUR_function(Number_True_Interactions))*weight");
      UQuarkContent = (TH1F*)gDirectory->Get(UQ2.c_str());
      gPad->Close();
      //
      string DQ1 = Form("D Quark Content >> DQC%i(10,0,10)",9);
      string DQ2 = Form("DQC%i",9);
      AnalysisChain.Draw(DQ1.c_str(),"(PUR_function(Number_True_Interactions))*weight");
      DQuarkContent = (TH1F*)gDirectory->Get(DQ2.c_str());
      gPad->Close();
      //
      string SQ1 = Form("S Quark Content >> SQC%i(10,0,10)",9);
      string SQ2 = Form("SQC%i",9);
      AnalysisChain.Draw(SQ1.c_str(),"(PUR_function(Number_True_Interactions))*weight");
      SQuarkContent = (TH1F*)gDirectory->Get(SQ2.c_str());
      gPad->Close();
      //
      string CQ1 = Form("C Quark Content >> CQC%i(10,0,10)",9);
      string CQ2 = Form("CQC%i",9);
      AnalysisChain.Draw(CQ1.c_str(),"(PUR_function(Number_True_Interactions))*weight");
      CQuarkContent = (TH1F*)gDirectory->Get(CQ2.c_str());
      gPad->Close();
      //
      string BQ1 = Form("B Quark Content >> BQC%i(10,0,10)",9);
      string BQ2 = Form("BQC%i",9);
      AnalysisChain.Draw(BQ1.c_str(),"(PUR_function(Number_True_Interactions))*weight");
      BQuarkContent = (TH1F*)gDirectory->Get(BQ2.c_str());
      gPad->Close();
      //B-tagging Working point	  
      string BTL1 = Form("Number_CSVLbtagged_jets >> CSVL%i(10,0,10)",9);
      string BTL2 = Form("CSVL%i",9);
      AnalysisChain.Draw(BTL1.c_str(),"(PUR_function(Number_True_Interactions))*weight");
      CSVLB = (TH1F*)gDirectory->Get(BTL2.c_str());
      gPad->Close();
      //	  
      string BTM1 = Form("Number_CSVMbtagged_jets >> CSVM%i(10,0,10)",9);
      string BTM2 = Form("CSVM%i",9);
      AnalysisChain.Draw(BTM1.c_str(),"(PUR_function(Number_True_Interactions))*weight");
      CSVMB = (TH1F*)gDirectory->Get(BTM2.c_str());
      gPad->Close();
      //	  
      string BTT1 = Form("Number_CSVTbtagged_jets >> CSVT%i(10,0,10)",9);
      string BTT2 = Form("CSVT%i",9);
      AnalysisChain.Draw(BTT1.c_str(),"(PUR_function(Number_True_Interactions))*weight");
      CSVTB = (TH1F*)gDirectory->Get(BTT2.c_str());
      gPad->Close();
      //	  
      string VT1 = Form("Vertices >> VTX%i(40,1,41)",9);
      string VT2 = Form("VTX%i",9);
      AnalysisChain.Draw(VT1.c_str(),"(PUR_function(Number_True_Interactions))*weight*(THT>630)"); //PUR_function(Number_True_Interactions)*(THT>630)");
      Vtcs = (TH1F*)gDirectory->Get(VT2.c_str());
      gPad->Close();
      //
      /*string CLB1 = Form("CorrectLB_tag >> CLooseB%i(6,0,6)",9);
      string CLB2 = Form("CLooseB%i",9);
      AnalysisChain.Draw(CLB1.c_str(),"(PUR_function(Number_True_Interactions))*weight");
      CLBT = (TH1F*)gDirectory->Get(CLB2.c_str());
      gPad->Close();
      //
      string CMB1 = Form("CorrectMB_tag >> CMedB%i(6,0,6)",9);
      string CMB2 = Form("CMedB%i",9);
      AnalysisChain.Draw(CMB1.c_str(),"(PUR_function(Number_True_Interactions))*weight");
      CMBT = (TH1F*)gDirectory->Get(CMB2.c_str());
      gPad->Close();
      //
      string CTB1 = Form("CorrectTB_tag >> CTightB%i(6,0,6)",9);
      string CTB2 = Form("CTightB%i",9);
      AnalysisChain.Draw(CTB1.c_str(),"(PUR_function(Number_True_Interactions))*weight");
      CTBT = (TH1F*)gDirectory->Get(CTB2.c_str());
      gPad->Close();
      //
      string FLLB1 = Form("FakeLB_tag_Light >> FLLooseB%i(6,0,6)",9);
      string FLLB2 = Form("FLLooseB%i",9);
      AnalysisChain.Draw(FLLB1.c_str(),"(PUR_function(Number_True_Interactions))*weight");
      FLLBT = (TH1F*)gDirectory->Get(FLLB2.c_str());
      gPad->Close();
      //
      string FLMB1 = Form("FakeMB_tag_Light >> FLMedB%i(6,0,6)",9);
      string FLMB2 = Form("FLMedB%i",9);
      AnalysisChain.Draw(FLMB1.c_str(),"(PUR_function(Number_True_Interactions))*weight");
      FLMBT = (TH1F*)gDirectory->Get(FLMB2.c_str());
      gPad->Close();
      //
      string FLTB1 = Form("FakeTB_tag_Light >> FLTightB%i(6,0,6)",9);
      string FLTB2 = Form("FLTightB%i",9);
      AnalysisChain.Draw(FLTB1.c_str(),"(PUR_function(Number_True_Interactions))*weight");
      FLTBT = (TH1F*)gDirectory->Get(FLTB2.c_str());
      gPad->Close();
      //
      string FCLB1 = Form("FakeLB_tag_C >> FCLooseB%i(6,0,6)",9);
      string FCLB2 = Form("FCLooseB%i",9);
      AnalysisChain.Draw(FCLB1.c_str(),"(PUR_function(Number_True_Interactions))*weight");
      FCLBT = (TH1F*)gDirectory->Get(FCLB2.c_str());
      gPad->Close();
      //
      string FCMB1 = Form("FakeMB_tag_C >> FCMedB%i(6,0,6)",9);
      string FCMB2 = Form("FCMedB%i",9);
      AnalysisChain.Draw(FCMB1.c_str(),"(PUR_function(Number_True_Interactions))*weight");
      FCMBT = (TH1F*)gDirectory->Get(FCMB2.c_str());
      gPad->Close();
      //
      string FCTB1 = Form("FakeTB_tag_C >> FCTightB%i(6,0,6)",9);
      string FCTB2 = Form("FCTightB%i",9);
      AnalysisChain.Draw(FCTB1.c_str(),"(PUR_function(Number_True_Interactions))*weight");
      FCTBT = (TH1F*)gDirectory->Get(FCTB2.c_str());
      gPad->Close();*/
      //
      string HTCSVM1 = Form("Number_CSVMbtagged_jets : THT >> HT_CSVM%i(65,300,1600,10)",9);
      string HTCSVM2 = Form("HT_CSVM%i",9);
      AnalysisChain.Draw(HTCSVM1.c_str(),"(PUR_function(Number_True_Interactions))*weight");
      HTCSVM = (TH2F*)gDirectory->Get(HTCSVM2.c_str());
      cout << "Correlation Factor: " << HTCSVM->GetCorrelationFactor() << endl;
      gPad->Close();
      //////////////////////////////////////
      /////////ABCD BKG Estimation//////////
      //////////////////////////////////////
      HM1 = Form("Reconstructed_Higgs.M() >> HMBE0%i(36,60,180)",9);
      HM2 = Form("HMBE0%i",9);
      AnalysisChain.Draw(HM1.c_str(),"Reconstructed_Higgs.Pt()<200 || Reconstructed_Top.Pt()<200");
      HiggsMassReversedHptTpt[0] = (TH1F*)gDirectory->Get(HM2.c_str());
      gPad->Close();
      //	  
      A1 = Form("Reconstructed_Tprime.M() >> TprimeMassBE0%i(60,400,1600)",9);
      A2 = Form("TprimeMassBE0%i",9);
      AnalysisChain.Draw(A1.c_str(),"Reconstructed_Higgs.Pt()<200 || Reconstructed_Top.Pt()<200");
      TprimeMassReversedHptTpt[0] = (TH1F*)gDirectory->Get(A2.c_str());
      gPad->Close();
      //	  
      HM1 = Form("Reconstructed_Higgs.M() >> HMBE1%i(36,60,180)",9);
      HM2 = Form("HMBE1%i",9);
      AnalysisChain.Draw(HM1.c_str(),"(Reconstructed_Higgs.Pt()<200 || Reconstructed_Top.Pt()<200) && (DeltaR_of_W_Higgs>=2.2 && DeltaR_of_W_Higgs<=3.5)");
      HiggsMassReversedHptTpt[1] = (TH1F*)gDirectory->Get(HM2.c_str());
      gPad->Close();
      //	  
      A1 = Form("Reconstructed_Tprime.M() >> TprimeMassBE1%i(60,400,1600)",9);
      A2 = Form("TprimeMassBE1%i",9);
      AnalysisChain.Draw(A1.c_str(),"(Reconstructed_Higgs.Pt()<200 || Reconstructed_Top.Pt()<200) && (DeltaR_of_W_Higgs>=2.2 && DeltaR_of_W_Higgs<=3.5)");
      TprimeMassReversedHptTpt[1] = (TH1F*)gDirectory->Get(A2.c_str());
      gPad->Close();
      //	  
      HM1 = Form("Reconstructed_Higgs.M() >> HMBE2%i(36,60,180)",9);
      HM2 = Form("HMBE2%i",9);
      AnalysisChain.Draw(HM1.c_str(),"(Reconstructed_Higgs.Pt()<200 || Reconstructed_Top.Pt()<200) && (DeltaR_of_W_Higgs>=2.2 && DeltaR_of_W_Higgs<=3.5) && (TMath::Abs(First_Higgs_Jet.Phi()-Second_Higgs_Jet.Phi())<=1.2 && TMath::Abs(Top_Jet.Phi()-Reconstructed_W.Phi())<=1.2)");
      HiggsMassReversedHptTpt[2] = (TH1F*)gDirectory->Get(HM2.c_str());
      gPad->Close();
      //	  
      A1 = Form("Reconstructed_Tprime.M() >> TprimeMassBE2%i(60,400,1600)",9);
      A2 = Form("TprimeMassBE2%i",9);
      AnalysisChain.Draw(A1.c_str(),"(Reconstructed_Higgs.Pt()<200 || Reconstructed_Top.Pt()<200) && (DeltaR_of_W_Higgs>=2.2 && DeltaR_of_W_Higgs<=3.5) && (TMath::Abs(First_Higgs_Jet.Phi()-Second_Higgs_Jet.Phi())<=1.2 && TMath::Abs(Top_Jet.Phi()-Reconstructed_W.Phi())<=1.2)");
      TprimeMassReversedHptTpt[2] = (TH1F*)gDirectory->Get(A2.c_str());
      gPad->Close();
      //	  
      HM1 = Form("Reconstructed_Higgs.M() >> HMBE3%i(36,60,180)",9);
      HM2 = Form("HMBE3%i",9);
      AnalysisChain.Draw(HM1.c_str(),"(Reconstructed_Higgs.Pt()<200 || Reconstructed_Top.Pt()<200) && (DeltaR_of_W_Higgs>=2.2 && DeltaR_of_W_Higgs<=3.5) && (TMath::Abs(First_Higgs_Jet.Phi()-Second_Higgs_Jet.Phi())<=1.2 && TMath::Abs(Top_Jet.Phi()-Reconstructed_W.Phi())<=1.2) && TMath::Abs(First_W_Jet.Phi()-Second_W_Jet.Phi())<=2.0");
      HiggsMassReversedHptTpt[3] = (TH1F*)gDirectory->Get(HM2.c_str());
      gPad->Close();
      //	  
      A1 = Form("Reconstructed_Tprime.M() >> TprimeMassBE3%i(60,400,1600)",9);
      A2 = Form("TprimeMassBE3%i",9);
      AnalysisChain.Draw(A1.c_str(),"(Reconstructed_Higgs.Pt()<200 || Reconstructed_Top.Pt()<200) && (DeltaR_of_W_Higgs>=2.2 && DeltaR_of_W_Higgs<=3.5) && (TMath::Abs(First_Higgs_Jet.Phi()-Second_Higgs_Jet.Phi())<=1.2 && TMath::Abs(Top_Jet.Phi()-Reconstructed_W.Phi())<=1.2) && TMath::Abs(First_W_Jet.Phi()-Second_W_Jet.Phi())<=2.0");
      TprimeMassReversedHptTpt[3] = (TH1F*)gDirectory->Get(A2.c_str());
      gPad->Close();
      //	  
      HM1 = Form("Reconstructed_Higgs.M() >> HMBE4%i(36,60,180)",9);
      HM2 = Form("HMBE4%i",9);
      AnalysisChain.Draw(HM1.c_str(),"(Reconstructed_Higgs.Pt()<200 || Reconstructed_Top.Pt()<200) && (DeltaR_of_W_Higgs>=2.2 && DeltaR_of_W_Higgs<=3.5) && (TMath::Abs(First_Higgs_Jet.Phi()-Second_Higgs_Jet.Phi())<=1.2 && TMath::Abs(Top_Jet.Phi()-Reconstructed_W.Phi())<=1.2) && TMath::Abs(First_W_Jet.Phi()-Second_W_Jet.Phi())<=2.0 && Relative_THT>=0.65");
      HiggsMassReversedHptTpt[4] = (TH1F*)gDirectory->Get(HM2.c_str());
      gPad->Close();
      //	  
      A1 = Form("Reconstructed_Tprime.M() >> TprimeMassBE4%i(60,400,1600)",9);
      A2 = Form("TprimeMassBE4%i",9);
      AnalysisChain.Draw(A1.c_str(),"(Reconstructed_Higgs.Pt()<200 || Reconstructed_Top.Pt()<200) && (DeltaR_of_W_Higgs>=2.2 && DeltaR_of_W_Higgs<=3.5) && (TMath::Abs(First_Higgs_Jet.Phi()-Second_Higgs_Jet.Phi())<=1.2 && TMath::Abs(Top_Jet.Phi()-Reconstructed_W.Phi())<=1.2) && TMath::Abs(First_W_Jet.Phi()-Second_W_Jet.Phi())<=2.0 && Relative_THT>=0.65");
      TprimeMassReversedHptTpt[4] = (TH1F*)gDirectory->Get(A2.c_str());
      //	  
      HM1 = Form("Reconstructed_Higgs.M() >> HMBE5%i(36,60,180)",9);
      HM2 = Form("HMBE5%i",9);
      AnalysisChain.Draw(HM1.c_str(),"(Reconstructed_Higgs.Pt()<200 || Reconstructed_Top.Pt()<200) && (DeltaR_of_W_Higgs>=2.2 && DeltaR_of_W_Higgs<=3.5) && (TMath::Abs(First_Higgs_Jet.Phi()-Second_Higgs_Jet.Phi())<=1.2 && TMath::Abs(Top_Jet.Phi()-Reconstructed_W.Phi())<=1.2) && TMath::Abs(First_W_Jet.Phi()-Second_W_Jet.Phi())<=2.0 && Relative_THT>=0.65 && (Relative_Mass>=0.3 && Relative_Mass<=0.5)");
      HiggsMassReversedHptTpt[5] = (TH1F*)gDirectory->Get(HM2.c_str());
      gPad->Close();
      //	  
      A1 = Form("Reconstructed_Tprime.M() >> TprimeMassBE5%i(60,400,1600)",9);
      A2 = Form("TprimeMassBE5%i",9);
      AnalysisChain.Draw(A1.c_str(),"(Reconstructed_Higgs.Pt()<200 || Reconstructed_Top.Pt()<200) && (DeltaR_of_W_Higgs>=2.2 && DeltaR_of_W_Higgs<=3.5) && (TMath::Abs(First_Higgs_Jet.Phi()-Second_Higgs_Jet.Phi())<=1.2 && TMath::Abs(Top_Jet.Phi()-Reconstructed_W.Phi())<=1.2) && TMath::Abs(First_W_Jet.Phi()-Second_W_Jet.Phi())<=2.0 && Relative_THT>=0.65 && (Relative_Mass>=0.3 && Relative_Mass<=0.5)");
      TprimeMassReversedHptTpt[5] = (TH1F*)gDirectory->Get(A2.c_str());
      gPad->Close();
      //
      string HPTTPT1 = Form("Reconstructed_Higgs.Pt():Reconstructed_Top.Pt() >> HPtTPt0%i(40,10,800,40,10,800)",9);
      string HPTTPT2 = Form("HPtTPt0%i",9);
      AnalysisChain.Draw(HPTTPT1.c_str());
      HptTpt[0] = (TH2F*)gDirectory->Get(HPTTPT2.c_str());
      cout << "Correlation Factor: " << HptTpt[0]->GetCorrelationFactor() << endl;
      gPad->Close();
      //
      HPTTPT1 = Form("Reconstructed_Higgs.Pt():Reconstructed_Top.Pt() >> HPtTPt1%i(40,10,800,40,10,800)",9);
      HPTTPT2 = Form("HPtTPt1%i",9);
      AnalysisChain.Draw(HPTTPT1.c_str(),"(DeltaR_of_W_Higgs>=2.2 && DeltaR_of_W_Higgs<=3.5)");
      HptTpt[1] = (TH2F*)gDirectory->Get(HPTTPT2.c_str());
      cout << "Correlation Factor: " << HptTpt[1]->GetCorrelationFactor() << endl;
      gPad->Close();
      //
      HPTTPT1 = Form("Reconstructed_Higgs.Pt():Reconstructed_Top.Pt() >> HPtTPt2%i(40,10,800,40,10,800)",9);
      HPTTPT2 = Form("HPtTPt2%i",9);
      AnalysisChain.Draw(HPTTPT1.c_str(),"(DeltaR_of_W_Higgs>=2.2 && DeltaR_of_W_Higgs<=3.5) && (TMath::Abs(First_Higgs_Jet.Phi()-Second_Higgs_Jet.Phi())<=1.2 && TMath::Abs(Top_Jet.Phi()-Reconstructed_W.Phi())<=1.2)");
      HptTpt[2] = (TH2F*)gDirectory->Get(HPTTPT2.c_str());
      cout << "Correlation Factor: " << HptTpt[2]->GetCorrelationFactor() << endl;
      gPad->Close();
      //
      HPTTPT1 = Form("Reconstructed_Higgs.Pt():Reconstructed_Top.Pt() >> HPtTPt3%i(40,10,800,40,10,800)",9);
      HPTTPT2 = Form("HPtTPt3%i",9);
      AnalysisChain.Draw(HPTTPT1.c_str(),"(DeltaR_of_W_Higgs>=2.2 && DeltaR_of_W_Higgs<=3.5) && (TMath::Abs(First_Higgs_Jet.Phi()-Second_Higgs_Jet.Phi())<=1.2 && TMath::Abs(Top_Jet.Phi()-Reconstructed_W.Phi())<=1.2) && TMath::Abs(First_W_Jet.Phi()-Second_W_Jet.Phi())<=2.0");
      HptTpt[3] = (TH2F*)gDirectory->Get(HPTTPT2.c_str());
      cout << "Correlation Factor: " << HptTpt[3]->GetCorrelationFactor() << endl;
      gPad->Close();
      //
      HPTTPT1 = Form("Reconstructed_Higgs.Pt():Reconstructed_Top.Pt() >> HPtTPt4%i(40,10,800,40,10,800)",9);
      HPTTPT2 = Form("HPtTPt4%i",9);
      AnalysisChain.Draw(HPTTPT1.c_str(),"(DeltaR_of_W_Higgs>=2.2 && DeltaR_of_W_Higgs<=3.5) && (TMath::Abs(First_Higgs_Jet.Phi()-Second_Higgs_Jet.Phi())<=1.2 && TMath::Abs(Top_Jet.Phi()-Reconstructed_W.Phi())<=1.2) && TMath::Abs(First_W_Jet.Phi()-Second_W_Jet.Phi())<=2.0 && Relative_THT>=0.65");
      HptTpt[4] = (TH2F*)gDirectory->Get(HPTTPT2.c_str());
      cout << "Correlation Factor: " << HptTpt[4]->GetCorrelationFactor() << endl;
      gPad->Close();
      //
      HPTTPT1 = Form("Reconstructed_Higgs.Pt():Reconstructed_Top.Pt() >> HPtTPt5%i(40,10,800,40,10,800)",9);
      HPTTPT2 = Form("HPtTPt5%i",9);
      AnalysisChain.Draw(HPTTPT1.c_str(),"(DeltaR_of_W_Higgs>=2.2 && DeltaR_of_W_Higgs<=3.5) && (TMath::Abs(First_Higgs_Jet.Phi()-Second_Higgs_Jet.Phi())<=1.2 && TMath::Abs(Top_Jet.Phi()-Reconstructed_W.Phi())<=1.2) && TMath::Abs(First_W_Jet.Phi()-Second_W_Jet.Phi())<=2.0 && Relative_THT>=0.65 && (Relative_Mass>=0.3 && Relative_Mass<=0.5)");
      HptTpt[5] = (TH2F*)gDirectory->Get(HPTTPT2.c_str());
      cout << "Correlation Factor: " << HptTpt[5]->GetCorrelationFactor() << endl;
      gPad->Close(); //(Reconstructed_Higgs.Pt()<200 || Reconstructed_Top.Pt()<200) && 
      //
      HM1 = Form("Reconstructed_Higgs.M() >> HM0%i(36,60,180)",9);
      HM2 = Form("HM0%i",9);
      AnalysisChain.Draw(HM1.c_str(),"Reconstructed_Higgs.Pt()>200 && Reconstructed_Top.Pt()>200");
      HiggsMassCuts[0] = (TH1F*)gDirectory->Get(HM2.c_str());
      gPad->Close();
      //	  
      HM1 = Form("Reconstructed_Higgs.M() >> HM1%i(36,60,180)",9);
      HM2 = Form("HM1%i",9);
      AnalysisChain.Draw(HM1.c_str(),"(Reconstructed_Higgs.Pt()>200 && Reconstructed_Top.Pt()>200) && (DeltaR_of_W_Higgs>=2.2 && DeltaR_of_W_Higgs<=3.5)");
      HiggsMassCuts[1] = (TH1F*)gDirectory->Get(HM2.c_str());
      gPad->Close();
      //	  
      HM1 = Form("Reconstructed_Higgs.M() >> HM2%i(36,60,180)",9);
      HM2 = Form("HM2%i",9);
      AnalysisChain.Draw(HM1.c_str(),"(Reconstructed_Higgs.Pt()>200 && Reconstructed_Top.Pt()>200) && (DeltaR_of_W_Higgs>=2.2 && DeltaR_of_W_Higgs<=3.5) && (TMath::Abs(First_Higgs_Jet.Phi()-Second_Higgs_Jet.Phi())<=1.2 && TMath::Abs(Top_Jet.Phi()-Reconstructed_W.Phi())<=1.2)");
      HiggsMassCuts[2] = (TH1F*)gDirectory->Get(HM2.c_str());
      gPad->Close();
      //	  
      HM1 = Form("Reconstructed_Higgs.M() >> HM3%i(36,60,180)",9);
      HM2 = Form("HM3%i",9);
      AnalysisChain.Draw(HM1.c_str(),"(Reconstructed_Higgs.Pt()>200 && Reconstructed_Top.Pt()>200) && (DeltaR_of_W_Higgs>=2.2 && DeltaR_of_W_Higgs<=3.5) && (TMath::Abs(First_Higgs_Jet.Phi()-Second_Higgs_Jet.Phi())<=1.2 && TMath::Abs(Top_Jet.Phi()-Reconstructed_W.Phi())<=1.2) && TMath::Abs(First_W_Jet.Phi()-Second_W_Jet.Phi())<=2.0");
      HiggsMassCuts[3] = (TH1F*)gDirectory->Get(HM2.c_str());
      gPad->Close();
      //
      HM1 = Form("Reconstructed_Higgs.M() >> HM4%i(36,60,180)",9);
      HM2 = Form("HM4%i",9);
      AnalysisChain.Draw(HM1.c_str(),"(Reconstructed_Higgs.Pt()>200 && Reconstructed_Top.Pt()>200) && (DeltaR_of_W_Higgs>=2.2 && DeltaR_of_W_Higgs<=3.5) && (TMath::Abs(First_Higgs_Jet.Phi()-Second_Higgs_Jet.Phi())<=1.2 && TMath::Abs(Top_Jet.Phi()-Reconstructed_W.Phi())<=1.2) && TMath::Abs(First_W_Jet.Phi()-Second_W_Jet.Phi())<=2.0 && Relative_THT>=0.65");
      HiggsMassCuts[4] = (TH1F*)gDirectory->Get(HM2.c_str());
      gPad->Close();
      //
      HM1 = Form("Reconstructed_Higgs.M() >> HM5%i(36,60,180)",9);
      HM2 = Form("HM5%i",9);
      AnalysisChain.Draw(HM1.c_str(),"(Reconstructed_Higgs.Pt()>200 && Reconstructed_Top.Pt()>200) && (DeltaR_of_W_Higgs>=2.2 && DeltaR_of_W_Higgs<=3.5) && (TMath::Abs(First_Higgs_Jet.Phi()-Second_Higgs_Jet.Phi())<=1.2 && TMath::Abs(Top_Jet.Phi()-Reconstructed_W.Phi())<=1.2) && TMath::Abs(First_W_Jet.Phi()-Second_W_Jet.Phi())<=2.0 && Relative_THT>=0.65 && (Relative_Mass>=0.3 && Relative_Mass<=0.5)");
      HiggsMassCuts[5] = (TH1F*)gDirectory->Get(HM2.c_str());
      gPad->Close();
      ///////////////////////
      //Sideband estimation//
      ///////////////////////
      string WFH1 = Form("W_From_Higgs.M() >> WMFH0%i(44,60.0,500)",9);
      string WFH2 = Form("WMFH0%i",9);
      AnalysisChain.Draw(WFH1.c_str(),"Reconstructed_Higgs.Pt()>200 && Reconstructed_Top.Pt()>200");
      WMassFromHiggs[0] = (TH1F*)gDirectory->Get(WFH2.c_str());
      gPad->Close();
      //
      WFH1 = Form("W_From_Higgs.M() >> WMFH1%i(44,60.0,500)",9);
      WFH2 = Form("WMFH1%i",9);
      AnalysisChain.Draw(WFH1.c_str(),"(Reconstructed_Higgs.Pt()>200 && Reconstructed_Top.Pt()>200) && (DeltaR_of_W_Higgs>=2.2 && DeltaR_of_W_Higgs<=3.5)");
      WMassFromHiggs[1] = (TH1F*)gDirectory->Get(WFH2.c_str());
      gPad->Close();
      //
      WFH1 = Form("W_From_Higgs.M() >> WMFH2%i(44,60.0,500)",9);
      WFH2 = Form("WMFH2%i",9);
      AnalysisChain.Draw(WFH1.c_str(),"(Reconstructed_Higgs.Pt()>200 && Reconstructed_Top.Pt()>200) && (DeltaR_of_W_Higgs>=2.2 && DeltaR_of_W_Higgs<=3.5) && (TMath::Abs(First_Higgs_Jet.Phi()-Second_Higgs_Jet.Phi())<=1.2 && TMath::Abs(Top_Jet.Phi()-Reconstructed_W.Phi())<=1.2)");
      WMassFromHiggs[2] = (TH1F*)gDirectory->Get(WFH2.c_str());
      gPad->Close();
      //
      WFH1 = Form("W_From_Higgs.M() >> WMFH3%i(44,60.0,500)",9);
      WFH2 = Form("WMFH3%i",9);
      AnalysisChain.Draw(WFH1.c_str(),"(Reconstructed_Higgs.Pt()>200 && Reconstructed_Top.Pt()>200) && (DeltaR_of_W_Higgs>=2.2 && DeltaR_of_W_Higgs<=3.5) && (TMath::Abs(First_Higgs_Jet.Phi()-Second_Higgs_Jet.Phi())<=1.2 && TMath::Abs(Top_Jet.Phi()-Reconstructed_W.Phi())<=1.2) && TMath::Abs(First_W_Jet.Phi()-Second_W_Jet.Phi())<=2.0");
      WMassFromHiggs[3] = (TH1F*)gDirectory->Get(WFH2.c_str());
      gPad->Close();
      //
      WFH1 = Form("W_From_Higgs.M() >> WMFH4%i(44,60.0,500)",9);
      WFH2 = Form("WMFH4%i",9);
      AnalysisChain.Draw(WFH1.c_str(),"(Reconstructed_Higgs.Pt()>200 && Reconstructed_Top.Pt()>200) && (DeltaR_of_W_Higgs>=2.2 && DeltaR_of_W_Higgs<=3.5) && (TMath::Abs(First_Higgs_Jet.Phi()-Second_Higgs_Jet.Phi())<=1.2 && TMath::Abs(Top_Jet.Phi()-Reconstructed_W.Phi())<=1.2) && TMath::Abs(First_W_Jet.Phi()-Second_W_Jet.Phi())<=2.0 && Relative_THT>=0.65");
      WMassFromHiggs[4] = (TH1F*)gDirectory->Get(WFH2.c_str());
      gPad->Close();
      //
      WFH1 = Form("W_From_Higgs.M() >> WMFH5%i(44,60.0,500)",9);
      WFH2 = Form("WMFH5%i",9);
      AnalysisChain.Draw(WFH1.c_str(),"(Reconstructed_Higgs.Pt()>200 && Reconstructed_Top.Pt()>200) && (DeltaR_of_W_Higgs>=2.2 && DeltaR_of_W_Higgs<=3.5) && (TMath::Abs(First_Higgs_Jet.Phi()-Second_Higgs_Jet.Phi())<=1.2 && TMath::Abs(Top_Jet.Phi()-Reconstructed_W.Phi())<=1.2) && TMath::Abs(First_W_Jet.Phi()-Second_W_Jet.Phi())<=2.0 && Relative_THT>=0.65 && (Relative_Mass>=0.3 && Relative_Mass<=0.5)");
      WMassFromHiggs[5] = (TH1F*)gDirectory->Get(WFH2.c_str());
      gPad->Close();
      //Chi2///////
      WFH1 = Form("W_From_Higgs_Chi2.M() >> WMFHC0%i(44,60.0,500)",9);
      WFH2 = Form("WMFHC0%i",9);
      AnalysisChain.Draw(WFH1.c_str(),"Reconstructed_Higgs.Pt()>200 && Reconstructed_Top.Pt()>200");
      WMassFromHiggsChi2[0] = (TH1F*)gDirectory->Get(WFH2.c_str());
      gPad->Close();
      //
      WFH1 = Form("W_From_Higgs_Chi2.M() >> WMFHC1%i(44,60.0,500)",9);
      WFH2 = Form("WMFHC1%i",9);
      AnalysisChain.Draw(WFH1.c_str(),"(Reconstructed_Higgs.Pt()>200 && Reconstructed_Top.Pt()>200) && (DeltaR_of_W_Higgs>=2.2 && DeltaR_of_W_Higgs<=3.5)");
      WMassFromHiggsChi2[1] = (TH1F*)gDirectory->Get(WFH2.c_str());
      gPad->Close();
      //
      WFH1 = Form("W_From_Higgs_Chi2.M() >> WMFHC2%i(44,60.0,500)",9);
      WFH2 = Form("WMFHC2%i",9);
      AnalysisChain.Draw(WFH1.c_str(),"(Reconstructed_Higgs.Pt()>200 && Reconstructed_Top.Pt()>200) && (DeltaR_of_W_Higgs>=2.2 && DeltaR_of_W_Higgs<=3.5) && (TMath::Abs(First_Higgs_Jet.Phi()-Second_Higgs_Jet.Phi())<=1.2 && TMath::Abs(Top_Jet.Phi()-Reconstructed_W.Phi())<=1.2)");
      WMassFromHiggsChi2[2] = (TH1F*)gDirectory->Get(WFH2.c_str());
      gPad->Close();
      //
      WFH1 = Form("W_From_Higgs_Chi2.M() >> WMFHC3%i(44,60.0,500)",9);
      WFH2 = Form("WMFHC3%i",9);
      AnalysisChain.Draw(WFH1.c_str(),"(Reconstructed_Higgs.Pt()>200 && Reconstructed_Top.Pt()>200) && (DeltaR_of_W_Higgs>=2.2 && DeltaR_of_W_Higgs<=3.5) && (TMath::Abs(First_Higgs_Jet.Phi()-Second_Higgs_Jet.Phi())<=1.2 && TMath::Abs(Top_Jet.Phi()-Reconstructed_W.Phi())<=1.2) && TMath::Abs(First_W_Jet.Phi()-Second_W_Jet.Phi())<=2.0");
      WMassFromHiggsChi2[3] = (TH1F*)gDirectory->Get(WFH2.c_str());
      gPad->Close();
      //
      WFH1 = Form("W_From_Higgs_Chi2.M() >> WMFHC4%i(44,60.0,500)",9);
      WFH2 = Form("WMFHC4%i",9);
      AnalysisChain.Draw(WFH1.c_str(),"(Reconstructed_Higgs.Pt()>200 && Reconstructed_Top.Pt()>200) && (DeltaR_of_W_Higgs>=2.2 && DeltaR_of_W_Higgs<=3.5) && (TMath::Abs(First_Higgs_Jet.Phi()-Second_Higgs_Jet.Phi())<=1.2 && TMath::Abs(Top_Jet.Phi()-Reconstructed_W.Phi())<=1.2) && TMath::Abs(First_W_Jet.Phi()-Second_W_Jet.Phi())<=2.0 && Relative_THT>=0.65");
      WMassFromHiggsChi2[4] = (TH1F*)gDirectory->Get(WFH2.c_str());
      gPad->Close();
      //
      WFH1 = Form("W_From_Higgs_Chi2.M() >> WMFHC5%i(44,60.0,500)",9);
      WFH2 = Form("WMFHC5%i",9);
      AnalysisChain.Draw(WFH1.c_str(),"(Reconstructed_Higgs.Pt()>200 && Reconstructed_Top.Pt()>200) && (DeltaR_of_W_Higgs>=2.2 && DeltaR_of_W_Higgs<=3.5) && (TMath::Abs(First_Higgs_Jet.Phi()-Second_Higgs_Jet.Phi())<=1.2 && TMath::Abs(Top_Jet.Phi()-Reconstructed_W.Phi())<=1.2) && TMath::Abs(First_W_Jet.Phi()-Second_W_Jet.Phi())<=2.0 && Relative_THT>=0.65 && (Relative_Mass>=0.3 && Relative_Mass<=0.5)");
      WMassFromHiggsChi2[5] = (TH1F*)gDirectory->Get(WFH2.c_str());
      gPad->Close();
      /////////////////////////////////////////////////////////////
      string TPNWH1 = Form("Reconstructed_Tprime.M() >> TPNWFH%i(60,400,1600)",9);
      string TPNWH2 = Form("TPNWFH%i",9);
      AnalysisChain.Draw(TPNWH1.c_str(),"(Reconstructed_Higgs.Pt()>200 && Reconstructed_Top.Pt()>200) && (DeltaR_of_W_Higgs>=2.2 && DeltaR_of_W_Higgs<=3.5) && (TMath::Abs(First_Higgs_Jet.Phi()-Second_Higgs_Jet.Phi())<=1.2 && TMath::Abs(Top_Jet.Phi()-Reconstructed_W.Phi())<=1.2) && TMath::Abs(First_W_Jet.Phi()-Second_W_Jet.Phi())<=2.0 && Relative_THT>=0.65 && (Relative_Mass>=0.3 && Relative_Mass<=0.5) && Number_of_Tops<2"); //(W_From_Higgs.M()<70 || W_From_Higgs.M()>110)");
      TprimeNotWFromHiggs = (TH1F*)gDirectory->Get(TPNWH2.c_str());
      gPad->Close();
      //
      string D1 = Form("ABs(First_W_Jet.Phi(),Second_W_Jet.Phi()):THT >> DPWJ_HT_%i(65,300,1600,60,0.5,3.5)",18);
      string D2 = Form("DPWJ_HT_%i",18);
      AnalysisChain.Draw(D1.c_str(),"(Top_From_Higgs_Chi2.M()<140 || Top_From_Higgs_Chi2.M()>230)"); //,"(Reconstructed_Higgs.Pt()<200 || Reconstructed_Top.Pt()<200) && (DeltaR_of_W_Higgs>=2.2 && DeltaR_of_W_Higgs<=3.5) && (TMath::Abs(First_Higgs_Jet.Phi()-Second_Higgs_Jet.Phi())<=1.2 && TMath::Abs(Top_Jet.Phi()-Reconstructed_W.Phi())<=1.2) && Relative_THT>=0.65 && (Relative_Mass>=0.3 && Relative_Mass<=0.5) && (Top_From_Higgs_Chi2.M()<140 || Top_From_Higgs_Chi2.M()>230) && (Reconstructed_Higgs.M()>=110 && Reconstructed_Higgs.M()<=140)");
      DPWJ_HT = (TH2F*)gDirectory->Get(D2.c_str());
      cout << "Correlation Factor DPWJ and HT: " << DPWJ_HT->GetCorrelationFactor() << " Integral->" << DPWJ_HT->Integral() << endl;
      gPad->Close();
      //
      D1 = Form("ABs(First_W_Jet.Phi(),Second_W_Jet.Phi()):THT >> DPWJ_HT_BE_%i(65,300,1600,60,0.5,3.5)",18);
      D2 = Form("DPWJ_HT_BE_%i",18);
      AnalysisChain.Draw(D1.c_str(),"(Top_From_Higgs_Chi2.M()>=140 && Top_From_Higgs_Chi2.M()<=230)"); //,"(Reconstructed_Higgs.Pt()<200 || Reconstructed_Top.Pt()<200) && (DeltaR_of_W_Higgs>=2.2 && DeltaR_of_W_Higgs<=3.5) && (TMath::Abs(First_Higgs_Jet.Phi()-Second_Higgs_Jet.Phi())<=1.2 && TMath::Abs(Top_Jet.Phi()-Reconstructed_W.Phi())<=1.2) && Relative_THT>=0.65 && (Relative_Mass>=0.3 && Relative_Mass<=0.5) && (Top_From_Higgs_Chi2.M()<140 || Top_From_Higgs_Chi2.M()>230) && (Reconstructed_Higgs.M()>=110 && Reconstructed_Higgs.M()<=140)");
      DPWJ_HT_BE = (TH2F*)gDirectory->Get(D2.c_str());
      cout << "Correlation Factor DPWJ and HT: " << DPWJ_HT->GetCorrelationFactor() << " Integral->" << DPWJ_HT->Integral() << endl;
      gPad->Close();
      //////////////////////////////////////////////////////////////////////////////
      //////////////////////////////////////////////////////////////////////////////
      //////////////////////////////////////////////////////////////////////////////
      A1 = Form("Reconstructed_Tprime3L.M() >> TprimeMass3L%i(60,400,1600)",9);
      A2 = Form("TprimeMass3L%i",9);
      AnalysisChain3L.Draw(A1.c_str());
      FiveJetsMass3L = (TH1F*)gDirectory->Get(A2.c_str());
      gPad->Close();
      /////////////////////
      ///////////////////// ADDTIONAL INFO ON JETS
      /////////////////////
      A1 = Form("JET1.Eta() >> jet1_eta%i(100,-5,5)",9);
      A2 = Form("jet1_eta%i",9);
      AnalysisChain.Draw(A1.c_str(),finalcut);
      JET1ETA = (TH1F*)gDirectory->Get(A2.c_str());
      gPad->Close();
      //
      A1 = Form("JET2.Eta() >> jet2_eta%i(100,-5,5)",9);
      A2 = Form("jet2_eta%i",9);
      AnalysisChain.Draw(A1.c_str(),finalcut);
      JET2ETA = (TH1F*)gDirectory->Get(A2.c_str());
      gPad->Close();
      //
      A1 = Form("JET3.Eta() >> jet3_eta%i(100,-5,5)",9);
      A2 = Form("jet3_eta%i",9);
      AnalysisChain.Draw(A1.c_str(),finalcut);
      JET3ETA = (TH1F*)gDirectory->Get(A2.c_str());
      gPad->Close();
      //
      A1 = Form("JET4.Eta() >> jet4_eta%i(100,-5,5)",9);
      A2 = Form("jet4_eta%i",9);
      AnalysisChain.Draw(A1.c_str(),finalcut);
      JET4ETA = (TH1F*)gDirectory->Get(A2.c_str());
      gPad->Close();
      //
      A1 = Form("JET5.Eta() >> jet5_eta%i(100,-5,5)",9);
      A2 = Form("jet5_eta%i",9);
      AnalysisChain.Draw(A1.c_str(),finalcut);
      JET5ETA = (TH1F*)gDirectory->Get(A2.c_str());
      gPad->Close();
      //
      A1 = Form("JET6.Eta() >> jet6_eta%i(100,-5,5)",9);
      A2 = Form("jet6_eta%i",9);
      AnalysisChain.Draw(A1.c_str(),finalcut);
      JET6ETA = (TH1F*)gDirectory->Get(A2.c_str());
      gPad->Close();
      //
      A1 = Form("JET1.Phi() >> jet1_phi%i(60,-3,3)",9);
      A2 = Form("jet1_phi%i",9);
      AnalysisChain.Draw(A1.c_str(),finalcut);
      JET1PHI = (TH1F*)gDirectory->Get(A2.c_str());
      gPad->Close();
      //
      A1 = Form("JET2.Phi() >> jet2_phi%i(60,-3,3)",9);
      A2 = Form("jet2_phi%i",9);
      AnalysisChain.Draw(A1.c_str(),finalcut);
      JET2PHI = (TH1F*)gDirectory->Get(A2.c_str());
      gPad->Close();
      //
      A1 = Form("JET3.Phi() >> jet3_phi%i(60,-3,3)",9);
      A2 = Form("jet3_phi%i",9);
      AnalysisChain.Draw(A1.c_str(),finalcut);
      JET3PHI = (TH1F*)gDirectory->Get(A2.c_str());
      gPad->Close();
      //
      A1 = Form("JET4.Phi() >> jet4_phi%i(60,-3,3)",9);
      A2 = Form("jet4_phi%i",9);
      AnalysisChain.Draw(A1.c_str(),finalcut);
      JET4PHI = (TH1F*)gDirectory->Get(A2.c_str());
      gPad->Close();
      //
      A1 = Form("JET5.Phi() >> jet5_phi%i(60,-3,3)",9);
      A2 = Form("jet5_phi%i",9);
      AnalysisChain.Draw(A1.c_str(),finalcut);
      JET5PHI = (TH1F*)gDirectory->Get(A2.c_str());
      gPad->Close();
      //
      A1 = Form("JET6.Phi() >> jet6_phi%i(60,-3,3)",9);
      A2 = Form("jet6_phi%i",9);
      AnalysisChain.Draw(A1.c_str(),finalcut);
      JET6PHI = (TH1F*)gDirectory->Get(A2.c_str());
      gPad->Close();
      //
      A1 = Form("Number_Jets >> Num_jets%i(20,0,20)",9);
      A2 = Form("Num_jets%i",9);
      AnalysisChain.Draw(A1.c_str(),finalcut);
      JETMULTI = (TH1F*)gDirectory->Get(A2.c_str());
      gPad->Close();
    }
    
  FiveJetsMass->Scale(Lumi*Xs/ProcessedEvents[9]);
  FiveJetsMass3L->Scale(Lumi*Xs/ProcessedEvents[9]);
  LeadingJetPT->Scale(Lumi*Xs/ProcessedEvents[9]);
  Leading2JetPT->Scale(Lumi*Xs/ProcessedEvents[9]);
  Leading3JetPT->Scale(Lumi*Xs/ProcessedEvents[9]);
  Leading4JetPT->Scale(Lumi*Xs/ProcessedEvents[9]);
  Leading5JetPT->Scale(Lumi*Xs/ProcessedEvents[9]);
  Leading6JetPT->Scale(Lumi*Xs/ProcessedEvents[9]);
  THT->Scale(Lumi*Xs/ProcessedEvents[9]);
  DRHjets->Scale(Lumi*Xs/ProcessedEvents[9]);
  DRWjets->Scale(Lumi*Xs/ProcessedEvents[9]);
  Hpt->Scale(Lumi*Xs/ProcessedEvents[9]);
  Tpt->Scale(Lumi*Xs/ProcessedEvents[9]);
  DRWH->Scale(Lumi*Xs/ProcessedEvents[9]);
  DPHjets->Scale(Lumi*Xs/ProcessedEvents[9]);
  DPWjets->Scale(Lumi*Xs/ProcessedEvents[9]);
  DPTjets->Scale(Lumi*Xs/ProcessedEvents[9]);
  HiggsMass->Scale(Lumi*Xs/ProcessedEvents[9]);
  RelHT->Scale(Lumi*Xs/ProcessedEvents[9]);
  DRTH->Scale(Lumi*Xs/ProcessedEvents[9]);
  PtNormalizedMass->Scale(Lumi*Xs/ProcessedEvents[9]);
  RelativeMass->Scale(Lumi*Xs/ProcessedEvents[9]);
  MotherPtNormalizedMass->Scale(Lumi*Xs/ProcessedEvents[9]);  
  NumberOfTops->Scale(Lumi*Xs/ProcessedEvents[9]);
  HiggsMassOverTopMass->Scale(Lumi*Xs/ProcessedEvents[9]);
  HiggsTopAsymmetry->Scale(Lumi*Xs/ProcessedEvents[9]);
  ThirdLooseBtag->Scale(Lumi*Xs/ProcessedEvents[9]);
  TopMass->Scale(Lumi*Xs/ProcessedEvents[9]);
  Chi2->Scale(Lumi*Xs/ProcessedEvents[9]);
  UQuarkContent->Scale(Lumi*Xs/ProcessedEvents[9]);
  DQuarkContent->Scale(Lumi*Xs/ProcessedEvents[9]);
  SQuarkContent->Scale(Lumi*Xs/ProcessedEvents[9]);
  CQuarkContent->Scale(Lumi*Xs/ProcessedEvents[9]);
  BQuarkContent->Scale(Lumi*Xs/ProcessedEvents[9]);
  CSVLB->Scale(Lumi*Xs/ProcessedEvents[9]);
  CSVMB->Scale(Lumi*Xs/ProcessedEvents[9]);
  CSVTB->Scale(Lumi*Xs/ProcessedEvents[9]);
  Vtcs->Scale(Lumi*Xs/ProcessedEvents[9]);	 
  for (int i=0; i<6; i++) HptTpt[i]->Scale(Lumi*Xs/ProcessedEvents[9]); 
  /*CLBT->Scale(Lumi*Xs/ProcessedEvents[9]); 
  CMBT->Scale(Lumi*Xs/ProcessedEvents[9]); 
  CTBT->Scale(Lumi*Xs/ProcessedEvents[9]);  
  FLLBT->Scale(Lumi*Xs/ProcessedEvents[9]); 
  FLMBT->Scale(Lumi*Xs/ProcessedEvents[9]); 
  FLTBT->Scale(Lumi*Xs/ProcessedEvents[9]); 
  FCLBT->Scale(Lumi*Xs/ProcessedEvents[9]); 
  FCMBT->Scale(Lumi*Xs/ProcessedEvents[9]); 
  FCTBT->Scale(Lumi*Xs/ProcessedEvents[9]); */
  HTCSVM->Scale(Lumi*Xs/ProcessedEvents[9]);
  for (int i=0; i<6; i++) HiggsMassReversedHptTpt[i]->Scale(Lumi*Xs/ProcessedEvents[9]);
  for (int i=0; i<6; i++) TprimeMassReversedHptTpt[i]->Scale(Lumi*Xs/ProcessedEvents[9]);    
  for (int i=0; i<6; i++) HiggsMassCuts[i]->Scale(Lumi*Xs/ProcessedEvents[9]);
  for (int i=0; i<6; i++) WMassFromHiggs[i]->Scale(Lumi*Xs/ProcessedEvents[9]);
  for (int i=0; i<6; i++) WMassFromHiggsChi2[i]->Scale(Lumi*Xs/ProcessedEvents[9]);
  TprimeNotWFromHiggs->Scale(Lumi*Xs/ProcessedEvents[9]);
  //
  FiveJetsMassBE->Scale(Lumi*Xs/ProcessedEvents[9]); FiveJetsMassBE->Scale((TopMassFromHiggs->Integral(TopMassFromHiggs->GetXaxis()->FindBin(0.0),TopMassFromHiggs->GetXaxis()->FindBin(140.0))+TopMassFromHiggs->Integral(TopMassFromHiggs->GetXaxis()->FindBin(230.0),TopMassFromHiggs->GetXaxis()->FindBin(10000.0)))/TopMassFromHiggs->Integral(TopMassFromHiggs->GetXaxis()->FindBin(140.0),TopMassFromHiggs->GetXaxis()->FindBin(230.0)));
  FiveJetsMassLC->Scale(Lumi*Xs/ProcessedEvents[9]);
  TopMassFromHiggs->Scale(Lumi*Xs/ProcessedEvents[9]);
  TopMassFromHiggsLC->Scale(Lumi*Xs/ProcessedEvents[9]);
  FiveJetsMassLCoverBE->Scale(Lumi*Xs/ProcessedEvents[9]); FiveJetsMassLCoverBE->Scale(FiveJetsMassBE->Integral()/FiveJetsMassLCoverBE->Integral()); FiveJetsMassLCoverBE->Divide(FiveJetsMassBE);
  //
  HMBE->Scale(Lumi*Xs/ProcessedEvents[9]);
  HMLC->Scale(Lumi*Xs/ProcessedEvents[9]);
  HMLCoverBE->Scale(Lumi*Xs/ProcessedEvents[9]); HMLCoverBE->Scale(HMBE->Integral()/HMLCoverBE->Integral()); HMLCoverBE->Divide(HMBE);
  FJLC_QCDBE->Scale(Lumi*Xs/ProcessedEvents[9]);
  FJLC_QCDLC->Scale(Lumi*Xs/ProcessedEvents[9]);
  FJLC_QCDLCoverBE->Scale(Lumi*Xs/ProcessedEvents[9]); FJLC_QCDLCoverBE->Scale(FJLC_QCDBE->Integral()/FJLC_QCDLCoverBE->Integral()); FJLC_QCDLCoverBE->Divide(FJLC_QCDBE);
  //
  FiveJetsMassBE_HHM->Scale(Lumi*Xs/ProcessedEvents[9]); //FiveJetsMassBE_HHM->Scale((TopMassFromHiggs_HHM->Integral(TopMassFromHiggs_HHM->GetXaxis()->FindBin(0.0),TopMassFromHiggs_HHM->GetXaxis()->FindBin(140.0))+TopMassFromHiggs_HHM->Integral(TopMassFromHiggs_HHM->GetXaxis()->FindBin(230.0),TopMassFromHiggs_HHM->GetXaxis()->FindBin(10000.0)))/TopMassFromHiggs_HHM->Integral(TopMassFromHiggs_HHM->GetXaxis()->FindBin(140.0),TopMassFromHiggs_HHM->GetXaxis()->FindBin(230.0)));
  FiveJetsMassLC_HHM->Scale(Lumi*Xs/ProcessedEvents[9]);
  TopMassFromHiggs_HHM->Scale(Lumi*Xs/ProcessedEvents[9]);
  TopMassFromHiggs_HHMoverFHM->Scale(Lumi*Xs/ProcessedEvents[9]); TopMassFromHiggs_HHMoverFHM->Scale(TopMassFromHiggs->Integral()/TopMassFromHiggs_HHMoverFHM->Integral()); TopMassFromHiggs_HHMoverFHM->Divide(TopMassFromHiggs);
  FiveJetsMassLCoverBE_HHM->Scale(Lumi*Xs/ProcessedEvents[9]); FiveJetsMassLCoverBE_HHM->Scale(FiveJetsMassBE_HHM->Integral()/FiveJetsMassLCoverBE_HHM->Integral()); FiveJetsMassLCoverBE_HHM->Divide(FiveJetsMassBE_HHM);
  FiveJetsMassBE_LHM->Scale(Lumi*Xs/ProcessedEvents[9]);
  FiveJetsMassLC_LHM->Scale(Lumi*Xs/ProcessedEvents[9]);
  FiveJetsMassLCoverBE_LHM->Scale(Lumi*Xs/ProcessedEvents[9]); FiveJetsMassLCoverBE_LHM->Scale(FiveJetsMassBE_LHM->Integral()/FiveJetsMassLCoverBE_LHM->Integral()); FiveJetsMassLCoverBE_LHM->Divide(FiveJetsMassBE_LHM);
  FiveJetsMassLC_LHMoverHHM->Scale(Lumi*Xs/ProcessedEvents[9]); FiveJetsMassLC_LHMoverHHM->Scale(FiveJetsMassLC_HHM->Integral()/FiveJetsMassLC_LHMoverHHM->Integral()); FiveJetsMassLC_LHMoverHHM->Divide(FiveJetsMassLC_HHM);
  //
  TM_LC->Scale(Lumi*Xs/ProcessedEvents[9]); TM_BE->Scale(Lumi*Xs/ProcessedEvents[9]); 
  TM_LCoverBE->Scale(Lumi*Xs/ProcessedEvents[9]); TM_LCoverBE->Scale(TM_BE->Integral()/TM_LCoverBE->Integral()); TM_LCoverBE->Divide(TM_BE);
  TM_LC_HMW->Scale(Lumi*Xs/ProcessedEvents[9]);
  //
  HM_W_BE->Scale(Lumi*Xs/ProcessedEvents[9]); HM_W_LC->Scale(Lumi*Xs/ProcessedEvents[9]);
  HM_W_LCoverBE->Scale(Lumi*Xs/ProcessedEvents[9]); HM_W_LCoverBE->Scale(HM_W_BE->Integral()/HM_W_LCoverBE->Integral()); HM_W_LCoverBE->Divide(HM_W_BE);
  FJM_W_BE->Scale(Lumi*Xs/ProcessedEvents[9]); FJM_W_LC->Scale(Lumi*Xs/ProcessedEvents[9]);
  FJM_W_LCoverBE->Scale(Lumi*Xs/ProcessedEvents[9]); FJM_W_LCoverBE->Scale(FJM_W_BE->Integral()/FJM_W_LCoverBE->Integral()); FJM_W_LCoverBE->Divide(FJM_W_BE);
  //
  FiveJetsMassLC_LRHT->Scale(Lumi*Xs/ProcessedEvents[9]);
  FiveJetsMassLC_HRHT->Scale(Lumi*Xs/ProcessedEvents[9]);
  FiveJetsMassLC_LRM->Scale(Lumi*Xs/ProcessedEvents[9]);
  FiveJetsMassLC_HRM->Scale(Lumi*Xs/ProcessedEvents[9]);
  FiveJetsMassLC_WRM->Scale(Lumi*Xs/ProcessedEvents[9]);
  //
  DPWJ_HT->Scale(Lumi*Xs/ProcessedEvents[9]);
  DPWJ_HT_BE->Scale(Lumi*Xs/ProcessedEvents[9]);
  //cout << "KKKKKKKKKKKKKKKKKKKKKGHGHGHGHGHGHGHG" << endl;
  JET1ETA->Scale(Lumi*Xs/ProcessedEvents[9]);
  JET2ETA->Scale(Lumi*Xs/ProcessedEvents[9]);
  JET3ETA->Scale(Lumi*Xs/ProcessedEvents[9]);
  JET4ETA->Scale(Lumi*Xs/ProcessedEvents[9]);
  JET5ETA->Scale(Lumi*Xs/ProcessedEvents[9]);
  JET6ETA->Scale(Lumi*Xs/ProcessedEvents[9]);
  JET1PHI->Scale(Lumi*Xs/ProcessedEvents[9]);
  JET2PHI->Scale(Lumi*Xs/ProcessedEvents[9]);
  JET3PHI->Scale(Lumi*Xs/ProcessedEvents[9]);
  JET4PHI->Scale(Lumi*Xs/ProcessedEvents[9]);
  JET5PHI->Scale(Lumi*Xs/ProcessedEvents[9]);
  JET6PHI->Scale(Lumi*Xs/ProcessedEvents[9]);
  JETMULTI->Scale(Lumi*Xs/ProcessedEvents[9]);
  if (SurvivalMarker) 
    {
      //Settings for TTbar
      FiveJetsMass->SetFillColor(kRed);
      FiveJetsMass->SetFillStyle(3345);
      FiveJetsMass3L->SetFillColor(kRed);
      FiveJetsMass3L->SetFillStyle(3345);
      LeadingJetPT->SetFillColor(kRed);
      LeadingJetPT->SetFillStyle(3345);
      Leading2JetPT->SetFillColor(kRed);
      Leading2JetPT->SetFillStyle(3345);
      Leading3JetPT->SetFillColor(kRed);
      Leading3JetPT->SetFillStyle(3345);
      Leading4JetPT->SetFillColor(kRed);
      Leading4JetPT->SetFillStyle(3345);
      Leading5JetPT->SetFillColor(kRed);
      Leading5JetPT->SetFillStyle(3345);
      Leading6JetPT->SetFillColor(kRed);
      Leading6JetPT->SetFillStyle(3345);
      THT->SetFillColor(kRed);
      THT->SetFillStyle(3345);
      DRHjets->SetFillColor(kRed);
      DRHjets->SetFillStyle(3345);
      DRWjets->SetFillColor(kRed);
      DRWjets->SetFillStyle(3345);
      Hpt->SetFillColor(kRed);
      Hpt->SetFillStyle(3345);
      Tpt->SetFillColor(kRed);
      Tpt->SetFillStyle(3345);
      DRWH->SetFillColor(kRed);
      DRWH->SetFillStyle(3345);
      DPHjets->SetFillColor(kRed);
      DPHjets->SetFillStyle(3345);
      DPWjets->SetFillColor(kRed);
      DPWjets->SetFillStyle(3345);
      DPTjets->SetFillColor(kRed);
      DPTjets->SetFillStyle(3345);
      HiggsMass->SetFillColor(kRed);
      HiggsMass->SetFillStyle(3345);
      RelHT->SetFillColor(kRed);
      RelHT->SetFillStyle(3345);
      DRTH->SetFillColor(kRed);
      DRTH->SetFillStyle(3345);
      PtNormalizedMass->SetFillColor(kRed);
      PtNormalizedMass->SetFillStyle(3345);
      RelativeMass->SetFillColor(kRed);
      RelativeMass->SetFillStyle(3345);
      MotherPtNormalizedMass->SetFillColor(kRed);
      MotherPtNormalizedMass->SetFillStyle(3345);
      NumberOfTops->SetFillColor(kRed);
      NumberOfTops->SetFillStyle(3345);
      HiggsMassOverTopMass->SetFillColor(kRed);
      HiggsMassOverTopMass->SetFillStyle(3345);
      HiggsTopAsymmetry->SetFillColor(kRed);
      HiggsTopAsymmetry->SetFillStyle(3345);
      ThirdLooseBtag->SetFillColor(kRed);
      ThirdLooseBtag->SetFillStyle(3345);
      TopMass->SetFillColor(kRed);
      TopMass->SetFillStyle(3345);
      Chi2->SetFillColor(kRed);
      Chi2->SetFillStyle(3345);
      UQuarkContent->SetFillColor(kRed);
      UQuarkContent->SetFillStyle(3345);
      DQuarkContent->SetFillColor(kRed);
      DQuarkContent->SetFillStyle(3345);
      SQuarkContent->SetFillColor(kRed);
      SQuarkContent->SetFillStyle(3345);
      CQuarkContent->SetFillColor(kRed);
      CQuarkContent->SetFillStyle(3345);
      BQuarkContent->SetFillColor(kRed);
      BQuarkContent->SetFillStyle(3345);
      CSVLB->SetFillColor(kRed);
      CSVLB->SetFillStyle(3345);
      CSVMB->SetFillColor(kRed);
      CSVMB->SetFillStyle(3345);
      CSVTB->SetFillColor(kRed);
      CSVTB->SetFillStyle(3345);
      Vtcs->SetFillColor(kRed);
      Vtcs->SetFillStyle(3345);
      /*CLBT->SetFillColor(kRed);
      CLBT->SetFillStyle(3345);
      CMBT->SetFillColor(kRed);
      CMBT->SetFillStyle(3345);
      CTBT->SetFillColor(kRed);
      CTBT->SetFillStyle(3345);
      FLLBT->SetFillColor(kRed);
      FLLBT->SetFillStyle(3345);
      FLMBT->SetFillColor(kRed);
      FLMBT->SetFillStyle(3345);
      FLTBT->SetFillColor(kRed);
      FLTBT->SetFillStyle(3345);
      FCLBT->SetFillColor(kRed);
      FCLBT->SetFillStyle(3345);
      FCMBT->SetFillColor(kRed);
      FCMBT->SetFillStyle(3345);
      FCTBT->SetFillColor(kRed);
      FCTBT->SetFillStyle(3345);*/
      for (int i=0; i<6; i++) {HiggsMassReversedHptTpt[i]->SetFillColor(kRed); HiggsMassReversedHptTpt[i]->SetFillStyle(3345);}
      for (int i=0; i<6; i++) {TprimeMassReversedHptTpt[i]->SetFillColor(kRed); TprimeMassReversedHptTpt[i]->SetFillStyle(3345);}
      for (int i=0; i<6; i++) {HiggsMassCuts[i]->SetFillColor(kRed); HiggsMassCuts[i]->SetFillStyle(3345);}
      for (int i=0; i<6; i++) {WMassFromHiggs[i]->SetFillColor(kRed); WMassFromHiggs[i]->SetFillStyle(3345);}
      for (int i=0; i<6; i++) {WMassFromHiggsChi2[i]->SetFillColor(kRed); WMassFromHiggsChi2[i]->SetFillStyle(3345);}
      TprimeNotWFromHiggs->SetFillColor(kRed); TprimeNotWFromHiggs->SetFillStyle(3345);
      FiveJetsMassLC->SetFillColor(kRed); FiveJetsMassLC->SetFillStyle(3345);
      FiveJetsMassBE->SetFillColor(kRed); FiveJetsMassBE->SetFillStyle(3345);
      TopMassFromHiggs->SetFillColor(kRed); TopMassFromHiggs->SetFillStyle(3345);
      TopMassFromHiggsLC->SetFillColor(kRed); TopMassFromHiggsLC->SetFillStyle(3345);
      HMBE->SetFillColor(kRed); HMBE->SetFillStyle(3345);
      HMLC->SetFillColor(kRed); HMLC->SetFillStyle(3345);
      FJLC_QCDBE->SetFillColor(kRed); FJLC_QCDBE->SetFillStyle(3345);
      FJLC_QCDLC->SetFillColor(kRed); FJLC_QCDLC->SetFillStyle(3345);
      //
      FiveJetsMassBE_HHM->SetFillColor(kRed); FiveJetsMassBE_HHM->SetFillStyle(3345);
      FiveJetsMassLC_HHM->SetFillColor(kRed); FiveJetsMassLC_HHM->SetFillStyle(3345);
      TopMassFromHiggs_HHM->SetFillColor(kRed); TopMassFromHiggs_HHM->SetFillStyle(3345);
      FiveJetsMassBE_LHM->SetFillColor(kRed); FiveJetsMassBE_LHM->SetFillStyle(3345); 
      FiveJetsMassLC_LHM->SetFillColor(kRed); FiveJetsMassLC_LHM->SetFillStyle(3345);
      TM_LC->SetFillColor(kRed); TM_LC->SetFillStyle(3345);
      TM_BE->SetFillColor(kRed); TM_BE->SetFillStyle(3345);
      TM_LC_HMW->SetFillColor(kRed); TM_LC_HMW->SetFillStyle(3345);
      //
      HM_W_BE->SetFillColor(kRed); HM_W_BE->SetFillStyle(3345);
      HM_W_LC->SetFillColor(kRed); HM_W_LC->SetFillStyle(3345);
      FJM_W_BE->SetFillColor(kRed); FJM_W_BE->SetFillStyle(3345);
      FJM_W_LC->SetFillColor(kRed); FJM_W_LC->SetFillStyle(3345);
      //
      FiveJetsMassLC_LRHT->SetFillColor(kRed); FiveJetsMassLC_LRHT->SetFillStyle(3345);
      FiveJetsMassLC_HRHT->SetFillColor(kRed); FiveJetsMassLC_HRHT->SetFillStyle(3345);
      FiveJetsMassLC_LRM->SetFillColor(kRed); FiveJetsMassLC_LRM->SetFillStyle(3345);
      FiveJetsMassLC_HRM->SetFillColor(kRed); FiveJetsMassLC_HRM->SetFillStyle(3345);
      FiveJetsMassLC_WRM->SetFillColor(kRed); FiveJetsMassLC_WRM->SetFillStyle(3345);
      //
      JET1ETA->SetFillColor(kRed); JET1ETA->SetFillStyle(3345);
      JET2ETA->SetFillColor(kRed); JET2ETA->SetFillStyle(3345);
      JET3ETA->SetFillColor(kRed); JET3ETA->SetFillStyle(3345);
      JET4ETA->SetFillColor(kRed); JET4ETA->SetFillStyle(3345);
      JET5ETA->SetFillColor(kRed); JET5ETA->SetFillStyle(3345);
      JET6ETA->SetFillColor(kRed); JET6ETA->SetFillStyle(3345);
      JET1PHI->SetFillColor(kRed); JET1PHI->SetFillStyle(3345);
      JET2PHI->SetFillColor(kRed); JET2PHI->SetFillStyle(3345);
      JET3PHI->SetFillColor(kRed); JET3PHI->SetFillStyle(3345);
      JET4PHI->SetFillColor(kRed); JET4PHI->SetFillStyle(3345);
      JET5PHI->SetFillColor(kRed); JET5PHI->SetFillStyle(3345);
      JET6PHI->SetFillColor(kRed); JET6PHI->SetFillStyle(3345);
      JETMULTI->SetFillColor(kRed); JETMULTI->SetFillStyle(3345);
      //   
      TFile f("TTJets.root", "RECREATE");
      FiveJetsMass->Write();
      FiveJetsMass3L->Write();
      LeadingJetPT->Write();
      Leading2JetPT->Write();
      Leading3JetPT->Write();
      Leading4JetPT->Write();
      Leading5JetPT->Write();
      Leading6JetPT->Write();
      THT->Write();
      DRHjets->Write();
      DRWjets->Write();
      Hpt->Write();
      Tpt->Write();
      DRWH->Write();
      DPHjets->Write();
      DPWjets->Write();
      DPTjets->Write();
      HiggsMass->Write();
      RelHT->Write();
      DRTH->Write();
      PtNormalizedMass->Write();
      RelativeMass->Write();
      MotherPtNormalizedMass->Write();
      NumberOfTops->Write();
      HiggsMassOverTopMass->Write();
      HiggsTopAsymmetry->Write();
      ThirdLooseBtag->Write();
      TopMass->Write();
      Chi2->Write();
      UQuarkContent->Write();
      DQuarkContent->Write();
      SQuarkContent->Write();
      CQuarkContent->Write();
      BQuarkContent->Write();
      CSVLB->Write();
      CSVMB->Write();
      CSVTB->Write();
      Vtcs->Write();
      /*CLBT->Write();
      CMBT->Write();
      CTBT->Write();
      FLLBT->Write();
      FLMBT->Write();
      FLTBT->Write();
      FCLBT->Write();
      FCMBT->Write();
      FCTBT->Write();*/
      for (int i=0; i<6; i++) HptTpt[i]->Write();
      HTCSVM->Write();
      for (int i=0; i<6; i++) HiggsMassReversedHptTpt[i]->Write();
      for (int i=0; i<6; i++) TprimeMassReversedHptTpt[i]->Write();
      for (int i=0; i<6; i++) HiggsMassCuts[i]->Write();
      for (int i=0; i<6; i++) WMassFromHiggs[i]->Write();
      for (int i=0; i<6; i++) 
	{
	  //Fitting
	  //TF1 *fit = new TF1("fitting",fit_function,60,100,5);
	  //TF1 *peak = new TF1("GB",GaussBreitW,60,100,3);
	  //peak->SetLineColor(4);
	  //Double_t par[5];
	  //fit->SetParameter(3,WMassFromHiggsChi2[i]->GetMean());
	  //fit->SetParameter(4,WMassFromHiggsChi2[i]->GetRMS());
	  //WMassFromHiggsChi2[i]->Fit("fitting");
	  //fit->GetParameters(par);
	  //peak->SetParameters(&par[2]);
	  //peak->Draw("same");
	  WMassFromHiggsChi2[i]->Write();
	}
      TprimeNotWFromHiggs->Write();
      FiveJetsMassLC->Write();
      FiveJetsMassBE->Write();
      TopMassFromHiggs->Write();
      TopMassFromHiggsLC->Write();
      FiveJetsMassLCoverBE->Write();
      HMBE->Write();
      HMLC->Write();
      HMLCoverBE->Write();
      FJLC_QCDBE->Write();
      FJLC_QCDLC->Write();
      FJLC_QCDLCoverBE->Write();
      //
      FiveJetsMassBE_HHM->Write(); FiveJetsMassLC_HHM->Write(); FiveJetsMassLCoverBE_HHM->Write(); TopMassFromHiggs_HHM->Write(); TopMassFromHiggs_HHMoverFHM->Write();
      TM_LC->Write(); TM_BE->Write(); TM_LCoverBE->Write();
      TM_LC_HMW->Write();
      //
      cout << "GHGHGHGHGHGHGHG" << endl;
      HM_W_BE->Write(); HM_W_LC->Write(); 
      HM_W_LCoverBE->Write();
      FJM_W_BE->Write(); FJM_W_LC->Write(); 
      FJM_W_LCoverBE->Write();
      FiveJetsMassBE_LHM->Write(); FiveJetsMassLC_LHM->Write();
      FiveJetsMassLCoverBE_LHM->Write();
      FiveJetsMassLC_LHMoverHHM->Write();
      //
      FiveJetsMassLC_LRHT->Write();
      FiveJetsMassLC_HRHT->Write();
      FiveJetsMassLC_LRM->Write();
      FiveJetsMassLC_HRM->Write();
      FiveJetsMassLC_WRM->Write();
      //
      DPWJ_HT->Write();
      DPWJ_HT_BE->Write();
      //
      JET1ETA->Write();
      JET2ETA->Write();
      JET3ETA->Write();
      JET4ETA->Write();
      JET5ETA->Write();
      JET6ETA->Write();
      JET1PHI->Write();
      JET2PHI->Write();
      JET3PHI->Write();
      JET4PHI->Write();
      JET5PHI->Write();
      JET6PHI->Write();
      JETMULTI->Write();
    }  

  double RegionA=0; double RegionB=0; double RegionC=0; double RegionD=0;
  for (int i=0; i<6; i++) 
    {
      RegionA=HptTpt[i]->Integral(HptTpt[i]->GetXaxis()->FindBin(200.0),HptTpt[i]->GetNbinsX(),HptTpt[i]->GetYaxis()->FindBin(200.0),HptTpt[i]->GetNbinsY());
      RegionB=HptTpt[i]->Integral(HptTpt[i]->GetXaxis()->FindBin(200.0),HptTpt[i]->GetNbinsX(),HptTpt[i]->GetYaxis()->FindBin(0.0),HptTpt[i]->GetYaxis()->FindBin(200.0));
      RegionC=HptTpt[i]->Integral(HptTpt[i]->GetXaxis()->FindBin(0.0),HptTpt[i]->GetXaxis()->FindBin(200.0),HptTpt[i]->GetYaxis()->FindBin(200.0),HptTpt[i]->GetNbinsY());
      RegionD=HptTpt[i]->Integral(HptTpt[i]->GetXaxis()->FindBin(0.0),HptTpt[i]->GetXaxis()->FindBin(200.0),HptTpt[i]->GetYaxis()->FindBin(0.0),HptTpt[i]->GetYaxis()->FindBin(200.0));

      cout << "ABCD Method Info for Hpt Top pt:" << endl;
      cout << "Number of events on region A (signal enriched)" <<  RegionA << endl;
      cout << "Number of events on region B " <<  RegionB << endl;
      cout << "Number of events on region C " <<  RegionC << endl;
      cout << "Number of events on region D " <<  RegionD << endl;
    }

  double Region2A=0; double Region2B=0; double Region2C=0; double Region2D=0;
  Region2A=HTCSVM->Integral(HTCSVM->GetXaxis()->FindBin(630.0),HTCSVM->GetNbinsX(),HTCSVM->GetYaxis()->FindBin(3.0),HTCSVM->GetNbinsY());
  Region2B=HTCSVM->Integral(HTCSVM->GetXaxis()->FindBin(630.0),HTCSVM->GetNbinsX(),HTCSVM->GetYaxis()->FindBin(0.0),HTCSVM->GetYaxis()->FindBin(3.0));
  Region2C=HTCSVM->Integral(HTCSVM->GetXaxis()->FindBin(0.0),HTCSVM->GetXaxis()->FindBin(630.0),HTCSVM->GetYaxis()->FindBin(3.0),HTCSVM->GetNbinsY());
  Region2D=HTCSVM->Integral(HTCSVM->GetXaxis()->FindBin(0.0),HTCSVM->GetXaxis()->FindBin(630.0),HTCSVM->GetYaxis()->FindBin(0.0),HTCSVM->GetYaxis()->FindBin(3.0));

  cout << "ABCD Method Info for Hpt Top pt:" << endl;
  cout << "Number of events on region A (signal enriched) " <<  Region2A << endl;
  cout << "Number of events on region B " <<  Region2B << endl;
  cout << "Number of events on region C " <<  Region2C << endl;
  cout << "Number of events on region D " <<  Region2D << endl;

  exit(0);
  
}
