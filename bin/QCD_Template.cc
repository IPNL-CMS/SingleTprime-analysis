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

const TString MainFolder = "file:/sps/cms/ruizalva/CMSSW_5_3_9_patch3/src/Extractors/PatExtractor/bin/WithTrigger/SUFFIX/";

const int NumberOfHistos=23;
const TString Histos[NumberOfHistos] = {"Cut_0", "Cut_1", "Cut_2", "Cut_3", "Cut_chi2", "Cut_4", "Cut_5", "Cut_6", "Cut_7", "Cut_8", "Cut_9", "Cut_10", "Cut_11", "Cut_12", "Cut_13", "Cut_14", "Cut_15", "Cut_16", "Cut_17", "Cut_18", "Cut_19", "Cut_20", "Cut_21"};

//Computing weights (Everything in pb)
//const float XS[NumberOfProcesses]={8.258, 33.72, 54.838, 2.82, 47, 10.7, 1.57, 25, 10.7, 234.0, 156293.3, 34138.15, 1759.549, 113.8791, 26.9921, 3.550036, 0.15}; 
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

//float Lumi=2000.;

float PUR_function(int TI) //function with input the Number of True Interactions
{

  if (TI>=PUBins) return 1;
  else return PU_weight[TI];

}

double ABs(double X, double Y)
{
  return TMath::Abs(X-Y);
}

void QCD()
{
  int ParcialTestMax=18; //Don't change!!!
  int ParcialTestMin=10; //Don't change!!!
  TH1F *FiveJetsMass[NumberOfProcesses];
  TH1F *FiveJetsMass3L[NumberOfProcesses];
  TH1F *FiveJetsMass3T[NumberOfProcesses];
  TH1F *FiveJetsMass3L_L[NumberOfProcesses];
  TH1F *FiveJetsMass3T_L[NumberOfProcesses];
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
  TH1F *CLBT[NumberOfProcesses];   
  TH1F *CMBT[NumberOfProcesses]; 
  TH1F *CTBT[NumberOfProcesses];  
  TH1F *FLLBT[NumberOfProcesses]; 
  TH1F *FLMBT[NumberOfProcesses]; 
  TH1F *FLTBT[NumberOfProcesses]; 
  TH1F *FCLBT[NumberOfProcesses]; 
  TH1F *FCMBT[NumberOfProcesses]; 
  TH1F *FCTBT[NumberOfProcesses];   
  TH2F *HptTpt[6][NumberOfProcesses]; 
  TH2F *HTCSVM[NumberOfProcesses];  
  TH2F *RelHTNTops[NumberOfProcesses];  
  TH1F *HiggsMassCuts[6][NumberOfProcesses];
  TH1F *HiggsMassReversedHptTpt[6][NumberOfProcesses];
  TH1F *TprimeMassReversedHptTpt[6][NumberOfProcesses];    
  TH1F *WMassFromHiggs[6][NumberOfProcesses];      
  TH1F *WMassFromHiggsChi2[6][NumberOfProcesses]; 
  
  TH1F *FiveJetsMassBE[NumberOfProcesses]; //5 jets mass inside the last cut
  TH1F *FiveJetsMassLC[NumberOfProcesses]; //5 jets mass outside last cut
  TH1F *TopMassFromHiggs[NumberOfProcesses];
  TH1F *FiveJetsMassLCoverBE[NumberOfProcesses]; //5 jets mass ratio plot in vs out the cut
  TH1F *HMBE[NumberOfProcesses];
  TH1F *HMLC[NumberOfProcesses];
  TH1F *HMLCoverBE[NumberOfProcesses];
  TH1F *FJLC_QCDBE[NumberOfProcesses];
  TH1F *FJLC_QCDLC[NumberOfProcesses];
  TH1F *FJLC_QCDLCoverBE[NumberOfProcesses];

  TH1F *TM_LC_HMW[NumberOfProcesses]; //Top mass after ttbar estimation cut inside the Higgs mass window

  TH2F *WM_TM[NumberOfProcesses]; //First W mass with first top mass
  TH2F *HT_TM[NumberOfProcesses]; //HT with first top mass
  TH2F *DPWJ_TM[NumberOfProcesses]; //DeltaR W jets with first top mass

  //RM or RHT sideband estimations
  TH1F *FiveJetsMassLC_LRHT[NumberOfProcesses]; TH1F *FiveJetsMassLC_HRHT[NumberOfProcesses];
  TH1F *FiveJetsMassLC_LRM[NumberOfProcesses]; TH1F *FiveJetsMassLC_HRM[NumberOfProcesses]; TH1F *FiveJetsMassLC_WRM[NumberOfProcesses];

  //HT estimation
  TH2F *DPWJ_HT[NumberOfProcesses];

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
      TChain AnalysisChain3T("stp3T");
      TH1F *ALLCuts[NumberOfHistos];
      cout << "XS for " << k <<  " sample is: " << XS[k] << " with " << ProcessedEvents[k] << " processed events. Implying a weight of " << 0.5*Lumi*XS[k]/ProcessedEvents[k] << " for " << Lumi << "pb-1" << endl;
      if (k==10)
	{
	  CutsChain.Add(MainFolder + "QCD_PT_120_170-v3_Full_analyzed.root");
	  AnalysisChain.Add(MainFolder + "QCD_PT_120_170-v3_Full_analyzed.root");
	  AnalysisChain3L.Add(MainFolder + "QCD_PT_120_170-v3_Full_analyzed.root");
	  AnalysisChain3T.Add(MainFolder + "QCD_PT_120_170-v3_Full_analyzed.root");
	}
      else if (k==11)
	{
	  CutsChain.Add(MainFolder + "QCD_PT_170_300_Full_analyzed.root");
	  AnalysisChain.Add(MainFolder + "QCD_PT_170_300_Full_analyzed.root");
	  AnalysisChain3L.Add(MainFolder + "QCD_PT_170_300_Full_analyzed.root");
	  AnalysisChain3T.Add(MainFolder + "QCD_PT_170_300_Full_analyzed.root");
	}
      else if (k==12)
	{
	  CutsChain.Add(MainFolder + "QCD_PT_300_470_Full_analyzed.root");
	  AnalysisChain.Add(MainFolder + "QCD_PT_300_470_Full_analyzed.root");
	  AnalysisChain3L.Add(MainFolder + "QCD_PT_300_470_Full_analyzed.root");
	  AnalysisChain3T.Add(MainFolder + "QCD_PT_300_470_Full_analyzed.root");
	}
      else if (k==13)
	{
	  CutsChain.Add(MainFolder + "QCD_PT_470_600_Full_analyzed.root");
	  AnalysisChain.Add(MainFolder + "QCD_PT_470_600_Full_analyzed.root");
	  AnalysisChain3L.Add(MainFolder + "QCD_PT_470_600_Full_analyzed.root");
	  AnalysisChain3T.Add(MainFolder + "QCD_PT_470_600_Full_analyzed.root");
	}
      else if (k==14)
	{
	  cout << "1 Marker" << endl; 
	  CutsChain.Add(MainFolder + "QCD_PT_600_800_Full_analyzed.root");
	  AnalysisChain.Add(MainFolder + "QCD_PT_600_800_Full_analyzed.root");
	  AnalysisChain3L.Add(MainFolder + "QCD_PT_600_800_Full_analyzed.root");
	  AnalysisChain3T.Add(MainFolder + "QCD_PT_600_800_Full_analyzed.root");
	  cout << "2 Marker" << endl;
	}
      else if (k==15) //continue;
	{
	  CutsChain.Add(MainFolder + "QCD_PT_800_1000_Full_analyzed.root");
	  AnalysisChain.Add(MainFolder + "QCD_PT_800_1000_Full_analyzed.root");
	  AnalysisChain3L.Add(MainFolder + "QCD_PT_800_1000_Full_analyzed.root");
	  AnalysisChain3T.Add(MainFolder + "QCD_PT_800_1000_Full_analyzed.root");
	  }
      else if (k==16)
        {
          CutsChain.Add(MainFolder + "QCD_HT_500_1000_Full_analyzed.root");
          AnalysisChain.Add(MainFolder + "QCD_HT_500_1000_Full_analyzed.root");
          AnalysisChain3L.Add(MainFolder + "QCD_HT_500_1000_Full_analyzed.root");
          AnalysisChain3T.Add(MainFolder + "QCD_HT_500_1000_Full_analyzed.root");
        }
      else if (k==17)
        {
          CutsChain.Add(MainFolder + "QCD_HT_1000_Inf_Full_analyzed.root");
          AnalysisChain.Add(MainFolder + "QCD_HT_1000_Inf_Full_analyzed.root");
          AnalysisChain3L.Add(MainFolder + "QCD_HT_1000_Inf_Full_analyzed.root");
          AnalysisChain3T.Add(MainFolder + "QCD_HT_1000_Inf_Full_analyzed.root");
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
	  AnalysisChain.Draw(A1.c_str(),finalcut);
	  FiveJetsMass[k] = (TH1F*)gDirectory->Get(A2.c_str());
	  gPad->Close();
	  SurvivalMarker[k]=true;
	  if (MarkerCounter==0) FirstGoodMarker=k;
	  MarkerCounter++;
	  LastGoodMarker=k;
	  //
	  A1 = Form("Reconstructed_Tprime.M() >> TprimeMassBkgE%i(60,400,1600)",k);
	  A2 = Form("TprimeMassBkgE%i",k);
	  AnalysisChain.Draw(A1.c_str(),"Top_From_Higgs_Chi2.M()>=140 && Top_From_Higgs_Chi2.M()<=230");
	  FiveJetsMassBE[k] = (TH1F*)gDirectory->Get(A2.c_str());
	  gPad->Close();
	  //
	  A1 = Form("Reconstructed_Tprime.M() >> TprimeMassLC%i(60,400,1600)",k);
	  A2 = Form("TprimeMassLC%i",k);
	  AnalysisChain.Draw(A1.c_str(),"Top_From_Higgs_Chi2.M()<140 || Top_From_Higgs_Chi2.M()>230");
	  FiveJetsMassLC[k] = (TH1F*)gDirectory->Get(A2.c_str());
	  gPad->Close();
	  //
	  A1 = Form("Reconstructed_Tprime.M() >> TprimeMassLCoverBE%i(60,400,1600)",k);
	  A2 = Form("TprimeMassLCoverBE%i",k);
	  AnalysisChain.Draw(A1.c_str(),"Top_From_Higgs_Chi2.M()<140 || Top_From_Higgs_Chi2.M()>230");
	  FiveJetsMassLCoverBE[k] = (TH1F*)gDirectory->Get(A2.c_str());
	  gPad->Close();
	  //
	  string HMBEstim1 = Form("Reconstructed_Higgs.M() >> HiggsMassBE%i(36,60,180)",k);
	  string HMBEstim2 = Form("HiggsMassBE%i",k);
	  AnalysisChain.Draw(HMBEstim1.c_str(),"Top_From_Higgs_Chi2.M()>=140 && Top_From_Higgs_Chi2.M()<=230");
	  HMBE[k] = (TH1F*)gDirectory->Get(HMBEstim2.c_str());
	  gPad->Close();
	  //
	  HMBEstim1 = Form("Reconstructed_Higgs.M() >> HiggsMassLC%i(36,60,180)",k);
	  HMBEstim2 = Form("HiggsMassLC%i",k);
	  AnalysisChain.Draw(HMBEstim1.c_str(),"Top_From_Higgs_Chi2.M()<140 || Top_From_Higgs_Chi2.M()>230");
	  HMLC[k] = (TH1F*)gDirectory->Get(HMBEstim2.c_str());
	  gPad->Close();
	  //
	  HMBEstim1 = Form("Reconstructed_Higgs.M() >> HiggsMassLCoverBE%i(36,60,180)",k);
	  HMBEstim2 = Form("HiggsMassLCoverBE%i",k);
	  AnalysisChain.Draw(HMBEstim1.c_str(),"Top_From_Higgs_Chi2.M()<140 || Top_From_Higgs_Chi2.M()>230");
	  HMLCoverBE[k] = (TH1F*)gDirectory->Get(HMBEstim2.c_str());
	  gPad->Close();
	  //
	  HMBEstim1 = Form("Reconstructed_Tprime.M() >> FJLCQCDBE%i(60,400,1600)",k);
	  HMBEstim2 = Form("FJLCQCDBE%i",k);
	  AnalysisChain.Draw(HMBEstim1.c_str(),"(Top_From_Higgs_Chi2.M()<140 || Top_From_Higgs_Chi2.M()>230) && (Reconstructed_Higgs.M()<110 || Reconstructed_Higgs.M()>140)");
	  FJLC_QCDBE[k] = (TH1F*)gDirectory->Get(HMBEstim2.c_str());
	  gPad->Close();
	  //
	  HMBEstim1 = Form("Reconstructed_Tprime.M() >> FJLCQCDLC%i(60,400,1600)",k);
	  HMBEstim2 = Form("FJLCQCDLC%i",k);
	  AnalysisChain.Draw(HMBEstim1.c_str(),"(Top_From_Higgs_Chi2.M()<140 || Top_From_Higgs_Chi2.M()>230) && (Reconstructed_Higgs.M()>=110 && Reconstructed_Higgs.M()<=140)");
	  FJLC_QCDLC[k] = (TH1F*)gDirectory->Get(HMBEstim2.c_str());
	  gPad->Close();
	  //
	  HMBEstim1 = Form("Reconstructed_Tprime.M() >> FJLCQCDLCoverBE%i(60,400,1600)",k);
	  HMBEstim2 = Form("FJLCQCDLCoverBE%i",k);
	  AnalysisChain.Draw(HMBEstim1.c_str(),"(Top_From_Higgs_Chi2.M()<140 || Top_From_Higgs_Chi2.M()>230) && (Reconstructed_Higgs.M()>=110 && Reconstructed_Higgs.M()<=140)");
	  FJLC_QCDLCoverBE[k] = (TH1F*)gDirectory->Get(HMBEstim2.c_str());
	  gPad->Close();
	  //////////
	  //////////
	  A1 = Form("Top_From_Higgs_Chi2.M() >> TopFromHiggsChi2Mass%i(95,50,1000)",k);
	  A2 = Form("TopFromHiggsChi2Mass%i",k);
	  AnalysisChain.Draw(A1.c_str(),"(PUR_function(Number_True_Interactions))*weight");
	  TopMassFromHiggs[k] = (TH1F*)gDirectory->Get(A2.c_str());
	  gPad->Close();
	  //////////
	  //////////
	  string TM1_LC = Form("Reconstructed_Top.M() >> TMass_LC_HMW%i(60,100,700)",k);
	  string TM2_LC = Form("TMass_LC_HMW%i",k);
	  AnalysisChain.Draw(TM1_LC.c_str(),"(Top_From_Higgs_Chi2.M()<140 || Top_From_Higgs_Chi2.M()>230) && (Reconstructed_Higgs.M()>=110 && Reconstructed_Higgs.M()<=140)");
	  TM_LC_HMW[k] = (TH1F*)gDirectory->Get(TM2_LC.c_str());
	  gPad->Close();
	  //////////
	  //////////	  
	  string ABCD1 = Form("Reconstructed_W.M():Reconstructed_Top.M() >> WM_TM_%i(62,80,700,48,20,500)",k);
	  string ABCD2 = Form("WM_TM_%i",k);
	  AnalysisChain.Draw(ABCD1.c_str(),"(Top_From_Higgs_Chi2.M()<140 || Top_From_Higgs_Chi2.M()>230) && (Reconstructed_Higgs.M()>=110 && Reconstructed_Higgs.M()<=140)");
	  WM_TM[k] = (TH2F*)gDirectory->Get(ABCD2.c_str());
	  cout << "Correlation Factor WM and TM: " << WM_TM[k]->GetCorrelationFactor() << endl;
	  gPad->Close();
	  //
	  string D1 = Form("Reconstructed_Top.M():THT >> HT_TM_%i(65,300,1600,62,80,700)",k);
	  string D2 = Form("HT_TM_%i",k);
	  AnalysisChain.Draw(D1.c_str(),"(Top_From_Higgs_Chi2.M()<140 || Top_From_Higgs_Chi2.M()>230) && (Reconstructed_Higgs.M()>=110 && Reconstructed_Higgs.M()<=140)");
	  HT_TM[k] = (TH2F*)gDirectory->Get(D2.c_str());
	  cout << "Correlation Factor HT and TM: " << HT_TM[k]->GetCorrelationFactor() << " Integral->" << HT_TM[k]->Integral() << endl;
	  gPad->Close();
	  //
	  D1 = Form("ABs(First_W_Jet.Phi(),Second_W_Jet.Phi()):Reconstructed_Top.M() >> DPWJ_TM_%i(62,80,700,60,0.5,3.5)",k);
	  D2 = Form("DPWJ_TM_%i",k);
	  AnalysisChain.Draw(D1.c_str()); //,"(Reconstructed_Higgs.Pt()<200 || Reconstructed_Top.Pt()<200) && (DeltaR_of_W_Higgs>=2.2 && DeltaR_of_W_Higgs<=3.5) && (TMath::Abs(First_Higgs_Jet.Phi()-Second_Higgs_Jet.Phi())<=1.2 && TMath::Abs(Top_Jet.Phi()-Reconstructed_W.Phi())<=1.2) && Relative_THT>=0.65 && (Relative_Mass>=0.3 && Relative_Mass<=0.5) && (Top_From_Higgs_Chi2.M()<140 || Top_From_Higgs_Chi2.M()>230) && (Reconstructed_Higgs.M()>=110 && Reconstructed_Higgs.M()<=140)");
	  DPWJ_TM[k] = (TH2F*)gDirectory->Get(D2.c_str());
	  cout << "Correlation Factor DPWJ and TM: " << DPWJ_TM[k]->GetCorrelationFactor() << " Integral->" << DPWJ_TM[k]->Integral() << endl;
	  gPad->Close();
	  //
	  D1 = Form("ABs(First_W_Jet.Phi(),Second_W_Jet.Phi()):THT >> DPWJ_HT_%i(65,300,1600,60,0.5,3.5)",k);
	  D2 = Form("DPWJ_HT_%i",k);
	  AnalysisChain.Draw(D1.c_str()); //,"(Reconstructed_Higgs.Pt()<200 || Reconstructed_Top.Pt()<200) && (DeltaR_of_W_Higgs>=2.2 && DeltaR_of_W_Higgs<=3.5) && (TMath::Abs(First_Higgs_Jet.Phi()-Second_Higgs_Jet.Phi())<=1.2 && TMath::Abs(Top_Jet.Phi()-Reconstructed_W.Phi())<=1.2) && Relative_THT>=0.65 && (Relative_Mass>=0.3 && Relative_Mass<=0.5) && (Top_From_Higgs_Chi2.M()<140 || Top_From_Higgs_Chi2.M()>230) && (Reconstructed_Higgs.M()>=110 && Reconstructed_Higgs.M()<=140)");
	  DPWJ_HT[k] = (TH2F*)gDirectory->Get(D2.c_str());
	  cout << "Correlation Factor DPWJ and HT: " << DPWJ_HT[k]->GetCorrelationFactor() << " Integral->" << DPWJ_HT[k]->Integral() << endl;
	  gPad->Close();
	  //////////DPWJ_HT
	  //////////
	  A1 = Form("Reconstructed_Tprime.M() >> TprimeMassLC_LRHT%i(60,400,1600)",k);
	  A2 = Form("TprimeMassLC_LRHT%i",k);
	  AnalysisChain.Draw(A1.c_str(),"(Reconstructed_Higgs.Pt()>200 && Reconstructed_Top.Pt()>200) && (DeltaR_of_W_Higgs>=2.2 && DeltaR_of_W_Higgs<=3.5) && (TMath::Abs(First_Higgs_Jet.Phi()-Second_Higgs_Jet.Phi())<=1.2 && TMath::Abs(Top_Jet.Phi()-Reconstructed_W.Phi())<=1.2) && TMath::Abs(First_W_Jet.Phi()-Second_W_Jet.Phi())<=2.0 && (Relative_Mass>=0.3 && Relative_Mass<=0.5) && (Reconstructed_Higgs.M()>=110 && Reconstructed_Higgs.M()<=140) && (Top_From_Higgs_Chi2.M()<140 || Top_From_Higgs_Chi2.M()>230) && Relative_THT<0.65");
	  FiveJetsMassLC_LRHT[k] = (TH1F*)gDirectory->Get(A2.c_str());
	  gPad->Close();
	  //
	  A1 = Form("Reconstructed_Tprime.M() >> TprimeMassLC_HRHT%i(60,400,1600)",k);
	  A2 = Form("TprimeMassLC_HRHT%i",k);
	  AnalysisChain.Draw(A1.c_str(),"(Reconstructed_Higgs.Pt()>200 && Reconstructed_Top.Pt()>200) && (DeltaR_of_W_Higgs>=2.2 && DeltaR_of_W_Higgs<=3.5) && (TMath::Abs(First_Higgs_Jet.Phi()-Second_Higgs_Jet.Phi())<=1.2 && TMath::Abs(Top_Jet.Phi()-Reconstructed_W.Phi())<=1.2) && TMath::Abs(First_W_Jet.Phi()-Second_W_Jet.Phi())<=2.0 && (Relative_Mass>=0.3 && Relative_Mass<=0.5) && (Reconstructed_Higgs.M()>=110 && Reconstructed_Higgs.M()<=140) && (Top_From_Higgs_Chi2.M()<140 || Top_From_Higgs_Chi2.M()>230) && Relative_THT>=0.65");
	  FiveJetsMassLC_HRHT[k] = (TH1F*)gDirectory->Get(A2.c_str());
	  gPad->Close();
	  //
	  A1 = Form("Reconstructed_Tprime.M() >> TprimeMassLC_LRM%i(60,400,1600)",k);
	  A2 = Form("TprimeMassLC_LRM%i",k);
	  AnalysisChain.Draw(A1.c_str(),"(Reconstructed_Higgs.Pt()>200 && Reconstructed_Top.Pt()>200) && (DeltaR_of_W_Higgs>=2.2 && DeltaR_of_W_Higgs<=3.5) && (TMath::Abs(First_Higgs_Jet.Phi()-Second_Higgs_Jet.Phi())<=1.2 && TMath::Abs(Top_Jet.Phi()-Reconstructed_W.Phi())<=1.2) && TMath::Abs(First_W_Jet.Phi()-Second_W_Jet.Phi())<=2.0 && Relative_THT>=0.65 && (Reconstructed_Higgs.M()>=110 && Reconstructed_Higgs.M()<=140) && (Top_From_Higgs_Chi2.M()<140 || Top_From_Higgs_Chi2.M()>230) && Relative_Mass<0.3");
	  FiveJetsMassLC_LRM[k] = (TH1F*)gDirectory->Get(A2.c_str());
	  gPad->Close();
	  //
	  A1 = Form("Reconstructed_Tprime.M() >> TprimeMassLC_HRM%i(60,400,1600)",k);
	  A2 = Form("TprimeMassLC_HRM%i",k);
	  AnalysisChain.Draw(A1.c_str(),"(Reconstructed_Higgs.Pt()>200 && Reconstructed_Top.Pt()>200) && (DeltaR_of_W_Higgs>=2.2 && DeltaR_of_W_Higgs<=3.5) && (TMath::Abs(First_Higgs_Jet.Phi()-Second_Higgs_Jet.Phi())<=1.2 && TMath::Abs(Top_Jet.Phi()-Reconstructed_W.Phi())<=1.2) && TMath::Abs(First_W_Jet.Phi()-Second_W_Jet.Phi())<=2.0 && Relative_THT>=0.65 && (Reconstructed_Higgs.M()>=110 && Reconstructed_Higgs.M()<=140) && (Top_From_Higgs_Chi2.M()<140 || Top_From_Higgs_Chi2.M()>230) && Relative_Mass>0.5");
	  FiveJetsMassLC_HRM[k] = (TH1F*)gDirectory->Get(A2.c_str());
	  gPad->Close();
	  //
	  A1 = Form("Reconstructed_Tprime.M() >> TprimeMassLC_WRM%i(60,400,1600)",k);
	  A2 = Form("TprimeMassLC_WRM%i",k);
	  AnalysisChain.Draw(A1.c_str(),"(Reconstructed_Higgs.Pt()>200 && Reconstructed_Top.Pt()>200) && (DeltaR_of_W_Higgs>=2.2 && DeltaR_of_W_Higgs<=3.5) && (TMath::Abs(First_Higgs_Jet.Phi()-Second_Higgs_Jet.Phi())<=1.2 && TMath::Abs(Top_Jet.Phi()-Reconstructed_W.Phi())<=1.2) && TMath::Abs(First_W_Jet.Phi()-Second_W_Jet.Phi())<=2.0 && Relative_THT>=0.65 && (Reconstructed_Higgs.M()>=110 && Reconstructed_Higgs.M()<=140) && (Top_From_Higgs_Chi2.M()<140 || Top_From_Higgs_Chi2.M()>230) && (Relative_Mass>=0.3 && Relative_Mass<=0.5)");
	  FiveJetsMassLC_WRM[k] = (TH1F*)gDirectory->Get(A2.c_str());
	  gPad->Close();
	  //////////
	  //////////
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
	  string THT1 = Form("THT >> THT%i(65,300,1600)",k); //*PUR_function(Number_True_Interactions)
	  string THT2 = Form("THT%i",k);
	  AnalysisChain.Draw(THT1.c_str(),"(PUR_function(Number_True_Interactions))*weight"); //PU_weight*(THT>630)"); 
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
	  string TLBT1 = Form("Number_of_Loose_and_non_med_b_tags >> TLBTag%i(10,0,10.)",k);
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
	  AnalysisChain.Draw(VT1.c_str(),"(PUR_function(Number_True_Interactions))*weight"); //(1./PUR_function(Number_True_Interactions))*(THT>630)"); //,"PU_weight");
	  Vtcs[k] = (TH1F*)gDirectory->Get(VT2.c_str());
	  gPad->Close();
	  //
	  /*string CLB1 = Form("CorrectLB_tag >> CLooseB%i(6,0,6)",k);
	  string CLB2 = Form("CLooseB%i",k);
	  AnalysisChain.Draw(CLB1.c_str(),"(PUR_function(Number_True_Interactions))*weight");
	  CLBT[k] = (TH1F*)gDirectory->Get(CLB2.c_str());
	  gPad->Close();
	  //
	  string CMB1 = Form("CorrectMB_tag >> CMedB%i(6,0,6)",k);
	  string CMB2 = Form("CMedB%i",k);
	  AnalysisChain.Draw(CMB1.c_str(),"(PUR_function(Number_True_Interactions))*weight");
	  CMBT[k] = (TH1F*)gDirectory->Get(CMB2.c_str());
	  gPad->Close();
	  //
	  string CTB1 = Form("CorrectTB_tag >> CTightB%i(6,0,6)",k);
	  string CTB2 = Form("CTightB%i",k);
	  AnalysisChain.Draw(CTB1.c_str(),"(PUR_function(Number_True_Interactions))*weight");
	  CTBT[k] = (TH1F*)gDirectory->Get(CTB2.c_str());
	  gPad->Close();
	  //
	  string FLLB1 = Form("FakeLB_tag_Light >> FLLooseB%i(6,0,6)",k);
	  string FLLB2 = Form("FLLooseB%i",k);
	  AnalysisChain.Draw(FLLB1.c_str(),"(PUR_function(Number_True_Interactions))*weight");
	  FLLBT[k] = (TH1F*)gDirectory->Get(FLLB2.c_str());
	  gPad->Close();
	  //
	  string FLMB1 = Form("FakeMB_tag_Light >> FLMedB%i(6,0,6)",k);
	  string FLMB2 = Form("FLMedB%i",k);
	  AnalysisChain.Draw(FLMB1.c_str(),"(PUR_function(Number_True_Interactions))*weight");
	  FLMBT[k] = (TH1F*)gDirectory->Get(FLMB2.c_str());
	  gPad->Close();
	  //
	  string FLTB1 = Form("FakeTB_tag_Light >> FLTightB%i(6,0,6)",k);
	  string FLTB2 = Form("FLTightB%i",k);
	  AnalysisChain.Draw(FLTB1.c_str(),"(PUR_function(Number_True_Interactions))*weight");
	  FLTBT[k] = (TH1F*)gDirectory->Get(FLTB2.c_str());
	  gPad->Close();
	  //
	  string FCLB1 = Form("FakeLB_tag_C >> FCLooseB%i(6,0,6)",k);
	  string FCLB2 = Form("FCLooseB%i",k);
	  AnalysisChain.Draw(FCLB1.c_str(),"(PUR_function(Number_True_Interactions))*weight");
	  FCLBT[k] = (TH1F*)gDirectory->Get(FCLB2.c_str());
	  gPad->Close();
	  //
	  string FCMB1 = Form("FakeMB_tag_C >> FCMedB%i(6,0,6)",k);
	  string FCMB2 = Form("FCMedB%i",k);
	  AnalysisChain.Draw(FCMB1.c_str(),"(PUR_function(Number_True_Interactions))*weight");
	  FCMBT[k] = (TH1F*)gDirectory->Get(FCMB2.c_str());
	  gPad->Close();
	  //
	  string FCTB1 = Form("FakeTB_tag_C >> FCTightB%i(6,0,6)",k);
	  string FCTB2 = Form("FCTightB%i",k);
	  AnalysisChain.Draw(FCTB1.c_str(),"(PUR_function(Number_True_Interactions))*weight");
	  FCTBT[k] = (TH1F*)gDirectory->Get(FCTB2.c_str());
	  gPad->Close();*/
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
	  //
	  string RelHTNTops1 = Form("Number_of_Tops : Relative_THT >> RHTNT%i(30,0,1,8,0,8)",k);
	  string RelHTNTops2 = Form("RHTNT%i",k);
	  AnalysisChain.Draw(RelHTNTops1.c_str(),"(PUR_function(Number_True_Interactions))*weight");
	  RelHTNTops[k] = (TH2F*)gDirectory->Get(RelHTNTops2.c_str());
	  cout << "Correlation Factor RelHTNtops: " << RelHTNTops[k]->GetCorrelationFactor() << endl;
	  gPad->Close();
	  //////////////////////////////////////
	  /////////ABCD BKG Estimation//////////
	  //////////////////////////////////////
	  HM1 = Form("Reconstructed_Higgs.M() >> HMBE0%i(36,60,180)",k);
	  HM2 = Form("HMBE0%i",k);
	  AnalysisChain.Draw(HM1.c_str(),"Reconstructed_Higgs.Pt()<200 || Reconstructed_Top.Pt()<200");
	  HiggsMassReversedHptTpt[0][k] = (TH1F*)gDirectory->Get(HM2.c_str());
	  gPad->Close();
	  //	  
	  A1 = Form("Reconstructed_Tprime.M() >> TprimeMassBE0%i(60,400,1600)",k);
	  A2 = Form("TprimeMassBE0%i",k);
	  AnalysisChain.Draw(A1.c_str(),"Reconstructed_Higgs.Pt()<200 || Reconstructed_Top.Pt()<200");
	  TprimeMassReversedHptTpt[0][k] = (TH1F*)gDirectory->Get(A2.c_str());
	  gPad->Close();
	  //	  
	  HM1 = Form("Reconstructed_Higgs.M() >> HMBE1%i(36,60,180)",k);
	  HM2 = Form("HMBE1%i",k);
	  AnalysisChain.Draw(HM1.c_str(),"(Reconstructed_Higgs.Pt()<200 || Reconstructed_Top.Pt()<200) && (DeltaR_of_W_Higgs>=2.2 && DeltaR_of_W_Higgs<=3.5)");
	  HiggsMassReversedHptTpt[1][k] = (TH1F*)gDirectory->Get(HM2.c_str());
	  gPad->Close();
	  //	  
	  A1 = Form("Reconstructed_Tprime.M() >> TprimeMassBE1%i(60,400,1600)",k);
	  A2 = Form("TprimeMassBE1%i",k);
	  AnalysisChain.Draw(A1.c_str(),"(Reconstructed_Higgs.Pt()<200 || Reconstructed_Top.Pt()<200) && (DeltaR_of_W_Higgs>=2.2 && DeltaR_of_W_Higgs<=3.5)");
	  TprimeMassReversedHptTpt[1][k] = (TH1F*)gDirectory->Get(A2.c_str());
	  gPad->Close();
	  //	  
	  HM1 = Form("Reconstructed_Higgs.M() >> HMBE2%i(36,60,180)",k);
	  HM2 = Form("HMBE2%i",k);
	  AnalysisChain.Draw(HM1.c_str(),"(Reconstructed_Higgs.Pt()<200 || Reconstructed_Top.Pt()<200) && (DeltaR_of_W_Higgs>=2.2 && DeltaR_of_W_Higgs<=3.5) && (TMath::Abs(First_Higgs_Jet.Phi()-Second_Higgs_Jet.Phi())<=1.2 && TMath::Abs(Top_Jet.Phi()-Reconstructed_W.Phi())<=1.2)");
	  HiggsMassReversedHptTpt[2][k] = (TH1F*)gDirectory->Get(HM2.c_str());
	  gPad->Close();
	  //	  
	  A1 = Form("Reconstructed_Tprime.M() >> TprimeMassBE2%i(60,400,1600)",k);
	  A2 = Form("TprimeMassBE2%i",k);
	  AnalysisChain.Draw(A1.c_str(),"(Reconstructed_Higgs.Pt()<200 || Reconstructed_Top.Pt()<200) && (DeltaR_of_W_Higgs>=2.2 && DeltaR_of_W_Higgs<=3.5) && (TMath::Abs(First_Higgs_Jet.Phi()-Second_Higgs_Jet.Phi())<=1.2 && TMath::Abs(Top_Jet.Phi()-Reconstructed_W.Phi())<=1.2)");
	  TprimeMassReversedHptTpt[2][k] = (TH1F*)gDirectory->Get(A2.c_str());
	  gPad->Close();
	  //	  
	  HM1 = Form("Reconstructed_Higgs.M() >> HMBE3%i(36,60,180)",k);
	  HM2 = Form("HMBE3%i",k);
	  AnalysisChain.Draw(HM1.c_str(),"(Reconstructed_Higgs.Pt()<200 || Reconstructed_Top.Pt()<200) && (DeltaR_of_W_Higgs>=2.2 && DeltaR_of_W_Higgs<=3.5) && (TMath::Abs(First_Higgs_Jet.Phi()-Second_Higgs_Jet.Phi())<=1.2 && TMath::Abs(Top_Jet.Phi()-Reconstructed_W.Phi())<=1.2) && TMath::Abs(First_W_Jet.Phi()-Second_W_Jet.Phi())<=2.0");
	  HiggsMassReversedHptTpt[3][k] = (TH1F*)gDirectory->Get(HM2.c_str());
	  gPad->Close();
	  //	  
	  A1 = Form("Reconstructed_Tprime.M() >> TprimeMassBE3%i(60,400,1600)",k);
	  A2 = Form("TprimeMassBE3%i",k);
	  AnalysisChain.Draw(A1.c_str(),"(Reconstructed_Higgs.Pt()<200 || Reconstructed_Top.Pt()<200) && (DeltaR_of_W_Higgs>=2.2 && DeltaR_of_W_Higgs<=3.5) && (TMath::Abs(First_Higgs_Jet.Phi()-Second_Higgs_Jet.Phi())<=1.2 && TMath::Abs(Top_Jet.Phi()-Reconstructed_W.Phi())<=1.2) && TMath::Abs(First_W_Jet.Phi()-Second_W_Jet.Phi())<=2.0");
	  TprimeMassReversedHptTpt[3][k] = (TH1F*)gDirectory->Get(A2.c_str());
	  gPad->Close();
	  //	  
	  HM1 = Form("Reconstructed_Higgs.M() >> HMBE4%i(36,60,180)",k);
	  HM2 = Form("HMBE4%i",k);
	  AnalysisChain.Draw(HM1.c_str(),"(Reconstructed_Higgs.Pt()<200 || Reconstructed_Top.Pt()<200) && (DeltaR_of_W_Higgs>=2.2 && DeltaR_of_W_Higgs<=3.5) && (TMath::Abs(First_Higgs_Jet.Phi()-Second_Higgs_Jet.Phi())<=1.2 && TMath::Abs(Top_Jet.Phi()-Reconstructed_W.Phi())<=1.2) && TMath::Abs(First_W_Jet.Phi()-Second_W_Jet.Phi())<=2.0 && Relative_THT>=0.65");
	  HiggsMassReversedHptTpt[4][k] = (TH1F*)gDirectory->Get(HM2.c_str());
	  gPad->Close();
	  //	  
	  A1 = Form("Reconstructed_Tprime.M() >> TprimeMassBE4%i(60,400,1600)",k);
	  A2 = Form("TprimeMassBE4%i",k);
	  AnalysisChain.Draw(A1.c_str(),"(Reconstructed_Higgs.Pt()<200 || Reconstructed_Top.Pt()<200) && (DeltaR_of_W_Higgs>=2.2 && DeltaR_of_W_Higgs<=3.5) && (TMath::Abs(First_Higgs_Jet.Phi()-Second_Higgs_Jet.Phi())<=1.2 && TMath::Abs(Top_Jet.Phi()-Reconstructed_W.Phi())<=1.2) && TMath::Abs(First_W_Jet.Phi()-Second_W_Jet.Phi())<=2.0 && Relative_THT>=0.65");
	  TprimeMassReversedHptTpt[4][k] = (TH1F*)gDirectory->Get(A2.c_str());
	  //	  
	  HM1 = Form("Reconstructed_Higgs.M() >> HMBE5%i(36,60,180)",k);
	  HM2 = Form("HMBE5%i",k);
	  AnalysisChain.Draw(HM1.c_str(),"(Reconstructed_Higgs.Pt()<200 || Reconstructed_Top.Pt()<200) && (DeltaR_of_W_Higgs>=2.2 && DeltaR_of_W_Higgs<=3.5) && (TMath::Abs(First_Higgs_Jet.Phi()-Second_Higgs_Jet.Phi())<=1.2 && TMath::Abs(Top_Jet.Phi()-Reconstructed_W.Phi())<=1.2) && TMath::Abs(First_W_Jet.Phi()-Second_W_Jet.Phi())<=2.0 && Relative_THT>=0.65 && (Relative_Mass>=0.3 && Relative_Mass<=0.5)");
	  HiggsMassReversedHptTpt[5][k] = (TH1F*)gDirectory->Get(HM2.c_str());
	  gPad->Close();
	  //	  
	  A1 = Form("Reconstructed_Tprime.M() >> TprimeMassBE5%i(60,400,1600)",k);
	  A2 = Form("TprimeMassBE5%i",k);
	  AnalysisChain.Draw(A1.c_str(),"(Reconstructed_Higgs.Pt()<200 || Reconstructed_Top.Pt()<200) && (DeltaR_of_W_Higgs>=2.2 && DeltaR_of_W_Higgs<=3.5) && (TMath::Abs(First_Higgs_Jet.Phi()-Second_Higgs_Jet.Phi())<=1.2 && TMath::Abs(Top_Jet.Phi()-Reconstructed_W.Phi())<=1.2) && TMath::Abs(First_W_Jet.Phi()-Second_W_Jet.Phi())<=2.0 && Relative_THT>=0.65 && (Relative_Mass>=0.3 && Relative_Mass<=0.5)");
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
	  //
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
	  ///////////////////////
	  //Sideband estimation//
	  ///////////////////////
	  string WFH1 = Form("W_From_Higgs.M() >> WMFH0%i(44,60.0,500)",k);
	  string WFH2 = Form("WMFH0%i",k);
	  AnalysisChain.Draw(WFH1.c_str(),"Reconstructed_Higgs.Pt()>200 && Reconstructed_Top.Pt()>200");
	  WMassFromHiggs[0][k] = (TH1F*)gDirectory->Get(WFH2.c_str());
	  gPad->Close();
	  //
	  WFH1 = Form("W_From_Higgs.M() >> WMFH1%i(44,60.0,500)",k);
	  WFH2 = Form("WMFH1%i",k);
	  AnalysisChain.Draw(WFH1.c_str(),"(Reconstructed_Higgs.Pt()>200 && Reconstructed_Top.Pt()>200) && (DeltaR_of_W_Higgs>=2.2 && DeltaR_of_W_Higgs<=3.5)");
	  WMassFromHiggs[1][k] = (TH1F*)gDirectory->Get(WFH2.c_str());
	  gPad->Close();
	  //
	  WFH1 = Form("W_From_Higgs.M() >> WMFH2%i(44,60.0,500)",k);
	  WFH2 = Form("WMFH2%i",k);
	  AnalysisChain.Draw(WFH1.c_str(),"(Reconstructed_Higgs.Pt()>200 && Reconstructed_Top.Pt()>200) && (DeltaR_of_W_Higgs>=2.2 && DeltaR_of_W_Higgs<=3.5) && (TMath::Abs(First_Higgs_Jet.Phi()-Second_Higgs_Jet.Phi())<=1.2 && TMath::Abs(Top_Jet.Phi()-Reconstructed_W.Phi())<=1.2)");
	  WMassFromHiggs[2][k] = (TH1F*)gDirectory->Get(WFH2.c_str());
	  gPad->Close();
	  //
	  WFH1 = Form("W_From_Higgs.M() >> WMFH3%i(44,60.0,500)",k);
	  WFH2 = Form("WMFH3%i",k);
	  AnalysisChain.Draw(WFH1.c_str(),"(Reconstructed_Higgs.Pt()>200 && Reconstructed_Top.Pt()>200) && (DeltaR_of_W_Higgs>=2.2 && DeltaR_of_W_Higgs<=3.5) && (TMath::Abs(First_Higgs_Jet.Phi()-Second_Higgs_Jet.Phi())<=1.2 && TMath::Abs(Top_Jet.Phi()-Reconstructed_W.Phi())<=1.2) && TMath::Abs(First_W_Jet.Phi()-Second_W_Jet.Phi())<=2.0");
	  WMassFromHiggs[3][k] = (TH1F*)gDirectory->Get(WFH2.c_str());
	  gPad->Close();
	  //
	  WFH1 = Form("W_From_Higgs.M() >> WMFH4%i(44,60.0,500)",k);
	  WFH2 = Form("WMFH4%i",k);
	  AnalysisChain.Draw(WFH1.c_str(),"(Reconstructed_Higgs.Pt()>200 && Reconstructed_Top.Pt()>200) && (DeltaR_of_W_Higgs>=2.2 && DeltaR_of_W_Higgs<=3.5) && (TMath::Abs(First_Higgs_Jet.Phi()-Second_Higgs_Jet.Phi())<=1.2 && TMath::Abs(Top_Jet.Phi()-Reconstructed_W.Phi())<=1.2) && TMath::Abs(First_W_Jet.Phi()-Second_W_Jet.Phi())<=2.0 && Relative_THT>=0.65");
	  WMassFromHiggs[4][k] = (TH1F*)gDirectory->Get(WFH2.c_str());
	  gPad->Close();
	  //
	  WFH1 = Form("W_From_Higgs.M() >> WMFH5%i(44,60.0,500)",k);
	  WFH2 = Form("WMFH5%i",k);
	  AnalysisChain.Draw(WFH1.c_str(),"(Reconstructed_Higgs.Pt()>200 && Reconstructed_Top.Pt()>200) && (DeltaR_of_W_Higgs>=2.2 && DeltaR_of_W_Higgs<=3.5) && (TMath::Abs(First_Higgs_Jet.Phi()-Second_Higgs_Jet.Phi())<=1.2 && TMath::Abs(Top_Jet.Phi()-Reconstructed_W.Phi())<=1.2) && TMath::Abs(First_W_Jet.Phi()-Second_W_Jet.Phi())<=2.0 && Relative_THT>=0.65 && (Relative_Mass>=0.3 && Relative_Mass<=0.5)");
	  WMassFromHiggs[5][k] = (TH1F*)gDirectory->Get(WFH2.c_str());
	  gPad->Close();
	  //////Chi2/////////////
	  WFH1 = Form("W_From_Higgs_Chi2.M() >> WMFHC0%i(44,60.0,500)",k);
	  WFH2 = Form("WMFHC0%i",k);
	  AnalysisChain.Draw(WFH1.c_str(),"Reconstructed_Higgs.Pt()>200 && Reconstructed_Top.Pt()>200");
	  WMassFromHiggsChi2[0][k] = (TH1F*)gDirectory->Get(WFH2.c_str());
	  gPad->Close();
	  //
	  WFH1 = Form("W_From_Higgs_Chi2.M() >> WMFHC1%i(44,60.0,500)",k);
	  WFH2 = Form("WMFHC1%i",k);
	  AnalysisChain.Draw(WFH1.c_str(),"(Reconstructed_Higgs.Pt()>200 && Reconstructed_Top.Pt()>200) && (DeltaR_of_W_Higgs>=2.2 && DeltaR_of_W_Higgs<=3.5)");
	  WMassFromHiggsChi2[1][k] = (TH1F*)gDirectory->Get(WFH2.c_str());
	  gPad->Close();
	  //
	  WFH1 = Form("W_From_Higgs_Chi2.M() >> WMFHC2%i(44,60.0,500)",k);
	  WFH2 = Form("WMFHC2%i",k);
	  AnalysisChain.Draw(WFH1.c_str(),"(Reconstructed_Higgs.Pt()>200 && Reconstructed_Top.Pt()>200) && (DeltaR_of_W_Higgs>=2.2 && DeltaR_of_W_Higgs<=3.5) && (TMath::Abs(First_Higgs_Jet.Phi()-Second_Higgs_Jet.Phi())<=1.2 && TMath::Abs(Top_Jet.Phi()-Reconstructed_W.Phi())<=1.2)");
	  WMassFromHiggsChi2[2][k] = (TH1F*)gDirectory->Get(WFH2.c_str());
	  gPad->Close();
	  //
	  WFH1 = Form("W_From_Higgs_Chi2.M() >> WMFHC3%i(44,60.0,500)",k);
	  WFH2 = Form("WMFHC3%i",k);
	  AnalysisChain.Draw(WFH1.c_str(),"(Reconstructed_Higgs.Pt()>200 && Reconstructed_Top.Pt()>200) && (DeltaR_of_W_Higgs>=2.2 && DeltaR_of_W_Higgs<=3.5) && (TMath::Abs(First_Higgs_Jet.Phi()-Second_Higgs_Jet.Phi())<=1.2 && TMath::Abs(Top_Jet.Phi()-Reconstructed_W.Phi())<=1.2) && TMath::Abs(First_W_Jet.Phi()-Second_W_Jet.Phi())<=2.0");
	  WMassFromHiggsChi2[3][k] = (TH1F*)gDirectory->Get(WFH2.c_str());
	  gPad->Close();
	  //
	  WFH1 = Form("W_From_Higgs_Chi2.M() >> WMFHC4%i(44,60.0,500)",k);
	  WFH2 = Form("WMFHC4%i",k);
	  AnalysisChain.Draw(WFH1.c_str(),"(Reconstructed_Higgs.Pt()>200 && Reconstructed_Top.Pt()>200) && (DeltaR_of_W_Higgs>=2.2 && DeltaR_of_W_Higgs<=3.5) && (TMath::Abs(First_Higgs_Jet.Phi()-Second_Higgs_Jet.Phi())<=1.2 && TMath::Abs(Top_Jet.Phi()-Reconstructed_W.Phi())<=1.2) && TMath::Abs(First_W_Jet.Phi()-Second_W_Jet.Phi())<=2.0 && Relative_THT>=0.65");
	  WMassFromHiggsChi2[4][k] = (TH1F*)gDirectory->Get(WFH2.c_str());
	  gPad->Close();
	  //
	  WFH1 = Form("W_From_Higgs_Chi2.M() >> WMFHC5%i(44,60.0,500)",k);
	  WFH2 = Form("WMFHC5%i",k);
	  AnalysisChain.Draw(WFH1.c_str(),"(Reconstructed_Higgs.Pt()>200 && Reconstructed_Top.Pt()>200) && (DeltaR_of_W_Higgs>=2.2 && DeltaR_of_W_Higgs<=3.5) && (TMath::Abs(First_Higgs_Jet.Phi()-Second_Higgs_Jet.Phi())<=1.2 && TMath::Abs(Top_Jet.Phi()-Reconstructed_W.Phi())<=1.2) && TMath::Abs(First_W_Jet.Phi()-Second_W_Jet.Phi())<=2.0 && Relative_THT>=0.65 && (Relative_Mass>=0.3 && Relative_Mass<=0.5)");
	  WMassFromHiggsChi2[5][k] = (TH1F*)gDirectory->Get(WFH2.c_str());
	  gPad->Close();
	  //////////////////////////////////////////////////////////////////////////////
	  //////////////////////////////////////////////////////////////////////////////
	  //////////////////////////////////////////////////////////////////////////////
	  A1 = Form("Reconstructed_Tprime3L.M() >> TprimeMass3L%i(60,400,1600)",k);
	  A2 = Form("TprimeMass3L%i",k);
	  AnalysisChain3L.Draw(A1.c_str(),"(Reconstructed_Top3L.M()>=230)");
	  FiveJetsMass3L[k] = (TH1F*)gDirectory->Get(A2.c_str());
	  gPad->Close();
	  //
	  A1 = Form("Reconstructed_Tprime3L.M() >> TprimeMass3L_L%i(60,400,1600)",k);
	  A2 = Form("TprimeMass3L_L%i",k);
	  AnalysisChain3L.Draw(A1.c_str(),"(Reconstructed_Top3L.M()<230)");
	  FiveJetsMass3L_L[k] = (TH1F*)gDirectory->Get(A2.c_str());
	  gPad->Close();
	  //	
	  /*A1 = Form("Reconstructed_Tprime3T.M() >> TprimeMass3T%i(60,400,1600)",k);
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
  cout << "HERE1" << endl;
  cout << FirstGoodMarker << LastGoodMarker << endl;
  //TprimeMass histo recollection
  string A=Form("TprimeMass%i",FirstGoodMarker);
  TH1F *QCDTprimeMass=(TH1F*)gDirectory->Get(A.c_str());  
  QCDTprimeMass->Scale(Lumi*XS[FirstGoodMarker]*0.5/ProcessedEvents[FirstGoodMarker]);
  //TprimeMass3L histo recollection
  A=Form("TprimeMass3L%i",FirstGoodMarker);
  TH1F *QCDTprimeMass3L=(TH1F*)gDirectory->Get(A.c_str());  
  QCDTprimeMass3L->Scale(Lumi*XS[FirstGoodMarker]*0.5/ProcessedEvents[FirstGoodMarker]);
  //TprimeMass3T histo recollection
  //A=Form("TprimeMass3T%i",FirstGoodMarker);
  //TH1F *QCDTprimeMass3T=(TH1F*)gDirectory->Get(A.c_str());  
  //QCDTprimeMass3T->Scale(Lumi*XS[FirstGoodMarker]*0.5/ProcessedEvents[FirstGoodMarker]);
  //TprimeMass3L_L histo recollection
  A=Form("TprimeMass3L_L%i",FirstGoodMarker);
  TH1F *QCDTprimeMass3L_L=(TH1F*)gDirectory->Get(A.c_str());  
  QCDTprimeMass3L_L->Scale(Lumi*XS[FirstGoodMarker]*0.5/ProcessedEvents[FirstGoodMarker]);
  //TprimeMass3T_L histo recollection
  //A=Form("TprimeMass3T_L%i",FirstGoodMarker);
  //TH1F *QCDTprimeMass3T_L=(TH1F*)gDirectory->Get(A.c_str());  
  //QCDTprimeMass3T_L->Scale(Lumi*XS[FirstGoodMarker]*0.5/ProcessedEvents[FirstGoodMarker]);
  //LeadingJetPT histo recollection
  A=Form("jet1_pt%i",FirstGoodMarker);
  TH1F *QCDTprimeMassLJPT=(TH1F*)gDirectory->Get(A.c_str()); 
  QCDTprimeMassLJPT->Scale(Lumi*XS[FirstGoodMarker]*0.5/ProcessedEvents[FirstGoodMarker]); 
  //LeadingJetPT histo recollection
  A=Form("jet2_pt%i",FirstGoodMarker);
  TH1F *QCDTprimeMassL2JPT=(TH1F*)gDirectory->Get(A.c_str());
  QCDTprimeMassL2JPT->Scale(Lumi*XS[FirstGoodMarker]*0.5/ProcessedEvents[FirstGoodMarker]);
  //LeadingJetPT histo recollection
  A=Form("jet3_pt%i",FirstGoodMarker);
  TH1F *QCDTprimeMassL3JPT=(TH1F*)gDirectory->Get(A.c_str());
  QCDTprimeMassL3JPT->Scale(Lumi*XS[FirstGoodMarker]*0.5/ProcessedEvents[FirstGoodMarker]);
  //LeadingJetPT histo recollection
  A=Form("jet4_pt%i",FirstGoodMarker);
  TH1F *QCDTprimeMassL4JPT=(TH1F*)gDirectory->Get(A.c_str());
  QCDTprimeMassL4JPT->Scale(Lumi*XS[FirstGoodMarker]*0.5/ProcessedEvents[FirstGoodMarker]);
  //LeadingJetPT histo recollection
  A=Form("jet5_pt%i",FirstGoodMarker);
  TH1F *QCDTprimeMassL5JPT=(TH1F*)gDirectory->Get(A.c_str());
  QCDTprimeMassL5JPT->Scale(Lumi*XS[FirstGoodMarker]*0.5/ProcessedEvents[FirstGoodMarker]);
  //LeadingJetPT histo recollection
  A=Form("jet6_pt%i",FirstGoodMarker);
  TH1F *QCDTprimeMassL6JPT=(TH1F*)gDirectory->Get(A.c_str());
  QCDTprimeMassL6JPT->Scale(Lumi*XS[FirstGoodMarker]*0.5/ProcessedEvents[FirstGoodMarker]);
  //THT histo recollection
  A=Form("THT%i",FirstGoodMarker);
  TH1F *QCDTprimeMassTHT=(TH1F*)gDirectory->Get(A.c_str());  
  QCDTprimeMassTHT->Scale(Lumi*XS[FirstGoodMarker]*0.5/ProcessedEvents[FirstGoodMarker]);
  //DR Higgs jets histo recollection
  A=Form("DeltaR_of_Higgs_Jets%i",FirstGoodMarker);
  TH1F *QCDTprimeMassDRHJ=(TH1F*)gDirectory->Get(A.c_str());
  QCDTprimeMassDRHJ->Scale(Lumi*XS[FirstGoodMarker]*0.5/ProcessedEvents[FirstGoodMarker]);
  //DR W jets histo recollection
  A=Form("DeltaR_of_W_Jets%i",FirstGoodMarker);
  TH1F *QCDTprimeMassDRWJ=(TH1F*)gDirectory->Get(A.c_str());
  QCDTprimeMassDRWJ->Scale(Lumi*XS[FirstGoodMarker]*0.5/ProcessedEvents[FirstGoodMarker]);
  //Higgs PT histo recollection
  A=Form("HPt%i",FirstGoodMarker);
  TH1F *QCDTprimeMassHPT=(TH1F*)gDirectory->Get(A.c_str());
  QCDTprimeMassHPT->Scale(Lumi*XS[FirstGoodMarker]*0.5/ProcessedEvents[FirstGoodMarker]);
  //Top PT histo recollection
  A=Form("TPt%i",FirstGoodMarker);
  TH1F *QCDTprimeMassTPT=(TH1F*)gDirectory->Get(A.c_str());
  QCDTprimeMassTPT->Scale(Lumi*XS[FirstGoodMarker]*0.5/ProcessedEvents[FirstGoodMarker]);
  //DR W Higgs histo recollection
  A=Form("DeltaR_of_W_Higgs%i",FirstGoodMarker);
  TH1F *QCDTprimeMassDRWH=(TH1F*)gDirectory->Get(A.c_str());
  QCDTprimeMassDRWH->Scale(Lumi*XS[FirstGoodMarker]*0.5/ProcessedEvents[FirstGoodMarker]);
  //DP Higgs jets histo recollection
  A=Form("DeltaPhi_of_Higgs_jets%i",FirstGoodMarker);
  TH1F *QCDTprimeMassDPHJ=(TH1F*)gDirectory->Get(A.c_str());
  QCDTprimeMassDPHJ->Scale(Lumi*XS[FirstGoodMarker]*0.5/ProcessedEvents[FirstGoodMarker]);
  //DP W jets histo recollection
  A=Form("DeltaPhi_of_W_jets%i",FirstGoodMarker);
  TH1F *QCDTprimeMassDPWJ=(TH1F*)gDirectory->Get(A.c_str());
  QCDTprimeMassDPWJ->Scale(Lumi*XS[FirstGoodMarker]*0.5/ProcessedEvents[FirstGoodMarker]);
  //DP Top jet histo recollection
  A=Form("DeltaPhi_of_T_jet%i",FirstGoodMarker);
  TH1F *QCDTprimeMassDPTJ=(TH1F*)gDirectory->Get(A.c_str());
  QCDTprimeMassDPTJ->Scale(Lumi*XS[FirstGoodMarker]*0.5/ProcessedEvents[FirstGoodMarker]);
  //Higgs Mass histo recollection
  A=Form("HM%i",FirstGoodMarker);
  TH1F *QCDTprimeMassHM=(TH1F*)gDirectory->Get(A.c_str());
  QCDTprimeMassHM->Scale(Lumi*XS[FirstGoodMarker]*0.5/ProcessedEvents[FirstGoodMarker]);
  //Relative HT histo recollection
  A=Form("RelHT%i",FirstGoodMarker);
  TH1F *QCDTprimeMassRelHT=(TH1F*)gDirectory->Get(A.c_str());
  QCDTprimeMassRelHT->Scale(Lumi*XS[FirstGoodMarker]*0.5/ProcessedEvents[FirstGoodMarker]);
  //DR Top Higgs histo recollection
  A=Form("DeltaR_of_Top_Higgs%i",FirstGoodMarker);
  TH1F *QCDTprimeMassDRTH=(TH1F*)gDirectory->Get(A.c_str());
  QCDTprimeMassDRTH->Scale(Lumi*XS[FirstGoodMarker]*0.5/ProcessedEvents[FirstGoodMarker]);
  //PT Normalized Mass histo recollection
  A=Form("PT_Normalized_Mass%i",FirstGoodMarker);
  TH1F *QCDTprimeMassPTNM=(TH1F*)gDirectory->Get(A.c_str());
  QCDTprimeMassPTNM->Scale(Lumi*XS[FirstGoodMarker]*0.5/ProcessedEvents[FirstGoodMarker]);
  //Relative Mass histo recollection
  A=Form("Relative_Mass%i",FirstGoodMarker);
  TH1F *QCDTprimeMassRM=(TH1F*)gDirectory->Get(A.c_str());
  QCDTprimeMassRM->Scale(Lumi*XS[FirstGoodMarker]*0.5/ProcessedEvents[FirstGoodMarker]);
  //Mother PT Normalized Mass histo recollection
  A=Form("Mother_PT_Normalized_Mass%i",FirstGoodMarker);
  TH1F *QCDTprimeMassMPTNM=(TH1F*)gDirectory->Get(A.c_str());
  QCDTprimeMassMPTNM->Scale(Lumi*XS[FirstGoodMarker]*0.5/ProcessedEvents[FirstGoodMarker]);
  //Number of Tops histo recollection
  A=Form("Number_of_Tops%i",FirstGoodMarker);
  TH1F *QCDTprimeMassNT=(TH1F*)gDirectory->Get(A.c_str());
  QCDTprimeMassNT->Scale(Lumi*XS[FirstGoodMarker]*0.5/ProcessedEvents[FirstGoodMarker]);
  //Higgs Mass over Top Mass histo recollection
  A=Form("HMoverTM%i",FirstGoodMarker);
  TH1F *QCDTprimeMassHMTM=(TH1F*)gDirectory->Get(A.c_str());
  QCDTprimeMassHMTM->Scale(Lumi*XS[FirstGoodMarker]*0.5/ProcessedEvents[FirstGoodMarker]);
  //Higgs Top Asym histo recollection
  A=Form("HTAsym%i",FirstGoodMarker);
  TH1F *QCDTprimeMassHTAsym=(TH1F*)gDirectory->Get(A.c_str());
  QCDTprimeMassHTAsym->Scale(Lumi*XS[FirstGoodMarker]*0.5/ProcessedEvents[FirstGoodMarker]);
  //Number of loose and non medium b-tag histo recollection
  A=Form("TLBTag%i",FirstGoodMarker);
  TH1F *QCDTprimeMassTLBT=(TH1F*)gDirectory->Get(A.c_str());
  QCDTprimeMassTLBT->Scale(Lumi*XS[FirstGoodMarker]*0.5/ProcessedEvents[FirstGoodMarker]);
  //Top Mass histo recollection
  A=Form("TMass%i",FirstGoodMarker);
  TH1F *QCDTprimeMassTM=(TH1F*)gDirectory->Get(A.c_str());
  QCDTprimeMassTM->Scale(Lumi*XS[FirstGoodMarker]*0.5/ProcessedEvents[FirstGoodMarker]);
  //Chi2 histo recollection
  A=Form("ChiSq%i",FirstGoodMarker);
  TH1F *QCDTprimeMassChi2=(TH1F*)gDirectory->Get(A.c_str());
  QCDTprimeMassChi2->Scale(Lumi*XS[FirstGoodMarker]*0.5/ProcessedEvents[FirstGoodMarker]);
  //U quark content histo recollection
  A=Form("UQC%i",FirstGoodMarker);
  TH1F *QCDTprimeMassUQC=(TH1F*)gDirectory->Get(A.c_str());
  QCDTprimeMassUQC->Scale(Lumi*XS[FirstGoodMarker]*0.5/ProcessedEvents[FirstGoodMarker]);
  //D quark content histo recollection
  A=Form("DQC%i",FirstGoodMarker);
  TH1F *QCDTprimeMassDQC=(TH1F*)gDirectory->Get(A.c_str());
  QCDTprimeMassDQC->Scale(Lumi*XS[FirstGoodMarker]*0.5/ProcessedEvents[FirstGoodMarker]);
  //S quark content histo recollection
  A=Form("SQC%i",FirstGoodMarker);
  TH1F *QCDTprimeMassSQC=(TH1F*)gDirectory->Get(A.c_str());
  QCDTprimeMassSQC->Scale(Lumi*XS[FirstGoodMarker]*0.5/ProcessedEvents[FirstGoodMarker]);
  //C quark content histo recollection
  A=Form("CQC%i",FirstGoodMarker);
  TH1F *QCDTprimeMassCQC=(TH1F*)gDirectory->Get(A.c_str());
  QCDTprimeMassCQC->Scale(Lumi*XS[FirstGoodMarker]*0.5/ProcessedEvents[FirstGoodMarker]);
  //B quark content histo recollection
  A=Form("BQC%i",FirstGoodMarker);
  TH1F *QCDTprimeMassBQC=(TH1F*)gDirectory->Get(A.c_str());
  QCDTprimeMassBQC->Scale(Lumi*XS[FirstGoodMarker]*0.5/ProcessedEvents[FirstGoodMarker]);
  //CSVL tagging
  A=Form("CSVL%i",FirstGoodMarker);
  TH1F *QCDTprimeMassCSVL=(TH1F*)gDirectory->Get(A.c_str());
  QCDTprimeMassCSVL->Scale(Lumi*XS[FirstGoodMarker]*0.5/ProcessedEvents[FirstGoodMarker]);
  //CSVM tagging
  A=Form("CSVM%i",FirstGoodMarker);
  TH1F *QCDTprimeMassCSVM=(TH1F*)gDirectory->Get(A.c_str());
  QCDTprimeMassCSVM->Scale(Lumi*XS[FirstGoodMarker]*0.5/ProcessedEvents[FirstGoodMarker]);
  //CSVT tagging
  A=Form("CSVT%i",FirstGoodMarker);
  TH1F *QCDTprimeMassCSVT=(TH1F*)gDirectory->Get(A.c_str());
  QCDTprimeMassCSVT->Scale(Lumi*XS[FirstGoodMarker]*0.5/ProcessedEvents[FirstGoodMarker]);
  //Vertices
  A=Form("VTX%i",FirstGoodMarker);
  TH1F *QCDTprimeMassVTX=(TH1F*)gDirectory->Get(A.c_str());
  QCDTprimeMassVTX->Scale(Lumi*XS[FirstGoodMarker]*0.5/ProcessedEvents[FirstGoodMarker]);
  for (int k=ParcialTestMin; k<ParcialTestMax; k++)
  //Correct Loose B-tag
  /*A=Form("CLooseB%i",FirstGoodMarker);
  TH1F *QCDTprimeMassCLBT=(TH1F*)gDirectory->Get(A.c_str());
  QCDTprimeMassCLBT->Scale(Lumi*XS[FirstGoodMarker]*0.5/ProcessedEvents[FirstGoodMarker]);
  //Correct Med B-tag
  A=Form("CMedB%i",FirstGoodMarker);
  TH1F *QCDTprimeMassCMBT=(TH1F*)gDirectory->Get(A.c_str());
  QCDTprimeMassCMBT->Scale(Lumi*XS[FirstGoodMarker]*0.5/ProcessedEvents[FirstGoodMarker]);
  //Correct Tight B-tag
  A=Form("CTightB%i",FirstGoodMarker);
  TH1F *QCDTprimeMassCTBT=(TH1F*)gDirectory->Get(A.c_str());
  QCDTprimeMassCTBT->Scale(Lumi*XS[FirstGoodMarker]*0.5/ProcessedEvents[FirstGoodMarker]);
  //Fake light Loose B-tag
  A=Form("FLLooseB%i",FirstGoodMarker);
  TH1F *QCDTprimeMassFLLBT=(TH1F*)gDirectory->Get(A.c_str());
  QCDTprimeMassFLLBT->Scale(Lumi*XS[FirstGoodMarker]*0.5/ProcessedEvents[FirstGoodMarker]);
  //Fake light Med B-tag
  A=Form("FLMedB%i",FirstGoodMarker);
  TH1F *QCDTprimeMassFLMBT=(TH1F*)gDirectory->Get(A.c_str());
  QCDTprimeMassFLMBT->Scale(Lumi*XS[FirstGoodMarker]*0.5/ProcessedEvents[FirstGoodMarker]);
  //Fake light Tight B-tag
  A=Form("FLTightB%i",FirstGoodMarker);
  TH1F *QCDTprimeMassFLTBT=(TH1F*)gDirectory->Get(A.c_str());
  QCDTprimeMassFLTBT->Scale(Lumi*XS[FirstGoodMarker]*0.5/ProcessedEvents[FirstGoodMarker]);
  //Fake C Loose B-tag
  A=Form("FCLooseB%i",FirstGoodMarker);
  TH1F *QCDTprimeMassFCLBT=(TH1F*)gDirectory->Get(A.c_str());
  QCDTprimeMassFCLBT->Scale(Lumi*XS[FirstGoodMarker]*0.5/ProcessedEvents[FirstGoodMarker]);
  //Fake C Med B-tag
  A=Form("FCMedB%i",FirstGoodMarker);
  TH1F *QCDTprimeMassFCMBT=(TH1F*)gDirectory->Get(A.c_str());
  QCDTprimeMassFCMBT->Scale(Lumi*XS[FirstGoodMarker]*0.5/ProcessedEvents[FirstGoodMarker]);
  //Fake C Tight B-tag
  A=Form("FCTightB%i",FirstGoodMarker);
  TH1F *QCDTprimeMassFCTBT=(TH1F*)gDirectory->Get(A.c_str());
  QCDTprimeMassFCTBT->Scale(Lumi*XS[FirstGoodMarker]*0.5/ProcessedEvents[FirstGoodMarker]);*/
  //HTCSVM
  A=Form("HT_CSVM%i",FirstGoodMarker);
  TH2F *QCD_HTCSVM=(TH2F*)gDirectory->Get(A.c_str());
  QCD_HTCSVM->Scale(Lumi*XS[FirstGoodMarker]*0.5/ProcessedEvents[FirstGoodMarker]);
  //RelHTNTops
  A=Form("RHTNT%i",FirstGoodMarker);
  TH2F *QCD_RHTNT=(TH2F*)gDirectory->Get(A.c_str());
  QCD_RHTNT->Scale(Lumi*XS[FirstGoodMarker]*0.5/ProcessedEvents[FirstGoodMarker]);
  ///////////BKG Estimation///////////////
  TH1F *QCDTprimeMassHMBE[6];
  TH1F *QCDTprimeMassBE[6];
  //Higgs Mass BE histo recollection
  A=Form("HMBE0%i",FirstGoodMarker);
  QCDTprimeMassHMBE[0]=(TH1F*)gDirectory->Get(A.c_str());
  QCDTprimeMassHMBE[0]->Scale(Lumi*XS[FirstGoodMarker]*0.5/ProcessedEvents[FirstGoodMarker]);
  //TprimeMass histo recollection
  A=Form("TprimeMassBE0%i",FirstGoodMarker);
  QCDTprimeMassBE[0]=(TH1F*)gDirectory->Get(A.c_str());
  QCDTprimeMassBE[0]->Scale(Lumi*XS[FirstGoodMarker]*0.5/ProcessedEvents[FirstGoodMarker]);
  //
  A=Form("HMBE1%i",FirstGoodMarker);
  QCDTprimeMassHMBE[1]=(TH1F*)gDirectory->Get(A.c_str());
  QCDTprimeMassHMBE[1]->Scale(Lumi*XS[FirstGoodMarker]*0.5/ProcessedEvents[FirstGoodMarker]);
  //
  A=Form("TprimeMassBE1%i",FirstGoodMarker);
  QCDTprimeMassBE[1]=(TH1F*)gDirectory->Get(A.c_str());
  QCDTprimeMassBE[1]->Scale(Lumi*XS[FirstGoodMarker]*0.5/ProcessedEvents[FirstGoodMarker]);
  //
  A=Form("HMBE2%i",FirstGoodMarker);
  QCDTprimeMassHMBE[2]=(TH1F*)gDirectory->Get(A.c_str());
  QCDTprimeMassHMBE[2]->Scale(Lumi*XS[FirstGoodMarker]*0.5/ProcessedEvents[FirstGoodMarker]);
  //
  A=Form("TprimeMassBE2%i",FirstGoodMarker);
  QCDTprimeMassBE[2]=(TH1F*)gDirectory->Get(A.c_str());
  QCDTprimeMassBE[2]->Scale(Lumi*XS[FirstGoodMarker]*0.5/ProcessedEvents[FirstGoodMarker]);
  //
  A=Form("HMBE3%i",FirstGoodMarker);
  QCDTprimeMassHMBE[3]=(TH1F*)gDirectory->Get(A.c_str());
  QCDTprimeMassHMBE[3]->Scale(Lumi*XS[FirstGoodMarker]*0.5/ProcessedEvents[FirstGoodMarker]);
  //
  A=Form("TprimeMassBE3%i",FirstGoodMarker);
  QCDTprimeMassBE[3]=(TH1F*)gDirectory->Get(A.c_str());
  QCDTprimeMassBE[3]->Scale(Lumi*XS[FirstGoodMarker]*0.5/ProcessedEvents[FirstGoodMarker]);
  //
  A=Form("HMBE4%i",FirstGoodMarker);
  QCDTprimeMassHMBE[4]=(TH1F*)gDirectory->Get(A.c_str());
  QCDTprimeMassHMBE[4]->Scale(Lumi*XS[FirstGoodMarker]*0.5/ProcessedEvents[FirstGoodMarker]);
  //
  A=Form("TprimeMassBE4%i",FirstGoodMarker);
  QCDTprimeMassBE[4]=(TH1F*)gDirectory->Get(A.c_str());
  QCDTprimeMassBE[4]->Scale(Lumi*XS[FirstGoodMarker]*0.5/ProcessedEvents[FirstGoodMarker]);
  //
  A=Form("HMBE5%i",FirstGoodMarker);
  QCDTprimeMassHMBE[5]=(TH1F*)gDirectory->Get(A.c_str());
  QCDTprimeMassHMBE[5]->Scale(Lumi*XS[FirstGoodMarker]*0.5/ProcessedEvents[FirstGoodMarker]);
  //
  A=Form("TprimeMassBE5%i",FirstGoodMarker);
  QCDTprimeMassBE[5]=(TH1F*)gDirectory->Get(A.c_str());
  QCDTprimeMassBE[5]->Scale(Lumi*XS[FirstGoodMarker]*0.5/ProcessedEvents[FirstGoodMarker]);  
  //HptToppt
  TH2F *QCD_HptTpt[6];
  A=Form("HPtTPt0%i",FirstGoodMarker);
  QCD_HptTpt[0]=(TH2F*)gDirectory->Get(A.c_str());
  QCD_HptTpt[0]->Scale(Lumi*XS[FirstGoodMarker]*0.5/ProcessedEvents[FirstGoodMarker]);
  //HptToppt
  A=Form("HPtTPt1%i",FirstGoodMarker);
  QCD_HptTpt[1]=(TH2F*)gDirectory->Get(A.c_str());
  QCD_HptTpt[1]->Scale(Lumi*XS[FirstGoodMarker]*0.5/ProcessedEvents[FirstGoodMarker]);
  //HptToppt
  A=Form("HPtTPt2%i",FirstGoodMarker);
  QCD_HptTpt[2]=(TH2F*)gDirectory->Get(A.c_str());
  QCD_HptTpt[2]->Scale(Lumi*XS[FirstGoodMarker]*0.5/ProcessedEvents[FirstGoodMarker]);
  //HptToppt
  A=Form("HPtTPt3%i",FirstGoodMarker);
  QCD_HptTpt[3]=(TH2F*)gDirectory->Get(A.c_str());
  QCD_HptTpt[3]->Scale(Lumi*XS[FirstGoodMarker]*0.5/ProcessedEvents[FirstGoodMarker]);
  //HptToppt
  A=Form("HPtTPt4%i",FirstGoodMarker);
  QCD_HptTpt[4]=(TH2F*)gDirectory->Get(A.c_str());
  QCD_HptTpt[4]->Scale(Lumi*XS[FirstGoodMarker]*0.5/ProcessedEvents[FirstGoodMarker]);
  //HptToppt
  A=Form("HPtTPt5%i",FirstGoodMarker);
  QCD_HptTpt[5]=(TH2F*)gDirectory->Get(A.c_str());
  QCD_HptTpt[5]->Scale(Lumi*XS[FirstGoodMarker]*0.5/ProcessedEvents[FirstGoodMarker]);
  TH1F *QCDTprimeMassHMCuts[6];
  //Higgs Mass histo recollection
  A=Form("HM0%i",FirstGoodMarker);
  QCDTprimeMassHMCuts[0]=(TH1F*)gDirectory->Get(A.c_str());
  QCDTprimeMassHMCuts[0]->Scale(Lumi*XS[FirstGoodMarker]*0.5/ProcessedEvents[FirstGoodMarker]);
  //
  A=Form("HM1%i",FirstGoodMarker);
  QCDTprimeMassHMCuts[1]=(TH1F*)gDirectory->Get(A.c_str());
  QCDTprimeMassHMCuts[1]->Scale(Lumi*XS[FirstGoodMarker]*0.5/ProcessedEvents[FirstGoodMarker]);
  //
  A=Form("HM2%i",FirstGoodMarker);
  QCDTprimeMassHMCuts[2]=(TH1F*)gDirectory->Get(A.c_str());
  QCDTprimeMassHMCuts[2]->Scale(Lumi*XS[FirstGoodMarker]*0.5/ProcessedEvents[FirstGoodMarker]);
  //
  A=Form("HM3%i",FirstGoodMarker);
  QCDTprimeMassHMCuts[3]=(TH1F*)gDirectory->Get(A.c_str());
  QCDTprimeMassHMCuts[3]->Scale(Lumi*XS[FirstGoodMarker]*0.5/ProcessedEvents[FirstGoodMarker]);
  //
  A=Form("HM4%i",FirstGoodMarker);
  QCDTprimeMassHMCuts[4]=(TH1F*)gDirectory->Get(A.c_str());
  QCDTprimeMassHMCuts[4]->Scale(Lumi*XS[FirstGoodMarker]*0.5/ProcessedEvents[FirstGoodMarker]);
  //
  A=Form("HM5%i",FirstGoodMarker);
  QCDTprimeMassHMCuts[5]=(TH1F*)gDirectory->Get(A.c_str()); 
  QCDTprimeMassHMCuts[5]->Scale(Lumi*XS[FirstGoodMarker]*0.5/ProcessedEvents[FirstGoodMarker]);
  //W Mass from Higgs histo recollection
  TH1F *QCDTprimeMassWMFHCuts[6];
  A=Form("WMFH0%i",FirstGoodMarker);
  QCDTprimeMassWMFHCuts[0]=(TH1F*)gDirectory->Get(A.c_str());
  QCDTprimeMassWMFHCuts[0]->Scale(Lumi*XS[FirstGoodMarker]*0.5/ProcessedEvents[FirstGoodMarker]); 
  //
  A=Form("WMFH1%i",FirstGoodMarker);
  QCDTprimeMassWMFHCuts[1]=(TH1F*)gDirectory->Get(A.c_str()); 
  QCDTprimeMassWMFHCuts[1]->Scale(Lumi*XS[FirstGoodMarker]*0.5/ProcessedEvents[FirstGoodMarker]);
  //
  A=Form("WMFH2%i",FirstGoodMarker);
  QCDTprimeMassWMFHCuts[2]=(TH1F*)gDirectory->Get(A.c_str()); 
  QCDTprimeMassWMFHCuts[2]->Scale(Lumi*XS[FirstGoodMarker]*0.5/ProcessedEvents[FirstGoodMarker]);
  //
  A=Form("WMFH3%i",FirstGoodMarker);
  QCDTprimeMassWMFHCuts[3]=(TH1F*)gDirectory->Get(A.c_str()); 
  QCDTprimeMassWMFHCuts[3]->Scale(Lumi*XS[FirstGoodMarker]*0.5/ProcessedEvents[FirstGoodMarker]);
  //
  A=Form("WMFH4%i",FirstGoodMarker);
  QCDTprimeMassWMFHCuts[4]=(TH1F*)gDirectory->Get(A.c_str()); 
  QCDTprimeMassWMFHCuts[4]->Scale(Lumi*XS[FirstGoodMarker]*0.5/ProcessedEvents[FirstGoodMarker]);
  //
  A=Form("WMFH5%i",FirstGoodMarker);
  QCDTprimeMassWMFHCuts[5]=(TH1F*)gDirectory->Get(A.c_str()); 
  QCDTprimeMassWMFHCuts[5]->Scale(Lumi*XS[FirstGoodMarker]*0.5/ProcessedEvents[FirstGoodMarker]);
  //  
  TH1F *QCDTprimeMassWMFHCCuts[6];
  A=Form("WMFHC0%i",FirstGoodMarker);
  QCDTprimeMassWMFHCCuts[0]=(TH1F*)gDirectory->Get(A.c_str()); 
  QCDTprimeMassWMFHCCuts[0]->Scale(Lumi*XS[FirstGoodMarker]*0.5/ProcessedEvents[FirstGoodMarker]);
  //
  A=Form("WMFHC1%i",FirstGoodMarker);
  QCDTprimeMassWMFHCCuts[1]=(TH1F*)gDirectory->Get(A.c_str()); 
  QCDTprimeMassWMFHCCuts[1]->Scale(Lumi*XS[FirstGoodMarker]*0.5/ProcessedEvents[FirstGoodMarker]);
  //
  A=Form("WMFHC2%i",FirstGoodMarker);
  QCDTprimeMassWMFHCCuts[2]=(TH1F*)gDirectory->Get(A.c_str());
  QCDTprimeMassWMFHCCuts[2]->Scale(Lumi*XS[FirstGoodMarker]*0.5/ProcessedEvents[FirstGoodMarker]); 
  //
  A=Form("WMFHC3%i",FirstGoodMarker);
  QCDTprimeMassWMFHCCuts[3]=(TH1F*)gDirectory->Get(A.c_str());
  QCDTprimeMassWMFHCCuts[3]->Scale(Lumi*XS[FirstGoodMarker]*0.5/ProcessedEvents[FirstGoodMarker]); 
  //
  A=Form("WMFHC4%i",FirstGoodMarker);
  QCDTprimeMassWMFHCCuts[4]=(TH1F*)gDirectory->Get(A.c_str()); 
  QCDTprimeMassWMFHCCuts[4]->Scale(Lumi*XS[FirstGoodMarker]*0.5/ProcessedEvents[FirstGoodMarker]);
  //
  A=Form("WMFHC5%i",FirstGoodMarker);
  QCDTprimeMassWMFHCCuts[5]=(TH1F*)gDirectory->Get(A.c_str()); 
  QCDTprimeMassWMFHCCuts[5]->Scale(Lumi*XS[FirstGoodMarker]*0.5/ProcessedEvents[FirstGoodMarker]);
  //TprimeMass BE histo recollection
  A=Form("TprimeMassBkgE%i",FirstGoodMarker);
  TH1F *QCDTprimeMassBETop=(TH1F*)gDirectory->Get(A.c_str());
  QCDTprimeMassBETop->Scale(Lumi*XS[FirstGoodMarker]*0.5/ProcessedEvents[FirstGoodMarker]);
  //TprimeMass LC histo recollection
  A=Form("TprimeMassLC%i",FirstGoodMarker);
  TH1F *QCDTprimeMassLC=(TH1F*)gDirectory->Get(A.c_str());
  QCDTprimeMassLC->Scale(Lumi*XS[FirstGoodMarker]*0.5/ProcessedEvents[FirstGoodMarker]);
  //TprimeMass LC over BE histo recollection
  A=Form("TprimeMassLCoverBE%i",FirstGoodMarker);
  TH1F *QCDTprimeMassLCoverBE=(TH1F*)gDirectory->Get(A.c_str());
  QCDTprimeMassLCoverBE->Scale(Lumi*XS[FirstGoodMarker]*0.5/ProcessedEvents[FirstGoodMarker]);
  //TopMassFromHiggs histo recollection
  A=Form("TopFromHiggsChi2Mass%i",FirstGoodMarker);
  TH1F *QCDTopMassFromHiggs=(TH1F*)gDirectory->Get(A.c_str());
  QCDTopMassFromHiggs->Scale(Lumi*XS[FirstGoodMarker]*0.5/ProcessedEvents[FirstGoodMarker]);
  
  //HM BE histo recollection
  A=Form("HiggsMassBE%i",FirstGoodMarker);
  TH1F *QCDHMBE=(TH1F*)gDirectory->Get(A.c_str());
  QCDHMBE->Scale(Lumi*XS[FirstGoodMarker]*0.5/ProcessedEvents[FirstGoodMarker]);
  //HM LC histo recollection
  A=Form("HiggsMassLC%i",FirstGoodMarker);
  TH1F *QCDHMLC=(TH1F*)gDirectory->Get(A.c_str());
  QCDHMLC->Scale(Lumi*XS[FirstGoodMarker]*0.5/ProcessedEvents[FirstGoodMarker]);
  //HM LC over BE histo recollection
  A=Form("HiggsMassLCoverBE%i",FirstGoodMarker);
  TH1F *QCDHMLCoverBE=(TH1F*)gDirectory->Get(A.c_str());
  QCDHMLCoverBE->Scale(Lumi*XS[FirstGoodMarker]*0.5/ProcessedEvents[FirstGoodMarker]);

  //TprimeMass LC_QCDBE histo recollection
  A=Form("FJLCQCDBE%i",FirstGoodMarker);
  TH1F *QCDFJLCQCDBE=(TH1F*)gDirectory->Get(A.c_str());
  QCDFJLCQCDBE->Scale(Lumi*XS[FirstGoodMarker]*0.5/ProcessedEvents[FirstGoodMarker]);
  //TprimeMass LC_QCDLC histo recollection
  A=Form("FJLCQCDLC%i",FirstGoodMarker);
  TH1F *QCDFJLCQCDLC=(TH1F*)gDirectory->Get(A.c_str());
  QCDFJLCQCDLC->Scale(Lumi*XS[FirstGoodMarker]*0.5/ProcessedEvents[FirstGoodMarker]);
  //TprimeMass LC_QCDLCoverBE histo recollection
  A=Form("FJLCQCDLCoverBE%i",FirstGoodMarker);
  TH1F *QCDFJLCQCDLCoverBE=(TH1F*)gDirectory->Get(A.c_str());
  QCDFJLCQCDLCoverBE->Scale(Lumi*XS[FirstGoodMarker]*0.5/ProcessedEvents[FirstGoodMarker]);

  //Top mass after ttbar estimation cut inside the Higgs mass window histo recollection
  A=Form("TMass_LC_HMW%i",FirstGoodMarker);
  TH1F *QCDTMass_LC_HMW=(TH1F*)gDirectory->Get(A.c_str());
  QCDTMass_LC_HMW->Scale(Lumi*XS[FirstGoodMarker]*0.5/ProcessedEvents[FirstGoodMarker]);
  
  //
  A=Form("WM_TM_%i",FirstGoodMarker);
  TH2F *QCD_WM_TM=(TH2F*)gDirectory->Get(A.c_str());
  QCD_WM_TM->Scale(Lumi*XS[FirstGoodMarker]*0.5/ProcessedEvents[FirstGoodMarker]);
  A=Form("HT_TM_%i",FirstGoodMarker);
  TH2F *QCD_HT_TM=(TH2F*)gDirectory->Get(A.c_str());
  QCD_HT_TM->Scale(Lumi*XS[FirstGoodMarker]*0.5/ProcessedEvents[FirstGoodMarker]);
  A=Form("DPWJ_TM_%i",FirstGoodMarker);
  TH2F *QCD_DPWJ_TM=(TH2F*)gDirectory->Get(A.c_str());
  QCD_DPWJ_TM->Scale(Lumi*XS[FirstGoodMarker]*0.5/ProcessedEvents[FirstGoodMarker]);
  A=Form("DPWJ_HT_%i",FirstGoodMarker);
  TH2F *QCD_DPWJ_HT=(TH2F*)gDirectory->Get(A.c_str());
  QCD_DPWJ_HT->Scale(Lumi*XS[FirstGoodMarker]*0.5/ProcessedEvents[FirstGoodMarker]);

  //
  A=Form("TprimeMassLC_LRHT%i",FirstGoodMarker);
  TH1F *QCDFJLCLRTH=(TH1F*)gDirectory->Get(A.c_str());
  QCDFJLCLRTH->Scale(Lumi*XS[FirstGoodMarker]*0.5/ProcessedEvents[FirstGoodMarker]);
  A=Form("TprimeMassLC_HRHT%i",FirstGoodMarker);
  TH1F *QCDFJLCHRTH=(TH1F*)gDirectory->Get(A.c_str());
  QCDFJLCHRTH->Scale(Lumi*XS[FirstGoodMarker]*0.5/ProcessedEvents[FirstGoodMarker]);
  A=Form("TprimeMassLC_LRM%i",FirstGoodMarker);
  TH1F *QCDFJLCLRM=(TH1F*)gDirectory->Get(A.c_str());
  QCDFJLCLRM->Scale(Lumi*XS[FirstGoodMarker]*0.5/ProcessedEvents[FirstGoodMarker]);
  A=Form("TprimeMassLC_HRM%i",FirstGoodMarker);
  TH1F *QCDFJLCHRM=(TH1F*)gDirectory->Get(A.c_str());
  QCDFJLCHRM->Scale(Lumi*XS[FirstGoodMarker]*0.5/ProcessedEvents[FirstGoodMarker]);
  A=Form("TprimeMassLC_WRM%i",FirstGoodMarker);
  TH1F *QCDFJLCWRM=(TH1F*)gDirectory->Get(A.c_str());
  QCDFJLCWRM->Scale(Lumi*XS[FirstGoodMarker]*0.5/ProcessedEvents[FirstGoodMarker]);
  
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
      if (k!=FirstGoodMarker)
	{
      FiveJetsMass[k]->Scale(Lumi*XS[k]*0.5/ProcessedEvents[k]);
      FiveJetsMass3L[k]->Scale(Lumi*XS[k]*0.5/ProcessedEvents[k]);
      //FiveJetsMass3T[k]->Scale(Lumi*XS[k]*0.5/ProcessedEvents[k]);
      FiveJetsMass3L_L[k]->Scale(Lumi*XS[k]*0.5/ProcessedEvents[k]);
      //FiveJetsMass3T_L[k]->Scale(Lumi*XS[k]*0.5/ProcessedEvents[k]);
      LeadingJetPT[k]->Scale(Lumi*XS[k]*0.5/ProcessedEvents[k]);
      Leading2JetPT[k]->Scale(Lumi*XS[k]*0.5/ProcessedEvents[k]);
      Leading3JetPT[k]->Scale(Lumi*XS[k]*0.5/ProcessedEvents[k]);
      Leading4JetPT[k]->Scale(Lumi*XS[k]*0.5/ProcessedEvents[k]);
      Leading5JetPT[k]->Scale(Lumi*XS[k]*0.5/ProcessedEvents[k]);
      Leading6JetPT[k]->Scale(Lumi*XS[k]*0.5/ProcessedEvents[k]);
      THT[k]->Scale(Lumi*XS[k]*0.5/ProcessedEvents[k]);
      DRHjets[k]->Scale(Lumi*XS[k]*0.5/ProcessedEvents[k]);
      DRWjets[k]->Scale(Lumi*XS[k]*0.5/ProcessedEvents[k]);
      Hpt[k]->Scale(Lumi*XS[k]*0.5/ProcessedEvents[k]);
      Tpt[k]->Scale(Lumi*XS[k]*0.5/ProcessedEvents[k]);
      DRWH[k]->Scale(Lumi*XS[k]*0.5/ProcessedEvents[k]);
      DPHjets[k]->Scale(Lumi*XS[k]*0.5/ProcessedEvents[k]);
      DPWjets[k]->Scale(Lumi*XS[k]*0.5/ProcessedEvents[k]);
      DPTjets[k]->Scale(Lumi*XS[k]*0.5/ProcessedEvents[k]);
      HiggsMass[k]->Scale(Lumi*XS[k]*0.5/ProcessedEvents[k]);
      RelHT[k]->Scale(Lumi*XS[k]*0.5/ProcessedEvents[k]);
      DRTH[k]->Scale(Lumi*XS[k]*0.5/ProcessedEvents[k]);
      PtNormalizedMass[k]->Scale(Lumi*XS[k]*0.5/ProcessedEvents[k]);
      RelativeMass[k]->Scale(Lumi*XS[k]*0.5/ProcessedEvents[k]);
      MotherPtNormalizedMass[k]->Scale(Lumi*XS[k]*0.5/ProcessedEvents[k]);
      NumberOfTops[k]->Scale(Lumi*XS[k]*0.5/ProcessedEvents[k]);
      HiggsMassOverTopMass[k]->Scale(Lumi*XS[k]*0.5/ProcessedEvents[k]);
      HiggsTopAsymmetry[k]->Scale(Lumi*XS[k]*0.5/ProcessedEvents[k]);
      ThirdLooseBtag[k]->Scale(Lumi*XS[k]*0.5/ProcessedEvents[k]);
      TopMass[k]->Scale(Lumi*XS[k]*0.5/ProcessedEvents[k]);
      Chi2[k]->Scale(Lumi*XS[k]*0.5/ProcessedEvents[k]);
      UQuarkContent[k]->Scale(Lumi*XS[k]*0.5/ProcessedEvents[k]);
      DQuarkContent[k]->Scale(Lumi*XS[k]*0.5/ProcessedEvents[k]);
      SQuarkContent[k]->Scale(Lumi*XS[k]*0.5/ProcessedEvents[k]);
      CQuarkContent[k]->Scale(Lumi*XS[k]*0.5/ProcessedEvents[k]);
      BQuarkContent[k]->Scale(Lumi*XS[k]*0.5/ProcessedEvents[k]);
      CSVLB[k]->Scale(Lumi*XS[k]*0.5/ProcessedEvents[k]);
      CSVMB[k]->Scale(Lumi*XS[k]*0.5/ProcessedEvents[k]);
      CSVTB[k]->Scale(Lumi*XS[k]*0.5/ProcessedEvents[k]);  
      Vtcs[k]->Scale(Lumi*XS[k]*0.5/ProcessedEvents[k]); 
      /*CLBT[k]->Scale(Lumi*XS[k]*0.5/ProcessedEvents[k]); 
      CMBT[k]->Scale(Lumi*XS[k]*0.5/ProcessedEvents[k]); 
      CTBT[k]->Scale(Lumi*XS[k]*0.5/ProcessedEvents[k]);  
      FLLBT[k]->Scale(Lumi*XS[k]*0.5/ProcessedEvents[k]); 
      FLMBT[k]->Scale(Lumi*XS[k]*0.5/ProcessedEvents[k]); 
      FLTBT[k]->Scale(Lumi*XS[k]*0.5/ProcessedEvents[k]); 
      FCLBT[k]->Scale(Lumi*XS[k]*0.5/ProcessedEvents[k]); 
      FCMBT[k]->Scale(Lumi*XS[k]*0.5/ProcessedEvents[k]); 
      FCTBT[k]->Scale(Lumi*XS[k]*0.5/ProcessedEvents[k]);*/ 
      for (int i=0; i<6; i++) HptTpt[i][k]->Scale(Lumi*XS[k]*0.5/ProcessedEvents[k]);
      HTCSVM[k]->Scale(Lumi*XS[k]*0.5/ProcessedEvents[k]);
      RelHTNTops[k]->Scale(Lumi*XS[k]*0.5/ProcessedEvents[k]);
      for (int i=0; i<6; i++) HiggsMassReversedHptTpt[i][k]->Scale(Lumi*XS[k]*0.5/ProcessedEvents[k]);
      for (int i=0; i<6; i++) TprimeMassReversedHptTpt[i][k]->Scale(Lumi*XS[k]*0.5/ProcessedEvents[k]);
      for (int i=0; i<6; i++) HiggsMassCuts[i][k]->Scale(Lumi*XS[k]*0.5/ProcessedEvents[k]);
      for (int i=0; i<6; i++) WMassFromHiggs[i][k]->Scale(Lumi*XS[k]*0.5/ProcessedEvents[k]);
      for (int i=0; i<6; i++) WMassFromHiggsChi2[i][k]->Scale(Lumi*XS[k]*0.5/ProcessedEvents[k]);

      FiveJetsMassBE[k]->Scale(Lumi*XS[k]*0.5/ProcessedEvents[k]); //FiveJetsMassBE[k]->Scale((TopMassFromHiggs[k]->Integral(TopMassFromHiggs[k]->GetXaxis()->FindBin(0.0),TopMassFromHiggs[k]->GetXaxis()->FindBin(140.0))+TopMassFromHiggs[k]->Integral(TopMassFromHiggs[k]->GetXaxis()->FindBin(230.0),TopMassFromHiggs[k]->GetXaxis()->FindBin(10000.0)))/TopMassFromHiggs[k]->Integral(TopMassFromHiggs[k]->GetXaxis()->FindBin(140.0),TopMassFromHiggs[k]->GetXaxis()->FindBin(230.0)));
      FiveJetsMassLC[k]->Scale(Lumi*XS[k]*0.5/ProcessedEvents[k]);
      TopMassFromHiggs[k]->Scale(Lumi*XS[k]*0.5/ProcessedEvents[k]);
      FiveJetsMassLCoverBE[k]->Scale(Lumi*XS[k]*0.5/ProcessedEvents[k]); //FiveJetsMassLCoverBE[k]->Divide(FiveJetsMassBE[k]);

      HMBE[k]->Scale(Lumi*XS[k]*0.5/ProcessedEvents[k]);
      HMLC[k]->Scale(Lumi*XS[k]*0.5/ProcessedEvents[k]);
      HMLCoverBE[k]->Scale(Lumi*XS[k]*0.5/ProcessedEvents[k]); //HMLCoverBE[k]->Divide(HMBE[k]);
      FJLC_QCDBE[k]->Scale(Lumi*XS[k]*0.5/ProcessedEvents[k]);
      FJLC_QCDLC[k]->Scale(Lumi*XS[k]*0.5/ProcessedEvents[k]);
      FJLC_QCDLCoverBE[k]->Scale(Lumi*XS[k]*0.5/ProcessedEvents[k]); //FJLC_QCDLCoverBE[k]->Divide(FJLC_QCDBE[k]);
      TM_LC_HMW[k]->Scale(Lumi*XS[k]*0.5/ProcessedEvents[k]);

      //
      WM_TM[k]->Scale(Lumi*XS[k]*0.5/ProcessedEvents[k]);
      HT_TM[k]->Scale(Lumi*XS[k]*0.5/ProcessedEvents[k]);
      DPWJ_TM[k]->Scale(Lumi*XS[k]*0.5/ProcessedEvents[k]);
      DPWJ_HT[k]->Scale(Lumi*XS[k]*0.5/ProcessedEvents[k]);
      FiveJetsMassLC_LRHT[k]->Scale(Lumi*XS[k]*0.5/ProcessedEvents[k]);
      FiveJetsMassLC_HRHT[k]->Scale(Lumi*XS[k]*0.5/ProcessedEvents[k]);
      FiveJetsMassLC_LRM[k]->Scale(Lumi*XS[k]*0.5/ProcessedEvents[k]);
      FiveJetsMassLC_HRM[k]->Scale(Lumi*XS[k]*0.5/ProcessedEvents[k]);
      FiveJetsMassLC_WRM[k]->Scale(Lumi*XS[k]*0.5/ProcessedEvents[k]);
      JET1ETA[k]->Scale(Lumi*XS[k]*0.5/ProcessedEvents[k]);
      JET2ETA[k]->Scale(Lumi*XS[k]*0.5/ProcessedEvents[k]);
      JET3ETA[k]->Scale(Lumi*XS[k]*0.5/ProcessedEvents[k]);
      JET4ETA[k]->Scale(Lumi*XS[k]*0.5/ProcessedEvents[k]);
      JET5ETA[k]->Scale(Lumi*XS[k]*0.5/ProcessedEvents[k]);
      JET6ETA[k]->Scale(Lumi*XS[k]*0.5/ProcessedEvents[k]);
      JET1PHI[k]->Scale(Lumi*XS[k]*0.5/ProcessedEvents[k]);
      JET2PHI[k]->Scale(Lumi*XS[k]*0.5/ProcessedEvents[k]);
      JET3PHI[k]->Scale(Lumi*XS[k]*0.5/ProcessedEvents[k]);
      JET4PHI[k]->Scale(Lumi*XS[k]*0.5/ProcessedEvents[k]);
      JET5PHI[k]->Scale(Lumi*XS[k]*0.5/ProcessedEvents[k]);
      JET6PHI[k]->Scale(Lumi*XS[k]*0.5/ProcessedEvents[k]);
      JETMULTI[k]->Scale(Lumi*XS[k]*0.5/ProcessedEvents[k]);

	}
      //Settings for QCD
      if (k>=ParcialTestMin && k<=ParcialTestMax)
	{
	  FiveJetsMass[k]->SetFillColor(kViolet); 
	  FiveJetsMass3L[k]->SetFillColor(kViolet);
	  //FiveJetsMass3T[k]->SetFillColor(kViolet);
	  FiveJetsMass3L_L[k]->SetFillColor(kViolet);
	  //FiveJetsMass3T_L[k]->SetFillColor(kViolet);
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
	  UQuarkContent[k]->SetFillColor(kViolet);
	  DQuarkContent[k]->SetFillColor(kViolet);
	  SQuarkContent[k]->SetFillColor(kViolet);
	  CQuarkContent[k]->SetFillColor(kViolet);
	  BQuarkContent[k]->SetFillColor(kViolet);
	  CSVLB[k]->SetFillColor(kViolet);
	  CSVMB[k]->SetFillColor(kViolet);
	  CSVTB[k]->SetFillColor(kViolet);
	  Vtcs[k]->SetFillColor(kViolet);	  
	  /*CLBT[k]->SetFillColor(kViolet); 
	  CMBT[k]->SetFillColor(kViolet); 
	  CTBT[k]->SetFillColor(kViolet);  
	  FLLBT[k]->SetFillColor(kViolet); 
	  FLMBT[k]->SetFillColor(kViolet); 
	  FLTBT[k]->SetFillColor(kViolet); 
	  FCLBT[k]->SetFillColor(kViolet); 
	  FCMBT[k]->SetFillColor(kViolet); 
	  FCTBT[k]->SetFillColor(kViolet); */
	  for (int i=0; i<6; i++) HiggsMassReversedHptTpt[i][k]->SetFillColor(kViolet);
	  for (int i=0; i<6; i++) TprimeMassReversedHptTpt[i][k]->SetFillColor(kViolet);
	  for (int i=0; i<6; i++) HiggsMassCuts[i][k]->SetFillColor(kViolet);
	  for (int i=0; i<6; i++) WMassFromHiggs[i][k]->SetFillColor(kViolet);	  
	  for (int i=0; i<6; i++) WMassFromHiggsChi2[i][k]->SetFillColor(kViolet);
	  FiveJetsMassBE[k]->SetFillColor(kViolet);
	  FiveJetsMassLC[k]->SetFillColor(kViolet);
	  TopMassFromHiggs[k]->SetFillColor(kViolet);
	  FiveJetsMassLCoverBE[k]->SetFillColor(kViolet);
	  HMBE[k]->SetFillColor(kViolet);
	  HMLC[k]->SetFillColor(kViolet);
	  HMLCoverBE[k]->SetFillColor(kViolet);
	  FJLC_QCDBE[k]->SetFillColor(kViolet);
	  FJLC_QCDLC[k]->SetFillColor(kViolet);
	  FJLC_QCDLCoverBE[k]->SetFillColor(kViolet);
	  TM_LC_HMW[k]->SetFillColor(kViolet);
	  FiveJetsMassLC_LRHT[k]->SetFillColor(kViolet);
	  FiveJetsMassLC_HRHT[k]->SetFillColor(kViolet);
	  FiveJetsMassLC_LRM[k]->SetFillColor(kViolet);
	  FiveJetsMassLC_HRM[k]->SetFillColor(kViolet);
	  FiveJetsMassLC_WRM[k]->SetFillColor(kViolet);
	  JET1ETA[k]->SetFillColor(kViolet);
	  JET2ETA[k]->SetFillColor(kViolet);
	  JET3ETA[k]->SetFillColor(kViolet);
	  JET4ETA[k]->SetFillColor(kViolet);
	  JET5ETA[k]->SetFillColor(kViolet);
	  JET6ETA[k]->SetFillColor(kViolet);
	  JET1PHI[k]->SetFillColor(kViolet);
	  JET2PHI[k]->SetFillColor(kViolet);
	  JET3PHI[k]->SetFillColor(kViolet);
	  JET4PHI[k]->SetFillColor(kViolet);
	  JET5PHI[k]->SetFillColor(kViolet);
	  JET6PHI[k]->SetFillColor(kViolet);
	  JETMULTI[k]->SetFillColor(kViolet);
	  if (k!=ParcialTestMin) 
	    {
	      QCDTprimeMass->Add(FiveJetsMass[k]);
	      QCDTprimeMass3L->Add(FiveJetsMass3L[k]);
	      //QCDTprimeMass3T->Add(FiveJetsMass3T[k]);
	      QCDTprimeMass3L_L->Add(FiveJetsMass3L_L[k]);
	      //QCDTprimeMass3T_L->Add(FiveJetsMass3T_L[k]);
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
	      QCDTprimeMassUQC->Add(UQuarkContent[k]);
	      QCDTprimeMassDQC->Add(DQuarkContent[k]);
	      QCDTprimeMassSQC->Add(SQuarkContent[k]);
	      QCDTprimeMassCQC->Add(CQuarkContent[k]);
	      QCDTprimeMassBQC->Add(BQuarkContent[k]);
	      QCDTprimeMassCSVL->Add(CSVLB[k]);
	      QCDTprimeMassCSVM->Add(CSVMB[k]);
	      QCDTprimeMassCSVT->Add(CSVTB[k]);
	      QCDTprimeMassVTX->Add(Vtcs[k]);
	      /*QCDTprimeMassCLBT->Add(CLBT[k]);
	      QCDTprimeMassCMBT->Add(CMBT[k]);
	      QCDTprimeMassCTBT->Add(CTBT[k]);
	      QCDTprimeMassFLLBT->Add(FLLBT[k]);
	      QCDTprimeMassFLMBT->Add(FLMBT[k]);
	      QCDTprimeMassFLTBT->Add(FLTBT[k]);
	      QCDTprimeMassFCLBT->Add(FCLBT[k]);
	      QCDTprimeMassFCMBT->Add(FCMBT[k]);
	      QCDTprimeMassFCTBT->Add(FCTBT[k]);*/
	      for (int i=0; i<6; i++) {QCD_HptTpt[i]->Add(HptTpt[i][k]); cout << "Adding HptTpt " << i << k << ", with integral: " << QCD_HptTpt[i]->Integral() << endl;}
	      QCD_HTCSVM->Add(HTCSVM[k]);
	      QCD_RHTNT->Add(RelHTNTops[k]);
	      for (int i=0; i<6; i++) QCDTprimeMassHMBE[i]->Add(HiggsMassReversedHptTpt[i][k]);
	      for (int i=0; i<6; i++) QCDTprimeMassBE[i]->Add(TprimeMassReversedHptTpt[i][k]);
	      for (int i=0; i<6; i++) QCDTprimeMassHMCuts[i]->Add(HiggsMassCuts[i][k]);
	      for (int i=0; i<6; i++) QCDTprimeMassWMFHCuts[i]->Add(WMassFromHiggs[i][k]);
	      for (int i=0; i<6; i++) QCDTprimeMassWMFHCCuts[i]->Add(WMassFromHiggsChi2[i][k]);
	      
	      QCDTprimeMassBETop->Add(FiveJetsMassBE[k]);
	      QCDTprimeMassLC->Add(FiveJetsMassLC[k]);
	      QCDTopMassFromHiggs->Add(TopMassFromHiggs[k]);
	      QCDTprimeMassLCoverBE->Add(FiveJetsMassLCoverBE[k]);

	      QCDHMBE->Add(HMBE[k]);
	      QCDHMLC->Add(HMLC[k]);
	      QCDHMLCoverBE->Add(HMLCoverBE[k]);

	      QCDFJLCQCDBE->Add(FJLC_QCDBE[k]);
	      QCDFJLCQCDLC->Add(FJLC_QCDLC[k]);
	      QCDFJLCQCDLCoverBE->Add(FJLC_QCDLCoverBE[k]);
	      QCDTMass_LC_HMW->Add(TM_LC_HMW[k]);
	      
	      QCD_WM_TM->Add(WM_TM[k]);
	      QCD_HT_TM->Add(HT_TM[k]);
	      QCD_DPWJ_TM->Add(DPWJ_TM[k]);
	      QCD_DPWJ_HT->Add(DPWJ_HT[k]);

	      QCDFJLCLRTH->Add(FiveJetsMassLC_LRHT[k]);
	      QCDFJLCHRTH->Add(FiveJetsMassLC_HRHT[k]);
	      QCDFJLCLRM->Add(FiveJetsMassLC_LRM[k]);
	      QCDFJLCHRM->Add(FiveJetsMassLC_HRM[k]);
	      QCDFJLCWRM->Add(FiveJetsMassLC_WRM[k]);
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
	      TFile f("QCD.root", "RECREATE"); 
	      QCDTprimeMass->Write();
	      QCDTprimeMass3L->Write();
	      //QCDTprimeMass3T->Write();
	      QCDTprimeMass3L_L->Write();
	      //QCDTprimeMass3T_L->Write();
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
	      QCDTprimeMassUQC->Write();
	      QCDTprimeMassDQC->Write();
	      QCDTprimeMassSQC->Write();
	      QCDTprimeMassCQC->Write();
	      QCDTprimeMassBQC->Write();
	      QCDTprimeMassCSVL->Write();
	      QCDTprimeMassCSVM->Write();
	      QCDTprimeMassCSVT->Write();
	      QCDTprimeMassVTX->Write();	      
	      /*QCDTprimeMassCLBT->Write();
	      QCDTprimeMassCMBT->Write();
	      QCDTprimeMassCTBT->Write();
	      QCDTprimeMassFLLBT->Write();
	      QCDTprimeMassFLMBT->Write();
	      QCDTprimeMassFLTBT->Write();
	      QCDTprimeMassFCLBT->Write();
	      QCDTprimeMassFCMBT->Write();
	      QCDTprimeMassFCTBT->Write();*/
	      for (int i=0; i<6; i++) QCD_HptTpt[i]->Write();
	      QCD_HTCSVM->Write();
	      QCD_RHTNT->Write();
	      for (int i=0; i<6; i++) QCDTprimeMassHMBE[i]->Write();
	      for (int i=0; i<6; i++) QCDTprimeMassBE[i]->Write();
	      for (int i=0; i<6; i++) QCDTprimeMassHMCuts[i]->Write();
	      for (int i=0; i<6; i++) QCDTprimeMassWMFHCuts[i]->Write();
	      for (int i=0; i<6; i++) QCDTprimeMassWMFHCCuts[i]->Write();
	      QCDTprimeMassBETop->Write();
	      QCDTprimeMassLC->Write();
	      QCDTopMassFromHiggs->Write();
	      QCDTprimeMassLCoverBE->Scale(QCDTprimeMassBETop->Integral()/QCDTprimeMassLCoverBE->Integral()); QCDTprimeMassLCoverBE->Divide(QCDTprimeMassBETop); QCDTprimeMassLCoverBE->Write();

	      QCDHMBE->Write();
	      QCDHMLC->Write();
	      QCDHMLCoverBE->Scale(QCDHMBE->Integral()/QCDHMLCoverBE->Integral()); QCDHMLCoverBE->Divide(QCDHMBE); QCDHMLCoverBE->Write();

	      QCDFJLCQCDBE->Write();
	      QCDFJLCQCDLC->Write();
	      QCDFJLCQCDLCoverBE->Scale(QCDFJLCQCDBE->Integral()/QCDFJLCQCDLCoverBE->Integral()); QCDFJLCQCDLCoverBE->Divide(QCDFJLCQCDBE); QCDFJLCQCDLCoverBE->Write();
	      
	      QCDTMass_LC_HMW->Write();

	      QCD_WM_TM->Write();
	      QCD_HT_TM->Write();
	      QCD_DPWJ_TM->Write();
	      QCD_DPWJ_HT->Write();

	      QCDFJLCLRTH->Write();
	      QCDFJLCHRTH->Write();
	      QCDFJLCLRM->Write();
	      QCDFJLCHRM->Write();
	      QCDFJLCWRM->Write();

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
	      
	      f.Close();
	    }
	}
    }
  
  double RegionA=0; double RegionB=0; double RegionC=0; double RegionD=0;
  TFile f("QCD_ABCD.root", "RECREATE");
  for (int i=0; i<6; i++) 
    {
  RegionA=QCD_HptTpt[i]->Integral(QCD_HptTpt[i]->GetXaxis()->FindBin(200.0),QCD_HptTpt[i]->GetNbinsX(),QCD_HptTpt[i]->GetYaxis()->FindBin(200.0),QCD_HptTpt[i]->GetNbinsY());
  RegionB=QCD_HptTpt[i]->Integral(QCD_HptTpt[i]->GetXaxis()->FindBin(200.0),QCD_HptTpt[i]->GetNbinsX(),QCD_HptTpt[i]->GetYaxis()->FindBin(0.0),QCD_HptTpt[i]->GetYaxis()->FindBin(200.0));
  RegionC=QCD_HptTpt[i]->Integral(QCD_HptTpt[i]->GetXaxis()->FindBin(0.0),QCD_HptTpt[i]->GetXaxis()->FindBin(200.0),QCD_HptTpt[i]->GetYaxis()->FindBin(200.0),QCD_HptTpt[i]->GetNbinsY());
  RegionD=QCD_HptTpt[i]->Integral(QCD_HptTpt[i]->GetXaxis()->FindBin(0.0),QCD_HptTpt[i]->GetXaxis()->FindBin(200.0),QCD_HptTpt[i]->GetYaxis()->FindBin(0.0),QCD_HptTpt[i]->GetYaxis()->FindBin(200.0));

  /*for (int j=0; j<6; j++)*/ cout << "Correlation Factor Hpt Tpt: " << QCD_HptTpt[i]->GetCorrelationFactor() << endl;

  cout << "ABCD Method Info for Hpt Top pt:" << endl;
  cout << "Number of events on region A (signal enriched) " <<  RegionA << endl;
  cout << "Number of events on region B " <<  RegionB << endl;
  cout << "Number of events on region C " <<  RegionC << endl;
  cout << "Number of events on region D " <<  RegionD << endl;
    
  //for (int i=0; i<6; i++) 
  //  {
      QCDTprimeMassHMBE[i]->Scale((RegionB*(RegionC/RegionD))/QCDTprimeMassHMBE[i]->Integral());
      QCDTprimeMassBE[i]->Scale((RegionB*(RegionC/RegionD))/QCDTprimeMassBE[i]->Integral());
      QCDTprimeMassHMBE[i]->Write();
      QCDTprimeMassBE[i]->Write();
    }
  f.Close();

  double Region2A=0; double Region2B=0; double Region2C=0; double Region2D=0;
  Region2A=QCD_HTCSVM->Integral(QCD_HTCSVM->GetXaxis()->FindBin(630.0),QCD_HTCSVM->GetNbinsX(),QCD_HTCSVM->GetYaxis()->FindBin(3.0),QCD_HTCSVM->GetNbinsY());
  Region2B=QCD_HTCSVM->Integral(QCD_HTCSVM->GetXaxis()->FindBin(630.0),QCD_HTCSVM->GetNbinsX(),QCD_HTCSVM->GetYaxis()->FindBin(0.0),QCD_HTCSVM->GetYaxis()->FindBin(3.0));
  Region2C=QCD_HTCSVM->Integral(QCD_HTCSVM->GetXaxis()->FindBin(0.0),QCD_HTCSVM->GetXaxis()->FindBin(630.0),QCD_HTCSVM->GetYaxis()->FindBin(3.0),QCD_HTCSVM->GetNbinsY());
  Region2D=QCD_HTCSVM->Integral(QCD_HTCSVM->GetXaxis()->FindBin(0.0),QCD_HTCSVM->GetXaxis()->FindBin(630.0),QCD_HTCSVM->GetYaxis()->FindBin(0.0),QCD_HTCSVM->GetYaxis()->FindBin(3.0));

  cout << "ABCD Method Info for HT and CSVM b-tags:" << endl;
  cout << "Number of events on region A (signal enriched) " <<  Region2A << endl;
  cout << "Number of events on region B " <<  Region2B << endl;
  cout << "Number of events on region C " <<  Region2C << endl;
  cout << "Number of events on region D " <<  Region2D << endl;  

  cout << "Correlation Factor RelHTNtops: " << QCD_RHTNT->GetCorrelationFactor() << endl;

  double Region3A=0; double Region3B=0; double Region3C=0; double Region3D=0;
  Region3A=QCD_RHTNT->Integral(QCD_RHTNT->GetXaxis()->FindBin(.65),QCD_RHTNT->GetNbinsX(),QCD_RHTNT->GetYaxis()->FindBin(0.0),QCD_RHTNT->GetYaxis()->FindBin(2.0));
  Region3B=QCD_RHTNT->Integral(QCD_RHTNT->GetXaxis()->FindBin(.65),QCD_RHTNT->GetNbinsX(),QCD_RHTNT->GetYaxis()->FindBin(2.0),QCD_RHTNT->GetNbinsY());
  Region3C=QCD_RHTNT->Integral(QCD_RHTNT->GetXaxis()->FindBin(0.0),QCD_RHTNT->GetXaxis()->FindBin(.65),QCD_RHTNT->GetYaxis()->FindBin(0.0),QCD_RHTNT->GetYaxis()->FindBin(2.0));
  Region3D=QCD_RHTNT->Integral(QCD_RHTNT->GetXaxis()->FindBin(0.0),QCD_RHTNT->GetXaxis()->FindBin(.65),QCD_RHTNT->GetYaxis()->FindBin(2.0),QCD_RHTNT->GetNbinsY());

  cout << "ABCD Method Info for RelHT and Ntops:" << endl;
  cout << "Number of events on region A (signal enriched) " <<  Region3A << endl;
  cout << "Number of events on region B " <<  Region3B << endl;
  cout << "Number of events on region C " <<  Region3C << endl;
  cout << "Number of events on region D " <<  Region3D << endl;

  cout << "Correlation Factor WM_TM: " << QCD_WM_TM->GetCorrelationFactor() << endl;

  double Region4A=0; double Region4B[2]={0,0}; double Region4C[4]={0,0,0,0}; double Region4D[2]={0,0};
  Region4A=QCD_WM_TM->Integral(QCD_WM_TM->GetXaxis()->FindBin(60.0),QCD_WM_TM->GetXaxis()->FindBin(120.0),QCD_WM_TM->GetYaxis()->FindBin(140.0),QCD_WM_TM->GetYaxis()->FindBin(230.0));

  Region4B[0]=QCD_WM_TM->Integral(QCD_WM_TM->GetXaxis()->FindBin(60.0),QCD_WM_TM->GetXaxis()->FindBin(120.0),QCD_WM_TM->GetYaxis()->FindBin(0.0),QCD_WM_TM->GetYaxis()->FindBin(140.0));
  Region4B[1]=QCD_WM_TM->Integral(QCD_WM_TM->GetXaxis()->FindBin(60.0),QCD_WM_TM->GetXaxis()->FindBin(120.0),QCD_WM_TM->GetYaxis()->FindBin(230.0),QCD_WM_TM->GetYaxis()->FindBin(1000.0));

  Region4C[0]=QCD_WM_TM->Integral(QCD_WM_TM->GetXaxis()->FindBin(120.0),QCD_WM_TM->GetXaxis()->FindBin(1000.0),QCD_WM_TM->GetYaxis()->FindBin(230.0),QCD_WM_TM->GetYaxis()->FindBin(1000.0));
  Region4C[1]=QCD_WM_TM->Integral(QCD_WM_TM->GetXaxis()->FindBin(0.0),QCD_WM_TM->GetXaxis()->FindBin(60.0),QCD_WM_TM->GetYaxis()->FindBin(0.0),QCD_WM_TM->GetYaxis()->FindBin(140.0));
  Region4C[2]=QCD_WM_TM->Integral(QCD_WM_TM->GetXaxis()->FindBin(120.0),QCD_WM_TM->GetXaxis()->FindBin(1000.0),QCD_WM_TM->GetYaxis()->FindBin(0.0),QCD_WM_TM->GetYaxis()->FindBin(140.0));
  Region4C[3]=QCD_WM_TM->Integral(QCD_WM_TM->GetXaxis()->FindBin(0.0),QCD_WM_TM->GetXaxis()->FindBin(60.0),QCD_WM_TM->GetYaxis()->FindBin(230.0),QCD_WM_TM->GetYaxis()->FindBin(1000.0));

  Region4D[0]=QCD_WM_TM->Integral(QCD_WM_TM->GetXaxis()->FindBin(0.0),QCD_WM_TM->GetXaxis()->FindBin(60.0),QCD_WM_TM->GetYaxis()->FindBin(140.0),QCD_WM_TM->GetYaxis()->FindBin(230.0));
  Region4D[0]=QCD_WM_TM->Integral(QCD_WM_TM->GetXaxis()->FindBin(120.0),QCD_WM_TM->GetXaxis()->FindBin(1000.0),QCD_WM_TM->GetYaxis()->FindBin(140.0),QCD_WM_TM->GetYaxis()->FindBin(230.0));

  cout << "ABCD Method Info for WM_TM (Veto approach):" << endl;
  cout << "Number of events on region A (signal enriched) " <<  Region4A << endl;
  cout << "Number of events on region B " <<  Region4B[0]+Region4B[1] << endl;
  cout << "Number of events on region C " <<  Region4C[0]+Region4C[1]+Region4C[2]+Region4C[3] << endl;
  cout << "Number of events on region D " <<  Region4D[0]+Region4D[1] << endl;

  exit(0);
  
}
