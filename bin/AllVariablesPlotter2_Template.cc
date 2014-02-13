#include <TH1.h>
#include <TH2F.h>
#include <TLegend.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <THStack.h>
#include <TPostScript.h>
#include <TPDF.h>
#include <TPaveText.h>
#include <TGraph.h>
#include <TMultiGraph.h>
#include <TChain.h>
#include <TBranch.h>
#include <TFile.h>

#include "TMath.h"
#include "TF1.h"

using namespace std;

const TString StorageDirPrefix="file:/afs/cern.ch/work/j/jruizalv/private/Analysis/CMSSW_5_3_9_patch3/src/Extractors/PatExtractor/bin/WithTrigger/SUFFIX/";
const int NOS=4;
//const TString Samples[NOS] = {"Diboson.root", "SingleTop.root", "TTJets.root", "Zjets.root", "Wjets.root", "QCD.root", "Signal.root"};
const TString Samples[NOS] = {"SingleTop.root", "TTJets.root", "QCD.root", "Signal.root"};

void AllVariablesPlotter2()
{
  TFile *CurrentFile[NOS];
  TH1F *TprimeHistos[NOS];
  TH1F *LeadingJetPT[NOS];
  TH1F *Leading2JetPT[NOS];
  TH1F *Leading3JetPT[NOS];
  TH1F *Leading4JetPT[NOS];
  TH1F *Leading5JetPT[NOS];
  TH1F *Leading6JetPT[NOS];
  TH1F *THT[NOS];
  TH1F *DRHjets[NOS];
  TH1F *DRWjets[NOS];
  TH1F *Hpt[NOS];
  TH1F *Tpt[NOS];
  TH1F *DRWH[NOS];
  TH1F *DPHjets[NOS];
  TH1F *DPWjets[NOS];
  TH1F *DPTjets[NOS];
  TH1F *HiggsMass[NOS];
  TH1F *RelHT[NOS];
  TH1F *DRTH[NOS];
  TH1F *PtNormalizedMass[NOS];
  TH1F *RelativeMass[NOS];
  TH1F *MotherPtNormalizedMass[NOS];
  TH1F *NumberOfTops[NOS];
  TH1F *HiggsMassOverTopMass[NOS];
  TH1F *HiggsTopAsymmetry[NOS];
  TH1F *ThirdLooseBtag[NOS];
  TH1F *TopMass[NOS];
  TH1F *Chi2[NOS];
  //TH2F *HptTpt[NOS];
  
  for (int i=0; i<NOS; i++)
    {
      CurrentFile[i] = new TFile(StorageDirPrefix + Samples[i], "READ");
      if ( CurrentFile[i]->IsOpen() ) printf( Samples[i] + " File opened successfully\n");
      if (i==0) 
	{
	  TprimeHistos[0]->SetDefaultSumw2(); TprimeHistos[0]= (TH1F*)gDirectory->Get("TprimeMass4");
	  LeadingJetPT[0]->SetDefaultSumw2(); LeadingJetPT[0]= (TH1F*)gDirectory->Get("jet1_pt4");
	  Leading2JetPT[0]->SetDefaultSumw2(); Leading2JetPT[0]= (TH1F*)gDirectory->Get("jet2_pt4");
	  Leading3JetPT[0]->SetDefaultSumw2(); Leading3JetPT[0]= (TH1F*)gDirectory->Get("jet3_pt4");
	  Leading4JetPT[0]->SetDefaultSumw2(); Leading4JetPT[0]= (TH1F*)gDirectory->Get("jet4_pt4");
	  Leading5JetPT[0]->SetDefaultSumw2(); Leading5JetPT[0]= (TH1F*)gDirectory->Get("jet5_pt4");
	  Leading6JetPT[0]->SetDefaultSumw2(); Leading6JetPT[0]= (TH1F*)gDirectory->Get("jet6_pt4");
	  THT[0]->SetDefaultSumw2(); THT[0]= (TH1F*)gDirectory->Get("THT4");
	  DRHjets[0]->SetDefaultSumw2(); DRHjets[0]= (TH1F*)gDirectory->Get("DeltaR_of_Higgs_Jets4");
	  DRWjets[0]->SetDefaultSumw2(); DRWjets[0]= (TH1F*)gDirectory->Get("DeltaR_of_W_Jets4");
	  Hpt[0]->SetDefaultSumw2(); Hpt[0]= (TH1F*)gDirectory->Get("HPt4");
	  Tpt[0]->SetDefaultSumw2(); Tpt[0]= (TH1F*)gDirectory->Get("TPt4");
	  DRWH[0]->SetDefaultSumw2(); DRWH[0]= (TH1F*)gDirectory->Get("DeltaR_of_W_Higgs4");
	  DPHjets[0]->SetDefaultSumw2(); DPHjets[0]= (TH1F*)gDirectory->Get("DeltaPhi_of_Higgs_jets4");
	  DPWjets[0]->SetDefaultSumw2(); DPWjets[0]= (TH1F*)gDirectory->Get("DeltaPhi_of_W_jets4");
	  DPTjets[0]->SetDefaultSumw2(); DPTjets[0]= (TH1F*)gDirectory->Get("DeltaPhi_of_T_jet4");
	  HiggsMass[0]->SetDefaultSumw2(); HiggsMass[0]= (TH1F*)gDirectory->Get("HM4");
	  RelHT[0]->SetDefaultSumw2(); RelHT[0]= (TH1F*)gDirectory->Get("RelHT4");
	  DRTH[0]->SetDefaultSumw2(); DRTH[0]= (TH1F*)gDirectory->Get("DeltaR_of_Top_Higgs4");
	  PtNormalizedMass[0]->SetDefaultSumw2(); PtNormalizedMass[0]= (TH1F*)gDirectory->Get("PT_Normalized_Mass4");
	  RelativeMass[0]->SetDefaultSumw2(); RelativeMass[0]= (TH1F*)gDirectory->Get("Relative_Mass4");
	  MotherPtNormalizedMass[0]->SetDefaultSumw2(); MotherPtNormalizedMass[0]= (TH1F*)gDirectory->Get("Mother_PT_Normalized_Mass4");
	  NumberOfTops[0]->SetDefaultSumw2();NumberOfTops[0]= (TH1F*)gDirectory->Get("Number_of_Tops4");
	  HiggsMassOverTopMass[0]->SetDefaultSumw2();HiggsMassOverTopMass[0]= (TH1F*)gDirectory->Get("HMoverTM4");
	  HiggsTopAsymmetry[0]->SetDefaultSumw2();HiggsTopAsymmetry[0]= (TH1F*)gDirectory->Get("HTAsym4");
	  ThirdLooseBtag[0]->SetDefaultSumw2();ThirdLooseBtag[0]= (TH1F*)gDirectory->Get("TLBTag4");
	  TopMass[0]->SetDefaultSumw2(); TopMass[0]= (TH1F*)gDirectory->Get("TMass4");
	  Chi2[0]->SetDefaultSumw2(); Chi2[0]= (TH1F*)gDirectory->Get("ChiSq4");
	}
      else if (i==1) 
	{
	  TprimeHistos[1]->SetDefaultSumw2(); TprimeHistos[1]= (TH1F*)gDirectory->Get("TprimeMass9");
	  LeadingJetPT[1]->SetDefaultSumw2(); LeadingJetPT[1]= (TH1F*)gDirectory->Get("jet1_pt9");
	  Leading2JetPT[1]->SetDefaultSumw2(); Leading2JetPT[1]= (TH1F*)gDirectory->Get("jet2_pt9");
	  Leading3JetPT[1]->SetDefaultSumw2(); Leading3JetPT[1]= (TH1F*)gDirectory->Get("jet3_pt9");
	  Leading4JetPT[1]->SetDefaultSumw2(); Leading4JetPT[1]= (TH1F*)gDirectory->Get("jet4_pt9");
	  Leading5JetPT[1]->SetDefaultSumw2(); Leading5JetPT[1]= (TH1F*)gDirectory->Get("jet5_pt9");
	  Leading6JetPT[1]->SetDefaultSumw2(); Leading6JetPT[1]= (TH1F*)gDirectory->Get("jet6_pt9");
	  THT[1]->SetDefaultSumw2(); THT[1]= (TH1F*)gDirectory->Get("THT9");
	  DRHjets[1]->SetDefaultSumw2(); DRHjets[1]= (TH1F*)gDirectory->Get("DeltaR_of_Higgs_Jets9");
	  DRWjets[1]->SetDefaultSumw2(); DRWjets[1]= (TH1F*)gDirectory->Get("DeltaR_of_W_Jets9");
	  Hpt[1]->SetDefaultSumw2(); Hpt[1]= (TH1F*)gDirectory->Get("HPt9");
	  Tpt[1]->SetDefaultSumw2(); Tpt[1]= (TH1F*)gDirectory->Get("TPt9");
	  DRWH[1]->SetDefaultSumw2(); DRWH[1]= (TH1F*)gDirectory->Get("DeltaR_of_W_Higgs9");
	  DPHjets[1]->SetDefaultSumw2(); DPHjets[1]= (TH1F*)gDirectory->Get("DeltaPhi_of_Higgs_jets9");
	  DPWjets[1]->SetDefaultSumw2(); DPWjets[1]= (TH1F*)gDirectory->Get("DeltaPhi_of_W_jets9");
	  DPTjets[1]->SetDefaultSumw2(); DPTjets[1]= (TH1F*)gDirectory->Get("DeltaPhi_of_T_jet9");
	  HiggsMass[1]->SetDefaultSumw2(); HiggsMass[1]= (TH1F*)gDirectory->Get("HM9");
	  RelHT[1]->SetDefaultSumw2(); RelHT[1]= (TH1F*)gDirectory->Get("RelHT9");
	  DRTH[1]->SetDefaultSumw2(); DRTH[1]= (TH1F*)gDirectory->Get("DeltaR_of_Top_Higgs9");
	  PtNormalizedMass[1]->SetDefaultSumw2(); PtNormalizedMass[1]= (TH1F*)gDirectory->Get("PT_Normalized_Mass9");
	  RelativeMass[1]->SetDefaultSumw2(); RelativeMass[1]= (TH1F*)gDirectory->Get("Relative_Mass9");
	  MotherPtNormalizedMass[1]->SetDefaultSumw2(); MotherPtNormalizedMass[1]= (TH1F*)gDirectory->Get("Mother_PT_Normalized_Mass9");
	  NumberOfTops[1]->SetDefaultSumw2();NumberOfTops[1]= (TH1F*)gDirectory->Get("Number_of_Tops9");
	  HiggsMassOverTopMass[1]->SetDefaultSumw2();HiggsMassOverTopMass[1]= (TH1F*)gDirectory->Get("HMoverTM9");
	  HiggsTopAsymmetry[1]->SetDefaultSumw2();HiggsTopAsymmetry[1]= (TH1F*)gDirectory->Get("HTAsym9");
	  ThirdLooseBtag[1]->SetDefaultSumw2();ThirdLooseBtag[1]= (TH1F*)gDirectory->Get("TLBTag9");
	  TopMass[1]->SetDefaultSumw2(); TopMass[1]= (TH1F*)gDirectory->Get("TMass9");
	  Chi2[1]->SetDefaultSumw2(); Chi2[1]= (TH1F*)gDirectory->Get("ChiSq9");
	}
      else if (i==2) 
	{
	  TprimeHistos[2]->SetDefaultSumw2(); TprimeHistos[2]= (TH1F*)gDirectory->Get("TprimeMass10");
	  LeadingJetPT[2]->SetDefaultSumw2(); LeadingJetPT[2]= (TH1F*)gDirectory->Get("jet1_pt10");
	  Leading2JetPT[2]->SetDefaultSumw2(); Leading2JetPT[2]= (TH1F*)gDirectory->Get("jet2_pt10");
	  Leading3JetPT[2]->SetDefaultSumw2(); Leading3JetPT[2]= (TH1F*)gDirectory->Get("jet3_pt10");
	  Leading4JetPT[2]->SetDefaultSumw2(); Leading4JetPT[2]= (TH1F*)gDirectory->Get("jet4_pt10");
	  Leading5JetPT[2]->SetDefaultSumw2(); Leading5JetPT[2]= (TH1F*)gDirectory->Get("jet5_pt10");
	  Leading6JetPT[2]->SetDefaultSumw2(); Leading6JetPT[2]= (TH1F*)gDirectory->Get("jet6_pt10");
	  THT[2]->SetDefaultSumw2(); THT[2]= (TH1F*)gDirectory->Get("THT10");
	  DRHjets[2]->SetDefaultSumw2(); DRHjets[2]= (TH1F*)gDirectory->Get("DeltaR_of_Higgs_Jets10");
	  DRWjets[2]->SetDefaultSumw2(); DRWjets[2]= (TH1F*)gDirectory->Get("DeltaR_of_W_Jets10");
	  Hpt[2]->SetDefaultSumw2(); Hpt[2]= (TH1F*)gDirectory->Get("HPt10");
	  Tpt[2]->SetDefaultSumw2(); Tpt[2]= (TH1F*)gDirectory->Get("TPt10");
	  DRWH[2]->SetDefaultSumw2(); DRWH[2]= (TH1F*)gDirectory->Get("DeltaR_of_W_Higgs10");
	  DPHjets[2]->SetDefaultSumw2(); DPHjets[2]= (TH1F*)gDirectory->Get("DeltaPhi_of_Higgs_jets10");
	  DPWjets[2]->SetDefaultSumw2(); DPWjets[2]= (TH1F*)gDirectory->Get("DeltaPhi_of_W_jets10");
	  DPTjets[2]->SetDefaultSumw2(); DPTjets[2]= (TH1F*)gDirectory->Get("DeltaPhi_of_T_jet10");
	  HiggsMass[2]->SetDefaultSumw2(); HiggsMass[2]= (TH1F*)gDirectory->Get("HM10");
	  RelHT[2]->SetDefaultSumw2(); RelHT[2]= (TH1F*)gDirectory->Get("RelHT10");
	  DRTH[2]->SetDefaultSumw2(); DRTH[2]= (TH1F*)gDirectory->Get("DeltaR_of_Top_Higgs10");
	  PtNormalizedMass[2]->SetDefaultSumw2(); PtNormalizedMass[2]= (TH1F*)gDirectory->Get("PT_Normalized_Mass10");
	  RelativeMass[2]->SetDefaultSumw2(); RelativeMass[2]= (TH1F*)gDirectory->Get("Relative_Mass10");
	  MotherPtNormalizedMass[2]->SetDefaultSumw2(); MotherPtNormalizedMass[2]= (TH1F*)gDirectory->Get("Mother_PT_Normalized_Mass10");
	  NumberOfTops[2]->SetDefaultSumw2();NumberOfTops[2]= (TH1F*)gDirectory->Get("Number_of_Tops10");
	  HiggsMassOverTopMass[2]->SetDefaultSumw2();HiggsMassOverTopMass[2]= (TH1F*)gDirectory->Get("HMoverTM10");
	  HiggsTopAsymmetry[2]->SetDefaultSumw2();HiggsTopAsymmetry[2]= (TH1F*)gDirectory->Get("HTAsym10");
	  ThirdLooseBtag[2]->SetDefaultSumw2();ThirdLooseBtag[2]= (TH1F*)gDirectory->Get("TLBTag10");
	  TopMass[2]->SetDefaultSumw2(); TopMass[2]= (TH1F*)gDirectory->Get("TMass10");
	  Chi2[2]->SetDefaultSumw2(); Chi2[2]= (TH1F*)gDirectory->Get("ChiSq10");
	}
      else if (i==3) 
	{
	  TprimeHistos[3]->SetDefaultSumw2(); TprimeHistos[3]= (TH1F*)gDirectory->Get("TprimeMass15");
	  LeadingJetPT[3]->SetDefaultSumw2(); LeadingJetPT[3]= (TH1F*)gDirectory->Get("jet1_pt15");
	  Leading2JetPT[3]->SetDefaultSumw2(); Leading2JetPT[3]= (TH1F*)gDirectory->Get("jet2_pt15");
	  Leading3JetPT[3]->SetDefaultSumw2(); Leading3JetPT[3]= (TH1F*)gDirectory->Get("jet3_pt15");
	  Leading4JetPT[3]->SetDefaultSumw2(); Leading4JetPT[3]= (TH1F*)gDirectory->Get("jet4_pt15");
	  Leading5JetPT[3]->SetDefaultSumw2(); Leading5JetPT[3]= (TH1F*)gDirectory->Get("jet5_pt15");
	  Leading6JetPT[3]->SetDefaultSumw2(); Leading6JetPT[3]= (TH1F*)gDirectory->Get("jet6_pt15");	  
	  THT[3]->SetDefaultSumw2(); THT[3]= (TH1F*)gDirectory->Get("THT15");
	  DRHjets[3]->SetDefaultSumw2(); DRHjets[3]= (TH1F*)gDirectory->Get("DeltaR_of_Higgs_Jets15");
	  DRWjets[3]->SetDefaultSumw2(); DRWjets[3]= (TH1F*)gDirectory->Get("DeltaR_of_W_Jets15");
	  Hpt[3]->SetDefaultSumw2(); Hpt[3]= (TH1F*)gDirectory->Get("HPt15");
	  Tpt[3]->SetDefaultSumw2(); Tpt[3]= (TH1F*)gDirectory->Get("TPt15");
	  DRWH[3]->SetDefaultSumw2(); DRWH[3]= (TH1F*)gDirectory->Get("DeltaR_of_W_Higgs15");
	  DPHjets[3]->SetDefaultSumw2(); DPHjets[3]= (TH1F*)gDirectory->Get("DeltaPhi_of_Higgs_jets15");
	  DPWjets[3]->SetDefaultSumw2(); DPWjets[3]= (TH1F*)gDirectory->Get("DeltaPhi_of_W_jets15");
	  DPTjets[3]->SetDefaultSumw2(); DPTjets[3]= (TH1F*)gDirectory->Get("DeltaPhi_of_T_jet15");
	  HiggsMass[3]->SetDefaultSumw2(); HiggsMass[3]= (TH1F*)gDirectory->Get("HM15");
	  RelHT[3]->SetDefaultSumw2(); RelHT[3]= (TH1F*)gDirectory->Get("RelHT15");
	  DRTH[3]->SetDefaultSumw2(); DRTH[3]= (TH1F*)gDirectory->Get("DeltaR_of_Top_Higgs15");
	  PtNormalizedMass[3]->SetDefaultSumw2(); PtNormalizedMass[3]= (TH1F*)gDirectory->Get("PT_Normalized_Mass15");
	  RelativeMass[3]->SetDefaultSumw2(); RelativeMass[3]= (TH1F*)gDirectory->Get("Relative_Mass15");
	  MotherPtNormalizedMass[3]->SetDefaultSumw2(); MotherPtNormalizedMass[3]= (TH1F*)gDirectory->Get("Mother_PT_Normalized_Mass15");
	  NumberOfTops[3]->SetDefaultSumw2();NumberOfTops[3]= (TH1F*)gDirectory->Get("Number_of_Tops15");
	  HiggsMassOverTopMass[3]->SetDefaultSumw2();HiggsMassOverTopMass[3]= (TH1F*)gDirectory->Get("HMoverTM15");
	  HiggsTopAsymmetry[3]->SetDefaultSumw2();HiggsTopAsymmetry[3]= (TH1F*)gDirectory->Get("HTAsym15");
	  ThirdLooseBtag[3]->SetDefaultSumw2();ThirdLooseBtag[3]= (TH1F*)gDirectory->Get("TLBTag15");
	  TopMass[3]->SetDefaultSumw2(); TopMass[3]= (TH1F*)gDirectory->Get("TMass15");	 
	  Chi2[3]->SetDefaultSumw2(); Chi2[3]= (TH1F*)gDirectory->Get("ChiSq15"); 
	}
      //CurrentFile[i]->Close();
    }

  //Calculation of Estimator
  double bkg[3]={0,0,0};
  double sig=0;

  sig=TprimeHistos[3]->Integral(TprimeHistos[3]->GetXaxis()->FindBin(710),TprimeHistos[3]->GetXaxis()->FindBin(750));
  for (int k=0; k<3; k++) bkg[k]=TprimeHistos[k]->Integral(TprimeHistos[k]->GetXaxis()->FindBin(710),TprimeHistos[k]->GetXaxis()->FindBin(750));
  
  cout << "Estimator under the peak!" << endl;
  cout << "Number of signal events: " << sig << ", Number of BKG events: " << bkg[0]+bkg[1]+bkg[2] << ", S/sqrt(S+B)=" << sig/sqrt(sig+bkg[0]+bkg[1]+bkg[2]) << endl;

  double bkgFR[3]={0,0,0};
  double sigFR=0;

  sigFR=TprimeHistos[3]->Integral();
  for (int k=0; k<3; k++) bkgFR[k]=TprimeHistos[k]->Integral();
  
  cout << "Estimator in full range!" << endl;
  cout << "Number of signal events: " << sigFR << ", Number of BKG events: " << bkgFR[0]+bkgFR[1]+bkgFR[2] << ", S/sqrt(S+B)=" << sigFR/sqrt(sigFR+bkgFR[0]+bkgFR[1]+bkgFR[2]) << endl;
  
  THStack *BKGandSignal = new THStack("BKGandSignal", "BKG and signal for 5 jets; M_{5j} GeV; Events");
  THStack *BKGLJPT = new THStack("BKGLJPT", "BKG for Leading Jet PT; pT(j_{1}) GeV; Events");
  THStack *BKGL2JPT = new THStack("BKGLJPT", "BKG for Leading 2 Jet PT; pT(j_{2}) GeV; Events");
  THStack *BKGL3JPT = new THStack("BKGLJPT", "BKG for Leading 3 Jet PT; pT(j_{3}) GeV; Events");
  THStack *BKGL4JPT = new THStack("BKGLJPT", "BKG for Leading 4 Jet PT; pT(j_{4}) GeV; Events");
  THStack *BKGL5JPT = new THStack("BKGLJPT", "BKG for Leading 5 Jet PT; pT(j_{5}) GeV; Events");
  THStack *BKGL6JPT = new THStack("BKGLJPT", "BKG for Leading 6 Jet PT; pT(j_{6}) GeV; Events");  
  THStack *BKGTHT = new THStack("BKGTHT", "BKG for HT; HT GeV; Events");
  THStack *BKGDRHjets= new THStack("BKGDRHjets", "BKG for DRHjets; #Delta R((jj)^{H}); Events");
  THStack *BKGDRWjets= new THStack("BKGDRWjets", "BKG for DRWjets; #Delta R((jj)^{W}); Events");
  THStack *BKGHpt= new THStack("BKGHpt", "BKG for Hpt; p_{T}(H) GeV; Events");
  THStack *BKGTpt= new THStack("BKGTpt", "BKG for Tpt; p_{T}(t) GeV; Events");
  THStack *BKGDRWH= new THStack("BKGDRWH", "BKG for DRWH; #Delta R(WH); Events");
  THStack *BKGDPHjets= new THStack("BKGDPHjets", "BKG for DPHjets; #Delta #phi((jj)^{H}); Events");
  THStack *BKGDPWjets= new THStack("BKGDPWjets", "BKG for DPWjets; #Delta #phi((jj)^{W}); Events");
  THStack *BKGDPTjets= new THStack("BKGDPTjets", "BKG for DPTjets; #Delta #phi(Wj^{t}); Events");
  THStack *BKGHiggsMass= new THStack("BKGHiggsMass", "BKG for HiggsMass; M_{H} GeV; Events");
  THStack *BKGRelHT= new THStack("BKGRelHT", "BKG for RelHT; (p_{T}(H)+p_{T}(t))/H_{T}; Events");
  THStack *BKGDRTH= new THStack("BKGDRTH", "BKG for DRTH; #Delta R(TH); Events");
  THStack *BKGPtNormalizedMass= new THStack("BKGPtNormalizedMass", "BKG for PtNormalizedMass; PTNM; Events");
  THStack *BKGRelativeMass= new THStack("BKGRelativeMass", "BKG for RelativeMass; RM; Events");
  THStack *BKGMotherPtNormalizedMass= new THStack("BKGMotherPtNormalizedMass", "BKG for MotherPtNormalizedMass; MPTNM; Events");
  THStack *BKGNumberOfTops= new THStack("BKGNumberOfTops", "BKG for Number of Tops; NTops; Events");
  THStack *BKGHMTM= new THStack("BKGHMTM", "BKG for Higgs Mass over Top Mass; M_{H}/M_{t}; Events");
  THStack *BKGHTAsym= new THStack("BKGHTAsym", "BKG for Higgs-Top Assymetry; (H_{PT}/M_{H}-t_{PT}/M_{t})/(H_{PT}/M_{H}+t_{PT}/M_{t}); Events");
  THStack *BKGTLBT= new THStack("BKGTLBT", "BKG for Third loose b-tag; n^{CSVLnotCSVM}_b; Events");
  THStack *BKGTopMass= new THStack("BKGTopMass", "BKG for TopMass; M_{t} GeV; Events");  
  THStack *BKGChi2= new THStack("BKGChi2", "BKG for Chi2; #Chi^{2}; Events");
  cout << "1 Marker" << endl;
  for (int i=0; i<NOS; i++)
    {
      //TprimeHistos[i]->GetXaxis()->SetRange(300,1600);
      //TprimeHistos[i]->GetYaxis()->SetRange(1,50000);
      BKGandSignal->Add(TprimeHistos[i]);
      if (i!=NOS-1)
	{
	  BKGLJPT->Add(LeadingJetPT[i]);
	  BKGL2JPT->Add(Leading2JetPT[i]);
	  BKGL3JPT->Add(Leading3JetPT[i]);
	  BKGL4JPT->Add(Leading4JetPT[i]);
	  BKGL5JPT->Add(Leading5JetPT[i]);
	  BKGL6JPT->Add(Leading6JetPT[i]);	  
	  BKGTHT->Add(THT[i]);
	  BKGDRHjets->Add(DRHjets[i]);
	  BKGDRWjets->Add(DRWjets[i]);
	  BKGHpt->Add(Hpt[i]);
	  BKGTpt->Add(Tpt[i]);
	  BKGDRWH->Add(DRWH[i]);
	  BKGDPHjets->Add(DPHjets[i]);
	  BKGDPWjets->Add(DPWjets[i]);
	  BKGDPTjets->Add(DPTjets[i]);
	  BKGHiggsMass->Add(HiggsMass[i]);
	  BKGRelHT->Add(RelHT[i]);
	  BKGDRTH->Add(DRTH[i]);
	  BKGPtNormalizedMass->Add(PtNormalizedMass[i]);
	  BKGRelativeMass->Add(RelativeMass[i]);
	  BKGMotherPtNormalizedMass->Add(MotherPtNormalizedMass[i]);
	  BKGNumberOfTops->Add(NumberOfTops[i]);
	  BKGHMTM->Add(HiggsMassOverTopMass[i]);
	  BKGHTAsym->Add(HiggsTopAsymmetry[i]);
	  BKGTLBT->Add(ThirdLooseBtag[i]);
	  BKGTopMass->Add(TopMass[i]);
	  BKGChi2->Add(Chi2[i]);
	}
    }
  cout << "2 Marker" << endl;

  //PostScript Plotting
  TLegend* BKGandSignallegend = new TLegend(0.75,0.65,0.90,0.9);
  BKGandSignallegend->AddEntry(TprimeHistos[2], "QCD", "f");
  BKGandSignallegend->AddEntry(TprimeHistos[1], "TTbar", "f");
  BKGandSignallegend->AddEntry(TprimeHistos[0], "SingleT", "f");
  BKGandSignallegend->AddEntry(TprimeHistos[3], "Signal", "f");
  char FileName[100];
  sprintf(FileName,"Single_Tprime_SUFFIX.pdf");
  //TPostScript *ps = new TPostScript(FileName,111);
  TPDF *ps = new TPDF(FileName,111);
  TCanvas *MyPlot = new TCanvas("MyPlot","Single t prime to top Higgs with backgrounds",600,800);
  TPaveText *p1;
  p1 = new TPaveText(3,0.5,8,3.5);
  char dateandtime[50];
  sprintf(dateandtime,"date: %s, time: %s",__DATE__,__TIME__);
  p1->AddText(dateandtime);
  //////////////
  //First Page//
  //////////////
  MyPlot->Clear();
  MyPlot->cd(1);
  ps->NewPage();
  gStyle->SetOptStat(0);//Remove the Stat Box
  //TprimeHistos[3]->Draw("hist");
  //TprimeHistos[3]->GetYaxis()->SetRangeUser(0.1,BKGandSignal->GetMaximum());
  BKGandSignal->Draw("hist");
  //BKGandSignal->GetYaxis()->SetRangeUser(0.1,BKGandSignal->GetMaximum());
  TprimeHistos[3]->Draw("histsame");
  gPad->SetLogy();
  //TprimeHistos[3]->Draw("histsame"); //Preserve double drawing of signal
  BKGandSignal->SetMinimum(0.1);
  BKGandSignallegend->Draw();
  gPad->RedrawAxis("");
  gPad->Update();
  MyPlot->Update();
  ///////////////
  //Second Page//
  ///////////////
  MyPlot->Clear();
  MyPlot->cd(1);
  ps->NewPage();
  BKGLJPT->Draw("hist");
  LeadingJetPT[3]->Draw("histsame");
  gPad->SetLogy();
  BKGLJPT->SetMinimum(0.1);
  BKGandSignallegend->Draw();
  gPad->RedrawAxis("");
  gPad->Update();
  MyPlot->Update();
  //
  MyPlot->Clear();
  MyPlot->cd(1);
  ps->NewPage();
  BKGL2JPT->Draw("hist");
  Leading2JetPT[3]->Draw("histsame");
  gPad->SetLogy();
  BKGL2JPT->SetMinimum(0.1);
  BKGandSignallegend->Draw();
  gPad->RedrawAxis("");
  gPad->Update();
  MyPlot->Update();
  //
  MyPlot->Clear();
  MyPlot->cd(1);
  ps->NewPage();
  BKGL3JPT->Draw("hist");
  Leading3JetPT[3]->Draw("histsame");
  gPad->SetLogy();
  BKGL3JPT->SetMinimum(0.1);
  BKGandSignallegend->Draw();
  gPad->RedrawAxis("");
  gPad->Update();
  MyPlot->Update();
  //
  MyPlot->Clear();
  MyPlot->cd(1);
  ps->NewPage();
  BKGL4JPT->Draw("hist");
  Leading4JetPT[3]->Draw("histsame");
  gPad->SetLogy();
  BKGL4JPT->SetMinimum(0.1);
  BKGandSignallegend->Draw();
  gPad->RedrawAxis("");
  gPad->Update();
  MyPlot->Update();
  //
  MyPlot->Clear();
  MyPlot->cd(1);
  ps->NewPage();
  BKGL5JPT->Draw("hist");
  Leading5JetPT[3]->Draw("histsame");
  gPad->SetLogy();
  BKGL5JPT->SetMinimum(0.1);
  BKGandSignallegend->Draw();
  gPad->RedrawAxis("");
  gPad->Update();
  MyPlot->Update();
  //
  MyPlot->Clear();
  MyPlot->cd(1);
  ps->NewPage();
  BKGL6JPT->Draw("hist");
  Leading6JetPT[3]->Draw("histsame");
  gPad->SetLogy();
  BKGL6JPT->SetMinimum(0.1);
  BKGandSignallegend->Draw();
  gPad->RedrawAxis("");
  gPad->Update();
  MyPlot->Update();
  //////////////
  //Third Page//
  //////////////
  MyPlot->Clear();
  MyPlot->cd(1);
  ps->NewPage();
  BKGTHT->Draw("hist");
  THT[3]->Draw("histsame");
  gPad->SetLogy();
  BKGTHT->SetMinimum(0.1);
  BKGandSignallegend->Draw();
  gPad->RedrawAxis("");
  gPad->Update();
  MyPlot->Update();
  ///////////////
  //Fourth Page//
  ///////////////
  MyPlot->Clear();
  MyPlot->cd(1);
  ps->NewPage();
  BKGDRHjets->Draw("hist");
  DRHjets[3]->Draw("histsame");
  gPad->SetLogy();
  BKGDRHjets->SetMinimum(0.1);
  BKGandSignallegend->Draw();
  gPad->RedrawAxis("");
  gPad->Update();
  MyPlot->Update();
  ////////////
  //5th Page//
  ////////////
  MyPlot->Clear();
  MyPlot->cd(1);
  ps->NewPage();
  BKGDRWjets->Draw("hist");
  DRWjets[3]->Draw("histsame");
  gPad->SetLogy();
  BKGDRWjets->SetMinimum(0.1);
  BKGandSignallegend->Draw();
  gPad->RedrawAxis("");
  gPad->Update();
  MyPlot->Update();
  ////////////
  //6th Page//
  ////////////
  MyPlot->Clear();
  MyPlot->cd(1);
  ps->NewPage();
  BKGHpt->Draw("hist");
  Hpt[3]->Draw("histsame");
  gPad->SetLogy();
  BKGHpt->SetMinimum(0.1);
  BKGandSignallegend->Draw();
  gPad->RedrawAxis("");
  gPad->Update();
  MyPlot->Update();
  ////////////
  //7th Page//
  ////////////
  MyPlot->Clear();
  MyPlot->cd(1);
  ps->NewPage();
  BKGTpt->Draw("hist");
  Tpt[3]->Draw("histsame");
  gPad->SetLogy();
  BKGTpt->SetMinimum(0.1);
  BKGandSignallegend->Draw();
  gPad->RedrawAxis("");
  gPad->Update();
  MyPlot->Update();
  ////////////
  //8th Page//
  ////////////
  MyPlot->Clear();
  MyPlot->cd(1);
  ps->NewPage();
  BKGDRWH->Draw("hist");
  DRWH[3]->Draw("histsame");
  gPad->SetLogy();
  BKGDRWH->SetMinimum(0.1);
  BKGandSignallegend->Draw();
  gPad->RedrawAxis("");
  gPad->Update();
  MyPlot->Update();
  ////////////
  //9th Page//
  ////////////
  MyPlot->Clear();
  MyPlot->cd(1);
  ps->NewPage();
  BKGDPHjets->Draw("hist");
  DPHjets[3]->Draw("histsame");
  gPad->SetLogy();
  BKGDPHjets->SetMinimum(0.1);
  BKGandSignallegend->Draw();
  gPad->RedrawAxis("");
  gPad->Update();
  MyPlot->Update();
  /////////////
  //10th Page//
  /////////////
  MyPlot->Clear();
  MyPlot->cd(1);
  ps->NewPage();
  BKGDPWjets->Draw("hist");
  DPWjets[3]->Draw("histsame");
  gPad->SetLogy();
  BKGDPWjets->SetMinimum(0.1);
  BKGandSignallegend->Draw();
  gPad->RedrawAxis("");
  gPad->Update();
  MyPlot->Update();
  /////////////
  //11th Page//
  /////////////
  MyPlot->Clear();
  MyPlot->cd(1);
  ps->NewPage();
  BKGDPTjets->Draw("hist");
  DPTjets[3]->Draw("histsame");
  gPad->SetLogy();
  BKGDPTjets->SetMinimum(0.1);
  BKGandSignallegend->Draw();
  gPad->RedrawAxis("");
  gPad->Update();
  MyPlot->Update();
  /////////////
  //12th Page//
  /////////////
  MyPlot->Clear();
  MyPlot->cd(1);
  ps->NewPage();
  BKGHiggsMass->Draw("hist");
  HiggsMass[3]->Draw("histsame");
  gPad->SetLogy();
  BKGHiggsMass->SetMinimum(0.1);
  BKGandSignallegend->Draw();
  gPad->RedrawAxis("");
  gPad->Update();
  MyPlot->Update();
  /////////////
  //13th Page//
  /////////////
  MyPlot->Clear();
  MyPlot->cd(1);
  ps->NewPage();
  BKGRelHT->Draw("hist");
  RelHT[3]->Draw("histsame");
  gPad->SetLogy();
  BKGRelHT->SetMinimum(0.1);
  BKGandSignallegend->Draw();
  gPad->RedrawAxis("");
  gPad->Update();
  MyPlot->Update();
  /////////////
  //14th Page//
  /////////////
  MyPlot->Clear();
  MyPlot->cd(1);
  ps->NewPage();
  BKGDRTH->Draw("hist");
  DRTH[3]->Draw("histsame");
  gPad->SetLogy();
  BKGDRTH->SetMinimum(0.1);
  BKGandSignallegend->Draw();
  gPad->RedrawAxis("");
  gPad->Update();
  MyPlot->Update();
  /////////////
  //15th Page//
  /////////////
  MyPlot->Clear();
  MyPlot->cd(1);
  ps->NewPage();
  BKGPtNormalizedMass->Draw("hist");
  PtNormalizedMass[3]->Draw("histsame");
  gPad->SetLogy();
  BKGPtNormalizedMass->SetMinimum(0.1);
  BKGandSignallegend->Draw();
  gPad->RedrawAxis("");
  gPad->Update();
  MyPlot->Update();
  /////////////
  //16th Page//
  /////////////
  MyPlot->Clear();
  MyPlot->cd(1);
  ps->NewPage();
  BKGRelativeMass->Draw("hist");
  RelativeMass[3]->Draw("histsame");
  gPad->SetLogy();
  BKGRelativeMass->SetMinimum(0.1);
  BKGandSignallegend->Draw();
  gPad->RedrawAxis("");
  gPad->Update();
  MyPlot->Update();
  /////////////
  //17th Page//
  /////////////
  MyPlot->Clear();
  MyPlot->cd(1);
  ps->NewPage();
  BKGMotherPtNormalizedMass->Draw("hist");
  MotherPtNormalizedMass[3]->Draw("histsame");
  gPad->SetLogy();
  BKGMotherPtNormalizedMass->SetMinimum(0.1);
  BKGandSignallegend->Draw();
  gPad->RedrawAxis("");
  gPad->Update();
  MyPlot->Update();
  //////////////
  //////////////
  //////////////
  MyPlot->Clear();
  MyPlot->cd(1);
  ps->NewPage();
  BKGNumberOfTops->Draw("hist");
  NumberOfTops[3]->Draw("histsame");
  gPad->SetLogy();
  BKGNumberOfTops->SetMinimum(0.1);
  BKGandSignallegend->Draw();
  gPad->RedrawAxis("");
  gPad->Update();
  MyPlot->Update();
  //////////////
  //////////////
  //////////////
  MyPlot->Clear();
  MyPlot->cd(1);
  ps->NewPage();
  BKGHMTM->Draw("hist");
  HiggsMassOverTopMass[3]->Draw("histsame");
  gPad->SetLogy();
  BKGHMTM->SetMinimum(0.1);
  BKGandSignallegend->Draw();
  gPad->RedrawAxis("");
  gPad->Update();
  MyPlot->Update();
  //////////////
  //////////////
  //////////////
  MyPlot->Clear();
  MyPlot->cd(1);
  ps->NewPage();
  BKGHTAsym->Draw("hist");
  HiggsTopAsymmetry[3]->Draw("histsame");
  gPad->SetLogy();
  BKGHTAsym->SetMinimum(0.1);
  BKGandSignallegend->Draw();
  gPad->RedrawAxis("");
  gPad->Update();
  MyPlot->Update();
  //////////////
  //////////////
  //////////////
  MyPlot->Clear();
  MyPlot->cd(1);
  ps->NewPage();
  BKGTLBT->Draw("hist");
  ThirdLooseBtag[3]->Draw("histsame");
  gPad->SetLogy();
  BKGTLBT->SetMinimum(0.1);
  BKGandSignallegend->Draw();
  gPad->RedrawAxis("");
  gPad->Update();
  MyPlot->Update();
  //////////////
  //////////////
  //////////////
  MyPlot->Clear();
  MyPlot->cd(1);
  ps->NewPage();
  BKGTopMass->Draw("hist");
  TopMass[3]->Draw("histsame");
  gPad->SetLogy();
  BKGTopMass->SetMinimum(0.1);
  BKGandSignallegend->Draw();
  gPad->RedrawAxis("");
  gPad->Update();
  MyPlot->Update();
  //////////////
  //////////////
  //////////////
  MyPlot->Clear();
  MyPlot->cd(1);
  ps->NewPage();
  BKGChi2->Draw("hist");
  Chi2[3]->Draw("histsame");
  gPad->SetLogy();
  BKGChi2->SetMinimum(0.1);
  BKGandSignallegend->Draw();
  gPad->RedrawAxis("");
  gPad->Update();
  MyPlot->Update();
  //////////////
  ps->Close();

  exit(0);
}
