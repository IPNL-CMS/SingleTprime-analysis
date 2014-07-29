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
#include "TMatrixD.h"
#include "TArrayD.h"
#include "TVectorD.h"
#include "TMatrixDLazy.h"
#include "TDecompLU.h"

using namespace std;

const TString StorageDirPrefixLoose="file:/afs/cern.ch/work/j/jruizalv/private/Analysis/CMSSW_5_3_9_patch3/src/Extractors/PatExtractor/bin/WithTrigger/SUFFIXLOOSE/";
const TString StorageDirPrefixMedium="file:/afs/cern.ch/work/j/jruizalv/private/Analysis/CMSSW_5_3_9_patch3/src/Extractors/PatExtractor/bin/WithTrigger/SUFFIXMEDIUM/";
const TString StorageDirPrefixTight="file:/afs/cern.ch/work/j/jruizalv/private/Analysis/CMSSW_5_3_9_patch3/src/Extractors/PatExtractor/bin/WithTrigger/SUFFIXTIGHT/";
const int NOS=5;
//const TString Samples[NOS] = {"Diboson.root", "SingleTop.root", "TTJets.root", "Zjets.root", "Wjets.root", "QCD.root", "Signal.root"};
const TString Samples[NOS] = {"Diboson.root", "SingleTop.root", "TTJets.root", "QCD.root", "Signal.root"};

int theNumber(TFile * theFile){
  //TH1F*histogram;
  for(int k=0;k<=16;k++){
    string name="TprimeMass";
    char kchar[20];sprintf(kchar,"%d",k);
    if((TH1F*)(theFile->Get((name+kchar).c_str()))) return k;
  }
  return -1;
}

void BKGEstimationTesterFull()
{
  //Reading files for loose sample
  TFile *CurrentFileL[NOS];
  TH1F *TprimeHistosL[NOS];
  TH1F *LeadingJetPTL[NOS];
  TH1F *Leading2JetPTL[NOS];
  TH1F *Leading3JetPTL[NOS];
  TH1F *Leading4JetPTL[NOS];
  TH1F *Leading5JetPTL[NOS];
  TH1F *Leading6JetPTL[NOS];
  TH1F *THTL[NOS];
  TH1F *DRHjetsL[NOS];
  TH1F *DRWjetsL[NOS];
  TH1F *HptL[NOS];
  TH1F *TptL[NOS];
  TH1F *DRWHL[NOS];
  TH1F *DPHjetsL[NOS];
  TH1F *DPWjetsL[NOS];
  TH1F *DPTjetsL[NOS];
  TH1F *HiggsMassL[NOS];
  TH1F *RelHTL[NOS];
  TH1F *DRTHL[NOS];
  TH1F *PtNormalizedMassL[NOS];
  TH1F *RelativeMassL[NOS];
  TH1F *MotherPtNormalizedMassL[NOS];
  TH1F *NumberOfTopsL[NOS];
  TH1F *HiggsMassOverTopMassL[NOS];
  TH1F *HiggsTopAsymmetryL[NOS];
  TH1F *ThirdLooseBtagL[NOS];
  TH1F *TopMassL[NOS];
  TH1F *Chi2L[NOS];

  for (int i=0; i<NOS; i++)
    {
      CurrentFileL[i] = new TFile(StorageDirPrefixLoose + Samples[i], "READ");
      if ( CurrentFileL[i]->IsOpen() ) printf( Samples[i] + " File opened successfully\n");      
      int NumberSuffix=theNumber(CurrentFileL[i]); cout << NumberSuffix << endl;
      string A=Form("TprimeMass%i",NumberSuffix);
      TprimeHistosL[i]->SetDefaultSumw2(); TprimeHistosL[i]= (TH1F*)gDirectory->Get(A.c_str());
      A=Form("jet1_pt%i",NumberSuffix);
      LeadingJetPTL[i]->SetDefaultSumw2(); LeadingJetPTL[i]= (TH1F*)gDirectory->Get(A.c_str());
      A=Form("jet2_pt%i",NumberSuffix);
      Leading2JetPTL[i]->SetDefaultSumw2(); Leading2JetPTL[i]= (TH1F*)gDirectory->Get(A.c_str());
      A=Form("jet3_pt%i",NumberSuffix);
      Leading3JetPTL[i]->SetDefaultSumw2(); Leading3JetPTL[i]= (TH1F*)gDirectory->Get(A.c_str());
      A=Form("jet4_pt%i",NumberSuffix);
      Leading4JetPTL[i]->SetDefaultSumw2(); Leading4JetPTL[i]= (TH1F*)gDirectory->Get(A.c_str());
      A=Form("jet5_p%i",NumberSuffix);
      Leading5JetPTL[i]->SetDefaultSumw2(); Leading5JetPTL[i]= (TH1F*)gDirectory->Get(A.c_str());
      A=Form("jet6_pt%i",NumberSuffix);
      Leading6JetPTL[i]->SetDefaultSumw2(); Leading6JetPTL[i]= (TH1F*)gDirectory->Get(A.c_str());
      A=Form("THT%i",NumberSuffix);
      THTL[i]->SetDefaultSumw2(); THTL[i]= (TH1F*)gDirectory->Get(A.c_str());
      A=Form("DeltaR_of_Higgs_Jets%i",NumberSuffix);
      DRHjetsL[i]->SetDefaultSumw2(); DRHjetsL[i]= (TH1F*)gDirectory->Get(A.c_str());
      A=Form("DeltaR_of_W_Jets%i",NumberSuffix);
      DRWjetsL[i]->SetDefaultSumw2(); DRWjetsL[i]= (TH1F*)gDirectory->Get(A.c_str());
      A=Form("HPt%i",NumberSuffix);
      HptL[i]->SetDefaultSumw2(); HptL[i]= (TH1F*)gDirectory->Get(A.c_str());
      A=Form("TPt%i",NumberSuffix);
      TptL[i]->SetDefaultSumw2(); TptL[i]= (TH1F*)gDirectory->Get(A.c_str());
      A=Form("DeltaR_of_W_Higgs%i",NumberSuffix);
      DRWHL[i]->SetDefaultSumw2(); DRWHL[i]= (TH1F*)gDirectory->Get(A.c_str());
      A=Form("DeltaPhi_of_Higgs_jets%i",NumberSuffix);
      DPHjetsL[i]->SetDefaultSumw2(); DPHjetsL[i]= (TH1F*)gDirectory->Get(A.c_str());
      A=Form("DeltaPhi_of_W_jets%i",NumberSuffix);
      DPWjetsL[i]->SetDefaultSumw2(); DPWjetsL[i]= (TH1F*)gDirectory->Get(A.c_str());
      A=Form("DeltaPhi_of_T_jet%i",NumberSuffix);
      DPTjetsL[i]->SetDefaultSumw2(); DPTjetsL[i]= (TH1F*)gDirectory->Get(A.c_str());
      A=Form("HM%i",NumberSuffix);
      HiggsMassL[i]->SetDefaultSumw2(); HiggsMassL[i]= (TH1F*)gDirectory->Get(A.c_str());
      A=Form("RelHT%i",NumberSuffix);
      RelHTL[i]->SetDefaultSumw2(); RelHTL[i]= (TH1F*)gDirectory->Get(A.c_str());
      A=Form("DeltaR_of_Top_Higgs%i",NumberSuffix);
      DRTHL[i]->SetDefaultSumw2(); DRTHL[i]= (TH1F*)gDirectory->Get(A.c_str());
      A=Form("PT_Normalized_Mass%i",NumberSuffix);
      PtNormalizedMassL[i]->SetDefaultSumw2(); PtNormalizedMassL[i]= (TH1F*)gDirectory->Get(A.c_str());
      A=Form("Relative_Mass%i",NumberSuffix);
      RelativeMassL[i]->SetDefaultSumw2(); RelativeMassL[i]= (TH1F*)gDirectory->Get(A.c_str());
      A=Form("Mother_PT_Normalized_Mass%i",NumberSuffix);
      MotherPtNormalizedMassL[i]->SetDefaultSumw2(); MotherPtNormalizedMassL[i]= (TH1F*)gDirectory->Get(A.c_str());
      A=Form("Number_of_Tops%i",NumberSuffix);
      NumberOfTopsL[i]->SetDefaultSumw2();NumberOfTopsL[i]= (TH1F*)gDirectory->Get(A.c_str());
      A=Form("HMoverTM%i",NumberSuffix);
      HiggsMassOverTopMassL[i]->SetDefaultSumw2();HiggsMassOverTopMassL[i]= (TH1F*)gDirectory->Get(A.c_str());
      A=Form("HTAsym%i",NumberSuffix);
      HiggsTopAsymmetryL[i]->SetDefaultSumw2();HiggsTopAsymmetryL[i]= (TH1F*)gDirectory->Get(A.c_str());
      A=Form("TLBTag%i",NumberSuffix);
      HiggsTopAsymmetryL[i]->SetDefaultSumw2();HiggsTopAsymmetryL[i]= (TH1F*)gDirectory->Get(A.c_str());
      A=Form("TMass%i",NumberSuffix);
      ThirdLooseBtagL[i]->SetDefaultSumw2();ThirdLooseBtagL[i]= (TH1F*)gDirectory->Get(A.c_str());
      A=Form("ChiSq%i",NumberSuffix);
      TopMassL[i]->SetDefaultSumw2(); TopMassL[i]= (TH1F*)gDirectory->Get(A.c_str());
      A=Form("UQC%i",NumberSuffix);
      Chi2L[i]->SetDefaultSumw2(); Chi2L[i]= (TH1F*)gDirectory->Get(A.c_str());
    }

  //Reading files for medium sample
  TFile *CurrentFileM[NOS];
  TH1F *TprimeHistosM[NOS];
  TH1F *LeadingJetPTM[NOS];
  TH1F *Leading2JetPTM[NOS];
  TH1F *Leading3JetPTM[NOS];
  TH1F *Leading4JetPTM[NOS];
  TH1F *Leading5JetPTM[NOS];
  TH1F *Leading6JetPTM[NOS];
  TH1F *THTM[NOS];
  TH1F *DRHjetsM[NOS];
  TH1F *DRWjetsM[NOS];
  TH1F *HptM[NOS];
  TH1F *TptM[NOS];
  TH1F *DRWHM[NOS];
  TH1F *DPHjetsM[NOS];
  TH1F *DPWjetsM[NOS];
  TH1F *DPTjetsM[NOS];
  TH1F *HiggsMassM[NOS];
  TH1F *RelHTM[NOS];
  TH1F *DRTHM[NOS];
  TH1F *PtNormalizedMassM[NOS];
  TH1F *RelativeMassM[NOS];
  TH1F *MotherPtNormalizedMassM[NOS];
  TH1F *NumberOfTopsM[NOS];
  TH1F *HiggsMassOverTopMassM[NOS];
  TH1F *HiggsTopAsymmetryM[NOS];
  TH1F *ThirdLooseBtagM[NOS];
  TH1F *TopMassM[NOS];
  TH1F *Chi2M[NOS];

  for (int i=0; i<NOS; i++)
    {
      CurrentFileM[i] = new TFile(StorageDirPrefixMedium + Samples[i], "READ");
      if ( CurrentFileM[i]->IsOpen() ) printf( Samples[i] + " File opened successfully\n");      
      int NumberSuffix=theNumber(CurrentFileM[i]); cout << NumberSuffix << endl;
      string A=Form("TprimeMass%i",NumberSuffix);
      TprimeHistosM[i]->SetDefaultSumw2(); TprimeHistosM[i]= (TH1F*)gDirectory->Get(A.c_str());
      A=Form("jet1_pt%i",NumberSuffix);
      LeadingJetPTM[i]->SetDefaultSumw2(); LeadingJetPTM[i]= (TH1F*)gDirectory->Get(A.c_str());
      A=Form("jet2_pt%i",NumberSuffix);
      Leading2JetPTM[i]->SetDefaultSumw2(); Leading2JetPTM[i]= (TH1F*)gDirectory->Get(A.c_str());
      A=Form("jet3_pt%i",NumberSuffix);
      Leading3JetPTM[i]->SetDefaultSumw2(); Leading3JetPTM[i]= (TH1F*)gDirectory->Get(A.c_str());
      A=Form("jet4_pt%i",NumberSuffix);
      Leading4JetPTM[i]->SetDefaultSumw2(); Leading4JetPTM[i]= (TH1F*)gDirectory->Get(A.c_str());
      A=Form("jet5_p%i",NumberSuffix);
      Leading5JetPTM[i]->SetDefaultSumw2(); Leading5JetPTM[i]= (TH1F*)gDirectory->Get(A.c_str());
      A=Form("jet6_pt%i",NumberSuffix);
      Leading6JetPTM[i]->SetDefaultSumw2(); Leading6JetPTM[i]= (TH1F*)gDirectory->Get(A.c_str());
      A=Form("THT%i",NumberSuffix);
      THTM[i]->SetDefaultSumw2(); THTM[i]= (TH1F*)gDirectory->Get(A.c_str());
      A=Form("DeltaR_of_Higgs_Jets%i",NumberSuffix);
      DRHjetsM[i]->SetDefaultSumw2(); DRHjetsM[i]= (TH1F*)gDirectory->Get(A.c_str());
      A=Form("DeltaR_of_W_Jets%i",NumberSuffix);
      DRWjetsM[i]->SetDefaultSumw2(); DRWjetsM[i]= (TH1F*)gDirectory->Get(A.c_str());
      A=Form("HPt%i",NumberSuffix);
      HptM[i]->SetDefaultSumw2(); HptM[i]= (TH1F*)gDirectory->Get(A.c_str());
      A=Form("TPt%i",NumberSuffix);
      TptM[i]->SetDefaultSumw2(); TptM[i]= (TH1F*)gDirectory->Get(A.c_str());
      A=Form("DeltaR_of_W_Higgs%i",NumberSuffix);
      DRWHM[i]->SetDefaultSumw2(); DRWHM[i]= (TH1F*)gDirectory->Get(A.c_str());
      A=Form("DeltaPhi_of_Higgs_jets%i",NumberSuffix);
      DPHjetsM[i]->SetDefaultSumw2(); DPHjetsM[i]= (TH1F*)gDirectory->Get(A.c_str());
      A=Form("DeltaPhi_of_W_jets%i",NumberSuffix);
      DPWjetsM[i]->SetDefaultSumw2(); DPWjetsM[i]= (TH1F*)gDirectory->Get(A.c_str());
      A=Form("DeltaPhi_of_T_jet%i",NumberSuffix);
      DPTjetsM[i]->SetDefaultSumw2(); DPTjetsM[i]= (TH1F*)gDirectory->Get(A.c_str());
      A=Form("HM%i",NumberSuffix);
      HiggsMassM[i]->SetDefaultSumw2(); HiggsMassM[i]= (TH1F*)gDirectory->Get(A.c_str());
      A=Form("RelHT%i",NumberSuffix);
      RelHTM[i]->SetDefaultSumw2(); RelHTM[i]= (TH1F*)gDirectory->Get(A.c_str());
      A=Form("DeltaR_of_Top_Higgs%i",NumberSuffix);
      DRTHM[i]->SetDefaultSumw2(); DRTHM[i]= (TH1F*)gDirectory->Get(A.c_str());
      A=Form("PT_Normalized_Mass%i",NumberSuffix);
      PtNormalizedMassM[i]->SetDefaultSumw2(); PtNormalizedMassM[i]= (TH1F*)gDirectory->Get(A.c_str());
      A=Form("Relative_Mass%i",NumberSuffix);
      RelativeMassM[i]->SetDefaultSumw2(); RelativeMassM[i]= (TH1F*)gDirectory->Get(A.c_str());
      A=Form("Mother_PT_Normalized_Mass%i",NumberSuffix);
      MotherPtNormalizedMassM[i]->SetDefaultSumw2(); MotherPtNormalizedMassM[i]= (TH1F*)gDirectory->Get(A.c_str());
      A=Form("Number_of_Tops%i",NumberSuffix);
      NumberOfTopsM[i]->SetDefaultSumw2();NumberOfTopsM[i]= (TH1F*)gDirectory->Get(A.c_str());
      A=Form("HMoverTM%i",NumberSuffix);
      HiggsMassOverTopMassM[i]->SetDefaultSumw2();HiggsMassOverTopMassM[i]= (TH1F*)gDirectory->Get(A.c_str());
      A=Form("HTAsym%i",NumberSuffix);
      HiggsTopAsymmetryM[i]->SetDefaultSumw2();HiggsTopAsymmetryM[i]= (TH1F*)gDirectory->Get(A.c_str());
      A=Form("TLBTag%i",NumberSuffix);
      HiggsTopAsymmetryM[i]->SetDefaultSumw2();HiggsTopAsymmetryM[i]= (TH1F*)gDirectory->Get(A.c_str());
      A=Form("TMass%i",NumberSuffix);
      ThirdLooseBtagM[i]->SetDefaultSumw2();ThirdLooseBtagM[i]= (TH1F*)gDirectory->Get(A.c_str());
      A=Form("ChiSq%i",NumberSuffix);
      TopMassM[i]->SetDefaultSumw2(); TopMassM[i]= (TH1F*)gDirectory->Get(A.c_str());
      A=Form("UQC%i",NumberSuffix);
      Chi2M[i]->SetDefaultSumw2(); Chi2M[i]= (TH1F*)gDirectory->Get(A.c_str());
    }
  
  //Reading files for tight sample
  TFile *CurrentFileT[NOS];
  TH1F *TprimeHistosT[NOS];
  TH1F *LeadingJetPTT[NOS];
  TH1F *Leading2JetPTT[NOS];
  TH1F *Leading3JetPTT[NOS];
  TH1F *Leading4JetPTT[NOS];
  TH1F *Leading5JetPTT[NOS];
  TH1F *Leading6JetPTT[NOS];
  TH1F *THTT[NOS];
  TH1F *DRHjetsT[NOS];
  TH1F *DRWjetsT[NOS];
  TH1F *HptT[NOS];
  TH1F *TptT[NOS];
  TH1F *DRWHT[NOS];
  TH1F *DPHjetsT[NOS];
  TH1F *DPWjetsT[NOS];
  TH1F *DPTjetsT[NOS];
  TH1F *HiggsMassT[NOS];
  TH1F *RelHTT[NOS];
  TH1F *DRTHT[NOS];
  TH1F *PtNormalizedMassT[NOS];
  TH1F *RelativeMassT[NOS];
  TH1F *MotherPtNormalizedMassT[NOS];
  TH1F *NumberOfTopsT[NOS];
  TH1F *HiggsMassOverTopMassT[NOS];
  TH1F *HiggsTopAsymmetryT[NOS];
  TH1F *ThirdLooseBtagT[NOS];
  TH1F *TopMassT[NOS];
  TH1F *Chi2T[NOS];

  for (int i=0; i<NOS; i++)
    {
      CurrentFileT[i] = new TFile(StorageDirPrefixTight + Samples[i], "READ");
      if ( CurrentFileT[i]->IsOpen() ) printf( Samples[i] + " File opened successfully\n");     
      int NumberSuffix=theNumber(CurrentFileT[i]); cout << NumberSuffix << endl;
      string A=Form("TprimeMass%i",NumberSuffix);
      TprimeHistosT[i]->SetDefaultSumw2(); TprimeHistosT[i]= (TH1F*)gDirectory->Get(A.c_str());
      A=Form("jet1_pt%i",NumberSuffix);
      LeadingJetPTT[i]->SetDefaultSumw2(); LeadingJetPTT[i]= (TH1F*)gDirectory->Get(A.c_str());
      A=Form("jet2_pt%i",NumberSuffix);
      Leading2JetPTT[i]->SetDefaultSumw2(); Leading2JetPTT[i]= (TH1F*)gDirectory->Get(A.c_str());
      A=Form("jet3_pt%i",NumberSuffix);
      Leading3JetPTT[i]->SetDefaultSumw2(); Leading3JetPTT[i]= (TH1F*)gDirectory->Get(A.c_str());
      A=Form("jet4_pt%i",NumberSuffix);
      Leading4JetPTT[i]->SetDefaultSumw2(); Leading4JetPTT[i]= (TH1F*)gDirectory->Get(A.c_str());
      A=Form("jet5_p%i",NumberSuffix);
      Leading5JetPTT[i]->SetDefaultSumw2(); Leading5JetPTT[i]= (TH1F*)gDirectory->Get(A.c_str());
      A=Form("jet6_pt%i",NumberSuffix);
      Leading6JetPTT[i]->SetDefaultSumw2(); Leading6JetPTT[i]= (TH1F*)gDirectory->Get(A.c_str());
      A=Form("THT%i",NumberSuffix);
      THTT[i]->SetDefaultSumw2(); THTT[i]= (TH1F*)gDirectory->Get(A.c_str());
      A=Form("DeltaR_of_Higgs_Jets%i",NumberSuffix);
      DRHjetsT[i]->SetDefaultSumw2(); DRHjetsT[i]= (TH1F*)gDirectory->Get(A.c_str());
      A=Form("DeltaR_of_W_Jets%i",NumberSuffix);
      DRWjetsT[i]->SetDefaultSumw2(); DRWjetsT[i]= (TH1F*)gDirectory->Get(A.c_str());
      A=Form("HPt%i",NumberSuffix);
      HptT[i]->SetDefaultSumw2(); HptT[i]= (TH1F*)gDirectory->Get(A.c_str());
      A=Form("TPt%i",NumberSuffix);
      TptT[i]->SetDefaultSumw2(); TptT[i]= (TH1F*)gDirectory->Get(A.c_str());
      A=Form("DeltaR_of_W_Higgs%i",NumberSuffix);
      DRWHT[i]->SetDefaultSumw2(); DRWHT[i]= (TH1F*)gDirectory->Get(A.c_str());
      A=Form("DeltaPhi_of_Higgs_jets%i",NumberSuffix);
      DPHjetsT[i]->SetDefaultSumw2(); DPHjetsT[i]= (TH1F*)gDirectory->Get(A.c_str());
      A=Form("DeltaPhi_of_W_jets%i",NumberSuffix);
      DPWjetsT[i]->SetDefaultSumw2(); DPWjetsT[i]= (TH1F*)gDirectory->Get(A.c_str());
      A=Form("DeltaPhi_of_T_jet%i",NumberSuffix);
      DPTjetsT[i]->SetDefaultSumw2(); DPTjetsT[i]= (TH1F*)gDirectory->Get(A.c_str());
      A=Form("HM%i",NumberSuffix);
      HiggsMassT[i]->SetDefaultSumw2(); HiggsMassT[i]= (TH1F*)gDirectory->Get(A.c_str());
      A=Form("RelHT%i",NumberSuffix);
      RelHTT[i]->SetDefaultSumw2(); RelHTT[i]= (TH1F*)gDirectory->Get(A.c_str());
      A=Form("DeltaR_of_Top_Higgs%i",NumberSuffix);
      DRTHT[i]->SetDefaultSumw2(); DRTHT[i]= (TH1F*)gDirectory->Get(A.c_str());
      A=Form("PT_Normalized_Mass%i",NumberSuffix);
      PtNormalizedMassT[i]->SetDefaultSumw2(); PtNormalizedMassT[i]= (TH1F*)gDirectory->Get(A.c_str());
      A=Form("Relative_Mass%i",NumberSuffix);
      RelativeMassT[i]->SetDefaultSumw2(); RelativeMassT[i]= (TH1F*)gDirectory->Get(A.c_str());
      A=Form("Mother_PT_Normalized_Mass%i",NumberSuffix);
      MotherPtNormalizedMassT[i]->SetDefaultSumw2(); MotherPtNormalizedMassT[i]= (TH1F*)gDirectory->Get(A.c_str());
      A=Form("Number_of_Tops%i",NumberSuffix);
      NumberOfTopsT[i]->SetDefaultSumw2();NumberOfTopsT[i]= (TH1F*)gDirectory->Get(A.c_str());
      A=Form("HMoverTM%i",NumberSuffix);
      HiggsMassOverTopMassT[i]->SetDefaultSumw2();HiggsMassOverTopMassT[i]= (TH1F*)gDirectory->Get(A.c_str());
      A=Form("HTAsym%i",NumberSuffix);
      HiggsTopAsymmetryT[i]->SetDefaultSumw2();HiggsTopAsymmetryT[i]= (TH1F*)gDirectory->Get(A.c_str());
      A=Form("TLBTag%i",NumberSuffix);
      HiggsTopAsymmetryT[i]->SetDefaultSumw2();HiggsTopAsymmetryT[i]= (TH1F*)gDirectory->Get(A.c_str());
      A=Form("TMass%i",NumberSuffix);
      ThirdLooseBtagT[i]->SetDefaultSumw2();ThirdLooseBtagT[i]= (TH1F*)gDirectory->Get(A.c_str());
      A=Form("ChiSq%i",NumberSuffix);
      TopMassT[i]->SetDefaultSumw2(); TopMassT[i]= (TH1F*)gDirectory->Get(A.c_str());
      A=Form("UQC%i",NumberSuffix);
      Chi2T[i]->SetDefaultSumw2(); Chi2T[i]= (TH1F*)gDirectory->Get(A.c_str());
    }

  //Range Definition for BKG estimation check
  int THTMax=1600; int THTMin=300; int THTstep=THTL[2]->GetBinWidth(1); int THTNumberOfBins=THTL[2]->GetNbinsX(); //(THTMax-THTMin)/THTstep;
  cout << THTMax << " " << THTMin << " " << THTstep << " " << THTNumberOfBins << THTL[2]->GetNbinsX() << endl;
  
  //Loose and tight number of entries per bin for each sample
  double NLTHT[THTNumberOfBins][NOS]; double NMTHT[THTNumberOfBins][NOS]; double NTTHT[THTNumberOfBins][NOS]; double NTotalLTHT[THTNumberOfBins]; double NTotalMTHT[THTNumberOfBins]; double NTotalTTHT[THTNumberOfBins];
  for (int j=0; j<THTNumberOfBins; j++) {NTotalLTHT[j]=0; NTotalMTHT[j]=0; NTotalTTHT[j]=0;}
  
  float FullIntegralT=0; float FullIntegralM=0; float FullIntegralL=0;

  for (int i=0; i<NOS; i++)
    {
      for (int j=0; j<THTNumberOfBins; j++)
	{
	  NLTHT[j][i]=0; NMTHT[j][i]=0; NTTHT[j][i]=0;
	  NLTHT[j][i]=THTL[i]->Integral(THTL[i]->GetXaxis()->FindBin(THTMin+j*THTstep),THTL[i]->GetXaxis()->FindBin(THTMin+(j+1)*THTstep));
	  NMTHT[j][i]=THTM[i]->Integral(THTM[i]->GetXaxis()->FindBin(THTMin+j*THTstep),THTM[i]->GetXaxis()->FindBin(THTMin+(j+1)*THTstep));
	  NTTHT[j][i]=THTT[i]->Integral(THTT[i]->GetXaxis()->FindBin(THTMin+j*THTstep),THTT[i]->GetXaxis()->FindBin(THTMin+(j+1)*THTstep));
	  NTotalLTHT[j]+=NLTHT[j][i];
	  NTotalMTHT[j]+=NMTHT[j][i];
	  NTotalTTHT[j]+=NTTHT[j][i];
	  //cout << i << j << " " << NLTHT[j][i] << " " << NMTHT[j][i] << " "<< NTTHT[j][i] << " " << NTotalLTHT[j] << " " << NTotalMTHT[j] << " " << NTotalTTHT[j] << endl;
	}
      
      FullIntegralL+=THTL[i]->Integral();
      FullIntegralM+=THTM[i]->Integral();
      FullIntegralT+=THTT[i]->Integral();
    }

  //Efficiencies from B-tag group
  //float ebL=0.8; float ebM=0.67; float ebT=0.42;
  //float ecL=0.6; float ecM=0.15; float ecT=0.01;
  //float elL=0.09; float elM=0.04; float elT=0.001;
  float ebL=0.9; float ebM=0.75; float ebT=0.65;
  float ecL=0.6; float ecM=0.4; float ecT=0.2;
  float elL=0.1; float elM=0.05; float elT=0.01;
  //float ebML=(ebM/ebL)*(0.953/0.987); float ecML=ecM/ecL; float elML=elM/elL; float errorebML=0.02; float errorecML=0.0026; float errorelML=0.0026;
  //float ebTL=(ebT/ebL)*(0.953/0.987); float ecTL=ecT/ecL; float elTL=elT/elL; float errorebTL=0.02; float errorecTL=0.0026; float errorelTL=0.0026;
  //float ebML=0.96591208; float ecML=0.15702889; float elML=0.00454222;
  //float ebTL=0.44658813; float ecTL=0.00042554; float elTL=0.000038543;

  //Matrix of efficiencies
  //Double_t col1[]={elL,elM,elT}; Double_t col2[]={ecL,ecM,ecT}; Double_t col3[]={ebL,ebM,ebT}; TVectorD C1; C1.Use(3,col1); TVectorD C2; C2.Use(3,col2); TVectorD C3; C3.Use(3,col3); 
  TMatrixD eff(3,3); //TMatrixD eff = THilbertMatrixD(3,3); //TMatrixDColumn(eff,0)=C1; TMatrixDColumn(eff,1)=C2; TMatrixDColumn(eff,2)=C3;
  //TArrayD ToFillEff(9); ToFillEff[0]=elL; ToFillEff[1]=ecL; ToFillEff[2]=ebL; ToFillEff[3]=elM; ToFillEff[4]=ecM; ToFillEff[5]=ebM; ToFillEff[6]=elT; ToFillEff[7]=ecT; ToFillEff[8]=ebT; 
  eff[0][0]=elL; eff[0][1]=ecL; eff[0][2]=ebL; eff[1][0]=elM; eff[1][1]=ecM; eff[1][2]=ebM; eff[2][0]=elT; eff[2][1]=ecT; eff[2][2]=ebT;
  //eff.SetMatrixArray(ToFillEff.GetArray());
  TMatrixD effInv = eff; effInv.Invert();

  //Printing matrices
  //eff.Print();
  //effInv.Print();
  //TMatrixD U2(eff,TMatrixD::kMult,effInv); U2.Print();
  TMatrixD eff2(eff,TMatrixD::kMult,eff);
  TMatrixD eff3(eff2,TMatrixD::kMult,eff);
  Double_t Deteff;
  TMatrixD eff3Inv = eff3; eff3Inv.Invert(&Deteff);
  eff3.Print();
  eff3Inv.Print();
  TDecompLU lu(3);
  lu.SetMatrix(eff3);
  lu.Decompose();
  TDecompLU lu2(3);
  lu2.SetMatrix(eff3);
  lu2.Decompose();
  TMatrixD U2(eff3,TMatrixD::kMult,eff3Inv); U2.Print();
  cout << "Element 00 analytically is: " << ((eff3[1][1]*eff3[2][2])-(eff3[1][2]*eff3[2][1]))/Deteff << " with the determinant: " << Deteff << endl; 

  //Vector of L-M-T
  //TMatrixD LMTBinned(3,1);
  TVectorD LMTBinned(3);
  //TMatrixD LMTOverall(3,1);
  TVectorD LMTOverall(3);

  //Signal and Background events from Loose-Tight ratio
  float NumberBTHT[THTNumberOfBins]; float NumberCTHT[THTNumberOfBins]; float NumberLTHT[THTNumberOfBins]; float ErrorNumberBTHT[THTNumberOfBins];  float ErrorNumberCTHT[THTNumberOfBins];float ErrorNumberLTHT[THTNumberOfBins];
  float NBTHT=0; float NCTHT=0; float NLiTHT=0;

  //Estimated Histogram
  TH1F *THTLEstimation= new TH1F("THTBKGEstimation","BKG Estimation for HT", THTNumberOfBins, THTMin, THTMax+THTstep);
  TH1F *THTCEstimation= new TH1F("THTCEstimation","C Estimation for HT", THTNumberOfBins, THTMin, THTMax+THTstep);
  TH1F *THTBEstimation= new TH1F("THTSignalEstimation","Signal Estimation for HT", THTNumberOfBins, THTMin, THTMax+THTstep);

  for (int j=0; j<THTNumberOfBins; j++)
    {
      NumberBTHT[j]=0; NumberCTHT[j]=0; NumberLTHT[j]=0;
      //LMTBinned[0][0]=NTotalLTHT[j]; LMTBinned[1][0]=NTotalMTHT[j]; LMTBinned[2][0]=NTotalTTHT[j];
      LMTBinned[0]=NTotalLTHT[j]; LMTBinned[1]=NTotalMTHT[j]; LMTBinned[2]=NTotalTTHT[j];
      //LMTBinned.Print();
      //TMatrixD Solution(3,1); Solution.Mult(eff3Inv,LMTBinned);
      //TMatrixD Solution(eff3Inv,TMatrixD::kMult,LMTBinned);
      //NumberBTHT[j]=Solution[0][0]; NumberCTHT[j]=Solution[1][0]; NumberLTHT[j]=Solution[2][0];
      //Bool_t ok;
      lu.Solve(LMTBinned);
      //NumberBTHT[j]=LMTBinned[0][0]; NumberCTHT[j]=LMTBinned[1][0]; NumberLTHT[j]=LMTBinned[2][0];
      NumberLTHT[j]=LMTBinned[0]; NumberCTHT[j]=LMTBinned[1]; NumberBTHT[j]=LMTBinned[2];
      cout << "In bin " << j << " Estimated light:" << NumberLTHT[j] << " Estimated C:" << NumberCTHT[j] << " Estimated B:" << NumberBTHT[j] << " Full loose estiamtion:" << NumberLTHT[j]+NumberCTHT[j]+NumberBTHT[j] << endl;
      cout << "In bin " << j << " Real loose:" << NTotalLTHT[j] << " Real medium:" << NTotalMTHT[j] << " Real tight:" << NTotalTTHT[j] << endl;
      ErrorNumberBTHT[j]=0; ErrorNumberCTHT[j]=0; ErrorNumberLTHT[j]=0;
      ErrorNumberLTHT[j]=0.01;
      ErrorNumberCTHT[j]=0.01;
      ErrorNumberBTHT[j]=0.01;
      THTLEstimation->SetBinContent(j,NumberLTHT[j]);
      THTLEstimation->SetBinError(j,ErrorNumberLTHT[j]);
      THTCEstimation->SetBinContent(j,NumberCTHT[j]);
      THTCEstimation->SetBinError(j,ErrorNumberCTHT[j]);
      THTBEstimation->SetBinContent(j,NumberBTHT[j]);
      THTBEstimation->SetBinError(j,ErrorNumberBTHT[j]);
    }
  
  //LMTOverall[0][0]=FullIntegralL; LMTOverall[1][0]=FullIntegralM; LMTOverall[2][0]=FullIntegralT;
  //TMatrixD SolutionO(3,1); SolutionO.Mult(eff3Inv,LMTOverall);
  LMTOverall[0]=FullIntegralL; LMTOverall[1]=FullIntegralM; LMTOverall[2]=FullIntegralT;

  //NLiTHT=SolutionO[0][0]; NCTHT=SolutionO[1][0]; NBTHT=SolutionO[2][0];
  //Bool_t ok;
  lu2.Solve(LMTOverall);
  NLiTHT=LMTOverall[0]; NCTHT=LMTOverall[1]; NBTHT=LMTOverall[2];

  //NLiTHT=(((ecTL-ebTL)/(ecML-ebML))*(FullIntegralM-ebML*FullIntegralL)-(FullIntegralT-ebTL*FullIntegralL))/((elML-ebML)*((ecTL-ebTL)/(ecML-ebML))-(elTL-ebTL));
  //NCTHT=(((elTL-ebTL)/(elML-ebML))*(FullIntegralM-ebML*FullIntegralL)-(FullIntegralT-ebTL*FullIntegralL))/((ecML-ebML)*((elTL-ebTL)/(elML-ebML))-(ecTL-ebTL));
  //NBTHT=FullIntegralL-NLiTHT-NCTHT;
  
  cout << "FOR BINNED ESTIMATION --> Full Integral Light Estimated: " << THTLEstimation->Integral() << " Full Integral C's Estimated: " << THTCEstimation->Integral() << " Full Integral B's Estimated: " << THTBEstimation->Integral() << " Full loose estimation:" << THTLEstimation->Integral()+THTCEstimation->Integral()+THTBEstimation->Integral() << endl;  
  cout << "FOR UNBINNED ESTIMATION --> Full Integral Light Estimated: " << NLiTHT << " Full Integral C's Estimated: " << NCTHT << " Full Integral B's Estimated: " << NBTHT << " Full loose estimation:" << NLiTHT+NCTHT+NBTHT << endl;
  cout << "Full Integral Loose: " << FullIntegralL << " Full Integral Medium: " << FullIntegralM << " Full Integral Tight: " << FullIntegralT << endl;
  
  //Printing Full Matrix Estimation
  cout << "-----Full Matrix Binned Estimation------" << endl;
  TMatrixD BinnedVector(3,1); BinnedVector[0][0]=THTLEstimation->Integral(); BinnedVector[1][0]=THTCEstimation->Integral(); BinnedVector[2][0]=THTBEstimation->Integral();
  TMatrixD BinnedResult(eff3,TMatrixD::kMult,BinnedVector); 
  cout << "  " << BinnedResult[0][0] << "    " << eff3[0][0]*BinnedVector[0][0] << " " << eff3[0][1]*BinnedVector[1][0] << " " << eff3[0][2]*BinnedVector[2][0] << "  " << endl;
  cout << "  " << BinnedResult[1][0] << "  = " << eff3[1][0]*BinnedVector[0][0] << " " << eff3[1][1]*BinnedVector[1][0] << " " << eff3[1][2]*BinnedVector[2][0] << "  " << endl;
  cout << "  " << BinnedResult[2][0] << "    " << eff3[2][0]*BinnedVector[0][0] << " " << eff3[2][1]*BinnedVector[1][0] << " " << eff3[2][2]*BinnedVector[2][0] << "  " << endl;
  
  cout << "-----Full Matrix Unbinned Estimation------" << endl;
  TMatrixD UnBinnedVector(3,1); UnBinnedVector[0][0]=LMTOverall[0]; UnBinnedVector[1][0]=LMTOverall[1]; UnBinnedVector[2][0]=LMTOverall[2];
  TMatrixD UnBinnedResult(eff3,TMatrixD::kMult,UnBinnedVector);
  cout << "  " << UnBinnedResult[0][0] << "    " << eff3[0][0]*LMTOverall[0] << " " << eff3[0][1]*LMTOverall[1] << " " << eff3[0][2]*LMTOverall[2] << "  " << endl;
  cout << "  " << UnBinnedResult[1][0] << "  = " << eff3[1][0]*LMTOverall[0] << " " << eff3[1][1]*LMTOverall[1] << " " << eff3[1][2]*LMTOverall[2] << "  " << endl;
  cout << "  " << UnBinnedResult[2][0] << "    " << eff3[2][0]*LMTOverall[0] << " " << eff3[2][1]*LMTOverall[1] << " " << eff3[2][2]*LMTOverall[2] << "  " << endl; 

  //cout << "Loose    " << ebL*NBTHT+ecL*NCTHT+elL*NLiTHT << " " << ebL*NBTHT << " " << ecL*NCTHT << " " << elL*NLiTHT << endl;
  //cout << "Medium   " << ebM*NBTHT+ecM*NCTHT+elM*NLiTHT << " " << ebM*NBTHT << " " << ecM*NCTHT << " " << elM*NLiTHT << endl;
  //cout << "Tight    " << ebT*NBTHT+ecT*NCTHT+elT*NLiTHT << " " << ebT*NBTHT << " " << ecT*NCTHT << " " << elT*NLiTHT << endl;

  //Stacks of backgrounds
  THStack *BKGTHTL = new THStack("BKGTHT", "BKG for HT; HT GeV; Events");
  for (int i=0; i<NOS-1; i++)
    {
      BKGTHTL->Add(THTL[i]);
    }

  THStack *SIGBKGTHTL = new THStack("SIGBKGTHT", "Sig and BKG for HT; HT GeV; Events");
  for (int i=0; i<NOS; i++)
    {
      SIGBKGTHTL->Add(THTL[i]);
    }

  TLegend* BKGlegend = new TLegend(0.75,0.65,0.90,0.9);
  BKGlegend->AddEntry(THTL[3], "QCD", "f");
  BKGlegend->AddEntry(THTL[2], "TTbar", "f");
  BKGlegend->AddEntry(THTL[1], "SingleT", "f");
  BKGlegend->AddEntry(THTL[0], "Diboson", "f");

  //Plotting!
  char FileName[100];
  sprintf(FileName,"TestingFullEstimationProcedure_SUFFIXLOOSE_SUFFIXMEDIUM_SUFFIXTIGHT.pdf");
  TPDF *ps = new TPDF(FileName,111);
  TCanvas *MyPlot = new TCanvas("MyPlot","Single t prime to top Higgs with backgrounds",600,800);
  TPad *pad1 = new TPad("pad1","pad1",0,0.3,1,1); //Line for ratio plot
  TPaveText *p1;
  p1 = new TPaveText(3,0.5,8,3.5);
  char dateandtime[50];
  sprintf(dateandtime,"date: %s, time: %s",__DATE__,__TIME__);
  p1->AddText(dateandtime);
  //////////////
  //First Page//
  //////////////
  MyPlot->Clear();
  pad1->SetBottomMargin(0); //Line for ratio plot
  pad1->Draw(); //Line for ratio plot
  pad1->cd(); //Line for ratio plot
  //MyPlot->cd(1);
  ps->NewPage();
  gStyle->SetOptStat(0);//Remove the Stat Box
  TH1 *h3=THTL[NOS-1]->DrawCopy("e hist"); //Line for ratio plot
  h3->SetMinimum(-100); //Line for ratio plot
  //THTL[NOS-1]->Draw("e hist");
  THTBEstimation->SetFillColor(kBlue);
  THTBEstimation->Draw("E1 SAME");
  MyPlot->cd(); //Line for ratio plot
  TPad *pad2 = new TPad("pad2","pad2",0,0,1,0.3); //Line for ratio plot
  pad2->SetTopMargin(0); //Line for ratio plot
  pad2->Draw(); //Line for ratio plot
  pad2->cd(); //Line for ratio plot
  THTL[NOS-1]->Sumw2(); //Line for ratio plot
  THTL[NOS-1]->SetStats(0); //Line for ratio plot
  THTL[NOS-1]->Divide(THTBEstimation); //Line for ratio plot
  THTL[NOS-1]->SetMarkerStyle(21); //Line for ratio plot
  THTL[NOS-1]->Draw("ep"); //Line for ratio plot
  MyPlot->cd(); //Line for ratio plot
  gPad->Update();
  MyPlot->Update();
  ///////////////
  //Second Page//
  ///////////////
  MyPlot->Clear();
  MyPlot->cd(1);
  ps->NewPage();
  gStyle->SetOptStat(0);//Remove the Stat Box
  BKGTHTL->Draw("e hist");
  THTLEstimation->SetFillColor(kBlue);
  THTLEstimation->Draw("E1 SAME");
  THTCEstimation->SetFillColor(kYellow);
  THTCEstimation->Draw("E1 SAME");
  BKGlegend->Draw();
  gPad->Update();
  MyPlot->Update();
  //////////////
  ps->Close();


  exit(0); 

}
