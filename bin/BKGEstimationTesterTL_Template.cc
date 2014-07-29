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

const TString StorageDirPrefixLoose="file:/afs/cern.ch/work/j/jruizalv/private/Analysis/CMSSW_5_3_9_patch3/src/Extractors/PatExtractor/bin/WithTrigger/SUFFIXLOOSE/";
const TString StorageDirPrefixTight="file:/afs/cern.ch/work/j/jruizalv/private/Analysis/CMSSW_5_3_9_patch3/src/Extractors/PatExtractor/bin/WithTrigger/SUFFIXTIGHT/";
const int NOS=5;
//const TString Samples[NOS] = {"Diboson.root", "SingleTop.root", "TTJets.root", "Zjets.root", "Wjets.root", "QCD.root", "Signal.root"};
const TString Samples[NOS] = {"Diboson.root", "SingleTop.root", "TTJets.root", "QCD.root", "Signal.root"};

int theNumber(TFile * theFile){
  //TH1F*histogram;
  for(int k=0;k<=20;k++){
    string name="TprimeMass";
    char kchar[20];sprintf(kchar,"%d",k);
    if((TH1F*)(theFile->Get((name+kchar).c_str()))) return k;
  }
  return -1;
}

void BKGEstimationTesterTL()
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
  double NLTHT[THTNumberOfBins][NOS]; double NTTHT[THTNumberOfBins][NOS]; double NTotalLTHT[THTNumberOfBins]; double NTotalTTHT[THTNumberOfBins];
  double SLTHT[THTNumberOfBins]; double STTHT[THTNumberOfBins]; double BKGLTHT[THTNumberOfBins]; double BKGTTHT[THTNumberOfBins];
  for (int j=0; j<THTNumberOfBins; j++) {NTotalLTHT[j]=0; NTotalTTHT[j]=0; SLTHT[j]=0; STTHT[j]=0; BKGLTHT[j]=0; BKGTTHT[j]=0;}
  
  float FullIntegralT=0; float FullIntegralL=0; float FullIntegralST=0; float FullIntegralSL=0; float FullIntegralBKGT=0; float FullIntegralBKGL=0;

  for (int i=0; i<NOS; i++)
    {
      for (int j=0; j<THTNumberOfBins; j++)
	{
	  NLTHT[j][i]=0; NTTHT[j][i]=0;
	  NLTHT[j][i]=THTL[i]->Integral(THTL[i]->GetXaxis()->FindBin(THTMin+j*THTstep),THTL[i]->GetXaxis()->FindBin(THTMin+(j+1)*THTstep));
	  NTTHT[j][i]=THTT[i]->Integral(THTT[i]->GetXaxis()->FindBin(THTMin+j*THTstep),THTT[i]->GetXaxis()->FindBin(THTMin+(j+1)*THTstep));
	  if (i==NOS-1) {SLTHT[j]=NLTHT[j][i]; STTHT[j]=NTTHT[j][i];}
	  else {BKGLTHT[j]+=NLTHT[j][i]; BKGTTHT[j]+=NTTHT[j][i];}
	  NTotalLTHT[j]+=NLTHT[j][i];
	  NTotalTTHT[j]+=NTTHT[j][i];
	  cout << i << j << " " << NLTHT[j][i] << " " << NTTHT[j][i] << " " << NTotalLTHT[j] << " " << NTotalTTHT[j] << endl;
	}
      if (i==NOS-1) {FullIntegralSL=THTL[i]->Integral(); FullIntegralST=THTT[i]->Integral();}
      else {FullIntegralBKGL+=THTL[i]->Integral(); FullIntegralBKGT+=THTT[i]->Integral();}
      FullIntegralL+=THTL[i]->Integral();
      FullIntegralT+=THTT[i]->Integral();
    }
  
  cout << "Full Integral Loose: " << FullIntegralL << " Full Integral Tight: " << FullIntegralT << endl;

  //Efficiencies from B-tag group
  //float ebL=0.8; float ebM=0.65; float ebT=0.46;
  //float ecL=0.6; float ecM=0.1; float ecT=0.01;
  //float elL=0.09; float elM=0.01; float elT=0.001;
  //float ebML=(ebM/ebL)*(0.953/0.987); float ecML=ecM/ecL; float elML=elM/elL; float errorebML=0.02; float errorecML=0.0026; float errorelML=0.0026;
  //float ebTL=(ebT/ebL)*(0.953/0.987); float ecTL=ecT/ecL; float elTL=elT/elL; float errorebTL=0.02; float errorecTL=0.0026; float errorelTL=0.0026;
  //float ebL=0.7; float elL=0.3; 
  //float ebM=0.3; float elM=0.04;
  //float eb=ebM/ebL; float el=elM/elL;
  //float eb=0.494; float el=0.109; //-------> For 3L 3M versions
  float eb=0.88; float el=0.288;

  //Signal and Background events from Loose-Tight ratio
  float NumberSTHT[THTNumberOfBins]; float NumberBKGTHT[THTNumberOfBins]; float ErrorNumberSTHT[THTNumberOfBins]; float ErrorNumberBKGTHT[THTNumberOfBins];

  //Estimated Histogram
  TH1F *THTBKGEstimation= new TH1F("THTBKGEstimation","BKG Estimation for HT", THTNumberOfBins, THTMin, THTMax+THTstep);
  TH1F *THTSignalEstimation= new TH1F("THTSignalEstimation","Signal Estimation for HT", THTNumberOfBins, THTMin, THTMax+THTstep);
  float NSTHT=0; float NBTHT=0;

  for (int j=0; j<THTNumberOfBins; j++)
    {
      NumberSTHT[j]=0; NumberBKGTHT[j]=0;
      NumberSTHT[j]=(NTotalTTHT[j]-el*NTotalLTHT[j])/(eb-el);
      NumberBKGTHT[j]=(eb*NTotalLTHT[j]-NTotalTTHT[j])/(eb-el);
      cout << "In bin " << j << " Estimated BKG:" << NumberBKGTHT[j] << " Estimated signal:" << NumberSTHT[j] << endl;
      cout << "Efficiency on signal: " << STTHT[j]/SLTHT[j] << ", Efficiency on bkgs: " << BKGTTHT[j]/BKGLTHT[j] << endl;
      ErrorNumberSTHT[j]=0; ErrorNumberBKGTHT[j]=0;
      ErrorNumberBKGTHT[j]=0.01; //NumberBKGTHT[j]*(((eb*NTotalLTHT[j]*((erroreb/eb)+(sqrt(NTotalLTHT[j])/NTotalLTHT[j]))+sqrt(NTotalTTHT[j]))/(NTotalTTHT[j]-eb*NTotalLTHT[j]))+((erroreb+errorel)/(el-eb)));
      ErrorNumberSTHT[j]=0.01; //NumberSTHT[j]*((((2*eb-el)*NTotalLTHT[j]*(((2*erroreb+errorel)/(2*eb-el))+(sqrt(NTotalLTHT[j])/NTotalLTHT[j]))+sqrt(NTotalTTHT[j]))/(NTotalLTHT[j]*(2*eb-el)-NTotalTTHT[j]))+((erroreb+errorel)/(eb-el)));
      THTBKGEstimation->SetBinContent(j+1,NumberBKGTHT[j]);
      THTBKGEstimation->SetBinError(j+1,ErrorNumberBKGTHT[j]);
      THTSignalEstimation->SetBinContent(j+1,NumberSTHT[j]);
      THTSignalEstimation->SetBinError(j+1,ErrorNumberSTHT[j]);
    }

  NSTHT=(FullIntegralT-el*FullIntegralL)/(eb-el);
  NBTHT=(eb*FullIntegralL-FullIntegralT)/(eb-el);

  cout << "FOR BINNED ESTIMATION --> Full Integral BKG Etstimated: " << THTBKGEstimation->Integral() << " Full Integral Signal Estimated: " << THTSignalEstimation->Integral() << endl;
  cout << "FOR UNBINNED ESTIMATION --> Full Integral BKG Etstimated: " << NBTHT << " Full Integral Signal Estimated: " << NSTHT << endl;

  cout << "Efficiency on signal full: " << FullIntegralST/FullIntegralSL << ", Efficiency on bkgs full: " << FullIntegralBKGT/FullIntegralBKGL << endl;

  cout << "Full Integral Loose: " << FullIntegralL << " Full Integral Tight: " << FullIntegralT << endl;

  cout << "-----Full Matrix Estimation------" << endl;
  cout << "Loose    " << NBTHT+NSTHT << " " << NSTHT << " " << NBTHT << endl;
  cout << "Tight    " << eb*NSTHT+el*NBTHT << " " << eb*NSTHT << " " << el*NBTHT << endl;

  cout << "-----Full Matrix Real------------" << endl;
  cout << "Loose    " << FullIntegralL << " " << FullIntegralSL << " " << FullIntegralBKGL << endl;
  cout << "Tight    " << FullIntegralT << " " << FullIntegralST << " " << FullIntegralBKGT << endl;

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
  sprintf(FileName,"TestingEstimationProcedure_SUFFIXLOOSE_SUFFIXTIGHT.pdf");
  TPDF *ps = new TPDF(FileName,111);
  TCanvas *MyPlot = new TCanvas("MyPlot","Single t prime to top Higgs with backgrounds",600,800);
  //TPad *pad1 = new TPad("pad1","pad1",0,0.3,1,1); //Line for ratio plot
  TPaveText *p1;
  p1 = new TPaveText(3,0.5,8,3.5);
  char dateandtime[50];
  sprintf(dateandtime,"date: %s, time: %s",__DATE__,__TIME__);
  p1->AddText(dateandtime);
  //////////////
  //First Page//
  //////////////
  MyPlot->Clear();
  //pad1->SetBottomMargin(0); //Line for ratio plot
  //pad1->Draw(); //Line for ratio plot
  //pad1->cd(); //Line for ratio plot
  MyPlot->cd(1);
  ps->NewPage();
  gStyle->SetOptStat(0);//Remove the Stat Box
  //TH1 *h3=THTL[NOS-1]->DrawCopy("e hist"); //Line for ratio plot
  //h3->SetMinimum(-100); //Line for ratio plot
  THTL[NOS-1]->Draw("e hist");
  THTSignalEstimation->SetFillColor(kBlue);
  THTSignalEstimation->Draw("E1 SAME");
  //MyPlot->cd(); //Line for ratio plot
  //TPad *pad2 = new TPad("pad2","pad2",0,0,1,0.3); //Line for ratio plot
  //pad2->SetTopMargin(0); //Line for ratio plot
  //pad2->Draw(); //Line for ratio plot
  //pad2->cd(); //Line for ratio plot
  //THTL[NOS-1]->Sumw2(); //Line for ratio plot
  //THTL[NOS-1]->SetStats(0); //Line for ratio plot
  //THTL[NOS-1]->Divide(THTBEstimation); //Line for ratio plot
  //THTL[NOS-1]->SetMarkerStyle(21); //Line for ratio plot
  //THTL[NOS-1]->Draw("ep"); //Line for ratio plot
  //MyPlot->cd(); //Line for ratio plot
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
  THTBKGEstimation->SetFillColor(kBlue);
  THTBKGEstimation->Draw("E1 SAME");
  BKGlegend->Draw();
  gPad->Update();
  MyPlot->Update();
  //////////////
  ps->Close();


  exit(0); 

}
