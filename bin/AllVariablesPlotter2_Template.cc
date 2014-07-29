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

const TString StorageDirPrefix="file:/sps/cms/ruizalva/CMSSW_5_3_9_patch3/src/Extractors/PatExtractor/bin/WithTrigger/SUFFIX/";
const int NOS=6;
//const TString Samples[NOS] = {"Diboson.root", "SingleTop.root", "TTJets.root", "Zjets.root", "Wjets.root", "QCD.root", "Signal.root"};
const TString Samples[NOS] = {"Diboson.root", "SingleTop.root", "TTJets.root", "DY.root", "QCD.root", "Signal.root"};
//const TString Samples[NOS] = {"Diboson.root", "SingleTop.root", "TTJets.root", "QCD.root", "Signal.root"};

int theNumber(TFile * theFile){
  //TH1F*histogram;
  for(int k=0;k<=19;k++){
    string name="TprimeMass";
    char kchar[20];sprintf(kchar,"%d",k);
    if((TH1F*)(theFile->Get((name+kchar).c_str()))) return k;
  }
  return -1;
}

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
  /*TH1F *Chi2[NOS];
  TH1F *UQuarkContent[NOS];
  TH1F *DQuarkContent[NOS];
  TH1F *SQuarkContent[NOS];
  TH1F *CQuarkContent[NOS];
  TH1F *BQuarkContent[NOS];
  TH1F *CSVLB[NOS];
  TH1F *CSVMB[NOS];
  TH1F *CSVTB[NOS];*/
  TH1F *VTX[NOS];
  //TH2F *HptTpt[NOS];
  /*TH1F *HiggsMassReversedHptTpt[6][NOS];
    TH1F *TprimeMassReversedHptTpt[6][NOS];*/
  
  for (int i=0; i<NOS; i++)
    {
      CurrentFile[i] = new TFile(StorageDirPrefix + Samples[i], "READ");
      if ( CurrentFile[i]->IsOpen() ) printf( Samples[i] + " File opened successfully\n");
      int NumberSuffix=theNumber(CurrentFile[i]);
      cout << NumberSuffix << endl;
      string A=Form("TprimeMass%i",NumberSuffix);
      //cout << A << endl;
      TprimeHistos[i]->SetDefaultSumw2(); TprimeHistos[i]= (TH1F*)gDirectory->Get(A.c_str());
      A=Form("jet1_pt%i",NumberSuffix);
      //cout << A << endl;
      LeadingJetPT[i]->SetDefaultSumw2(); LeadingJetPT[i]= (TH1F*)gDirectory->Get(A.c_str());
      A=Form("jet2_pt%i",NumberSuffix);
      //cout << A << endl;
      Leading2JetPT[i]->SetDefaultSumw2(); Leading2JetPT[i]= (TH1F*)gDirectory->Get(A.c_str());
      A=Form("jet3_pt%i",NumberSuffix);
      //cout << A << endl;
      Leading3JetPT[i]->SetDefaultSumw2(); Leading3JetPT[i]= (TH1F*)gDirectory->Get(A.c_str());
      A=Form("jet4_pt%i",NumberSuffix);
      //cout << A << endl;
      Leading4JetPT[i]->SetDefaultSumw2(); Leading4JetPT[i]= (TH1F*)gDirectory->Get(A.c_str());
      A=Form("jet5_pt%i",NumberSuffix);
      //cout << A << endl;
      Leading5JetPT[i]->SetDefaultSumw2(); Leading5JetPT[i]= (TH1F*)gDirectory->Get(A.c_str());
      A=Form("jet6_pt%i",NumberSuffix);
      //cout << A << endl;
      Leading6JetPT[i]->SetDefaultSumw2(); Leading6JetPT[i]= (TH1F*)gDirectory->Get(A.c_str());
      A=Form("THT%i",NumberSuffix);
      //cout << A << endl;
      THT[i]->SetDefaultSumw2(); THT[i]= (TH1F*)gDirectory->Get(A.c_str());
      A=Form("DeltaR_of_Higgs_Jets%i",NumberSuffix);
      //cout << A << endl;
      DRHjets[i]->SetDefaultSumw2(); DRHjets[i]= (TH1F*)gDirectory->Get(A.c_str());
      A=Form("DeltaR_of_W_Jets%i",NumberSuffix);
      //cout << A << endl;
      DRWjets[i]->SetDefaultSumw2(); DRWjets[i]= (TH1F*)gDirectory->Get(A.c_str());
      A=Form("HPt%i",NumberSuffix);
      //cout << A << endl;
      Hpt[i]->SetDefaultSumw2(); Hpt[i]= (TH1F*)gDirectory->Get(A.c_str());
      A=Form("TPt%i",NumberSuffix);
      //cout << A << endl;
      Tpt[i]->SetDefaultSumw2(); Tpt[i]= (TH1F*)gDirectory->Get(A.c_str());
      A=Form("DeltaR_of_W_Higgs%i",NumberSuffix);
      //cout << A << endl;
      DRWH[i]->SetDefaultSumw2(); DRWH[i]= (TH1F*)gDirectory->Get(A.c_str());
      A=Form("DeltaPhi_of_Higgs_jets%i",NumberSuffix);
      //cout << A << endl;
      DPHjets[i]->SetDefaultSumw2(); DPHjets[i]= (TH1F*)gDirectory->Get(A.c_str());
      A=Form("DeltaPhi_of_W_jets%i",NumberSuffix);
      //cout << A << endl;
      DPWjets[i]->SetDefaultSumw2(); DPWjets[i]= (TH1F*)gDirectory->Get(A.c_str());
      A=Form("DeltaPhi_of_T_jet%i",NumberSuffix);
      //cout << A << endl;
      DPTjets[i]->SetDefaultSumw2(); DPTjets[i]= (TH1F*)gDirectory->Get(A.c_str());
      A=Form("HM%i",NumberSuffix);
      //cout << A << endl;
      HiggsMass[i]->SetDefaultSumw2(); HiggsMass[i]= (TH1F*)gDirectory->Get(A.c_str());
      A=Form("RelHT%i",NumberSuffix);
      //cout << A << endl;
      RelHT[i]->SetDefaultSumw2(); RelHT[i]= (TH1F*)gDirectory->Get(A.c_str());
      A=Form("DeltaR_of_Top_Higgs%i",NumberSuffix);
      //cout << A << endl;
      DRTH[i]->SetDefaultSumw2(); DRTH[i]= (TH1F*)gDirectory->Get(A.c_str());
      A=Form("PT_Normalized_Mass%i",NumberSuffix);
      //cout << A << endl;
      PtNormalizedMass[i]->SetDefaultSumw2(); PtNormalizedMass[i]= (TH1F*)gDirectory->Get(A.c_str());
      A=Form("Relative_Mass%i",NumberSuffix);
      //cout << A << endl;
      RelativeMass[i]->SetDefaultSumw2(); RelativeMass[i]= (TH1F*)gDirectory->Get(A.c_str());
      A=Form("Mother_PT_Normalized_Mass%i",NumberSuffix);
      //cout << A << endl;
      MotherPtNormalizedMass[i]->SetDefaultSumw2(); MotherPtNormalizedMass[i]= (TH1F*)gDirectory->Get(A.c_str());
      A=Form("Number_of_Tops%i",NumberSuffix);
      //cout << A << endl;
      NumberOfTops[i]->SetDefaultSumw2();NumberOfTops[i]= (TH1F*)gDirectory->Get(A.c_str());
      A=Form("HMoverTM%i",NumberSuffix);
      //cout << A << endl;
      HiggsMassOverTopMass[i]->SetDefaultSumw2();HiggsMassOverTopMass[i]= (TH1F*)gDirectory->Get(A.c_str());
      A=Form("HTAsym%i",NumberSuffix);
      //cout << A << endl;
      HiggsTopAsymmetry[i]->SetDefaultSumw2();HiggsTopAsymmetry[i]= (TH1F*)gDirectory->Get(A.c_str());
      A=Form("TLBTa%i",NumberSuffix);
      //cout << A << endl;
      ThirdLooseBtag[i]->SetDefaultSumw2();ThirdLooseBtag[i]= (TH1F*)gDirectory->Get(A.c_str());
      A=Form("TMass%i",NumberSuffix);
      //cout << A << endl;
      TopMass[i]->SetDefaultSumw2(); TopMass[i]= (TH1F*)gDirectory->Get(A.c_str());
      /*A=Form("UQC%i",NumberSuffix);
      //cout << A << endl;
      Chi2[i]->SetDefaultSumw2(); Chi2[i]= (TH1F*)gDirectory->Get(A.c_str());
      A=Form("DQC%i",NumberSuffix);
      //cout << A << endl;
      DQuarkContent[i]->SetDefaultSumw2(); DQuarkContent[i]= (TH1F*)gDirectory->Get(A.c_str());
      A=Form("SQC%i",NumberSuffix);
      //cout << A << endl;
      SQuarkContent[i]->SetDefaultSumw2(); SQuarkContent[i]= (TH1F*)gDirectory->Get(A.c_str());
      A=Form("CQC%i",NumberSuffix);
      //cout << A << endl;
      CQuarkContent[i]->SetDefaultSumw2(); CQuarkContent[i]= (TH1F*)gDirectory->Get(A.c_str());
      A=Form("BQC%i",NumberSuffix);
      //cout << A << endl;
      BQuarkContent[i]->SetDefaultSumw2(); BQuarkContent[i]= (TH1F*)gDirectory->Get(A.c_str());
      A=Form("CSVL%i",NumberSuffix);
      //cout << A << endl;
      CSVLB[i]->SetDefaultSumw2(); CSVLB[i]= (TH1F*)gDirectory->Get(A.c_str());
      A=Form("CSVM%i",NumberSuffix);
      //cout << A << endl;
      CSVMB[i]->SetDefaultSumw2(); CSVMB[i]= (TH1F*)gDirectory->Get(A.c_str());
      A=Form("CSVT%i",NumberSuffix);
      //cout << A << endl;
      CSVTB[i]->SetDefaultSumw2(); CSVTB[i]= (TH1F*)gDirectory->Get(A.c_str());*/
      A=Form("VTX%i",NumberSuffix);
      //cout << A << endl;
      VTX[i]->SetDefaultSumw2(); VTX[i]= (TH1F*)gDirectory->Get(A.c_str());
      /*A=Form("HMBE0%i",NumberSuffix);
      HiggsMassReversedHptTpt[0][i]->SetDefaultSumw2(); HiggsMassReversedHptTpt[0][i]= (TH1F*)gDirectory->Get(A.c_str());
      A=Form("HMBE1%i",NumberSuffix);
      HiggsMassReversedHptTpt[1][i]->SetDefaultSumw2(); HiggsMassReversedHptTpt[1][i]= (TH1F*)gDirectory->Get(A.c_str());
      A=Form("HMBE2%i",NumberSuffix);
      HiggsMassReversedHptTpt[2][i]->SetDefaultSumw2(); HiggsMassReversedHptTpt[2][i]= (TH1F*)gDirectory->Get(A.c_str());
      A=Form("HMBE3%i",NumberSuffix);
      HiggsMassReversedHptTpt[3][i]->SetDefaultSumw2(); HiggsMassReversedHptTpt[3][i]= (TH1F*)gDirectory->Get(A.c_str());
      A=Form("HMBE4%i",NumberSuffix);
      HiggsMassReversedHptTpt[4][i]->SetDefaultSumw2(); HiggsMassReversedHptTpt[4][i]= (TH1F*)gDirectory->Get(A.c_str());
      A=Form("HMBE5%i",NumberSuffix);
      HiggsMassReversedHptTpt[5][i]->SetDefaultSumw2(); HiggsMassReversedHptTpt[5][i]= (TH1F*)gDirectory->Get(A.c_str());
      A=Form("TprimeMassBE0%i",NumberSuffix);
      TprimeMassReversedHptTpt[0][i]->SetDefaultSumw2(); TprimeMassReversedHptTpt[0][i]= (TH1F*)gDirectory->Get(A.c_str());
      A=Form("TprimeMassBE1%i",NumberSuffix);
      TprimeMassReversedHptTpt[1][i]->SetDefaultSumw2(); TprimeMassReversedHptTpt[1][i]= (TH1F*)gDirectory->Get(A.c_str());
      A=Form("TprimeMassBE2%i",NumberSuffix);
      TprimeMassReversedHptTpt[2][i]->SetDefaultSumw2(); TprimeMassReversedHptTpt[2][i]= (TH1F*)gDirectory->Get(A.c_str());
      A=Form("TprimeMassBE3%i",NumberSuffix);
      TprimeMassReversedHptTpt[3][i]->SetDefaultSumw2(); TprimeMassReversedHptTpt[3][i]= (TH1F*)gDirectory->Get(A.c_str());
      A=Form("TprimeMassBE4%i",NumberSuffix);
      TprimeMassReversedHptTpt[4][i]->SetDefaultSumw2(); TprimeMassReversedHptTpt[4][i]= (TH1F*)gDirectory->Get(A.c_str());
      A=Form("TprimeMassBE5%i",NumberSuffix);
      TprimeMassReversedHptTpt[5][i]->SetDefaultSumw2(); TprimeMassReversedHptTpt[5][i]= (TH1F*)gDirectory->Get(A.c_str());*/
      //CurrentFile[i]->Close();
    }

  //Calculation of Estimator
  double bkg[11][NOS-1]; //={0,0,0};
  double sig=0;
  float massStep=50;
  float IntegrationSigma=20;

  sig=TprimeHistos[NOS-1]->Integral(TprimeHistos[NOS-1]->GetXaxis()->FindBin(710),TprimeHistos[NOS-1]->GetXaxis()->FindBin(750));
  for (int i=0; i<10; i++)
    {
      for (int k=0; k<NOS-1; k++)
	{
	  bkg[i][k]=0;
	  bkg[i][k]=TprimeHistos[k]->Integral(TprimeHistos[k]->GetXaxis()->FindBin(550+(i*massStep)-IntegrationSigma),TprimeHistos[k]->GetXaxis()->FindBin(550+(i*massStep)+IntegrationSigma));
	}
    }
  for (int k=0; k<NOS-1; k++)
    {
      bkg[10][k]=0;
      bkg[10][k]=TprimeHistos[k]->Integral(TprimeHistos[k]->GetXaxis()->FindBin(734-IntegrationSigma),TprimeHistos[k]->GetXaxis()->FindBin(734+IntegrationSigma));
    }

  cout << "Estimator under the peak!" << endl;
  cout << "Number of signal events: " << sig << ", Number of BKG events: " << bkg[10][0]+bkg[10][1]+bkg[10][2] << ", S/sqrt(S+B)=" << sig/sqrt(sig+bkg[10][0]+bkg[10][1]+bkg[10][2]) << endl;
  double estimatedSig[10]={0,0,0,0,0,0,0,0,0,0};
  double principalXs=150;
  double Xsections[10]={280,225,179,144,120,98,84,70,61,52};
  double LostEfficiency[4]={0.9*0.9*0.9*0.9,0.9*0.9*0.9,0.9*0.9,0.9};
  for (int i=0; i<10; i++)
    {
      if (Xsections[i]<principalXs) estimatedSig[i]=sig*Xsections[i]/principalXs;
      else estimatedSig[i]=sig*LostEfficiency[i]*Xsections[i]/principalXs;
      cout << "Number of signal events: " << estimatedSig[i] <<  ", and Number of BKG events for " << i << " mass point: " << bkg[i][0]+bkg[i][1]+bkg[i][2] << ", S/sqrt(S+B)=" << estimatedSig[i]/sqrt(estimatedSig[i]+bkg[i][0]+bkg[i][1]+bkg[i][2]) << endl;
    }  
  
  double bkgFR[3]={0,0,0};
  double sigFR=0;

  sigFR=TprimeHistos[NOS-1]->Integral();
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
  /*THStack *BKGChi2= new THStack("BKGChi2", "BKG for Chi2; #Chi^{2}; Events"); 
  THStack *BKGUQC= new THStack("BKGUQC", "BKG for UQC; N_{u}; Events");
  THStack *BKGDQC= new THStack("BKGDQC", "BKG for DQC; N_{d}; Events");
  THStack *BKGSQC= new THStack("BKGSQC", "BKG for SQC; N_{s}; Events");
  THStack *BKGCQC= new THStack("BKGCQC", "BKG for CQC; N_{c}; Events");
  THStack *BKGBQC= new THStack("BKGBQC", "BKG for BQC; N_{b}; Events");
  THStack *BKGCSVL= new THStack("BKGCSVL", "BKG for CSVL; N_{b}^{CSVL}; Events");
  THStack *BKGCSVM= new THStack("BKGCSVM", "BKG for CSVM; N_{b}^{CSVM}; Events");
  THStack *BKGCSVT= new THStack("BKGCSVT", "BKG for CSVT; N_{b}^{CSVT}; Events");*/
  THStack *BKGVTX= new THStack("BKGVTX", "BKG for VTX; N_{vtx}; Events");
  /*THStack *BKGHMBE0= new THStack("BKGHMBE0", "BKG for HMBE0; M_{H}; Events");
  THStack *BKGHMBE1= new THStack("BKGHMBE1", "BKG for HMBE1; M_{H}; Events");
  THStack *BKGHMBE2= new THStack("BKGHMBE2", "BKG for HMBE2; M_{H}; Events");
  THStack *BKGHMBE3= new THStack("BKGHMBE3", "BKG for HMBE3; M_{H}; Events");
  THStack *BKGHMBE4= new THStack("BKGHMBE4", "BKG for HMBE4; M_{H}; Events");
  THStack *BKGHMBE5= new THStack("BKGHMBE5", "BKG for HMBE5; M_{H}; Events");
  THStack *BKGTprimeMassBE0= new THStack("BKGTprimeMassBE0", "BKG for TprimeMassBE0; M_{5j}; Events");
  THStack *BKGTprimeMassBE1= new THStack("BKGTprimeMassBE1", "BKG for TprimeMassBE1; M_{5j}; Events");
  THStack *BKGTprimeMassBE2= new THStack("BKGTprimeMassBE2", "BKG for TprimeMassBE2; M_{5j}; Events");
  THStack *BKGTprimeMassBE3= new THStack("BKGTprimeMassBE3", "BKG for TprimeMassBE3; M_{5j}; Events");
  THStack *BKGTprimeMassBE4= new THStack("BKGTprimeMassBE4", "BKG for TprimeMassBE4; M_{5j}; Events");
  THStack *BKGTprimeMassBE5= new THStack("BKGTprimeMassBE5", "BKG for TprimeMassBE5; M_{5j}; Events");*/
  
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
	  /*BKGChi2->Add(Chi2[i]);
	  BKGUQC->Add(UQuarkContent[i]);
	  BKGDQC->Add(DQuarkContent[i]);
	  BKGSQC->Add(SQuarkContent[i]);
	  BKGCQC->Add(CQuarkContent[i]);
	  BKGBQC->Add(BQuarkContent[i]);
	  BKGCSVL->Add(CSVLB[i]);
	  BKGCSVM->Add(CSVMB[i]);
	  BKGCSVT->Add(CSVTB[i]);*/
	  BKGVTX->Add(VTX[i]);
	  /*BKGHMBE0->Add(HiggsMassReversedHptTpt[0][i]);
	  BKGHMBE1->Add(HiggsMassReversedHptTpt[1][i]);
	  BKGHMBE2->Add(HiggsMassReversedHptTpt[2][i]);
	  BKGHMBE3->Add(HiggsMassReversedHptTpt[3][i]);
	  BKGHMBE4->Add(HiggsMassReversedHptTpt[4][i]);
	  BKGHMBE5->Add(HiggsMassReversedHptTpt[5][i]);
	  BKGTprimeMassBE0->Add(TprimeMassReversedHptTpt[0][i]);
	  BKGTprimeMassBE1->Add(TprimeMassReversedHptTpt[1][i]);
	  BKGTprimeMassBE2->Add(TprimeMassReversedHptTpt[2][i]);
	  BKGTprimeMassBE3->Add(TprimeMassReversedHptTpt[3][i]);
	  BKGTprimeMassBE4->Add(TprimeMassReversedHptTpt[4][i]);
	  BKGTprimeMassBE5->Add(TprimeMassReversedHptTpt[5][i]);*/
	  
	}
    }
  cout << "2 Marker" << endl;

  //TH1F *FakeDY; FakeDY->SetFillColor(kBlue);

  //PostScript Plotting
  cout << "MARKER" << endl;
  TLegend* BKGandSignallegend = new TLegend(0.75,0.65,0.90,0.9);
  TLegend* BKGlegend = new TLegend(0.75,0.65,0.90,0.9);
  if (NOS==5) 
    {
      BKGandSignallegend->AddEntry(TprimeHistos[3], "QCD", "f");
      //BKGandSignallegend->AddEntry(FakeDY, "Zjets", "f");      
      BKGandSignallegend->AddEntry(TprimeHistos[2], "TTbar", "f");
      BKGandSignallegend->AddEntry(TprimeHistos[1], "SingleT", "f");
      BKGandSignallegend->AddEntry(TprimeHistos[0], "Diboson", "f");
      BKGandSignallegend->AddEntry(TprimeHistos[4], "Signal", "f");
      BKGlegend->AddEntry(TprimeHistos[3], "QCD", "f");
      //BKGlegend->AddEntry(FakeDY, "Zjets", "f");      
      BKGlegend->AddEntry(TprimeHistos[2], "TTbar", "f");
      BKGlegend->AddEntry(TprimeHistos[1], "SingleT", "f");
      BKGlegend->AddEntry(TprimeHistos[0], "Diboson", "f");
    }
  else 
    {
      BKGandSignallegend->AddEntry(TprimeHistos[4], "QCD", "f");
      BKGandSignallegend->AddEntry(TprimeHistos[3], "Zjets", "f");
      BKGandSignallegend->AddEntry(TprimeHistos[2], "TTbar", "f");
      BKGandSignallegend->AddEntry(TprimeHistos[1], "SingleT", "f");
      BKGandSignallegend->AddEntry(TprimeHistos[0], "Diboson", "f");
      BKGandSignallegend->AddEntry(TprimeHistos[5], "Signal", "f");
      BKGlegend->AddEntry(TprimeHistos[4], "QCD", "f");
      BKGlegend->AddEntry(TprimeHistos[3], "Zjets", "f");
      BKGlegend->AddEntry(TprimeHistos[2], "TTbar", "f");
      BKGlegend->AddEntry(TprimeHistos[1], "SingleT", "f");
      BKGlegend->AddEntry(TprimeHistos[0], "Diboson", "f");
    }
  cout << "MARKER" << endl;
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

  cout << "3 Marker" << endl;
  //////////////
  //First Page//
  //////////////
  MyPlot->Clear();
  MyPlot->cd(1);
  ps->NewPage();
  gStyle->SetOptStat(0);//Remove the Stat Box
  //TprimeHistos[NOS-1]->Draw("hist");
  //TprimeHistos[NOS-1]->GetYaxis()->SetRangeUser(0.1,BKGandSignal->GetMaximum());
  BKGandSignal->Draw("hist");
  //BKGandSignal->GetYaxis()->SetRangeUser(0.1,BKGandSignal->GetMaximum());
  TprimeHistos[NOS-1]->Draw("histsame");
  gPad->SetLogy();
  //TprimeHistos[NOS-1]->Draw("histsame"); //Preserve double drawing of signal
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
  LeadingJetPT[NOS-1]->Draw("histsame");
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
  Leading2JetPT[NOS-1]->Draw("histsame");
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
  Leading3JetPT[NOS-1]->Draw("histsame");
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
  Leading4JetPT[NOS-1]->Draw("histsame");
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
  Leading5JetPT[NOS-1]->Draw("histsame");
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
  Leading6JetPT[NOS-1]->Draw("histsame");
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
  THT[NOS-1]->Draw("histsame");
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
  DRHjets[NOS-1]->Draw("histsame");
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
  DRWjets[NOS-1]->Draw("histsame");
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
  Hpt[NOS-1]->Draw("histsame");
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
  Tpt[NOS-1]->Draw("histsame");
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
  DRWH[NOS-1]->Draw("histsame");
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
  DPHjets[NOS-1]->Draw("histsame");
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
  DPWjets[NOS-1]->Draw("histsame");
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
  DPTjets[NOS-1]->Draw("histsame");
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
  HiggsMass[NOS-1]->Draw("histsame");
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
  RelHT[NOS-1]->Draw("histsame");
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
  DRTH[NOS-1]->Draw("histsame");
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
  PtNormalizedMass[NOS-1]->Draw("histsame");
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
  RelativeMass[NOS-1]->Draw("histsame");
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
  MotherPtNormalizedMass[NOS-1]->Draw("histsame");
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
  NumberOfTops[NOS-1]->Draw("histsame");
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
  HiggsMassOverTopMass[NOS-1]->Draw("histsame");
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
  HiggsTopAsymmetry[NOS-1]->Draw("histsame");
  //gPad->SetLogy();
  BKGHTAsym->SetMinimum(0.1);
  BKGandSignallegend->Draw();
  gPad->RedrawAxis("");
  gPad->Update();
  MyPlot->Update();
  //////////////
  //////////////
  //////////////
  /*MyPlot->Clear();
  MyPlot->cd(1);
  ps->NewPage();
  BKGTLBT->Draw("hist");
  ThirdLooseBtag[NOS-1]->Draw("histsame");
  //gPad->SetLogy();
  //BKGTLBT->SetMinimum(0.1);
  BKGandSignallegend->Draw();
  gPad->RedrawAxis("");
  gPad->Update();
  MyPlot->Update();
  cout << "Here 1" << endl;*/
  //////////////
  //////////////
  //////////////
  MyPlot->Clear();
  MyPlot->cd(1);
  gPad->Clear();
  ps->NewPage();
  BKGTopMass->Draw("hist");
  TopMass[NOS-1]->Draw("histsame");
  //gPad->SetLogy();
  //BKGTopMass->SetMinimum(0.1);
  BKGandSignallegend->Draw();
  gPad->RedrawAxis("");
  gPad->Update();
  MyPlot->Update();
  cout << "Here 1" << endl;  
  //////////////
  //////////////
  //////////////
  /*MyPlot->Clear();
  MyPlot->cd(1);
  ps->NewPage();
  BKGChi2->Draw("hist");
  Chi2[NOS-1]->Draw("histsame");
  //gPad->SetLogy();
  BKGChi2->SetMinimum(0.1);
  BKGandSignallegend->Draw();
  gPad->RedrawAxis("");
  gPad->Update();
  MyPlot->Update();*/
  //////////////
  //////////////
  //////////////
  /*MyPlot->Clear();
  MyPlot->cd(1);
  ps->NewPage();
  BKGUQC->Draw("hist");
  UQuarkContent[NOS-1]->Draw("histsame");
  //gPad->SetLogy();
  BKGUQC->SetMinimum(0.1);
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
  BKGDQC->Draw("hist");
  DQuarkContent[NOS-1]->Draw("histsame");
  //gPad->SetLogy();
  BKGDQC->SetMinimum(0.1);
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
  BKGSQC->Draw("hist");
  SQuarkContent[NOS-1]->Draw("histsame");
  //gPad->SetLogy();
  BKGSQC->SetMinimum(0.1);
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
  BKGCQC->Draw("hist");
  CQuarkContent[NOS-1]->Draw("histsame");
  //gPad->SetLogy();
  BKGCQC->SetMinimum(0.1);
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
  BKGBQC->Draw("hist");
  BQuarkContent[NOS-1]->Draw("histsame");
  //gPad->SetLogy();
  BKGBQC->SetMinimum(0.1);
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
  BKGCSVL->Draw("hist");
  CSVLB[NOS-1]->Draw("histsame");
  //gPad->SetLogy();
  BKGCSVL->SetMinimum(0.1);
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
  BKGCSVM->Draw("hist");
  CSVMB[NOS-1]->Draw("histsame");
  //gPad->SetLogy();
  BKGCSVM->SetMinimum(0.1);
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
  BKGCSVT->Draw("hist");
  CSVTB[NOS-1]->Draw("histsame");
  //gPad->SetLogy();
  BKGCSVT->SetMinimum(0.1);
  BKGandSignallegend->Draw();
  gPad->RedrawAxis("");
  gPad->Update();
  MyPlot->Update();*/
  //////////////
  //////////////
  //////////////
  MyPlot->Clear();
  MyPlot->cd(1);
  ps->NewPage();
  BKGVTX->Draw("hist");
  VTX[NOS-1]->Draw("histsame");
  //gPad->SetLogy();
  BKGVTX->SetMinimum(0.1);
  BKGandSignallegend->Draw();
  gPad->RedrawAxis("");
  gPad->Update();
  MyPlot->Update();
  //////////////
  //////////////
  //////////////
  /*MyPlot->Clear();
  MyPlot->cd(1);
  ps->NewPage();
  BKGHMBE0->Draw("hist");
  //gPad->SetLogy();
  BKGHMBE0->SetMinimum(0.1);
  BKGlegend->Draw();
  gPad->RedrawAxis("");
  gPad->Update();
  MyPlot->Update();
  //////////////
  //////////////
  //////////////
  MyPlot->Clear();
  MyPlot->cd(1);
  ps->NewPage();
  BKGHMBE1->Draw("hist");
  //gPad->SetLogy();
  BKGHMBE1->SetMinimum(0.1);
  BKGlegend->Draw();
  gPad->RedrawAxis("");
  gPad->Update();
  MyPlot->Update();
  //////////////
  //////////////
  //////////////
  MyPlot->Clear();
  MyPlot->cd(1);
  ps->NewPage();
  BKGHMBE2->Draw("hist");
  //gPad->SetLogy();
  BKGHMBE2->SetMinimum(0.1);
  BKGlegend->Draw();
  gPad->RedrawAxis("");
  gPad->Update();
  MyPlot->Update();
  //////////////
  //////////////
  //////////////
  MyPlot->Clear();
  MyPlot->cd(1);
  ps->NewPage();
  BKGHMBE3->Draw("hist");
  //gPad->SetLogy();
  BKGHMBE3->SetMinimum(0.1);
  BKGlegend->Draw();
  gPad->RedrawAxis("");
  gPad->Update();
  MyPlot->Update();
  //////////////
  //////////////
  //////////////
  MyPlot->Clear();
  MyPlot->cd(1);
  ps->NewPage();
  BKGHMBE4->Draw("hist");
  //gPad->SetLogy();
  BKGHMBE4->SetMinimum(0.1);
  BKGlegend->Draw();
  gPad->RedrawAxis("");
  gPad->Update();
  MyPlot->Update();
  //////////////
  //////////////
  //////////////
  MyPlot->Clear();
  MyPlot->cd(1);
  ps->NewPage();
  BKGHMBE5->Draw("hist");
  //gPad->SetLogy();
  BKGHMBE5->SetMinimum(0.1);
  BKGlegend->Draw();
  gPad->RedrawAxis("");
  gPad->Update();
  MyPlot->Update();
  //////////////
  //////////////
  //////////////
  MyPlot->Clear();
  MyPlot->cd(1);
  ps->NewPage();
  BKGTprimeMassBE0->Draw("hist");
  //gPad->SetLogy();
  BKGTprimeMassBE0->SetMinimum(0.1);
  BKGlegend->Draw();
  gPad->RedrawAxis("");
  gPad->Update();
  MyPlot->Update();
  //////////////
  //////////////
  //////////////
  MyPlot->Clear();
  MyPlot->cd(1);
  ps->NewPage();
  BKGTprimeMassBE1->Draw("hist");
  //gPad->SetLogy();
  BKGTprimeMassBE1->SetMinimum(0.1);
  BKGlegend->Draw();
  gPad->RedrawAxis("");
  gPad->Update();
  MyPlot->Update();
  //////////////
  //////////////
  //////////////
  MyPlot->Clear();
  MyPlot->cd(1);
  ps->NewPage();
  BKGTprimeMassBE2->Draw("hist");
  //gPad->SetLogy();
  BKGTprimeMassBE2->SetMinimum(0.1);
  BKGlegend->Draw();
  gPad->RedrawAxis("");
  gPad->Update();
  MyPlot->Update();
  //////////////
  //////////////
  //////////////
  MyPlot->Clear();
  MyPlot->cd(1);
  ps->NewPage();
  BKGTprimeMassBE3->Draw("hist");
  //gPad->SetLogy();
  BKGTprimeMassBE3->SetMinimum(0.1);
  BKGlegend->Draw();
  gPad->RedrawAxis("");
  gPad->Update();
  MyPlot->Update();
  //////////////
  //////////////
  //////////////
  MyPlot->Clear();
  MyPlot->cd(1);
  ps->NewPage();
  BKGTprimeMassBE4->Draw("hist");
  //gPad->SetLogy();
  BKGTprimeMassBE4->SetMinimum(0.1);
  BKGlegend->Draw();
  gPad->RedrawAxis("");
  gPad->Update();
  MyPlot->Update();
  //////////////
  //////////////
  //////////////
  MyPlot->Clear();
  MyPlot->cd(1);
  ps->NewPage();
  BKGTprimeMassBE5->Draw("hist");
  //gPad->SetLogy();
  BKGTprimeMassBE5->SetMinimum(0.1);
  BKGlegend->Draw();
  gPad->RedrawAxis("");
  gPad->Update();
  MyPlot->Update();*/
  //////////////
  ps->Close();

  exit(0);
}
