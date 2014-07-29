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
#include "PU_Reweighting.h"

using namespace std;

const TString StorageDirPrefixMC="file:/sps/cms/ruizalva/CMSSW_5_3_9_patch3/src/Extractors/PatExtractor/bin/WithTrigger/SUFFIXMC/";
const TString StorageDirPrefixDATA="file:/sps/cms/ruizalva/CMSSW_5_3_9_patch3/src/Extractors/PatExtractor/bin/WithTrigger/SUFFIXDATA/";
const int NOS=6;
//const TString Samples[NOS] = {"Diboson.root", "SingleTop.root", "TTJets.root", "Zjets.root", "Wjets.root", "QCD.root", "Signal.root"};
const TString Samples[NOS] = {"Diboson.root", "SingleTop.root", "TTJets.root", "DY.root", "QCD.root", "Signal.root"};
const int NOD=1;
const TString DATASamples[NOD] = {"Data.root"};

int theNumber(TFile * theFile){
  //TH1F*histogram;
  for(int k=0;k<=20;k++){
    string name="TprimeMass";
    char kchar[20];sprintf(kchar,"%d",k);
    if((TH1F*)(theFile->Get((name+kchar).c_str()))) return k;
  }
  return -1;
}

//const float ActualLumi = 397+16.888+119.575+273.700+215.538+155.678+24.291; // inverse pb
const float ActualLumi = 19694.471;

void DataMCComparison()
{
  int ReBin=2;

  double THTStackIntegral=0;
  double VtcsStackIntegral=0;

  //Reading files for MC samples
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
  TH1F *VtcsL[NOS];
  TH1F *TopFromHiggsL[NOS];
  TH1F *JET1ETAL[NOS];
  TH1F *JET2ETAL[NOS];
  TH1F *JET3ETAL[NOS]; 
  TH1F *JET4ETAL[NOS];
  TH1F *JET5ETAL[NOS];
  TH1F *JET6ETAL[NOS];
  TH1F *JET1PHIL[NOS];   
  TH1F *JET2PHIL[NOS]; 
  TH1F *JET3PHIL[NOS]; 
  TH1F *JET4PHIL[NOS];
  TH1F *JET5PHIL[NOS];
  TH1F *JET6PHIL[NOS];
  TH1F *JETMULTIL[NOS];

  for (int i=0; i<NOS; i++)
    {
      CurrentFileL[i] = new TFile(StorageDirPrefixMC + Samples[i], "READ");
      if ( CurrentFileL[i]->IsOpen() ) printf( Samples[i] + " File opened successfully\n");      
      int NumberSuffix=theNumber(CurrentFileL[i]);
      ///////////////M5J//////////////
      string A=Form("TprimeMass%i",NumberSuffix);
      TprimeHistosL[i]->SetDefaultSumw2(); TprimeHistosL[i]= (TH1F*)gDirectory->Get(A.c_str()); TprimeHistosL[i]->Scale(ActualLumi/Lumi); TprimeHistosL[i]->Rebin(ReBin);
      //if (i==0) {AllBKGSTPM->SetDefaultSumw2(); AllBKGSTPM= (TH1F*)TprimeHistosL[i]->Clone("ALLBKGSTPM");}
      ///////////////LJPT//////////////
      A=Form("jet1_pt%i",NumberSuffix);
      LeadingJetPTL[i]->SetDefaultSumw2(); LeadingJetPTL[i]= (TH1F*)gDirectory->Get(A.c_str()); /*THTStackIntegral+=THTL[i]->Integral();*/ LeadingJetPTL[i]->Scale(ActualLumi/Lumi); LeadingJetPTL[i]->Rebin(ReBin);
      //if (i==0) {AllBKGSLJPT->SetDefaultSumw2(); AllBKGSLJPT= (TH1F*)LeadingJetPTL[i]->Clone("ALLBKGSLJPT");}
      ///////////////SLJPT//////////////
      A=Form("jet2_pt%i",NumberSuffix);
      Leading2JetPTL[i]->SetDefaultSumw2(); Leading2JetPTL[i]= (TH1F*)gDirectory->Get(A.c_str());
      ///////////////SSLJPT//////////////
      A=Form("jet3_pt%i",NumberSuffix);
      Leading3JetPTL[i]->SetDefaultSumw2(); Leading3JetPTL[i]= (TH1F*)gDirectory->Get(A.c_str());
      ///////////////SSSLJPT//////////////
      A=Form("jet4_pt%i",NumberSuffix);
      Leading4JetPTL[i]->SetDefaultSumw2(); Leading4JetPTL[i]= (TH1F*)gDirectory->Get(A.c_str());
      ///////////////SSSSLJPT//////////////
      A=Form("jet5_pt%i",NumberSuffix);
      Leading5JetPTL[i]->SetDefaultSumw2(); Leading5JetPTL[i]= (TH1F*)gDirectory->Get(A.c_str());
      ///////////////SSSSSLJPT//////////////
      A=Form("jet6_pt%i",NumberSuffix);
      Leading6JetPTL[i]->SetDefaultSumw2(); Leading6JetPTL[i]= (TH1F*)gDirectory->Get(A.c_str());
      ///////////////THT//////////////
      A=Form("THT%i",NumberSuffix);
      THTL[i]->SetDefaultSumw2(); THTL[i]= (TH1F*)gDirectory->Get(A.c_str()); THTStackIntegral+=THTL[i]->Integral(); THTL[i]->Scale(ActualLumi/Lumi); THTL[i]->Rebin(ReBin);
      //if (i==0) {AllBKGSTHT->SetDefaultSumw2(); AllBKGSTHT= (TH1F*)THTL[i]->Clone("ALLBKGSTHT");} //(TH1F*)gDirectory->Get(A.c_str());}
      ///////////////DRHJ//////////////
      A=Form("DeltaR_of_Higgs_Jets%i",NumberSuffix);
      DRHjetsL[i]->SetDefaultSumw2(); DRHjetsL[i]= (TH1F*)gDirectory->Get(A.c_str()); DRHjetsL[i]->Scale(ActualLumi/Lumi); DRHjetsL[i]->Rebin(ReBin);
      //if (i==0) {AllBKGSDRHJ->SetDefaultSumw2(); AllBKGSDRHJ= (TH1F*)DRHjetsL[i]->Clone("ALLBKGSDRHJ");}
      ///////////////DRWJ//////////////
      A=Form("DeltaR_of_W_Jets%i",NumberSuffix);
      DRWjetsL[i]->SetDefaultSumw2(); DRWjetsL[i]= (TH1F*)gDirectory->Get(A.c_str()); DRWjetsL[i]->Rebin(ReBin);
      ///////////////HPT//////////////
      A=Form("HPt%i",NumberSuffix);
      HptL[i]->SetDefaultSumw2(); HptL[i]= (TH1F*)gDirectory->Get(A.c_str()); HptL[i]->Scale(ActualLumi/Lumi); HptL[i]->Rebin(ReBin);
      ///////////////TPT//////////////
      A=Form("TPt%i",NumberSuffix);
      TptL[i]->SetDefaultSumw2(); TptL[i]= (TH1F*)gDirectory->Get(A.c_str()); TptL[i]->Scale(ActualLumi/Lumi); TptL[i]->Rebin(ReBin);
      ///////////////DRWH//////////////
      A=Form("DeltaR_of_W_Higgs%i",NumberSuffix);
      DRWHL[i]->SetDefaultSumw2(); DRWHL[i]= (TH1F*)gDirectory->Get(A.c_str()); DRWHL[i]->Scale(ActualLumi/Lumi); DRWHL[i]->Rebin(ReBin);
      ///////////////DPHJ//////////////
      A=Form("DeltaPhi_of_Higgs_jets%i",NumberSuffix);
      DPHjetsL[i]->SetDefaultSumw2(); DPHjetsL[i]= (TH1F*)gDirectory->Get(A.c_str()); DPHjetsL[i]->Scale(ActualLumi/Lumi); DPHjetsL[i]->Rebin(ReBin);
      ///////////////DPWJ//////////////
      A=Form("DeltaPhi_of_W_jets%i",NumberSuffix);
      DPWjetsL[i]->SetDefaultSumw2(); DPWjetsL[i]= (TH1F*)gDirectory->Get(A.c_str()); DPWjetsL[i]->Scale(ActualLumi/Lumi); DPWjetsL[i]->Rebin(ReBin);
      ///////////////DPTJ//////////////
      A=Form("DeltaPhi_of_T_jet%i",NumberSuffix);
      DPTjetsL[i]->SetDefaultSumw2(); DPTjetsL[i]= (TH1F*)gDirectory->Get(A.c_str()); DPTjetsL[i]->Scale(ActualLumi/Lumi); DPTjetsL[i]->Rebin(ReBin);
      ///////////////HM//////////////
      A=Form("HM%i",NumberSuffix);
      HiggsMassL[i]->SetDefaultSumw2(); HiggsMassL[i]= (TH1F*)gDirectory->Get(A.c_str()); HiggsMassL[i]->Scale(ActualLumi/Lumi); HiggsMassL[i]->Rebin(ReBin);
      ///////////////RelHT//////////////
      A=Form("RelHT%i",NumberSuffix);
      RelHTL[i]->SetDefaultSumw2(); RelHTL[i]= (TH1F*)gDirectory->Get(A.c_str()); RelHTL[i]->Scale(ActualLumi/Lumi); RelHTL[i]->Rebin(ReBin);
      ///////////////DRTH//////////////
      A=Form("DeltaR_of_Top_Higgs%i",NumberSuffix);
      DRTHL[i]->SetDefaultSumw2(); DRTHL[i]= (TH1F*)gDirectory->Get(A.c_str());
      /////////////////////////////
      A=Form("PT_Normalized_Mass%i",NumberSuffix);
      PtNormalizedMassL[i]->SetDefaultSumw2(); PtNormalizedMassL[i]= (TH1F*)gDirectory->Get(A.c_str());
      ///////////////RM//////////////
      A=Form("Relative_Mass%i",NumberSuffix);
      RelativeMassL[i]->SetDefaultSumw2(); RelativeMassL[i]= (TH1F*)gDirectory->Get(A.c_str()); RelativeMassL[i]->Scale(ActualLumi/Lumi); RelativeMassL[i]->Rebin(ReBin);
      /////////////////////////////
      A=Form("Mother_PT_Normalized_Mass%i",NumberSuffix);
      MotherPtNormalizedMassL[i]->SetDefaultSumw2(); MotherPtNormalizedMassL[i]= (TH1F*)gDirectory->Get(A.c_str());
      /////////////////////////////
      A=Form("Number_of_Tops%i",NumberSuffix);
      NumberOfTopsL[i]->SetDefaultSumw2();NumberOfTopsL[i]= (TH1F*)gDirectory->Get(A.c_str());
      /////////////////////////////
      A=Form("HMoverTM%i",NumberSuffix);
      HiggsMassOverTopMassL[i]->SetDefaultSumw2();HiggsMassOverTopMassL[i]= (TH1F*)gDirectory->Get(A.c_str());
      /////////////////////////////
      A=Form("HTAsym%i",NumberSuffix);
      HiggsTopAsymmetryL[i]->SetDefaultSumw2();HiggsTopAsymmetryL[i]= (TH1F*)gDirectory->Get(A.c_str());
      /////////////////////////////
      A=Form("TLBTag%i",NumberSuffix);
      HiggsTopAsymmetryL[i]->SetDefaultSumw2();HiggsTopAsymmetryL[i]= (TH1F*)gDirectory->Get(A.c_str());
      ///////////////TopM//////////////
      A=Form("TMass%i",NumberSuffix);
      TopMassL[i]->SetDefaultSumw2();TopMassL[i]= (TH1F*)gDirectory->Get(A.c_str()); TopMassL[i]->Scale(ActualLumi/Lumi); TopMassL[i]->Rebin(ReBin);
      /////////////////////////////
      A=Form("ChiSq%i",NumberSuffix);
      //TopMassL[i]->SetDefaultSumw2(); TopMassL[i]= (TH1F*)gDirectory->Get(A.c_str());
      /////////////////////////////
      A=Form("UQC%i",NumberSuffix);
      Chi2L[i]->SetDefaultSumw2(); Chi2L[i]= (TH1F*)gDirectory->Get(A.c_str());
      ///////////////VTX//////////////
      A=Form("VTX%i",NumberSuffix);
      VtcsL[i]->SetDefaultSumw2(); VtcsL[i]= (TH1F*)gDirectory->Get(A.c_str()); VtcsStackIntegral+=VtcsL[i]->Integral(); VtcsL[i]->Scale(ActualLumi/Lumi);  VtcsL[i]->Rebin(ReBin);
      //if (i==0) {AllBKGSVtcs->SetDefaultSumw2(); AllBKGSVtcs= (TH1F*)VtcsL[i]->Clone("ALLBKGSVtcs");}
      ///////////////VTX//////////////
      A=Form("TopFromHiggsChi2Mass%i",NumberSuffix);
      TopFromHiggsL[i]->SetDefaultSumw2(); TopFromHiggsL[i]= (TH1F*)gDirectory->Get(A.c_str()); TopFromHiggsL[i]->Scale(ActualLumi/Lumi);  TopFromHiggsL[i]->Rebin(ReBin);
      ////////////////////////////////ADDITIONAL INFO OF JETS
      A=Form("jet1_eta%i",NumberSuffix);
      JET1ETAL[i]->SetDefaultSumw2(); JET1ETAL[i]= (TH1F*)gDirectory->Get(A.c_str()); JET1ETAL[i]->Scale(ActualLumi/Lumi);  JET1ETAL[i]->Rebin(ReBin);
      A=Form("jet2_eta%i",NumberSuffix);
      JET2ETAL[i]->SetDefaultSumw2(); JET2ETAL[i]= (TH1F*)gDirectory->Get(A.c_str()); JET2ETAL[i]->Scale(ActualLumi/Lumi);  JET2ETAL[i]->Rebin(ReBin);
      A=Form("jet3_eta%i",NumberSuffix);
      JET3ETAL[i]->SetDefaultSumw2(); JET3ETAL[i]= (TH1F*)gDirectory->Get(A.c_str()); JET3ETAL[i]->Scale(ActualLumi/Lumi);  JET3ETAL[i]->Rebin(ReBin);
      A=Form("jet4_eta%i",NumberSuffix);
      JET4ETAL[i]->SetDefaultSumw2(); JET4ETAL[i]= (TH1F*)gDirectory->Get(A.c_str()); JET4ETAL[i]->Scale(ActualLumi/Lumi);  JET4ETAL[i]->Rebin(ReBin);
      A=Form("jet5_eta%i",NumberSuffix);
      JET5ETAL[i]->SetDefaultSumw2(); JET5ETAL[i]= (TH1F*)gDirectory->Get(A.c_str()); JET5ETAL[i]->Scale(ActualLumi/Lumi);  JET5ETAL[i]->Rebin(ReBin);
      A=Form("jet6_eta%i",NumberSuffix);
      JET6ETAL[i]->SetDefaultSumw2(); JET6ETAL[i]= (TH1F*)gDirectory->Get(A.c_str()); JET6ETAL[i]->Scale(ActualLumi/Lumi);  JET6ETAL[i]->Rebin(ReBin);
      A=Form("jet1_phi%i",NumberSuffix);
      JET1PHIL[i]->SetDefaultSumw2(); JET1PHIL[i]= (TH1F*)gDirectory->Get(A.c_str()); JET1PHIL[i]->Scale(ActualLumi/Lumi);  JET1PHIL[i]->Rebin(ReBin);
      A=Form("jet2_phi%i",NumberSuffix);
      JET2PHIL[i]->SetDefaultSumw2(); JET2PHIL[i]= (TH1F*)gDirectory->Get(A.c_str()); JET2PHIL[i]->Scale(ActualLumi/Lumi);  JET2PHIL[i]->Rebin(ReBin);
      A=Form("jet3_phi%i",NumberSuffix);
      JET3PHIL[i]->SetDefaultSumw2(); JET3PHIL[i]= (TH1F*)gDirectory->Get(A.c_str()); JET3PHIL[i]->Scale(ActualLumi/Lumi);  JET3PHIL[i]->Rebin(ReBin);
      A=Form("jet4_phi%i",NumberSuffix);
      JET4PHIL[i]->SetDefaultSumw2(); JET4PHIL[i]= (TH1F*)gDirectory->Get(A.c_str()); JET4PHIL[i]->Scale(ActualLumi/Lumi);  JET4PHIL[i]->Rebin(ReBin);
      A=Form("jet5_phi%i",NumberSuffix);
      JET5PHIL[i]->SetDefaultSumw2(); JET5PHIL[i]= (TH1F*)gDirectory->Get(A.c_str()); JET5PHIL[i]->Scale(ActualLumi/Lumi);  JET5PHIL[i]->Rebin(ReBin);
      A=Form("jet6_phi%i",NumberSuffix);
      JET6PHIL[i]->SetDefaultSumw2(); JET6PHIL[i]= (TH1F*)gDirectory->Get(A.c_str()); JET6PHIL[i]->Scale(ActualLumi/Lumi);  JET6PHIL[i]->Rebin(ReBin);
      A=Form("Num_jets%i",NumberSuffix);
      JETMULTIL[i]->SetDefaultSumw2(); JETMULTIL[i]= (TH1F*)gDirectory->Get(A.c_str()); JETMULTIL[i]->Scale(ActualLumi/Lumi);  JETMULTIL[i]->Rebin(ReBin);
    }
  
  /////////////////////////////////
  //Reading files for data sample//
  /////////////////////////////////
  TFile *CurrentFileT[NOD];
  TH1F *TprimeHistosT[NOD];
  TH1F *LeadingJetPTT[NOD];
  TH1F *Leading2JetPTT[NOD];
  TH1F *Leading3JetPTT[NOD];
  TH1F *Leading4JetPTT[NOD];
  TH1F *Leading5JetPTT[NOD];
  TH1F *Leading6JetPTT[NOD];
  TH1F *THTT[NOD];
  TH1F *DRHjetsT[NOD];
  TH1F *DRWjetsT[NOD];
  TH1F *HptT[NOD];
  TH1F *TptT[NOD];
  TH1F *DRWHT[NOD];
  TH1F *DPHjetsT[NOD];
  TH1F *DPWjetsT[NOD];
  TH1F *DPTjetsT[NOD];
  TH1F *HiggsMassT[NOD];
  TH1F *RelHTT[NOD];
  TH1F *DRTHT[NOD];
  TH1F *PtNormalizedMassT[NOD];
  TH1F *RelativeMassT[NOD];
  TH1F *MotherPtNormalizedMassT[NOD];
  TH1F *NumberOfTopsT[NOD];
  TH1F *HiggsMassOverTopMassT[NOD];
  TH1F *HiggsTopAsymmetryT[NOD];
  TH1F *ThirdLooseBtagT[NOD];
  TH1F *TopMassT[NOD];
  TH1F *Chi2T[NOD];
  TH1F *VtcsT[NOD]; 
  TH1F *TopFromHiggsT[NOD];
  TH1F *JET1ETAT[NOD];
  TH1F *JET2ETAT[NOD];
  TH1F *JET3ETAT[NOD]; 
  TH1F *JET4ETAT[NOD];
  TH1F *JET5ETAT[NOD];
  TH1F *JET6ETAT[NOD];
  TH1F *JET1PHIT[NOD];   
  TH1F *JET2PHIT[NOD]; 
  TH1F *JET3PHIT[NOD]; 
  TH1F *JET4PHIT[NOD];
  TH1F *JET5PHIT[NOD];
  TH1F *JET6PHIT[NOD];
  TH1F *JETMULTIT[NOD];

  for (int i=0; i<NOD; i++)
    {
      //cout << StorageDirPrefixDATA + DATASamples[i] << endl;
      CurrentFileT[i] = new TFile(StorageDirPrefixDATA + DATASamples[i], "READ");
      cout << StorageDirPrefixDATA + DATASamples[i] << endl;
      if ( CurrentFileT[i]->IsOpen() ) printf( DATASamples[i] + " File opened successfully\n");  
      ///////////////////////////// 
      string A=Form("TprimeMass");
      TprimeHistosT[i]->SetDefaultSumw2(); TprimeHistosT[i]= (TH1F*)gDirectory->Get(A.c_str()); TprimeHistosT[i]->Rebin(ReBin);
      /////////////////////////////
      A=Form("jet1_pt");
      LeadingJetPTT[i]->SetDefaultSumw2(); LeadingJetPTT[i]= (TH1F*)gDirectory->Get(A.c_str()); LeadingJetPTT[i]->Rebin(ReBin);
      /////////////////////////////
      A=Form("jet2_pt");
      Leading2JetPTT[i]->SetDefaultSumw2(); Leading2JetPTT[i]= (TH1F*)gDirectory->Get(A.c_str());
      /////////////////////////////
      A=Form("jet3_pt");
      Leading3JetPTT[i]->SetDefaultSumw2(); Leading3JetPTT[i]= (TH1F*)gDirectory->Get(A.c_str());
      /////////////////////////////
      A=Form("jet4_pt");
      Leading4JetPTT[i]->SetDefaultSumw2(); Leading4JetPTT[i]= (TH1F*)gDirectory->Get(A.c_str());
      /////////////////////////////
      A=Form("jet5_p");
      Leading5JetPTT[i]->SetDefaultSumw2(); Leading5JetPTT[i]= (TH1F*)gDirectory->Get(A.c_str());
      /////////////////////////////
      A=Form("jet6_pt");
      Leading6JetPTT[i]->SetDefaultSumw2(); Leading6JetPTT[i]= (TH1F*)gDirectory->Get(A.c_str());
      /////////////////////////////
      A=Form("THT");
      THTT[i]->SetDefaultSumw2(); THTT[i]= (TH1F*)gDirectory->Get(A.c_str()); THTT[i]->Rebin(ReBin);
      /////////////////////////////
      A=Form("DeltaR_of_Higgs_Jets");
      DRHjetsT[i]->SetDefaultSumw2(); DRHjetsT[i]= (TH1F*)gDirectory->Get(A.c_str()); DRHjetsT[i]->Rebin(ReBin); //DRHjetsT[i]->GetXaxis()->SetRangeUser(0.0,3.0);
      /////////////////////////////
      A=Form("DeltaR_of_W_Jets");
      DRWjetsT[i]->SetDefaultSumw2(); DRWjetsT[i]= (TH1F*)gDirectory->Get(A.c_str()); DRWjetsT[i]->Rebin(ReBin); //DRWjetsT[i]->GetXaxis()->SetRangeUser(0.0,3.0);
      /////////////////////////////
      A=Form("HPt");
      HptT[i]->SetDefaultSumw2(); HptT[i]= (TH1F*)gDirectory->Get(A.c_str()); HptT[i]->Rebin(ReBin);
      /////////////////////////////
      A=Form("TPt");
      TptT[i]->SetDefaultSumw2(); TptT[i]= (TH1F*)gDirectory->Get(A.c_str()); TptT[i]->Rebin(ReBin);
      /////////////////////////////
      A=Form("DeltaR_of_W_Higgs");
      DRWHT[i]->SetDefaultSumw2(); DRWHT[i]= (TH1F*)gDirectory->Get(A.c_str()); DRWHT[i]->Rebin(ReBin); //DRWHT[i]->GetXaxis()->SetRangeUser(0.0,3.0);
      /////////////////////////////
      A=Form("DeltaPhi_of_Higgs_jets");
      DPHjetsT[i]->SetDefaultSumw2(); DPHjetsT[i]= (TH1F*)gDirectory->Get(A.c_str()); DPHjetsT[i]->Rebin(ReBin); //DPHjetsT[i]->GetXaxis()->SetRangeUser(0.0,3.0);
      /////////////////////////////
      A=Form("DeltaPhi_of_W_jets");
      DPWjetsT[i]->SetDefaultSumw2(); DPWjetsT[i]= (TH1F*)gDirectory->Get(A.c_str()); DPWjetsT[i]->Rebin(ReBin); //DPWjetsT[i]->GetXaxis()->SetRangeUser(0.0,3.0);
      /////////////////////////////
      A=Form("DeltaPhi_of_T_jet");
      DPTjetsT[i]->SetDefaultSumw2(); DPTjetsT[i]= (TH1F*)gDirectory->Get(A.c_str()); DPTjetsT[i]->Rebin(ReBin); //DPTjetsT[i]->GetXaxis()->SetRangeUser(0.0,3.0);
      /////////////////////////////
      A=Form("HM");
      HiggsMassT[i]->SetDefaultSumw2(); HiggsMassT[i]= (TH1F*)gDirectory->Get(A.c_str()); HiggsMassT[i]->Rebin(ReBin);
      /////////////////////////////
      A=Form("RelHT");
      RelHTT[i]->SetDefaultSumw2(); RelHTT[i]= (TH1F*)gDirectory->Get(A.c_str()); RelHTT[i]->Rebin(ReBin);
      /////////////////////////////
      A=Form("DeltaR_of_Top_Higgs");
      DRTHT[i]->SetDefaultSumw2(); DRTHT[i]= (TH1F*)gDirectory->Get(A.c_str());
      /////////////////////////////
      A=Form("PT_Normalized_Mass");
      PtNormalizedMassT[i]->SetDefaultSumw2(); PtNormalizedMassT[i]= (TH1F*)gDirectory->Get(A.c_str());
      /////////////////////////////
      A=Form("Relative_Mass");
      RelativeMassT[i]->SetDefaultSumw2(); RelativeMassT[i]= (TH1F*)gDirectory->Get(A.c_str()); RelativeMassT[i]->Rebin(ReBin);
      /////////////////////////////
      A=Form("Mother_PT_Normalized_Mass");
      MotherPtNormalizedMassT[i]->SetDefaultSumw2(); MotherPtNormalizedMassT[i]= (TH1F*)gDirectory->Get(A.c_str());
      /////////////////////////////
      A=Form("Number_of_Tops");
      NumberOfTopsT[i]->SetDefaultSumw2();NumberOfTopsT[i]= (TH1F*)gDirectory->Get(A.c_str());
      /////////////////////////////
      A=Form("HMoverTM");
      HiggsMassOverTopMassT[i]->SetDefaultSumw2();HiggsMassOverTopMassT[i]= (TH1F*)gDirectory->Get(A.c_str());
      /////////////////////////////
      A=Form("HTAsym");
      HiggsTopAsymmetryT[i]->SetDefaultSumw2();HiggsTopAsymmetryT[i]= (TH1F*)gDirectory->Get(A.c_str());
      /////////////////////////////
      A=Form("TLBTag");
      HiggsTopAsymmetryT[i]->SetDefaultSumw2();HiggsTopAsymmetryT[i]= (TH1F*)gDirectory->Get(A.c_str());
      /////////////////////////////
      A=Form("TMass");
      TopMassT[i]->SetDefaultSumw2();TopMassT[i]= (TH1F*)gDirectory->Get(A.c_str()); TopMassT[i]->Rebin(ReBin);
      /////////////////////////////
      A=Form("ChiSq");
      //TopMassT[i]->SetDefaultSumw2(); TopMassT[i]= (TH1F*)gDirectory->Get(A.c_str());
      /////////////////////////////
      A=Form("UQC");
      Chi2T[i]->SetDefaultSumw2(); Chi2T[i]= (TH1F*)gDirectory->Get(A.c_str());
      /////////////////////////////
      A=Form("VTX");
      VtcsT[i]->SetDefaultSumw2(); VtcsT[i]= (TH1F*)gDirectory->Get(A.c_str()); VtcsT[i]->Rebin(ReBin);
      ///////////////VTX//////////////
      A=Form("TopFromHiggsChi2Mass");
      TopFromHiggsT[i]->SetDefaultSumw2(); TopFromHiggsT[i]= (TH1F*)gDirectory->Get(A.c_str()); TopFromHiggsT[i]->Rebin(ReBin);
      ////////////////////////////////ADDITIONAL INFO OF JETS
      A=Form("jet1_eta");
      JET1ETAT[i]->SetDefaultSumw2(); JET1ETAT[i]= (TH1F*)gDirectory->Get(A.c_str()); JET1ETAT[i]->Rebin(ReBin);
      A=Form("jet2_eta");
      JET2ETAT[i]->SetDefaultSumw2(); JET2ETAT[i]= (TH1F*)gDirectory->Get(A.c_str()); JET2ETAT[i]->Rebin(ReBin);
      A=Form("jet3_eta");
      JET3ETAT[i]->SetDefaultSumw2(); JET3ETAT[i]= (TH1F*)gDirectory->Get(A.c_str()); JET3ETAT[i]->Rebin(ReBin);
      A=Form("jet4_eta");
      JET4ETAT[i]->SetDefaultSumw2(); JET4ETAT[i]= (TH1F*)gDirectory->Get(A.c_str()); JET4ETAT[i]->Rebin(ReBin);
      A=Form("jet5_eta");
      JET5ETAT[i]->SetDefaultSumw2(); JET5ETAT[i]= (TH1F*)gDirectory->Get(A.c_str()); JET5ETAT[i]->Rebin(ReBin);
      A=Form("jet6_eta");
      JET6ETAT[i]->SetDefaultSumw2(); JET6ETAT[i]= (TH1F*)gDirectory->Get(A.c_str()); JET6ETAT[i]->Rebin(ReBin);
      A=Form("jet1_phi");
      JET1PHIT[i]->SetDefaultSumw2(); JET1PHIT[i]= (TH1F*)gDirectory->Get(A.c_str()); JET1PHIT[i]->Rebin(ReBin);
      A=Form("jet2_phi");
      JET2PHIT[i]->SetDefaultSumw2(); JET2PHIT[i]= (TH1F*)gDirectory->Get(A.c_str()); JET2PHIT[i]->Rebin(ReBin);
      A=Form("jet3_phi");
      JET3PHIT[i]->SetDefaultSumw2(); JET3PHIT[i]= (TH1F*)gDirectory->Get(A.c_str()); JET3PHIT[i]->Rebin(ReBin);
      A=Form("jet4_phi");
      JET4PHIT[i]->SetDefaultSumw2(); JET4PHIT[i]= (TH1F*)gDirectory->Get(A.c_str()); JET4PHIT[i]->Rebin(ReBin);
      A=Form("jet5_phi");
      JET5PHIT[i]->SetDefaultSumw2(); JET5PHIT[i]= (TH1F*)gDirectory->Get(A.c_str()); JET5PHIT[i]->Rebin(ReBin);
      A=Form("jet6_phi");
      JET6PHIT[i]->SetDefaultSumw2(); JET6PHIT[i]= (TH1F*)gDirectory->Get(A.c_str()); JET6PHIT[i]->Rebin(ReBin);
      A=Form("Num_jets");
      JETMULTIT[i]->SetDefaultSumw2(); JETMULTIT[i]= (TH1F*)gDirectory->Get(A.c_str()); JETMULTIT[i]->Rebin(ReBin);
    }

  //Stacks of backgrounds
  THStack *BKGLJPTL = new THStack("BKGLJPT", "BKG for LJPT; Leading j pT GeV; Events/20 GeV");
  THStack *BKGTHTL = new THStack("BKGTHT", "BKG for HT; HT GeV; Events/20 GeV");
  THStack *BKGVtcsL = new THStack("BKGVtcsL", "BKG for Vtcs; N_{v}; Events");
  THStack *BKGDRHJL = new THStack("BKGDRHJ", "BKG for DRHJ; #Delta R_{bb}^{H}; Events/0.1");
  THStack *BKGTPML = new THStack("BKGTPM", "BKG for TPM; M_{5j} GeV; Events/20 GeV");
  THStack *BKGDRWJ = new THStack("BKGDRWJ", "BKG for DRWJ; #Delta R_{jj}^{W}; Events/0.1");
  THStack *BKGHPT = new THStack("BKGHPT", "BKG for HPT; p_{T}^{H} GeV; Events/20 GeV");
  THStack *BKGTPT = new THStack("BKGTPT", "BKG for TPT; p_{T}^{t} GeV; Events/20 GeV");
  THStack *BKGDRWH = new THStack("BKGDRWH", "BKG for DRWH; #Delta R_{WH}; Events/0.1");
  THStack *BKGDPHJ = new THStack("BKGDPHJ", "BKG for DPHJ; #Delta #phi_{bb}^{H}; Events/0.1");
  THStack *BKGDPWJ = new THStack("BKGDPWJ", "BKG for DPWJ; #Delta #phi_{jj}^{W}; Events/0.1");
  THStack *BKGDPTJ = new THStack("BKGDPTJ", "BKG for DPTJ; #Delta #phi_{Wj}^{t}; Events/0.1");
  THStack *BKGHM = new THStack("BKGHM", "BKG for HM; M_{H} GeV; Events/20 GeV");
  THStack *BKGRHT = new THStack("BKGRHT", "BKG for RHT; (M_{H}+M_{t})/HT ; Events/0.1");
  THStack *BKGRM = new THStack("BKGRM", "BKG for RM; (M_{H}+M_{t})/M_{j}^{5} ; Events/0.1");
  THStack *BKGTM = new THStack("BKGTM", "BKG for TM; M_{t} GeV; Events/20 GeV");
  THStack *BKGTFH = new THStack("BKGTFH", "BKG for TFHM; M_{t^{2nd}} GeV; Events/20 GeV");
  THStack *BKGJET1ETA = new THStack("BKGJET1ETA", "BKG for JET1ETA; #eta_{j_{1}}; Events/0.1"); 
  THStack *BKGJET2ETA = new THStack("BKGJET2ETA", "BKG for JET2ETA; #eta_{j_{2}}; Events/0.1");  
  THStack *BKGJET3ETA = new THStack("BKGJET3ETA", "BKG for JET3ETA; #eta_{j_{3}}; Events/0.1"); 
  THStack *BKGJET4ETA = new THStack("BKGJET4ETA", "BKG for JET4ETA; #eta_{j_{4}}; Events/0.1"); 
  THStack *BKGJET5ETA = new THStack("BKGJET5ETA", "BKG for JET5ETA; #eta_{j_{5}}; Events/0.1"); 
  THStack *BKGJET6ETA = new THStack("BKGJET6ETA", "BKG for JET6ETA; #eta_{j_{6}}; Events/0.1"); 
  THStack *BKGJET1PHI = new THStack("BKGJET1PHI", "BKG for JET1PHI; #phi_{j_{1}}; Events/0.1"); 
  THStack *BKGJET2PHI = new THStack("BKGJET2PHI", "BKG for JET2PHI; #phi_{j_{2}}; Events/0.1"); 
  THStack *BKGJET3PHI = new THStack("BKGJET3PHI", "BKG for JET3PHI; #phi_{j_{3}}; Events/0.1"); 
  THStack *BKGJET4PHI = new THStack("BKGJET4PHI", "BKG for JET4PHI; #phi_{j_{4}}; Events/0.1"); 
  THStack *BKGJET5PHI = new THStack("BKGJET5PHI", "BKG for JET5PHI; #phi_{j_{5}}; Events/0.1"); 
  THStack *BKGJET6PHI = new THStack("BKGJET6PHI", "BKG for JET6PHI; #phi_{j_{6}}; Events/0.1"); 
  THStack *BKGJETMULTI = new THStack("BKGJETMULTI", "BKG for JETMULTI; n_{j}; Events"); 
  /*THStack *BKG = new THStack("BKG", "BKG for ;  GeV; Events/20 GeV");
  THStack *BKG = new THStack("BKG", "BKG for ;  GeV; Events/20 GeV");
  THStack *BKG = new THStack("BKG", "BKG for ;  GeV; Events/20 GeV");*/
  
  cout << "Lumi for data: " << ActualLumi << endl;
  cout << "Factor between data and MC in THT: " <<  THTT[0]->Integral()/THTStackIntegral << endl;
  cout << "Factor between data and MC in Vtcs: " <<  VtcsT[0]->Integral()/VtcsStackIntegral << endl;
  cout << "Factor for actual lumi: " << ActualLumi/Lumi << endl;
  cout << "Factors rate: " << (ActualLumi/Lumi)/(THTT[0]->Integral()/THTStackIntegral) << endl;
  //cout << "Factor between data and Diboson MC in THT: " <<  THTT[0]->Integral()/THTL[0]->Integral() << endl;
  //cout << "Factor between data and Singletop MC in THT: " <<  THTT[0]->Integral()/THTL[1]->Integral() << endl;
  //cout << "Factor between data and TTjets MC in THT: " <<  THTT[0]->Integral()/THTL[2]->Integral() << endl;
  //cout << "Factor between data and DY MC in THT: " <<  THTT[0]->Integral()/THTL[3]->Integral() << endl;
  //cout << "Factor between data and QCD MC in THT: " <<  THTT[0]->Integral()/THTL[4]->Integral() << endl;
  
  //AllBKGSDRWJ
  TH1F *AllBKGSLJPT;
  TH1F *AllBKGSTHT;
  TH1F *AllBKGSVtcs;
  TH1F *AllBKGSDRHJ;
  TH1F *AllBKGSTPM; 
  TH1F *AllBKGSDRWJ;
  TH1F *AllBKGSHpt;
  TH1F *AllBKGSTpt;
  TH1F *AllBKGSDRWH;
  TH1F *AllBKGSDPHJ;
  TH1F *AllBKGSDPWJ;
  TH1F *AllBKGSDPTJ;
  TH1F *AllBKGSHM;
  TH1F *AllBKGSRelHT;
  TH1F *AllBKGSRM;
  TH1F *AllBKGSTM;
  TH1F *AllBKGSTFH;  
  TH1F *AllBKGSJET1ETA;  
  TH1F *AllBKGSJET2ETA;  
  TH1F *AllBKGSJET3ETA;  
  TH1F *AllBKGSJET4ETA;  
  TH1F *AllBKGSJET5ETA;  
  TH1F *AllBKGSJET6ETA;  
  TH1F *AllBKGSJET1PHI;   
  TH1F *AllBKGSJET2PHI;  
  TH1F *AllBKGSJET3PHI;  
  TH1F *AllBKGSJET4PHI;  
  TH1F *AllBKGSJET5PHI;  
  TH1F *AllBKGSJET6PHI;  
  TH1F *AllBKGSJETMULTI; 

  for (int i=0; i<NOS-1; i++)
    {
      //THTL[i]->Scale(THTT[0]->Integral()/THTStackIntegral);
      //VtcsL[i]->Scale(VtcsT[0]->Integral()/VtcsStackIntegral);
      //THTL[i]->Scale(1./THTL[i]->Integral()); VtcsL[i]->Scale(1./VtcsL[i]->Integral());
      //THTL[i]->Scale(THTT[i]->Integral()/(THTL[i]->Integral()/THTStackIntegral));
      //THTL[i]->Scale(ActualLumi/Lumi); // *(THTT[0]->Integral()/THTStackIntegral));
      /*DRHjetsL[i]->GetXaxis()->SetRangeUser(0.0,3.0);
      DRWjetsL[i]->GetXaxis()->SetRangeUser(0.0,3.0);
      DPHjetsL[i]->GetXaxis()->SetRangeUser(0.0,3.0);
      DPWjetsL[i]->GetXaxis()->SetRangeUser(0.0,3.0);
      DPTjetsL[i]->GetXaxis()->SetRangeUser(0.0,3.0); 
      DRWHL[i]->GetXaxis()->SetRangeUser(0.0,3.0);
      if (i==NOS-2)
	{
	  DRHjetsL[i+1]->GetXaxis()->SetRangeUser(0.0,3.0);
	  DRWjetsL[i+1]->GetXaxis()->SetRangeUser(0.0,3.0);
	  DPHjetsL[i+1]->GetXaxis()->SetRangeUser(0.0,3.0);
	  DPWjetsL[i+1]->GetXaxis()->SetRangeUser(0.0,3.0);
	  DPTjetsL[i+1]->GetXaxis()->SetRangeUser(0.0,3.0); 
	  DRWHL[i+1]->GetXaxis()->SetRangeUser(0.0,3.0);
	  }*/
      if (i==0)
	{
	  AllBKGSTPM->SetDefaultSumw2(); AllBKGSTPM= (TH1F*)TprimeHistosL[i]->Clone("ALLBKGSTPM");
	  AllBKGSLJPT->SetDefaultSumw2(); AllBKGSLJPT= (TH1F*)LeadingJetPTL[i]->Clone("ALLBKGSLJPT");
	  AllBKGSTHT->SetDefaultSumw2(); AllBKGSTHT= (TH1F*)THTL[i]->Clone("ALLBKGSTHT");
	  AllBKGSDRHJ->SetDefaultSumw2(); AllBKGSDRHJ= (TH1F*)DRHjetsL[i]->Clone("ALLBKGSDRHJ");
	  AllBKGSVtcs->SetDefaultSumw2(); AllBKGSVtcs= (TH1F*)VtcsL[i]->Clone("ALLBKGSVtcs");
	  AllBKGSDRWJ->SetDefaultSumw2(); AllBKGSDRWJ= (TH1F*)DRWjetsL[i]->Clone("AllBKGSDRWJ");
	  AllBKGSHpt->SetDefaultSumw2(); AllBKGSHpt= (TH1F*)HptL[i]->Clone("AllBKGSHpt");
	  AllBKGSTpt->SetDefaultSumw2(); AllBKGSTpt= (TH1F*)TptL[i]->Clone("AllBKGSTpt");
	  AllBKGSDRWH->SetDefaultSumw2(); AllBKGSDRWH= (TH1F*)DRWHL[i]->Clone("AllBKGSDRWH");
	  AllBKGSDPHJ->SetDefaultSumw2(); AllBKGSDPHJ= (TH1F*)DPHjetsL[i]->Clone("AllBKGSDPHJ");
	  AllBKGSDPWJ->SetDefaultSumw2(); AllBKGSDPWJ= (TH1F*)DPWjetsL[i]->Clone("AllBKGSDPWJ");
	  AllBKGSDPTJ->SetDefaultSumw2(); AllBKGSDPTJ= (TH1F*)DPTjetsL[i]->Clone("AllBKGSDPTJ");
	  AllBKGSHM->SetDefaultSumw2(); AllBKGSHM= (TH1F*)HiggsMassL[i]->Clone("AllBKGSHM");
	  AllBKGSRelHT->SetDefaultSumw2(); AllBKGSRelHT= (TH1F*)RelHTL[i]->Clone("AllBKGSRelHT");
	  AllBKGSRM->SetDefaultSumw2(); AllBKGSRM= (TH1F*)RelativeMassL[i]->Clone("AllBKGSRM");
	  AllBKGSTM->SetDefaultSumw2(); AllBKGSTM= (TH1F*)TopMassL[i]->Clone("AllBKGSTM");
	  AllBKGSTFH->SetDefaultSumw2(); AllBKGSTFH= (TH1F*)TopFromHiggsL[i]->Clone("AllBKGSTFH");
	  AllBKGSJET1ETA->SetDefaultSumw2(); AllBKGSJET1ETA= (TH1F*)JET1ETAL[i]->Clone("AllBKGSJET1ETA");
	  AllBKGSJET2ETA->SetDefaultSumw2(); AllBKGSJET2ETA= (TH1F*)JET2ETAL[i]->Clone("AllBKGSJET2ETA");
	  AllBKGSJET3ETA->SetDefaultSumw2(); AllBKGSJET3ETA= (TH1F*)JET3ETAL[i]->Clone("AllBKGSJET3ETA");
	  AllBKGSJET4ETA->SetDefaultSumw2(); AllBKGSJET4ETA= (TH1F*)JET4ETAL[i]->Clone("AllBKGSJET4ETA");
	  AllBKGSJET5ETA->SetDefaultSumw2(); AllBKGSJET5ETA= (TH1F*)JET5ETAL[i]->Clone("AllBKGSJET5ETA");
	  AllBKGSJET6ETA->SetDefaultSumw2(); AllBKGSJET6ETA= (TH1F*)JET6ETAL[i]->Clone("AllBKGSJET6ETA");
	  AllBKGSJET1PHI->SetDefaultSumw2(); AllBKGSJET1PHI= (TH1F*)JET1PHIL[i]->Clone("AllBKGSJET1PHI");
	  AllBKGSJET2PHI->SetDefaultSumw2(); AllBKGSJET2PHI= (TH1F*)JET2PHIL[i]->Clone("AllBKGSJET2PHI");
	  AllBKGSJET3PHI->SetDefaultSumw2(); AllBKGSJET3PHI= (TH1F*)JET3PHIL[i]->Clone("AllBKGSJET3PHI");
	  AllBKGSJET4PHI->SetDefaultSumw2(); AllBKGSJET4PHI= (TH1F*)JET4PHIL[i]->Clone("AllBKGSJET4PHI");
	  AllBKGSJET5PHI->SetDefaultSumw2(); AllBKGSJET5PHI= (TH1F*)JET5PHIL[i]->Clone("AllBKGSJET5PHI");
	  AllBKGSJET6PHI->SetDefaultSumw2(); AllBKGSJET6PHI= (TH1F*)JET6PHIL[i]->Clone("AllBKGSJET6PHI");
	  AllBKGSJETMULTI->SetDefaultSumw2(); AllBKGSJETMULTI= (TH1F*)JETMULTIL[i]->Clone("AllBKGSJETMULTI");
	}
      else 
	{
	  AllBKGSTHT->Sumw2(); AllBKGSTHT->Add(THTL[i]); 
	  AllBKGSVtcs->Sumw2(); AllBKGSVtcs->Add(VtcsL[i]); 
	  AllBKGSLJPT->Sumw2(); AllBKGSLJPT->Add(LeadingJetPTL[i]); 
	  AllBKGSDRHJ->Sumw2(); AllBKGSDRHJ->Add(DRHjetsL[i]); 
	  AllBKGSTPM->Sumw2(); AllBKGSTPM->Add(TprimeHistosL[i]);
	  AllBKGSDRWJ->Sumw2(); AllBKGSDRWJ->Add(DRWjetsL[i]);
	  AllBKGSHpt->Sumw2(); AllBKGSHpt->Add(HptL[i]);
	  AllBKGSTpt->Sumw2(); AllBKGSTpt->Add(TptL[i]);
	  AllBKGSDRWH->Sumw2(); AllBKGSDRWH->Add(DRWHL[i]);
	  AllBKGSDPHJ->Sumw2(); AllBKGSDPHJ->Add(DPHjetsL[i]);
	  AllBKGSDPWJ->Sumw2(); AllBKGSDPWJ->Add(DPWjetsL[i]);
	  AllBKGSDPTJ->Sumw2(); AllBKGSDPTJ->Add(DPTjetsL[i]);
	  AllBKGSHM->Sumw2(); AllBKGSHM->Add(HiggsMassL[i]);
	  AllBKGSRelHT->Sumw2(); AllBKGSRelHT->Add(RelHTL[i]);
	  AllBKGSRM->Sumw2(); AllBKGSRM->Add(RelativeMassL[i]);
	  AllBKGSTM->Sumw2(); AllBKGSTM->Add(TopMassL[i]);
	  AllBKGSTFH->Sumw2(); AllBKGSTFH->Add(TopFromHiggsL[i]);
	  AllBKGSJET1ETA->Sumw2(); AllBKGSJET1ETA->Add(JET1ETAL[i]);
	  AllBKGSJET2ETA->Sumw2(); AllBKGSJET2ETA->Add(JET2ETAL[i]);
	  AllBKGSJET3ETA->Sumw2(); AllBKGSJET3ETA->Add(JET3ETAL[i]);
	  AllBKGSJET4ETA->Sumw2(); AllBKGSJET4ETA->Add(JET4ETAL[i]);
	  AllBKGSJET5ETA->Sumw2(); AllBKGSJET5ETA->Add(JET5ETAL[i]);
	  AllBKGSJET6ETA->Sumw2(); AllBKGSJET6ETA->Add(JET6ETAL[i]);
	  AllBKGSJET1PHI->Sumw2(); AllBKGSJET1PHI->Add(JET1PHIL[i]);
	  AllBKGSJET2PHI->Sumw2(); AllBKGSJET2PHI->Add(JET2PHIL[i]);
	  AllBKGSJET3PHI->Sumw2(); AllBKGSJET3PHI->Add(JET3PHIL[i]);
	  AllBKGSJET4PHI->Sumw2(); AllBKGSJET4PHI->Add(JET4PHIL[i]);
	  AllBKGSJET5PHI->Sumw2(); AllBKGSJET5PHI->Add(JET5PHIL[i]);
	  AllBKGSJET6PHI->Sumw2(); AllBKGSJET6PHI->Add(JET6PHIL[i]);
	  AllBKGSJETMULTI->Sumw2(); AllBKGSJETMULTI->Add(JETMULTIL[i]);
	  
	}
      BKGLJPTL->Add(LeadingJetPTL[i]);
      BKGTHTL->Add(THTL[i]);
      BKGVtcsL->Add(VtcsL[i]);
      BKGDRHJL->Add(DRHjetsL[i]);
      BKGTPML->Add(TprimeHistosL[i]);
      BKGDRWJ->Add(DRWjetsL[i]);
      BKGHPT->Add(HptL[i]);
      BKGTPT->Add(TptL[i]);
      BKGDRWH->Add(DRWHL[i]);
      BKGDPHJ->Add(DPHjetsL[i]);
      BKGDPWJ->Add(DPWjetsL[i]);
      BKGDPTJ->Add(DPTjetsL[i]);
      BKGHM->Add(HiggsMassL[i]);
      BKGRHT->Add(RelHTL[i]);
      BKGRM->Add(RelativeMassL[i]);
      BKGTM->Add(TopMassL[i]);
      BKGTFH->Add(TopFromHiggsL[i]);
      BKGJET1ETA->Add(JET1ETAL[i]);
      BKGJET2ETA->Add(JET2ETAL[i]);
      BKGJET3ETA->Add(JET3ETAL[i]);
      BKGJET4ETA->Add(JET4ETAL[i]);
      BKGJET5ETA->Add(JET5ETAL[i]);
      BKGJET6ETA->Add(JET6ETAL[i]);
      BKGJET1PHI->Add(JET1PHIL[i]);
      BKGJET2PHI->Add(JET2PHIL[i]);
      BKGJET3PHI->Add(JET3PHIL[i]);
      BKGJET4PHI->Add(JET4PHIL[i]);
      BKGJET5PHI->Add(JET5PHIL[i]);
      BKGJET6PHI->Add(JET6PHIL[i]);
      BKGJETMULTI->Add(JETMULTIL[i]);
      
    }

  /*BKGDRHJL->GetXaxis()->SetRangeUser(0.0,3.0);
  BKGDRWJ->GetXaxis()->SetRangeUser(0.0,3.0);
  BKGDRWH->GetXaxis()->SetRangeUser(0.0,3.0);
  BKGDPHJ->GetXaxis()->SetRangeUser(0.0,3.0);
  BKGDPWJ->GetXaxis()->SetRangeUser(0.0,3.0);
  BKGDPTJ->GetXaxis()->SetRangeUser(0.0,3.0);*/

  //cout << THTT[0]->Integral()/AllBKGSTHT->Integral() << endl;
  //AllBKGSTHT->Scale(ActualLumi/Lumi);
  //AllBKGSTHT->Scale(THTT[0]->Integral()/AllBKGSTHT->Integral());
  //AllBKGSVtcs->Scale(VtcsT[0]->Integral()/AllBKGSVtcs->Integral());

  cout << "HERE1" << endl;

  //Stack of data
  THStack *DATATHT = new THStack("DATATHT", "DATA for HT; HT GeV; Events");
  THStack *DATAVtcs = new THStack("DATAVtcs", "DATA for Vtcs; n_{v}; Events");
  for (int i=0; i<NOD; i++)
    {
      DATATHT->Add(THTT[i]);
      DATAVtcs->Add(VtcsT[i]);
    }

  //THStack *SIGBKGTHTL = new THStack("SIGBKGTHT", "Sig and BKG for HT; HT GeV; Events");
  //for (int i=0; i<NOS; i++)
  //  {
  //    SIGBKGTHTL->Add(THTL[i]);
  //  }

  TLegend* BKGlegend = new TLegend(0.75,0.65,0.90,0.9);
  BKGlegend->AddEntry(THTT[0], "Data", "pe");
  BKGlegend->AddEntry(THTL[4], "QCD", "f");
  BKGlegend->AddEntry(THTL[3], "Zjets", "f");
  BKGlegend->AddEntry(THTL[2], "TTbar", "f");
  BKGlegend->AddEntry(THTL[1], "SingleT", "f");
  BKGlegend->AddEntry(THTL[0], "Diboson", "f");
  BKGlegend->AddEntry(THTL[NOS-1], "MT'=734 GeV", "f");
  TLegend* BKGlegend2 = new TLegend(0.75,0.65,0.90,0.9);
  BKGlegend2->AddEntry(THTL[4], "QCD", "f");
  BKGlegend2->AddEntry(THTL[3], "Zjets", "f");
  BKGlegend2->AddEntry(THTL[2], "TTbar", "f");
  BKGlegend2->AddEntry(THTL[1], "SingleT", "f");
  BKGlegend2->AddEntry(THTL[0], "Diboson", "f");
  BKGlegend2->AddEntry(THTL[NOS-1], "MT'=734 GeV", "f");
  cout << "HERE1" << endl;

  //Plotting!
  char FileName[100];
  sprintf(FileName,"DataMCComparison_SUFFIXMC_SUFFIXDATA.pdf");
  TPDF *ps = new TPDF(FileName,111);
  TCanvas *MyPlot = new TCanvas("MyPlot","Data MC Comparison",600,800);
  TPad *pad1 = new TPad("pad1","pad1",0,0.3,1,1); //Line for ratio plot
  TPaveText *p1;
  p1 = new TPaveText(3,0.5,8,3.5);
  char dateandtime[50];
  sprintf(dateandtime,"date: %s, time: %s",__DATE__,__TIME__);
  p1->AddText(dateandtime);
  cout << "HERE1" << endl;
  //////////////
  //First Page//
  //////////////
  MyPlot->Clear();
  MyPlot->cd(1);
  ps->NewPage();
  gStyle->SetOptStat(0);//Remove the Stat Box
  //AllBKGSTHT->Draw("e hist");
  BKGTHTL->Draw("hist");
  BKGTHTL->SetMinimum(1);
  //DATATHT->SetFillColor(kBlack);
  THTT[0]->Draw("E1 same");
  //DATATHT->Draw("E1 same");
  gPad->SetLogy();
  //BKGTHTL->SetMinimum(0.1);
  BKGlegend->Draw();
  gPad->Update();
  MyPlot->Update();
  cout << "HERE1" << endl;
  /*ps->NewPage();
  gStyle->SetOptStat(0);//Remove the Stat Box
  //TH1 *h3=AllBKGSTHT->DrawCopy("e hist"); //Line for ratio plot
  //h3->SetMinimum(-100); //Line for ratio plot
  AllBKGSTHT->Draw("e hist");
  THTT[0]->Draw("E1 same");
  MyPlot->cd(); //Line for ratio plot
  h3->Sumw2(); //Line for ratio plot
  h3->SetStats(0); //Line for ratio plot
  h3->Divide(THTT[0]); //Line for ratio plot
  h3->SetMarkerStyle(21); //Line for ratio plot
  h3->Draw("ep"); //Line for ratio plot
  MyPlot->cd(); //Line for ratio plot
  gPad->Update();
  MyPlot->Update();*/
  cout << "HERE1" << endl;
  ///////////////
  //Second Page//
  ///////////////
  MyPlot->Clear(); 
  MyPlot->cd(1); 
  ps->NewPage();
  pad1->SetBottomMargin(0); //Line for ratio plot
  pad1->Draw(); //Line for ratio plot
  pad1->cd(); //Line for ratio plot
  AllBKGSTHT->DrawCopy("e0 hist");
  AllBKGSTHT->SetMinimum(0.5);
  AllBKGSTHT->SetMaximum(1.5);
  BKGTHTL->Draw("hist same");
  BKGTHTL->SetMinimum(1);
  THTT[0]->Draw("E1 same");
  THTL[NOS-1]->Draw("hist same");
  //gPad->SetLogy();
  BKGlegend->Draw();
  MyPlot->cd();
  TPad *pad2 = new TPad("pad2","pad2",0,0,1,0.3); //Line for ratio plot
  pad2->SetTopMargin(0); //Line for ratio plot
  pad2->Draw(); //Line for ratio plot
  pad2->cd(); //Line for ratio plot
  AllBKGSTHT->Sumw2();
  AllBKGSTHT->SetStats(0);
  AllBKGSTHT->SetTitle(";HT GeV;MC/Data");
  AllBKGSTHT->Divide(THTT[0]);
  AllBKGSTHT->SetMarkerStyle(21);
  AllBKGSTHT->Draw("ep");
  TLine *L1 = new TLine(630,1.3,1600,1.3);
  TLine *L2 = new TLine(630,1.0,1600,1.0);
  TLine *L3 = new TLine(630,0.7,1600,0.7);
  L1->Draw(); L2->Draw(); L3->Draw(); 
  MyPlot->cd();
  MyPlot->Update();
  //////////////
  //Third Page//
  //////////////
  MyPlot->Clear();
  MyPlot->cd(1);
  ps->NewPage();
  TPad *pad3 = new TPad("pad3","pad3",0,0.3,1,1); //Line for ratio plot  
  pad3->SetBottomMargin(0); //Line for ratio plot
  pad3->Draw(); //Line for ratio plot
  pad3->cd(); //Line for ratio plot

  //gStyle->SetOptStat(0);//Remove the Stat Box
  //AllBKGSVtcs->Draw("e hist");
  //BKGVtcsL->Draw("e hist");
  //DATATHT->SetFillColor(kBlack);

  AllBKGSVtcs->DrawCopy("e0 hist");
  AllBKGSVtcs->SetMinimum(0.1);
  AllBKGSVtcs->SetMaximum(1.9);  
  BKGVtcsL->Draw("hist");
  BKGVtcsL->SetMinimum(1);
  VtcsT[0]->Draw("E1 same");
  //gPad->SetLogy();
  BKGlegend->Draw();
  MyPlot->cd();
  TPad *pad4 = new TPad("pad4","pad4",0,0,1,0.3); //Line for ratio plot  
  pad4->SetTopMargin(0); //Line for ratio plot
  pad4->Draw(); //Line for ratio plot
  pad4->cd(); //Line for ratio plot
  //DATAVtcs->Draw("E1 same");
  //BKGTHTL->SetMinimum(0.1);
  AllBKGSVtcs->Sumw2();
  AllBKGSVtcs->SetStats(0);
  AllBKGSVtcs->SetTitle(";n_{vtcs};MC/Data");
  AllBKGSVtcs->Divide(VtcsT[0]);
  AllBKGSVtcs->SetMarkerStyle(21);
  AllBKGSVtcs->Draw("ep");
  //TLine *L1 = new TLine(630,1.3,1600,1.3);
  //TLine *L2 = new TLine(630,1.0,1600,1.0);
  //TLine *L3 = new TLine(630,0.7,1600,0.7);
  L1->DrawLine(0,1.3,50,1.3); L2->DrawLine(0,1,50,1); L3->DrawLine(0,0.7,50,0.7); 
  //BKGlegend->Draw();
  MyPlot->cd();  
  MyPlot->Update();
  cout << "HERE1" << endl;
  ///////////////
  //Fourth Page//
  ///////////////
  MyPlot->Clear();
  MyPlot->cd(1);
  ps->NewPage();
  TPad *pad5 = new TPad("pad3","pad3",0,0.3,1,1); //Line for ratio plot  
  pad5->SetBottomMargin(0); //Line for ratio plot
  pad5->Draw(); //Line for ratio plot
  pad5->cd(); //Line for ratio plot
  AllBKGSLJPT->DrawCopy("e0 hist");
  AllBKGSLJPT->SetMinimum(0.1);
  AllBKGSLJPT->SetMaximum(1.9);  
  BKGLJPTL->Draw("hist");
  BKGLJPTL->SetMinimum(1);
  LeadingJetPTT[0]->Draw("E1 same");
  LeadingJetPTL[NOS-1]->Draw("hist same");
  gPad->SetLogy();
  BKGlegend->Draw();
  MyPlot->cd();
  TPad *pad6 = new TPad("pad4","pad4",0,0,1,0.3); //Line for ratio plot  
  pad6->SetTopMargin(0); //Line for ratio plot
  pad6->Draw(); //Line for ratio plot
  pad6->cd(); //Line for ratio plot
  AllBKGSLJPT->Sumw2();
  AllBKGSLJPT->SetStats(0);
  AllBKGSLJPT->SetTitle(";Leading jet p_{T} GeV;MC/Data");
  AllBKGSLJPT->Divide(LeadingJetPTT[0]);
  AllBKGSLJPT->SetMarkerStyle(21);
  AllBKGSLJPT->Draw("ep");
  L1->DrawLine(150,1.3,1000,1.3); L2->DrawLine(150,1,1000,1); L3->DrawLine(150,0.7,1000,0.7); 
  //BKGlegend->Draw();
  MyPlot->cd();  
  MyPlot->Update();
  cout << "HERE1" << endl;
  ///////////////
  ///////////////
  ///////////////
  MyPlot->Clear();
  MyPlot->cd(1);
  ps->NewPage();
  TPad *pad7 = new TPad("pad3","pad3",0,0.3,1,1); //Line for ratio plot  
  pad7->SetBottomMargin(0); //Line for ratio plot
  pad7->Draw(); //Line for ratio plot
  pad7->cd(); //Line for ratio plot
  AllBKGSDRHJ->DrawCopy("e0 hist");
  AllBKGSDRHJ->SetMinimum(0.1);
  AllBKGSDRHJ->SetMaximum(1.9);  
  //BKGDRHJL->Draw("hist");
  //BKGDRHJL->SetMinimum(1);
  DRHjetsT[0]->Draw("E1 same");
  DRHjetsL[NOS-1]->Draw("hist same");
  gPad->SetLogy();
  BKGlegend->Draw();
  MyPlot->cd();
  TPad *pad8 = new TPad("pad4","pad4",0,0,1,0.3); //Line for ratio plot  
  pad8->SetTopMargin(0); //Line for ratio plot
  pad8->Draw(); //Line for ratio plot
  pad8->cd(); //Line for ratio plot
  AllBKGSDRHJ->Sumw2();
  AllBKGSDRHJ->SetStats(0);
  AllBKGSDRHJ->SetTitle(";#Delta R(bb);MC/Data");
  AllBKGSDRHJ->Divide(DRHjetsT[0]);
  AllBKGSDRHJ->SetMarkerStyle(21);
  AllBKGSDRHJ->Draw("ep");
  L1->DrawLine(0.5,1.3,1.2,1.3); L2->DrawLine(0.5,1,1.2,1); L3->DrawLine(0.5,0.7,1.2,0.7); 
  //BKGlegend->Draw();
  MyPlot->cd();  
  MyPlot->Update();
  cout << "HERE1" << endl;
  ///////////////
  ///////////////
  ///////////////
  MyPlot->Clear();
  MyPlot->cd(1);
  ps->NewPage();
  TPad *pad9 = new TPad("pad3","pad3",0,0.3,1,1); //Line for ratio plot  
  pad9->SetBottomMargin(0); //Line for ratio plot
  pad9->Draw(); //Line for ratio plot
  pad9->cd(); //Line for ratio plot
  AllBKGSTPM->DrawCopy("e0 hist");
  AllBKGSTPM->SetMinimum(0.1);
  AllBKGSTPM->SetMaximum(1.9);  
  BKGTPML->Add(TprimeHistosL[NOS-1]);
  BKGTPML->Draw("hist");
  //BKGTPML->SetMinimum(1);
  //TprimeHistosT[0]->Draw("E1 same");
  TprimeHistosL[NOS-1]->Draw("hist same");
  //gPad->SetLogy();
  BKGlegend2->Draw();
  MyPlot->cd();
  MyPlot->Update();
  /*TPad *pad10 = new TPad("pad4","pad4",0,0,1,0.3); //Line for ratio plot  
  pad10->SetTopMargin(0); //Line for ratio plot
  pad10->Draw(); //Line for ratio plot
  pad10->cd(); //Line for ratio plot
  AllBKGSTPM->Sumw2();
  AllBKGSTPM->SetStats(0);
  AllBKGSTPM->SetTitle(";M(5j) GeV; MC/Data");
  AllBKGSTPM->Divide(TprimeHistosT[0]);
  AllBKGSTPM->SetMarkerStyle(21);
  AllBKGSTPM->Draw("ep");
  L1->DrawLine(400,1.3,1600,1.3); L2->DrawLine(400,1,1600,1); L3->DrawLine(400,0.7,1600,0.7); 
  //BKGlegend->Draw();
  MyPlot->cd();  
  MyPlot->Update();
  cout << "HERE1" << endl; */ 
  ///////////////
  ///////////////
  ///////////////
  MyPlot->Clear();
  MyPlot->cd(1);
  ps->NewPage();
  TPad *pad11 = new TPad("pad3","pad3",0,0.3,1,1); //Line for ratio plot  
  pad11->SetBottomMargin(0); //Line for ratio plot
  pad11->Draw(); //Line for ratio plot
  pad11->cd(); //Line for ratio plot
  AllBKGSDRWJ->DrawCopy("e0 hist");
  AllBKGSDRWJ->SetMinimum(0.1);
  AllBKGSDRWJ->SetMaximum(1.9);  
  BKGDRWJ->Draw("hist");
  BKGDRWJ->SetMinimum(1);
  DRWjetsT[0]->Draw("E1 same");
  DRWjetsL[NOS-1]->Draw("hist same");
  gPad->SetLogy();
  BKGlegend->Draw();
  MyPlot->cd();
  TPad *pad12 = new TPad("pad4","pad4",0,0,1,0.3); //Line for ratio plot  
  pad12->SetTopMargin(0); //Line for ratio plot
  pad12->Draw(); //Line for ratio plot
  pad12->cd(); //Line for ratio plot
  AllBKGSDRWJ->Sumw2();
  AllBKGSDRWJ->SetStats(0);
  AllBKGSDRWJ->SetTitle(";#Delta R(jj);MC/Data");
  AllBKGSDRWJ->Divide(DRWjetsT[0]);
  AllBKGSDRWJ->SetMarkerStyle(21);
  AllBKGSDRWJ->Draw("ep");
  L1->DrawLine(0.5,1.3,3.0,1.3); L2->DrawLine(0.5,1,3.0,1); L3->DrawLine(0.5,0.7,3.0,0.7); 
  //BKGlegend->Draw();
  MyPlot->cd();  
  MyPlot->Update();
  ///////////////
  ///////////////
  ///////////////
  MyPlot->Clear();
  MyPlot->cd(1);
  ps->NewPage();
  TPad *pad13 = new TPad("pad3","pad3",0,0.3,1,1); //Line for ratio plot  
  pad13->SetBottomMargin(0); //Line for ratio plot
  pad13->Draw(); //Line for ratio plot
  pad13->cd(); //Line for ratio plot
  AllBKGSHpt->DrawCopy("e0 hist");
  AllBKGSHpt->SetMinimum(0.1);
  AllBKGSHpt->SetMaximum(1.9);  
  BKGHPT->Draw("hist");
  BKGHPT->SetMinimum(1);
  HptT[0]->Draw("E1 same");
  HptL[NOS-1]->Draw("hist same");
  //gPad->SetLogy();
  BKGlegend->Draw();
  MyPlot->cd();
  TPad *pad14 = new TPad("pad4","pad4",0,0,1,0.3); //Line for ratio plot  
  pad14->SetTopMargin(0); //Line for ratio plot
  pad14->Draw(); //Line for ratio plot
  pad14->cd(); //Line for ratio plot
  AllBKGSHpt->Sumw2();
  AllBKGSHpt->SetStats(0);
  AllBKGSHpt->SetTitle(";p_{T}^{H} GeV;MC/Data");
  AllBKGSHpt->Divide(HptT[0]);
  AllBKGSHpt->SetMarkerStyle(21);
  AllBKGSHpt->Draw("ep");
  L1->DrawLine(100,1.3,800,1.3); L2->DrawLine(100,1,800,1); L3->DrawLine(100,0.7,800,0.7); 
  //BKGlegend->Draw();
  MyPlot->cd();  
  MyPlot->Update();
  ///////////////
  ///////////////
  ///////////////
  MyPlot->Clear();
  MyPlot->cd(1);
  ps->NewPage();
  TPad *pad15 = new TPad("pad3","pad3",0,0.3,1,1); //Line for ratio plot  
  pad15->SetBottomMargin(0); //Line for ratio plot
  pad15->Draw(); //Line for ratio plot
  pad15->cd(); //Line for ratio plot
  AllBKGSTpt->DrawCopy("e0 hist");
  AllBKGSTpt->SetMinimum(0.1);
  AllBKGSTpt->SetMaximum(1.9);  
  BKGTPT->Draw("hist");
  BKGTPT->SetMinimum(1);
  TptT[0]->Draw("E1 same");
  TptL[NOS-1]->Draw("hist same");
  //gPad->SetLogy();
  BKGlegend->Draw();
  MyPlot->cd();
  TPad *pad16 = new TPad("pad4","pad4",0,0,1,0.3); //Line for ratio plot  
  pad16->SetTopMargin(0); //Line for ratio plot
  pad16->Draw(); //Line for ratio plot
  pad16->cd(); //Line for ratio plot
  AllBKGSTpt->Sumw2();
  AllBKGSTpt->SetStats(0);
  AllBKGSTpt->SetTitle(";p_{T}^{t} GeV;MC/Data");
  AllBKGSTpt->Divide(TptT[0]);
  AllBKGSTpt->SetMarkerStyle(21);
  AllBKGSTpt->Draw("ep");
  L1->DrawLine(100,1.3,800,1.3); L2->DrawLine(100,1,800,1); L3->DrawLine(100,0.7,800,0.7); 
  //BKGlegend->Draw();
  MyPlot->cd();  
  MyPlot->Update();
  ///////////////
  ///////////////
  ///////////////
  MyPlot->Clear();
  MyPlot->cd(1);
  ps->NewPage();
  TPad *pad17 = new TPad("pad3","pad3",0,0.3,1,1); //Line for ratio plot  
  pad17->SetBottomMargin(0); //Line for ratio plot
  pad17->Draw(); //Line for ratio plot
  pad17->cd(); //Line for ratio plot
  AllBKGSDRWH->DrawCopy("e0 hist");
  AllBKGSDRWH->SetMinimum(0.1);
  AllBKGSDRWH->SetMaximum(1.9);  
  BKGDRWH->Draw("hist");
  BKGDRWH->SetMinimum(1);
  DRWHT[0]->Draw("E1 same");
  DRWHL[NOS-1]->Draw("hist same");
  gPad->SetLogy();
  BKGlegend->Draw();
  MyPlot->cd();
  TPad *pad18 = new TPad("pad4","pad4",0,0,1,0.3); //Line for ratio plot  
  pad18->SetTopMargin(0); //Line for ratio plot
  pad18->Draw(); //Line for ratio plot
  pad18->cd(); //Line for ratio plot
  AllBKGSDRWH->Sumw2();
  AllBKGSDRWH->SetStats(0);
  AllBKGSDRWH->SetTitle(";#Delta R(WH);MC/Data");
  AllBKGSDRWH->Divide(DRWHT[0]);
  AllBKGSDRWH->SetMarkerStyle(21);
  AllBKGSDRWH->Draw("ep");
  L1->DrawLine(2.0,1.3,4.0,1.3); L2->DrawLine(2.0,1,4.0,1); L3->DrawLine(2.0,0.7,4.0,0.7); 
  //BKGlegend->Draw();
  MyPlot->cd();  
  MyPlot->Update();
  ///////////////
  ///////////////
  ///////////////
  MyPlot->Clear();
  MyPlot->cd(1);
  ps->NewPage();
  TPad *pad19 = new TPad("pad3","pad3",0,0.3,1,1); //Line for ratio plot  
  pad19->SetBottomMargin(0); //Line for ratio plot
  pad19->Draw(); //Line for ratio plot
  pad19->cd(); //Line for ratio plot
  AllBKGSDPHJ->DrawCopy("e0 hist");
  AllBKGSDPHJ->SetMinimum(0.1);
  AllBKGSDPHJ->SetMaximum(1.9);  
  BKGDPHJ->Draw("hist");
  BKGDPHJ->SetMinimum(1);
  DPHjetsT[0]->Draw("E1 same");
  DPHjetsL[NOS-1]->Draw("hist same");
  gPad->SetLogy();
  BKGlegend->Draw();
  MyPlot->cd();
  TPad *pad20 = new TPad("pad4","pad4",0,0,1,0.3); //Line for ratio plot  
  pad20->SetTopMargin(0); //Line for ratio plot
  pad20->Draw(); //Line for ratio plot
  pad20->cd(); //Line for ratio plot
  AllBKGSDPHJ->Sumw2();
  AllBKGSDPHJ->SetStats(0);
  AllBKGSDPHJ->SetTitle(";#Delta #phi(jj)^{H};MC/Data");
  AllBKGSDPHJ->Divide(DPHjetsT[0]);
  AllBKGSDPHJ->SetMarkerStyle(21);
  AllBKGSDPHJ->Draw("ep");
  L1->DrawLine(0.0,1.3,3.0,1.3); L2->DrawLine(0.0,1,3.0,1); L3->DrawLine(0.0,0.7,3.0,0.7); 
  //BKGlegend->Draw();
  MyPlot->cd();  
  MyPlot->Update();
  ///////////////
  ///////////////
  ///////////////
  MyPlot->Clear();
  MyPlot->cd(1);
  ps->NewPage();
  TPad *pad21 = new TPad("pad3","pad3",0,0.3,1,1); //Line for ratio plot  
  pad21->SetBottomMargin(0); //Line for ratio plot
  pad21->Draw(); //Line for ratio plot
  pad21->cd(); //Line for ratio plot
  AllBKGSDPWJ->DrawCopy("e0 hist");
  AllBKGSDPWJ->SetMinimum(0.1);
  AllBKGSDPWJ->SetMaximum(1.9);  
  BKGDPWJ->Draw("hist");
  BKGDPWJ->SetMinimum(1);
  DPWjetsT[0]->Draw("E1 same");
  DPWjetsL[NOS-1]->Draw("hist same");
  gPad->SetLogy();
  BKGlegend->Draw();
  MyPlot->cd();
  TPad *pad22 = new TPad("pad4","pad4",0,0,1,0.3); //Line for ratio plot  
  pad22->SetTopMargin(0); //Line for ratio plot
  pad22->Draw(); //Line for ratio plot
  pad22->cd(); //Line for ratio plot
  AllBKGSDPWJ->Sumw2();
  AllBKGSDPWJ->SetStats(0);
  AllBKGSDPWJ->SetTitle(";#Delta #phi(jj)^{W};MC/Data");
  AllBKGSDPWJ->Divide(DPWjetsT[0]);
  AllBKGSDPWJ->SetMarkerStyle(21);
  AllBKGSDPWJ->Draw("ep");
  L1->DrawLine(0.0,1.3,3.0,1.3); L2->DrawLine(0.0,1,3.0,1); L3->DrawLine(0.0,0.7,3.0,0.7); 
  //BKGlegend->Draw();
  MyPlot->cd();  
  MyPlot->Update();
  ///////////////
  ///////////////
  ///////////////
  MyPlot->Clear();
  MyPlot->cd(1);
  ps->NewPage();
  TPad *pad23 = new TPad("pad3","pad3",0,0.3,1,1); //Line for ratio plot  
  pad23->SetBottomMargin(0); //Line for ratio plot
  pad23->Draw(); //Line for ratio plot
  pad23->cd(); //Line for ratio plot
  AllBKGSDPTJ->DrawCopy("e0 hist");
  AllBKGSDPTJ->SetMinimum(0.1);
  AllBKGSDPTJ->SetMaximum(1.9);  
  BKGDPTJ->Draw("hist");
  BKGDPTJ->SetMinimum(1);
  DPTjetsT[0]->Draw("E1 same");
  DPTjetsL[NOS-1]->Draw("hist same");
  gPad->SetLogy();
  BKGlegend->Draw();
  MyPlot->cd();
  TPad *pad24 = new TPad("pad4","pad4",0,0,1,0.3); //Line for ratio plot  
  pad24->SetTopMargin(0); //Line for ratio plot
  pad24->Draw(); //Line for ratio plot
  pad24->cd(); //Line for ratio plot
  AllBKGSDPTJ->Sumw2();
  AllBKGSDPTJ->SetStats(0);
  AllBKGSDPTJ->SetTitle(";#Delta #phi(Wj)^{t};MC/Data");
  AllBKGSDPTJ->Divide(DPWjetsT[0]);
  AllBKGSDPTJ->SetMarkerStyle(21);
  AllBKGSDPTJ->Draw("ep");
  L1->DrawLine(0.0,1.3,3.0,1.3); L2->DrawLine(0.0,1,3.0,1); L3->DrawLine(0.0,0.7,3.0,0.7); 
  //BKGlegend->Draw();
  MyPlot->cd();  
  MyPlot->Update();
  ///////////////
  ///////////////
  ///////////////
  MyPlot->Clear();
  MyPlot->cd(1);
  ps->NewPage();
  TPad *pad25 = new TPad("pad3","pad3",0,0.3,1,1); //Line for ratio plot  
  pad25->SetBottomMargin(0); //Line for ratio plot
  pad25->Draw(); //Line for ratio plot
  pad25->cd(); //Line for ratio plot
  AllBKGSHM->DrawCopy("e0 hist");
  AllBKGSHM->SetMinimum(0.1);
  AllBKGSHM->SetMaximum(1.9);  
  BKGHM->Draw("hist");
  BKGHM->SetMinimum(1);
  HiggsMassT[0]->Draw("E1 same");
  HiggsMassL[NOS-1]->Draw("hist same");
  //gPad->SetLogy();
  BKGlegend->Draw();
  MyPlot->cd();
  TPad *pad26 = new TPad("pad4","pad4",0,0,1,0.3); //Line for ratio plot  
  pad26->SetTopMargin(0); //Line for ratio plot
  pad26->Draw(); //Line for ratio plot
  pad26->cd(); //Line for ratio plot
  AllBKGSHM->Sumw2();
  AllBKGSHM->SetStats(0);
  AllBKGSHM->SetTitle(";M_{H} GeV;MC/Data");
  AllBKGSHM->Divide(HiggsMassT[0]);
  AllBKGSHM->SetMarkerStyle(21);
  AllBKGSHM->Draw("ep");
  L1->DrawLine(60,1.3,180,1.3); L2->DrawLine(60,1,180,1); L3->DrawLine(60,0.7,180,0.7); 
  //BKGlegend->Draw();
  MyPlot->cd();  
  MyPlot->Update();
  ///////////////
  ///////////////
  ///////////////
  MyPlot->Clear();
  MyPlot->cd(1);
  ps->NewPage();
  TPad *pad27 = new TPad("pad3","pad3",0,0.3,1,1); //Line for ratio plot  
  pad27->SetBottomMargin(0); //Line for ratio plot
  pad27->Draw(); //Line for ratio plot
  pad27->cd(); //Line for ratio plot
  AllBKGSRelHT->DrawCopy("e0 hist");
  AllBKGSRelHT->SetMinimum(0.1);
  AllBKGSRelHT->SetMaximum(1.9);  
  BKGRHT->Draw("hist");
  BKGRHT->SetMinimum(1);
  RelHTT[0]->Draw("E1 same");
  RelHTL[NOS-1]->Draw("hist same");
  gPad->SetLogy();
  BKGlegend->Draw();
  MyPlot->cd();
  TPad *pad28 = new TPad("pad4","pad4",0,0,1,0.3); //Line for ratio plot  
  pad28->SetTopMargin(0); //Line for ratio plot
  pad28->Draw(); //Line for ratio plot
  pad28->cd(); //Line for ratio plot
  AllBKGSRelHT->Sumw2();
  AllBKGSRelHT->SetStats(0);
  AllBKGSRelHT->SetTitle(";(M_{H}+M_{t})/H_{T};MC/Data");
  AllBKGSRelHT->Divide(RelHTT[0]);
  AllBKGSRelHT->SetMarkerStyle(21);
  AllBKGSRelHT->Draw("ep");
  L1->DrawLine(0.5,1.3,1.0,1.3); L2->DrawLine(0.5,1,1.0,1); L3->DrawLine(0.5,0.7,1.0,0.7); 
  //BKGlegend->Draw();
  MyPlot->cd();  
  MyPlot->Update();
  ///////////////
  ///////////////
  ///////////////
  MyPlot->Clear();
  MyPlot->cd(1);
  ps->NewPage();
  TPad *pad29 = new TPad("pad3","pad3",0,0.3,1,1); //Line for ratio plot  
  pad29->SetBottomMargin(0); //Line for ratio plot
  pad29->Draw(); //Line for ratio plot
  pad29->cd(); //Line for ratio plot
  AllBKGSRM->DrawCopy("e0 hist");
  AllBKGSRM->SetMinimum(0.1);
  AllBKGSRM->SetMaximum(1.9);  
  BKGRM->Draw("hist");
  BKGRM->SetMinimum(1);
  RelativeMassT[0]->Draw("E1 same");
  RelativeMassL[NOS-1]->Draw("hist same");
  gPad->SetLogy();
  BKGlegend->Draw();
  MyPlot->cd();
  TPad *pad30 = new TPad("pad4","pad4",0,0,1,0.3); //Line for ratio plot  
  pad30->SetTopMargin(0); //Line for ratio plot
  pad30->Draw(); //Line for ratio plot
  pad30->cd(); //Line for ratio plot
  AllBKGSRM->Sumw2();
  AllBKGSRM->SetStats(0);
  AllBKGSRM->SetTitle(";(M_{H}+M_{t})/M_{j}^{5};MC/Data");
  AllBKGSRM->Divide(RelativeMassT[0]);
  AllBKGSRM->SetMarkerStyle(21);
  AllBKGSRM->Draw("ep");
  L1->DrawLine(0.2,1.3,0.6,1.3); L2->DrawLine(0.2,1,0.6,1); L3->DrawLine(0.2,0.7,0.6,0.7); 
  //BKGlegend->Draw();
  MyPlot->cd();  
  MyPlot->Update();
  ///////////////
  ///////////////
  ///////////////
  MyPlot->Clear();
  MyPlot->cd(1);
  ps->NewPage();
  TPad *pad31 = new TPad("pad3","pad3",0,0.3,1,1); //Line for ratio plot  
  pad31->SetBottomMargin(0); //Line for ratio plot
  pad31->Draw(); //Line for ratio plot
  pad31->cd(); //Line for ratio plot
  AllBKGSTM->DrawCopy("e0 hist");
  AllBKGSTM->SetMinimum(0.1);
  AllBKGSTM->SetMaximum(1.9);  
  BKGTM->Draw("hist");
  BKGTM->SetMinimum(1);
  TopMassT[0]->Draw("E1 same");
  TopMassL[NOS-1]->Draw("hist same");
  //gPad->SetLogy();
  BKGlegend->Draw();
  MyPlot->cd();
  TPad *pad32 = new TPad("pad4","pad4",0,0,1,0.3); //Line for ratio plot  
  pad32->SetTopMargin(0); //Line for ratio plot
  pad32->Draw(); //Line for ratio plot
  pad32->cd(); //Line for ratio plot
  AllBKGSTM->Sumw2();
  AllBKGSTM->SetStats(0);
  AllBKGSTM->SetTitle(";M_{t} GeV;MC/Data");
  AllBKGSTM->Divide(TopMassT[0]);
  AllBKGSTM->SetMarkerStyle(21);
  AllBKGSTM->Draw("ep");
  L1->DrawLine(100,1.3,300,1.3); L2->DrawLine(100,1,300,1); L3->DrawLine(100,0.7,300,0.7); 
  //BKGlegend->Draw();
  MyPlot->cd();  
  MyPlot->Update();
  ///////////////
  ///////////////
  ///////////////
  MyPlot->Clear();
  MyPlot->cd(1);
  ps->NewPage();
  TPad *pad33 = new TPad("pad3","pad3",0,0.3,1,1); //Line for ratio plot  
  pad33->SetBottomMargin(0); //Line for ratio plot
  pad33->Draw(); //Line for ratio plot
  pad33->cd(); //Line for ratio plot
  AllBKGSTFH->DrawCopy("e0 hist");
  AllBKGSTFH->SetMinimum(0.1);
  AllBKGSTFH->SetMaximum(1.9);  
  BKGTFH->Draw("hist");
  BKGTFH->SetMinimum(1);
  TopFromHiggsT[0]->Draw("E1 same");
  TopFromHiggsL[NOS-1]->Draw("hist same");
  //gPad->SetLogy();
  BKGlegend->Draw();
  MyPlot->cd();
  TPad *pad34 = new TPad("pad4","pad4",0,0,1,0.3); //Line for ratio plot  
  pad34->SetTopMargin(0); //Line for ratio plot
  pad34->Draw(); //Line for ratio plot
  pad34->cd(); //Line for ratio plot
  AllBKGSTFH->Sumw2();
  AllBKGSTFH->SetStats(0);
  AllBKGSTFH->SetTitle(";M_{t^{2nd}} GeV;MC/Data");
  AllBKGSTFH->Divide(TopFromHiggsT[0]);
  AllBKGSTFH->SetMarkerStyle(21);
  AllBKGSTFH->Draw("ep");
  L1->DrawLine(100,1.3,800,1.3); L2->DrawLine(100,1,800,1); L3->DrawLine(100,0.7,800,0.7); 
  //BKGlegend->Draw();
  MyPlot->cd();  
  MyPlot->Update();
  ///////////////
  ///////////////
  ///////////////
  MyPlot->Clear();
  MyPlot->cd(1);
  ps->NewPage();
  TPad *pad35 = new TPad("pad3","pad3",0,0.3,1,1); //Line for ratio plot  
  pad35->SetBottomMargin(0); //Line for ratio plot
  pad35->Draw(); //Line for ratio plot
  pad35->cd(); //Line for ratio plot
  AllBKGSJET1ETA->DrawCopy("e0 hist");
  AllBKGSJET1ETA->SetMinimum(0.1);
  AllBKGSJET1ETA->SetMaximum(1.9);
  BKGJET1ETA->Draw("hist");
  BKGJET1ETA->SetMinimum(1);
  JET1ETAT[0]->Draw("E1 same");
  JET1ETAL[NOS-1]->Draw("hist same");
  //gPad->SetLogy();
  BKGlegend->Draw();
  MyPlot->cd();
  TPad *pad36 = new TPad("pad4","pad4",0,0,1,0.3); //Line for ratio plot  
  pad36->SetTopMargin(0); //Line for ratio plot
  pad36->Draw(); //Line for ratio plot
  pad36->cd(); //Line for ratio plot
  AllBKGSJET1ETA->Sumw2();
  AllBKGSJET1ETA->SetStats(0);
  AllBKGSJET1ETA->SetTitle(";M_{t^{2nd}} GeV;MC/Data");
  AllBKGSJET1ETA->Divide(JET1ETAT[0]);
  AllBKGSJET1ETA->SetMarkerStyle(21);
  AllBKGSJET1ETA->Draw("ep");
  L1->DrawLine(-5,1.3,5,1.3); L2->DrawLine(-5,1,5,1); L3->DrawLine(-5,0.7,5,0.7); 
  //BKGlegend->Draw();
  MyPlot->cd();
  MyPlot->Update();
  ///////////////
  ///////////////
  ///////////////
  MyPlot->Clear();
  MyPlot->cd(1);
  ps->NewPage();
  TPad *pad37 = new TPad("pad3","pad3",0,0.3,1,1); //Line for ratio plot  
  pad37->SetBottomMargin(0); //Line for ratio plot
  pad37->Draw(); //Line for ratio plot
  pad37->cd(); //Line for ratio plot
  AllBKGSJET2ETA->DrawCopy("e0 hist");
  AllBKGSJET2ETA->SetMinimum(0.1);
  AllBKGSJET2ETA->SetMaximum(1.9);
  BKGJET2ETA->Draw("hist");
  BKGJET2ETA->SetMinimum(1);
  JET2ETAT[0]->Draw("E1 same");
  JET2ETAL[NOS-1]->Draw("hist same");
  //gPad->SetLogy();
  BKGlegend->Draw();
  MyPlot->cd();
  TPad *pad38 = new TPad("pad4","pad4",0,0,1,0.3); //Line for ratio plot  
  pad38->SetTopMargin(0); //Line for ratio plot
  pad38->Draw(); //Line for ratio plot
  pad38->cd(); //Line for ratio plot
  AllBKGSJET2ETA->Sumw2();
  AllBKGSJET2ETA->SetStats(0);
  AllBKGSJET2ETA->SetTitle(";M_{t^{2nd}} GeV;MC/Data");
  AllBKGSJET2ETA->Divide(JET2ETAT[0]);
  AllBKGSJET2ETA->SetMarkerStyle(21);
  AllBKGSJET2ETA->Draw("ep");
  L1->DrawLine(-5,1.3,5,1.3); L2->DrawLine(-5,1,5,1); L3->DrawLine(-5,0.7,5,0.7); 
  //BKGlegend->Draw();
  MyPlot->cd();
  MyPlot->Update();
  ///////////////
  ///////////////
  ///////////////
  MyPlot->Clear();
  MyPlot->cd(1);
  ps->NewPage();
  TPad *pad39 = new TPad("pad3","pad3",0,0.3,1,1); //Line for ratio plot  
  pad39->SetBottomMargin(0); //Line for ratio plot
  pad39->Draw(); //Line for ratio plot
  pad39->cd(); //Line for ratio plot
  AllBKGSJET3ETA->DrawCopy("e0 hist");
  AllBKGSJET3ETA->SetMinimum(0.1);
  AllBKGSJET3ETA->SetMaximum(1.9);
  BKGJET3ETA->Draw("hist");
  BKGJET3ETA->SetMinimum(1);
  JET3ETAT[0]->Draw("E1 same");
  JET3ETAL[NOS-1]->Draw("hist same");
  //gPad->SetLogy();
  BKGlegend->Draw();
  MyPlot->cd();
  TPad *pad40 = new TPad("pad4","pad4",0,0,1,0.3); //Line for ratio plot  
  pad40->SetTopMargin(0); //Line for ratio plot
  pad40->Draw(); //Line for ratio plot
  pad40->cd(); //Line for ratio plot
  AllBKGSJET3ETA->Sumw2();
  AllBKGSJET3ETA->SetStats(0);
  AllBKGSJET3ETA->SetTitle(";M_{t^{2nd}} GeV;MC/Data");
  AllBKGSJET3ETA->Divide(JET3ETAT[0]);
  AllBKGSJET3ETA->SetMarkerStyle(21);
  AllBKGSJET3ETA->Draw("ep");
  L1->DrawLine(-5,1.3,5,1.3); L2->DrawLine(-5,1,5,1); L3->DrawLine(-5,0.7,5,0.7); 
  //BKGlegend->Draw();
  MyPlot->cd();
  MyPlot->Update();
  ///////////////
  ///////////////
  ///////////////
  MyPlot->Clear();
  MyPlot->cd(1);
  ps->NewPage();
  TPad *pad41 = new TPad("pad3","pad3",0,0.3,1,1); //Line for ratio plot  
  pad41->SetBottomMargin(0); //Line for ratio plot
  pad41->Draw(); //Line for ratio plot
  pad41->cd(); //Line for ratio plot
  AllBKGSJET4ETA->DrawCopy("e0 hist");
  AllBKGSJET4ETA->SetMinimum(0.1);
  AllBKGSJET4ETA->SetMaximum(1.9);
  BKGJET4ETA->Draw("hist");
  BKGJET4ETA->SetMinimum(1);
  JET4ETAT[0]->Draw("E1 same");
  JET4ETAL[NOS-1]->Draw("hist same");
  //gPad->SetLogy();
  BKGlegend->Draw();
  MyPlot->cd();
  TPad *pad42 = new TPad("pad4","pad4",0,0,1,0.3); //Line for ratio plot  
  pad42->SetTopMargin(0); //Line for ratio plot
  pad42->Draw(); //Line for ratio plot
  pad42->cd(); //Line for ratio plot
  AllBKGSJET4ETA->Sumw2();
  AllBKGSJET4ETA->SetStats(0);
  AllBKGSJET4ETA->SetTitle(";M_{t^{2nd}} GeV;MC/Data");
  AllBKGSJET4ETA->Divide(JET4ETAT[0]);
  AllBKGSJET4ETA->SetMarkerStyle(21);
  AllBKGSJET4ETA->Draw("ep");
  L1->DrawLine(-5,1.3,5,1.3); L2->DrawLine(-5,1,5,1); L3->DrawLine(-5,0.7,5,0.7); 
  //BKGlegend->Draw();
  MyPlot->cd();
  MyPlot->Update();
  ///////////////
  ///////////////
  ///////////////
  MyPlot->Clear();
  MyPlot->cd(1);
  ps->NewPage();
  TPad *pad43 = new TPad("pad3","pad3",0,0.3,1,1); //Line for ratio plot  
  pad43->SetBottomMargin(0); //Line for ratio plot
  pad43->Draw(); //Line for ratio plot
  pad43->cd(); //Line for ratio plot
  AllBKGSJET5ETA->DrawCopy("e0 hist");
  AllBKGSJET5ETA->SetMinimum(0.1);
  AllBKGSJET5ETA->SetMaximum(1.9);
  BKGJET5ETA->Draw("hist");
  BKGJET5ETA->SetMinimum(1);
  JET5ETAT[0]->Draw("E1 same");
  JET5ETAL[NOS-1]->Draw("hist same");
  //gPad->SetLogy();
  BKGlegend->Draw();
  MyPlot->cd();
  TPad *pad44 = new TPad("pad4","pad4",0,0,1,0.3); //Line for ratio plot  
  pad44->SetTopMargin(0); //Line for ratio plot
  pad44->Draw(); //Line for ratio plot
  pad44->cd(); //Line for ratio plot
  AllBKGSJET5ETA->Sumw2();
  AllBKGSJET5ETA->SetStats(0);
  AllBKGSJET5ETA->SetTitle(";M_{t^{2nd}} GeV;MC/Data");
  AllBKGSJET5ETA->Divide(JET5ETAT[0]);
  AllBKGSJET5ETA->SetMarkerStyle(21);
  AllBKGSJET5ETA->Draw("ep");
  L1->DrawLine(-5,1.3,5,1.3); L2->DrawLine(-5,1,5,1); L3->DrawLine(-5,0.7,5,0.7); 
  //BKGlegend->Draw();
  MyPlot->cd();
  MyPlot->Update();
  ///////////////
  ///////////////
  ///////////////
  MyPlot->Clear();
  MyPlot->cd(1);
  ps->NewPage();
  TPad *pad45 = new TPad("pad3","pad3",0,0.3,1,1); //Line for ratio plot  
  pad45->SetBottomMargin(0); //Line for ratio plot
  pad45->Draw(); //Line for ratio plot
  pad45->cd(); //Line for ratio plot
  AllBKGSJET6ETA->DrawCopy("e0 hist");
  AllBKGSJET6ETA->SetMinimum(0.1);
  AllBKGSJET6ETA->SetMaximum(1.9);
  BKGJET6ETA->Draw("hist");
  BKGJET6ETA->SetMinimum(1);
  JET6ETAT[0]->Draw("E1 same");
  JET6ETAL[NOS-1]->Draw("hist same");
  //gPad->SetLogy();
  BKGlegend->Draw();
  MyPlot->cd();
  TPad *pad46 = new TPad("pad4","pad4",0,0,1,0.3); //Line for ratio plot  
  pad46->SetTopMargin(0); //Line for ratio plot
  pad46->Draw(); //Line for ratio plot
  pad46->cd(); //Line for ratio plot
  AllBKGSJET6ETA->Sumw2();
  AllBKGSJET6ETA->SetStats(0);
  AllBKGSJET6ETA->SetTitle(";M_{t^{2nd}} GeV;MC/Data");
  AllBKGSJET6ETA->Divide(JET6ETAT[0]);
  AllBKGSJET6ETA->SetMarkerStyle(21);
  AllBKGSJET6ETA->Draw("ep");
  L1->DrawLine(-5,1.3,5,1.3); L2->DrawLine(-5,1,5,1); L3->DrawLine(-5,0.7,5,0.7); 
  //BKGlegend->Draw();
  MyPlot->cd();
  MyPlot->Update();
  ///////////////
  ///////////////
  ///////////////
  MyPlot->Clear();
  MyPlot->cd(1);
  ps->NewPage();
  TPad *pad47 = new TPad("pad3","pad3",0,0.3,1,1); //Line for ratio plot  
  pad47->SetBottomMargin(0); //Line for ratio plot
  pad47->Draw(); //Line for ratio plot
  pad47->cd(); //Line for ratio plot
  AllBKGSJET1PHI->DrawCopy("e0 hist");
  AllBKGSJET1PHI->SetMinimum(0.1);
  AllBKGSJET1PHI->SetMaximum(1.9);
  BKGJET1PHI->Draw("hist");
  BKGJET1PHI->SetMinimum(1);
  JET1PHIT[0]->Draw("E1 same");
  JET1PHIL[NOS-1]->Draw("hist same");
  //gPad->SetLogy();
  BKGlegend->Draw();
  MyPlot->cd();
  TPad *pad48 = new TPad("pad4","pad4",0,0,1,0.3); //Line for ratio plot  
  pad48->SetTopMargin(0); //Line for ratio plot
  pad48->Draw(); //Line for ratio plot
  pad48->cd(); //Line for ratio plot
  AllBKGSJET1PHI->Sumw2();
  AllBKGSJET1PHI->SetStats(0);
  AllBKGSJET1PHI->SetTitle(";M_{t^{2nd}} GeV;MC/Data");
  AllBKGSJET1PHI->Divide(JET1PHIT[0]);
  AllBKGSJET1PHI->SetMarkerStyle(21);
  AllBKGSJET1PHI->Draw("ep");
  L1->DrawLine(-5,1.3,5,1.3); L2->DrawLine(-5,1,5,1); L3->DrawLine(-5,0.7,5,0.7); 
  //BKGlegend->Draw();
  MyPlot->cd();
  MyPlot->Update();
  ///////////////
  ///////////////
  ///////////////
  MyPlot->Clear();
  MyPlot->cd(1);
  ps->NewPage();
  TPad *pad49 = new TPad("pad3","pad3",0,0.3,1,1); //Line for ratio plot  
  pad49->SetBottomMargin(0); //Line for ratio plot
  pad49->Draw(); //Line for ratio plot
  pad49->cd(); //Line for ratio plot
  AllBKGSJET2PHI->DrawCopy("e0 hist");
  AllBKGSJET2PHI->SetMinimum(0.1);
  AllBKGSJET2PHI->SetMaximum(1.9);
  BKGJET2PHI->Draw("hist");
  BKGJET2PHI->SetMinimum(1);
  JET2PHIT[0]->Draw("E1 same");
  JET2PHIL[NOS-1]->Draw("hist same");
  //gPad->SetLogy();
  BKGlegend->Draw();
  MyPlot->cd();
  TPad *pad50 = new TPad("pad4","pad4",0,0,1,0.3); //Line for ratio plot  
  pad50->SetTopMargin(0); //Line for ratio plot
  pad50->Draw(); //Line for ratio plot
  pad50->cd(); //Line for ratio plot
  AllBKGSJET2PHI->Sumw2();
  AllBKGSJET2PHI->SetStats(0);
  AllBKGSJET2PHI->SetTitle(";M_{t^{2nd}} GeV;MC/Data");
  AllBKGSJET2PHI->Divide(JET2PHIT[0]);
  AllBKGSJET2PHI->SetMarkerStyle(21);
  AllBKGSJET2PHI->Draw("ep");
  L1->DrawLine(-5,1.3,5,1.3); L2->DrawLine(-5,1,5,1); L3->DrawLine(-5,0.7,5,0.7); 
  //BKGlegend->Draw();
  MyPlot->cd();
  MyPlot->Update();
  ///////////////
  ///////////////
  ///////////////
  MyPlot->Clear();
  MyPlot->cd(1);
  ps->NewPage();
  TPad *pad51 = new TPad("pad3","pad3",0,0.3,1,1); //Line for ratio plot  
  pad51->SetBottomMargin(0); //Line for ratio plot
  pad51->Draw(); //Line for ratio plot
  pad51->cd(); //Line for ratio plot
  AllBKGSJET3PHI->DrawCopy("e0 hist");
  AllBKGSJET3PHI->SetMinimum(0.1);
  AllBKGSJET3PHI->SetMaximum(1.9);
  BKGJET3PHI->Draw("hist");
  BKGJET3PHI->SetMinimum(1);
  JET3PHIT[0]->Draw("E1 same");
  JET3PHIL[NOS-1]->Draw("hist same");
  //gPad->SetLogy();
  BKGlegend->Draw();
  MyPlot->cd();
  TPad *pad52 = new TPad("pad4","pad4",0,0,1,0.3); //Line for ratio plot  
  pad52->SetTopMargin(0); //Line for ratio plot
  pad52->Draw(); //Line for ratio plot
  pad52->cd(); //Line for ratio plot
  AllBKGSJET3PHI->Sumw2();
  AllBKGSJET3PHI->SetStats(0);
  AllBKGSJET3PHI->SetTitle(";M_{t^{2nd}} GeV;MC/Data");
  AllBKGSJET3PHI->Divide(JET3PHIT[0]);
  AllBKGSJET3PHI->SetMarkerStyle(21);
  AllBKGSJET3PHI->Draw("ep");
  L1->DrawLine(-5,1.3,5,1.3); L2->DrawLine(-5,1,5,1); L3->DrawLine(-5,0.7,5,0.7); 
  //BKGlegend->Draw();
  MyPlot->cd();
  MyPlot->Update();
  ///////////////
  ///////////////
  ///////////////
  MyPlot->Clear();
  MyPlot->cd(1);
  ps->NewPage();
  TPad *pad53 = new TPad("pad3","pad3",0,0.3,1,1); //Line for ratio plot  
  pad53->SetBottomMargin(0); //Line for ratio plot
  pad53->Draw(); //Line for ratio plot
  pad53->cd(); //Line for ratio plot
  AllBKGSJET4PHI->DrawCopy("e0 hist");
  AllBKGSJET4PHI->SetMinimum(0.1);
  AllBKGSJET4PHI->SetMaximum(1.9);
  BKGJET4PHI->Draw("hist");
  BKGJET4PHI->SetMinimum(1);
  JET4PHIT[0]->Draw("E1 same");
  JET4PHIL[NOS-1]->Draw("hist same");
  //gPad->SetLogy();
  BKGlegend->Draw();
  MyPlot->cd();
  TPad *pad54 = new TPad("pad4","pad4",0,0,1,0.3); //Line for ratio plot  
  pad54->SetTopMargin(0); //Line for ratio plot
  pad54->Draw(); //Line for ratio plot
  pad54->cd(); //Line for ratio plot
  AllBKGSJET4PHI->Sumw2();
  AllBKGSJET4PHI->SetStats(0);
  AllBKGSJET4PHI->SetTitle(";M_{t^{2nd}} GeV;MC/Data");
  AllBKGSJET4PHI->Divide(JET4PHIT[0]);
  AllBKGSJET4PHI->SetMarkerStyle(21);
  AllBKGSJET4PHI->Draw("ep");
  L1->DrawLine(-5,1.3,5,1.3); L2->DrawLine(-5,1,5,1); L3->DrawLine(-5,0.7,5,0.7); 
  //BKGlegend->Draw();
  MyPlot->cd();
  MyPlot->Update();
  ///////////////
  ///////////////
  ///////////////
  MyPlot->Clear();
  MyPlot->cd(1);
  ps->NewPage();
  TPad *pad55 = new TPad("pad3","pad3",0,0.3,1,1); //Line for ratio plot  
  pad55->SetBottomMargin(0); //Line for ratio plot
  pad55->Draw(); //Line for ratio plot
  pad55->cd(); //Line for ratio plot
  AllBKGSJET5PHI->DrawCopy("e0 hist");
  AllBKGSJET5PHI->SetMinimum(0.1);
  AllBKGSJET5PHI->SetMaximum(1.9);
  BKGJET5PHI->Draw("hist");
  BKGJET5PHI->SetMinimum(1);
  JET5PHIT[0]->Draw("E1 same");
  JET5PHIL[NOS-1]->Draw("hist same");
  //gPad->SetLogy();
  BKGlegend->Draw();
  MyPlot->cd();
  TPad *pad56 = new TPad("pad4","pad4",0,0,1,0.3); //Line for ratio plot  
  pad56->SetTopMargin(0); //Line for ratio plot
  pad56->Draw(); //Line for ratio plot
  pad56->cd(); //Line for ratio plot
  AllBKGSJET5PHI->Sumw2();
  AllBKGSJET5PHI->SetStats(0);
  AllBKGSJET5PHI->SetTitle(";M_{t^{2nd}} GeV;MC/Data");
  AllBKGSJET5PHI->Divide(JET5PHIT[0]);
  AllBKGSJET5PHI->SetMarkerStyle(21);
  AllBKGSJET5PHI->Draw("ep");
  L1->DrawLine(-5,1.3,5,1.3); L2->DrawLine(-5,1,5,1); L3->DrawLine(-5,0.7,5,0.7); 
  //BKGlegend->Draw();
  MyPlot->cd();
  MyPlot->Update();
  ///////////////
  ///////////////
  ///////////////
  MyPlot->Clear();
  MyPlot->cd(1);
  ps->NewPage();
  TPad *pad57 = new TPad("pad3","pad3",0,0.3,1,1); //Line for ratio plot  
  pad57->SetBottomMargin(0); //Line for ratio plot
  pad57->Draw(); //Line for ratio plot
  pad57->cd(); //Line for ratio plot
  AllBKGSJET6PHI->DrawCopy("e0 hist");
  AllBKGSJET6PHI->SetMinimum(0.1);
  AllBKGSJET6PHI->SetMaximum(1.9);
  BKGJET6PHI->Draw("hist");
  BKGJET6PHI->SetMinimum(1);
  JET6PHIT[0]->Draw("E1 same");
  JET6PHIL[NOS-1]->Draw("hist same");
  //gPad->SetLogy();
  BKGlegend->Draw();
  MyPlot->cd();
  TPad *pad58 = new TPad("pad4","pad4",0,0,1,0.3); //Line for ratio plot  
  pad58->SetTopMargin(0); //Line for ratio plot
  pad58->Draw(); //Line for ratio plot
  pad58->cd(); //Line for ratio plot
  AllBKGSJET6PHI->Sumw2();
  AllBKGSJET6PHI->SetStats(0);
  AllBKGSJET6PHI->SetTitle(";M_{t^{2nd}} GeV;MC/Data");
  AllBKGSJET6PHI->Divide(JET6PHIT[0]);
  AllBKGSJET6PHI->SetMarkerStyle(21);
  AllBKGSJET6PHI->Draw("ep");
  L1->DrawLine(-5,1.3,5,1.3); L2->DrawLine(-5,1,5,1); L3->DrawLine(-5,0.7,5,0.7); 
  //BKGlegend->Draw();
  MyPlot->cd();
  MyPlot->Update();
  ///////////////
  ///////////////
  ///////////////
  MyPlot->Clear();
  MyPlot->cd(1);
  ps->NewPage();
  TPad *pad59 = new TPad("pad3","pad3",0,0.3,1,1); //Line for ratio plot  
  pad59->SetBottomMargin(0); //Line for ratio plot
  pad59->Draw(); //Line for ratio plot
  pad59->cd(); //Line for ratio plot
  AllBKGSJETMULTI->DrawCopy("e0 hist");
  AllBKGSJETMULTI->SetMinimum(0.1);
  AllBKGSJETMULTI->SetMaximum(1.9);
  BKGJETMULTI->Draw("hist");
  BKGJETMULTI->SetMinimum(1);
  JETMULTIT[0]->Draw("E1 same");
  JETMULTIL[NOS-1]->Draw("hist same");
  //gPad->SetLogy();
  BKGlegend->Draw();
  MyPlot->cd();
  TPad *pad60 = new TPad("pad4","pad4",0,0,1,0.3); //Line for ratio plot  
  pad60->SetTopMargin(0); //Line for ratio plot
  pad60->Draw(); //Line for ratio plot
  pad60->cd(); //Line for ratio plot
  AllBKGSJETMULTI->Sumw2();
  AllBKGSJETMULTI->SetStats(0);
  AllBKGSJETMULTI->SetTitle(";M_{t^{2nd}} GeV;MC/Data");
  AllBKGSJETMULTI->Divide(JETMULTIT[0]);
  AllBKGSJETMULTI->SetMarkerStyle(21);
  AllBKGSJETMULTI->Draw("ep");
  L1->DrawLine(0,1.3,20,1.3); L2->DrawLine(0,1,20,1); L3->DrawLine(0,0.7,20,0.7); 
  //BKGlegend->Draw();
  MyPlot->cd();
  MyPlot->Update();
  //////////////
  ps->Close();


  exit(0); 

}
