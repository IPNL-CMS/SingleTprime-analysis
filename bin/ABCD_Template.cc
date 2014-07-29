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
const int NOS=6;
//const TString Samples[NOS] = {"Diboson.root", "SingleTop.root", "TTJets.root", "Zjets.root", "Wjets.root", "QCD.root", "Signal.root"};
const TString Samples[NOS] = {"Diboson.root", "SingleTop.root", "TTJets.root", "DY.root", "QCD.root", "Signal.root"};
//const TString Samples[NOS] = {"Diboson.root", "SingleTop.root", "TTJets.root", "QCD.root", "Signal.root"};

int theNumber(TFile * theFile){
  //TH1F*histogram;
  for(int k=0;k<=21;k++){
    string name="TprimeMass";
    char kchar[20];sprintf(kchar,"%d",k);
    if((TH1F*)(theFile->Get((name+kchar).c_str()))) return k;
  }
  return -1;
}

void ABCD()
{
  TFile *CurrentFile[NOS];
  TH2F *HptTpt[6][NOS];
  TH1F *HiggsMassReversedHptTpt[6][NOS-1];
  TH1F *TprimeMassReversedHptTpt[6][NOS-1];
  TH1F *HiggsMassCuts[6][NOS-1];
  int Index[NOS];
  
  for (int i=0; i<NOS; i++)
    {
      CurrentFile[i] = new TFile(StorageDirPrefix + Samples[i], "READ");
      if ( CurrentFile[i]->IsOpen() ) printf( Samples[i] + " File opened successfully\n");
      int NumberSuffix=theNumber(CurrentFile[i]);
      cout << NumberSuffix << endl;
      Index[i]=NumberSuffix;
      if (i!=NOS-1)
	{
	  string A=Form("HMBE0%i",NumberSuffix);
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
	  TprimeMassReversedHptTpt[5][i]->SetDefaultSumw2(); TprimeMassReversedHptTpt[5][i]= (TH1F*)gDirectory->Get(A.c_str());
	  A=Form("HM0%i",NumberSuffix);
	  HiggsMassCuts[0][i]->SetDefaultSumw2(); HiggsMassCuts[0][i]= (TH1F*)gDirectory->Get(A.c_str());
	  A=Form("HM1%i",NumberSuffix);
	  HiggsMassCuts[1][i]->SetDefaultSumw2(); HiggsMassCuts[1][i]= (TH1F*)gDirectory->Get(A.c_str());
	  A=Form("HM2%i",NumberSuffix);
	  HiggsMassCuts[2][i]->SetDefaultSumw2(); HiggsMassCuts[2][i]= (TH1F*)gDirectory->Get(A.c_str());
	  A=Form("HM3%i",NumberSuffix);
	  HiggsMassCuts[3][i]->SetDefaultSumw2(); HiggsMassCuts[3][i]= (TH1F*)gDirectory->Get(A.c_str());
	  A=Form("HM4%i",NumberSuffix);
	  HiggsMassCuts[4][i]->SetDefaultSumw2(); HiggsMassCuts[4][i]= (TH1F*)gDirectory->Get(A.c_str());
	  A=Form("HM5%i",NumberSuffix);
	  HiggsMassCuts[5][i]->SetDefaultSumw2(); HiggsMassCuts[5][i]= (TH1F*)gDirectory->Get(A.c_str());
	}
      //Proving ACBD method
      string A=Form("HPtTPt0%i",NumberSuffix);
      HptTpt[0][i]= (TH2F*)gDirectory->Get(A.c_str()); //HptTpt[0][i]->SetDefaultSumw2(); 
      A=Form("HPtTPt1%i",NumberSuffix);
      HptTpt[1][i]= (TH2F*)gDirectory->Get(A.c_str()); //HptTpt[1][i]->SetDefaultSumw2(); 
      A=Form("HPtTPt2%i",NumberSuffix);
      HptTpt[2][i]= (TH2F*)gDirectory->Get(A.c_str()); //HptTpt[2][i]->SetDefaultSumw2(); 
      A=Form("HPtTPt3%i",NumberSuffix);
      HptTpt[3][i]= (TH2F*)gDirectory->Get(A.c_str()); //HptTpt[3][i]->SetDefaultSumw2(); 
      A=Form("HPtTPt4%i",NumberSuffix);
      HptTpt[4][i]= (TH2F*)gDirectory->Get(A.c_str()); //HptTpt[4][i]->SetDefaultSumw2(); 
      A=Form("HPtTPt5%i",NumberSuffix);
      HptTpt[5][i]= (TH2F*)gDirectory->Get(A.c_str()); //HptTpt[5][i]->SetDefaultSumw2();
      //CurrentFile[i]->Close();
    }
  
  
  for (int i=0; i<NOS; i++) {for (int j=0; j<6; j++) {cout << "Correlation factor on Hpt Tpt on cut before adding up->" << i << j << " " << HptTpt[j][i]->GetCorrelationFactor() << endl;}}

  THStack *BKGHMBE0= new THStack("BKGHMBE0", "BKG for HMBE0; M_{H}; Events");
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
  THStack *BKGTprimeMassBE5= new THStack("BKGTprimeMassBE5", "BKG for TprimeMassBE5; M_{5j}; Events");
  string A=Form("HPtTPt0%i",Index[0]);
  TH2F *FullHPtTPt0 = (TH2F*)gDirectory->Get(A.c_str()); //HptTpt[0][0]->Clone();
  //gPad->Close();
  A=Form("HPtTPt1%i",Index[0]);
  TH2F *FullHPtTPt1 = (TH2F*)gDirectory->Get(A.c_str()); 
  //gPad->Close();
  A=Form("HPtTPt2%i",Index[0]);
  TH2F *FullHPtTPt2 = (TH2F*)gDirectory->Get(A.c_str()); 
  //gPad->Close();
  A=Form("HPtTPt3%i",Index[0]);
  TH2F *FullHPtTPt3 = (TH2F*)gDirectory->Get(A.c_str()); 
  //gPad->Close();
  A=Form("HPtTPt4%i",Index[0]);
  TH2F *FullHPtTPt4 = (TH2F*)gDirectory->Get(A.c_str()); 
  //gPad->Close();
  A=Form("HPtTPt5%i",Index[0]);
  TH2F *FullHPtTPt5 = (TH2F*)gDirectory->Get(A.c_str()); 
  //gPad->Close();
  A=Form("HMBE0%i",Index[0]);
  TH1F *BKGHMBEFull0 = (TH1F*)gDirectory->Get(A.c_str()); 
  //gPad->Close();
  A=Form("HMBE1%i",Index[0]);
  TH1F *BKGHMBEFull1 = (TH1F*)gDirectory->Get(A.c_str()); 
  //gPad->Close();
  A=Form("HMBE2%i",Index[0]);
  TH1F *BKGHMBEFull2 = (TH1F*)gDirectory->Get(A.c_str()); 
  //gPad->Close();
  A=Form("HMBE3%i",Index[0]);
  TH1F *BKGHMBEFull3 = (TH1F*)gDirectory->Get(A.c_str()); 
  //gPad->Close();
  A=Form("HMBE4%i",Index[0]);
  TH1F *BKGHMBEFull4 = (TH1F*)gDirectory->Get(A.c_str()); 
  //gPad->Close();
  A=Form("HMBE5%i",Index[0]);
  TH1F *BKGHMBEFull5 = (TH1F*)gDirectory->Get(A.c_str()); 
  //gPad->Close();  
  cout << "1 Marker" << endl;
  for (int i=0; i<NOS; i++)
    {
      if (i!=0)
	{
	  HptTpt[0][0]->Add(HptTpt[0][i]);
	  HptTpt[1][0]->Add(HptTpt[1][i]);
	  HptTpt[2][0]->Add(HptTpt[2][i]);
	  HptTpt[3][0]->Add(HptTpt[3][i]);
	  HptTpt[4][0]->Add(HptTpt[4][i]);
	  HptTpt[5][0]->Add(HptTpt[5][i]);
	  //FullHPtTPt0->Add(HptTpt[0][i]);
	  //FullHPtTPt1->Add(HptTpt[1][i]);
	  //FullHPtTPt2->Add(HptTpt[2][i]);
	  //FullHPtTPt3->Add(HptTpt[3][i]);
	  //FullHPtTPt4->Add(HptTpt[4][i]);
	  //FullHPtTPt5->Add(HptTpt[5][i]);
	}
      if (i!=NOS-1)
	{
	  BKGHMBE0->Add(HiggsMassReversedHptTpt[0][i]);
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
	  BKGTprimeMassBE5->Add(TprimeMassReversedHptTpt[5][i]);
	  if (i!=0)
	    {
	      HiggsMassReversedHptTpt[0][0]->Add(HiggsMassReversedHptTpt[0][i]);
	      HiggsMassReversedHptTpt[1][0]->Add(HiggsMassReversedHptTpt[1][i]);
	      HiggsMassReversedHptTpt[2][0]->Add(HiggsMassReversedHptTpt[2][i]);
	      HiggsMassReversedHptTpt[3][0]->Add(HiggsMassReversedHptTpt[3][i]);
	      HiggsMassReversedHptTpt[4][0]->Add(HiggsMassReversedHptTpt[4][i]);
	      HiggsMassReversedHptTpt[5][0]->Add(HiggsMassReversedHptTpt[5][i]);
	      //BKGHMBEFull0->Add(HiggsMassReversedHptTpt[0][i]);
	      //BKGHMBEFull1->Add(HiggsMassReversedHptTpt[1][i]);
	      //BKGHMBEFull2->Add(HiggsMassReversedHptTpt[2][i]);
	      //BKGHMBEFull3->Add(HiggsMassReversedHptTpt[3][i]);
	      //BKGHMBEFull4->Add(HiggsMassReversedHptTpt[4][i]);
	      //BKGHMBEFull5->Add(HiggsMassReversedHptTpt[5][i]);
	      HiggsMassCuts[0][0]->Add(HiggsMassCuts[0][i]);
	      HiggsMassCuts[1][0]->Add(HiggsMassCuts[1][i]);
	      HiggsMassCuts[2][0]->Add(HiggsMassCuts[2][i]);
	      HiggsMassCuts[3][0]->Add(HiggsMassCuts[3][i]);
	      HiggsMassCuts[4][0]->Add(HiggsMassCuts[4][i]);
	      HiggsMassCuts[5][0]->Add(HiggsMassCuts[5][i]);
	    }
	}
    }
  cout << "2 Marker" << endl;

  //ABCD Regions
  double RegionA[6]={0,0,0,0,0,0}; double RegionB[6]={0,0,0,0,0,0}; double RegionC[6]={0,0,0,0,0,0}; double RegionD[6]={0,0,0,0,0,0};
  cout << "Correlation factor on Hpt Tpt on cut 6: " << HptTpt[0][0]->GetCorrelationFactor() << endl;
  RegionA[0]=HptTpt[0][0]->Integral(HptTpt[0][0]->GetXaxis()->FindBin(200.0),HptTpt[0][0]->GetNbinsX(),HptTpt[0][0]->GetYaxis()->FindBin(200.0),HptTpt[0][0]->GetNbinsY());
  RegionB[0]=HptTpt[0][0]->Integral(HptTpt[0][0]->GetXaxis()->FindBin(200.0),HptTpt[0][0]->GetNbinsX(),HptTpt[0][0]->GetYaxis()->FindBin(0.0),HptTpt[0][0]->GetYaxis()->FindBin(200.0));
  RegionC[0]=HptTpt[0][0]->Integral(HptTpt[0][0]->GetXaxis()->FindBin(0.0),HptTpt[0][0]->GetXaxis()->FindBin(200.0),HptTpt[0][0]->GetYaxis()->FindBin(200.0),HptTpt[0][0]->GetNbinsY());
  RegionD[0]=HptTpt[0][0]->Integral(HptTpt[0][0]->GetXaxis()->FindBin(0.0),HptTpt[0][0]->GetXaxis()->FindBin(200.0),HptTpt[0][0]->GetYaxis()->FindBin(0.0),HptTpt[0][0]->GetYaxis()->FindBin(200.0));
  cout << "Correlation factor on Hpt Tpt on cut 7: " << HptTpt[1][0]->GetCorrelationFactor() << endl;
  RegionA[1]=HptTpt[1][0]->Integral(HptTpt[1][0]->GetXaxis()->FindBin(200.0),HptTpt[1][0]->GetNbinsX(),HptTpt[1][0]->GetYaxis()->FindBin(200.0),HptTpt[1][0]->GetNbinsY());
  RegionB[1]=HptTpt[1][0]->Integral(HptTpt[1][0]->GetXaxis()->FindBin(200.0),HptTpt[1][0]->GetNbinsX(),HptTpt[1][0]->GetYaxis()->FindBin(0.0),HptTpt[1][0]->GetYaxis()->FindBin(200.0));
  RegionC[1]=HptTpt[1][0]->Integral(HptTpt[1][0]->GetXaxis()->FindBin(0.0),HptTpt[1][0]->GetXaxis()->FindBin(200.0),HptTpt[1][0]->GetYaxis()->FindBin(200.0),HptTpt[1][0]->GetNbinsY());
  RegionD[1]=HptTpt[1][0]->Integral(HptTpt[1][0]->GetXaxis()->FindBin(0.0),HptTpt[1][0]->GetXaxis()->FindBin(200.0),HptTpt[1][0]->GetYaxis()->FindBin(0.0),HptTpt[1][0]->GetYaxis()->FindBin(200.0));
  cout << "Correlation factor on Hpt Tpt on cut 8: " << HptTpt[2][0]->GetCorrelationFactor() << endl;
  RegionA[2]=HptTpt[2][0]->Integral(HptTpt[2][0]->GetXaxis()->FindBin(200.0),HptTpt[2][0]->GetNbinsX(),HptTpt[2][0]->GetYaxis()->FindBin(200.0),HptTpt[2][0]->GetNbinsY());
  RegionB[2]=HptTpt[2][0]->Integral(HptTpt[2][0]->GetXaxis()->FindBin(200.0),HptTpt[2][0]->GetNbinsX(),HptTpt[2][0]->GetYaxis()->FindBin(0.0),HptTpt[2][0]->GetYaxis()->FindBin(200.0));
  RegionC[2]=HptTpt[2][0]->Integral(HptTpt[2][0]->GetXaxis()->FindBin(0.0),HptTpt[2][0]->GetXaxis()->FindBin(200.0),HptTpt[2][0]->GetYaxis()->FindBin(200.0),HptTpt[2][0]->GetNbinsY());
  RegionD[2]=HptTpt[2][0]->Integral(HptTpt[2][0]->GetXaxis()->FindBin(0.0),HptTpt[2][0]->GetXaxis()->FindBin(200.0),HptTpt[2][0]->GetYaxis()->FindBin(0.0),HptTpt[2][0]->GetYaxis()->FindBin(200.0));
  cout << "Correlation factor on Hpt Tpt on cut 9: " << HptTpt[3][0]->GetCorrelationFactor() << endl;
  RegionA[3]=HptTpt[3][0]->Integral(HptTpt[3][0]->GetXaxis()->FindBin(200.0),HptTpt[3][0]->GetNbinsX(),HptTpt[3][0]->GetYaxis()->FindBin(200.0),HptTpt[3][0]->GetNbinsY());
  RegionB[3]=HptTpt[3][0]->Integral(HptTpt[3][0]->GetXaxis()->FindBin(200.0),HptTpt[3][0]->GetNbinsX(),HptTpt[3][0]->GetYaxis()->FindBin(0.0),HptTpt[3][0]->GetYaxis()->FindBin(200.0));
  RegionC[3]=HptTpt[3][0]->Integral(HptTpt[3][0]->GetXaxis()->FindBin(0.0),HptTpt[3][0]->GetXaxis()->FindBin(200.0),HptTpt[3][0]->GetYaxis()->FindBin(200.0),HptTpt[3][0]->GetNbinsY());
  RegionD[3]=HptTpt[3][0]->Integral(HptTpt[3][0]->GetXaxis()->FindBin(0.0),HptTpt[3][0]->GetXaxis()->FindBin(200.0),HptTpt[3][0]->GetYaxis()->FindBin(0.0),HptTpt[3][0]->GetYaxis()->FindBin(200.0));
  cout << "Correlation factor on Hpt Tpt on cut 10: " << HptTpt[4][0]->GetCorrelationFactor() << endl;
  RegionA[4]=HptTpt[4][0]->Integral(HptTpt[4][0]->GetXaxis()->FindBin(200.0),HptTpt[4][0]->GetNbinsX(),HptTpt[4][0]->GetYaxis()->FindBin(200.0),HptTpt[4][0]->GetNbinsY());
  RegionB[4]=HptTpt[4][0]->Integral(HptTpt[4][0]->GetXaxis()->FindBin(200.0),HptTpt[4][0]->GetNbinsX(),HptTpt[4][0]->GetYaxis()->FindBin(0.0),HptTpt[4][0]->GetYaxis()->FindBin(200.0));
  RegionC[4]=HptTpt[4][0]->Integral(HptTpt[4][0]->GetXaxis()->FindBin(0.0),HptTpt[4][0]->GetXaxis()->FindBin(200.0),HptTpt[4][0]->GetYaxis()->FindBin(200.0),HptTpt[4][0]->GetNbinsY());
  RegionD[4]=HptTpt[4][0]->Integral(HptTpt[4][0]->GetXaxis()->FindBin(0.0),HptTpt[4][0]->GetXaxis()->FindBin(200.0),HptTpt[4][0]->GetYaxis()->FindBin(0.0),HptTpt[4][0]->GetYaxis()->FindBin(200.0));
  cout << "Correlation factor on Hpt Tpt on cut 11: " << HptTpt[5][0]->GetCorrelationFactor() << endl;
  RegionA[5]=HptTpt[5][0]->Integral(HptTpt[5][0]->GetXaxis()->FindBin(200.0),HptTpt[5][0]->GetNbinsX(),HptTpt[5][0]->GetYaxis()->FindBin(200.0),HptTpt[5][0]->GetNbinsY());
  RegionB[5]=HptTpt[5][0]->Integral(HptTpt[5][0]->GetXaxis()->FindBin(200.0),HptTpt[5][0]->GetNbinsX(),HptTpt[5][0]->GetYaxis()->FindBin(0.0),HptTpt[5][0]->GetYaxis()->FindBin(200.0));
  RegionC[5]=HptTpt[5][0]->Integral(HptTpt[5][0]->GetXaxis()->FindBin(0.0),HptTpt[5][0]->GetXaxis()->FindBin(200.0),HptTpt[5][0]->GetYaxis()->FindBin(200.0),HptTpt[5][0]->GetNbinsY());
  RegionD[5]=HptTpt[5][0]->Integral(HptTpt[5][0]->GetXaxis()->FindBin(0.0),HptTpt[5][0]->GetXaxis()->FindBin(200.0),HptTpt[5][0]->GetYaxis()->FindBin(0.0),HptTpt[5][0]->GetYaxis()->FindBin(200.0));

  for (int j=0; j<6; j++) cout << "Overall estimation N bkgs up to " << j << ": " <<  RegionB[j]*RegionC[j]/RegionD[j] << ", and Region A: " << RegionA[j] << ", and sum up of all regions: " << RegionA[j]+RegionB[j]+RegionC[j]+RegionD[j] << endl; 
  //for (int j=0; j<6; j++) cout << "Overall real N bkgs up to " << j << ": " <<  HiggsMassCuts[j][0]->Integral() << endl;
  
  /*
  //TH1F *FakeDY; FakeDY->SetFillColor(kBlue);

  //PostScript Plotting
  cout << "MARKER" << endl;
  TLegend* BKGandSignallegend = new TLegend(0.75,0.65,0.90,0.9);
  TLegend* BKGlegend = new TLegend(0.75,0.65,0.90,0.9);
  if (NOS==5) 
    {
      BKGandSignallegend->AddEntry(HiggsMassReversedHptTpt[0][3], "QCD", "f");
      //BKGandSignallegend->AddEntry(FakeDY, "Zjets", "f");      
      BKGandSignallegend->AddEntry(HiggsMassReversedHptTpt[0][2], "TTbar", "f");
      BKGandSignallegend->AddEntry(HiggsMassReversedHptTpt[0][1], "SingleT", "f");
      BKGandSignallegend->AddEntry(HiggsMassReversedHptTpt[0][0], "Diboson", "f");
      BKGandSignallegend->AddEntry(HiggsMassReversedHptTpt[0][4], "Signal", "f");
      BKGlegend->AddEntry(HiggsMassReversedHptTpt[0][3], "QCD", "f");
      //BKGlegend->AddEntry(FakeDY, "Zjets", "f");      
      BKGlegend->AddEntry(HiggsMassReversedHptTpt[0][2], "TTbar", "f");
      BKGlegend->AddEntry(HiggsMassReversedHptTpt[0][1], "SingleT", "f");
      BKGlegend->AddEntry(HiggsMassReversedHptTpt[0][0], "Diboson", "f");
    }
  else 
    {
      BKGandSignallegend->AddEntry(HiggsMassReversedHptTpt[0][4], "QCD", "f");
      BKGandSignallegend->AddEntry(HiggsMassReversedHptTpt[0][3], "Zjets", "f");
      BKGandSignallegend->AddEntry(HiggsMassReversedHptTpt[0][2], "TTbar", "f");
      BKGandSignallegend->AddEntry(HiggsMassReversedHptTpt[0][1], "SingleT", "f");
      BKGandSignallegend->AddEntry(HiggsMassReversedHptTpt[0][0], "Diboson", "f");
      BKGandSignallegend->AddEntry(HiggsMassReversedHptTpt[0][5], "Signal", "f");
      BKGlegend->AddEntry(HiggsMassReversedHptTpt[0][4], "QCD", "f");
      BKGlegend->AddEntry(HiggsMassReversedHptTpt[0][3], "Zjets", "f");
      BKGlegend->AddEntry(HiggsMassReversedHptTpt[0][2], "TTbar", "f");
      BKGlegend->AddEntry(HiggsMassReversedHptTpt[0][1], "SingleT", "f");
      BKGlegend->AddEntry(HiggsMassReversedHptTpt[0][0], "Diboson", "f");
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
  MyPlot->Update();
  //////////////
  //////////////
  //////////////
  MyPlot->Clear();
  MyPlot->cd(1);
  ps->NewPage();
  FullHPtTPt0->Draw("colz");
  gPad->Update();
  MyPlot->Update();
  //////////////
  //////////////
  //////////////
  MyPlot->Clear();
  MyPlot->cd(1);
  ps->NewPage();
  FullHPtTPt1->Draw("colz");
  gPad->Update();
  MyPlot->Update();
  //////////////
  //////////////
  //////////////
  MyPlot->Clear();
  MyPlot->cd(1);
  ps->NewPage();
  FullHPtTPt2->Draw("colz");
  gPad->Update();
  MyPlot->Update();
  //////////////
  //////////////
  //////////////
  MyPlot->Clear();
  MyPlot->cd(1);
  ps->NewPage();
  FullHPtTPt3->Draw("colz");
  gPad->Update();
  MyPlot->Update();
  //////////////
  //////////////
  //////////////
  MyPlot->Clear();
  MyPlot->cd(1);
  ps->NewPage();
  FullHPtTPt4->Draw("colz");
  gPad->Update();
  MyPlot->Update();
  //////////////
  //////////////
  //////////////
  MyPlot->Clear();
  MyPlot->cd(1);
  ps->NewPage();
  FullHPtTPt5->Draw("colz");
  gPad->Update();
  MyPlot->Update();
  //////////////  
  ps->Close();*/

  exit(0);
}
