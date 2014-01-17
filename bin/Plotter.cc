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
#include <TFile.h>

#include "TMath.h"
#include "TF1.h"

using namespace std;

const TString StorageDirPrefix="file:/home/cms/jruizalv/work/CMSSW_5_3_9_patch2/src/Extractors/PatExtractor/bin/";
const int NOS=5;
//const TString Samples[NOS] = {"Diboson.root", "SingleTop.root", "TTJets.root", "Zjets.root", "Wjets.root", "QCD.root", "Signal.root"};
const TString Samples[NOS] = {"Diboson.root", "SingleTop.root", "TTJets.root", "QCD.root", "Signal.root"};

void Plotter()
{
  TFile *CurrentFile[NOS];
  TH1F *TprimeHistos[NOS];
  
  for (int i=0; i<NOS; i++)
    {
      CurrentFile[i] = new TFile(StorageDirPrefix + Samples[i], "READ");
      if ( CurrentFile[i]->IsOpen() ) printf( Samples[i] + " File opened successfully\n");
      if (i==0) {TprimeHistos[0]->SetDefaultSumw2(); TprimeHistos[0]= (TH1F*)gDirectory->Get("TprimeMass0");}
      else if (i==1) {TprimeHistos[1]->SetDefaultSumw2(); TprimeHistos[1]= (TH1F*)gDirectory->Get("TprimeMass3");}
      else if (i==2) {TprimeHistos[2]->SetDefaultSumw2(); TprimeHistos[2]= (TH1F*)gDirectory->Get("TprimeMass9");}
      else if (i==3) {TprimeHistos[3]->SetDefaultSumw2(); TprimeHistos[3]= (TH1F*)gDirectory->Get("TprimeMass10");}
      else if (i==4) {TprimeHistos[4]->SetDefaultSumw2(); TprimeHistos[4]= (TH1F*)gDirectory->Get("TprimeMass15");}
      //else if (i==3) {TprimeHistos[3]->SetDefaultSumw2(); TprimeHistos[3]= (TH1F*)gDirectory->Get("TprimeMass9");}
      //else if (i==4) {TprimeHistos[4]->SetDefaultSumw2(); TprimeHistos[4]= (TH1F*)gDirectory->Get("TprimeMass10");}
      //else if (i==5) {TprimeHistos[5]->SetDefaultSumw2(); TprimeHistos[5]= (TH1F*)gDirectory->Get("TprimeMass11");}
      //else if (i==6) {TprimeHistos[6]->SetDefaultSumw2(); TprimeHistos[6]= (TH1F*)gDirectory->Get("TprimeMass12");}
    }
  
  THStack *BKGandSignal = new THStack("BKGandSignal", "BKG and signal for 5 jets; M_{5j} GeV; Events");
  for (int i=0; i<NOS; i++)
    {
      
      //float BinSize=TprimeHistos[i]->GetBinWidth(1);
      //int RebinSet=TMath::Nint(10./BinSize);
      //int NumberOfBinsToSet=TprimeHistos[2]->GetNbinsX();
      //int RebinSet=TMath::Nint(10./NumberOfBinsToSet);
      //cout << RebinSet << " " << BinSize << endl;
      //TprimeHistos[i]->Rebin(RebinSet);
      //TprimeHistos[i]->SetAxisRange(50,1600,"X");
      //TprimeHistos[i]->SetBit(TH1::kCanRebin);
      //if (RebinSet<3) 
      BKGandSignal->Add(TprimeHistos[i]);
    }
  
  TCanvas *MyPlot = new TCanvas("MyPlot","Single t prime to top Higgs with backgrounds",600,600);
  MyPlot->Clear();
  MyPlot->cd(1);
  BKGandSignal->Draw("hist");
  //BKGandSignallegend->Draw();
  TprimeHistos[4]->Draw("histsame");
  //TprimeHistos[6]->Draw("histsame");
  gPad->SetLogy();
  //BKGandSignal->SetMinimum(FiveJetsMass[0]->GetMinimum()+1.);  
  gPad->RedrawAxis("");
  gPad->Update();
  MyPlot->Update();

  //Calculation of Estimator
  double bkg[4]={0,0,0,0};
  double sig=0;

  sig=TprimeHistos[4]->Integral(TprimeHistos[4]->GetXaxis()->FindBin(710),TprimeHistos[4]->GetXaxis()->FindBin(750));
  for (int k=0; k<4; k++) bkg[k]=TprimeHistos[k]->Integral(TprimeHistos[k]->GetXaxis()->FindBin(710),TprimeHistos[k]->GetXaxis()->FindBin(750));
  
  cout << "Estimator under the peak!" << endl;
  cout << "Number of signal events: " << sig << ", Number of BKG events: " << bkg[0]+bkg[1]+bkg[2]+bkg[3] << ", S/sqrt(S+B)=" << sig/sqrt(sig+bkg[0]+bkg[1]+bkg[2]+bkg[3]) << endl;

  double bkgFR[4]={0,0,0,0};
  double sigFR=0;

  sigFR=TprimeHistos[4]->Integral();
  for (int k=0; k<4; k++) bkgFR[k]=TprimeHistos[k]->Integral();
  
  cout << "Estimator in full range!" << endl;
  cout << "Number of signal events: " << sigFR << ", Number of BKG events: " << bkgFR[0]+bkgFR[1]+bkgFR[2]+bkgFR[3] << ", S/sqrt(S+B)=" << sigFR/sqrt(sigFR+bkgFR[0]+bkgFR[1]+bkgFR[2]+bkgFR[3]) << endl;

}
