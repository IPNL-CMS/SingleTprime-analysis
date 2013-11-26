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
#include "FilesNames.C"

using namespace std;

// Global Parameters
// Location of root files

const TString StorageDirPrefix[NumberOfProcesses]={"file:/gridgroup/cms/jruizalv/Extracted_MC/ZZ/", "file:/gridgroup/cms/jruizalv/Extracted_MC/T-s/", "file:/gridgroup/cms/jruizalv/Extracted_MC/T-t/","file:/gridgroup/cms/jruizalv/Extracted_MC/T-tw/", "file:/gridgroup/cms/jruizalv/Extracted_MC/Tbar-s/", "file:/gridgroup/cms/jruizalv/Extracted_MC/Tbar-t/", "file:/gridgroup/cms/jruizalv/Extracted_MC/Tbar-tw/", "file:/gridgroup/cms/jruizalv/Extracted_MC/TTJets/", "file:/gridgroup/cms/jruizalv/Extracted_MC/DYToCC/", "file:/gridgroup/cms/jruizalv/Extracted_MC/DYToBB/", "file:/gridgroup/cms/jruizalv/Extracted_MC/Wjets_VBF/", "file:/gridgroup/cms/jruizalv/Extracted_MC/QCD_300_CSVM/", "file:/home/cms/jruizalv/work/CMSSW_5_3_9_patch2/src/Extractors/PatExtractor/test/"};

//QCD_PT_300_470/

const int NumberOfHistos=14;
const TString Histos[NumberOfHistos] = {"Cut_0", "Cut_1", "Cut_2", "Cut_3", "Cut_4", "Cut_5", "Cut_6", "Cut_7", "Cut_8", "Cut_9", "Cut_10", "Cut_11", "Cut_12", "Cut_13"};

//Computing weights (Everything in pb)
const float XS[NumberOfProcesses]={8.258, 2.82, 47, 10.7, 1.57, 25, 10.7, 234.0, 3060.099, 3840.86, 1759.549, 1205.0, 0.15}; // 113.8791, 26.9921, 3.550036, 54.838, 33.72, 1205.0}; //Signal,T_s,T_t,T_tW,Tbar_s,Tbar_t,Tbar_tW,TTbar,QCD_300,DYToBB,DYToCC,ZZ,    QCD_470,QCD_600,QCD_800,WW,WZ,Wjets_VBF
/*float T_s_xs=2.82;
float T_t_xs=47;
float T_tW_xs=10.7;
float Tbar_s_xs=1.57;
float Tbar_t_xs=25;
float Tbar_tW_xs=10.7;
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

float Lumi=20000;

void GeneralAnalyzer()
{
  int ParcialTestMax=NumberOfProcesses;
  int ParcialTestMin=0;
  TH1F *FiveJetsMass[NumberOfProcesses];
  int EntriePerSample[NumberOfProcesses];
  bool SurvivalMarker[NumberOfProcesses];
  //for (int k=0; k<NumberOfProcesses; k++)
  for (int k=ParcialTestMin; k<ParcialTestMax; k++)
    {
      TChain CutsChain("cuts");
      TChain AnalysisChain("stp");
      TH1F *ALLCuts[NumberOfHistos];
      for (int i=0; i<NumberOfSamples[k]; i++)
	{
	  if (k==NumberOfProcesses-1)
	    {
	      CutsChain.Add(StorageDirPrefix[k] + SamplesSignal[i]);
	      AnalysisChain.Add(StorageDirPrefix[k] + SamplesSignal[i]);
	    }
	  else if (k==1)
	    {
	      CutsChain.Add(StorageDirPrefix[k] + SamplesT_s[i]);
	      AnalysisChain.Add(StorageDirPrefix[k] + SamplesT_s[i]);
	    }
	  else if (k==2)
	    {
	      CutsChain.Add(StorageDirPrefix[k] + SamplesT_t[i]);
	      AnalysisChain.Add(StorageDirPrefix[k] + SamplesT_t[i]);
	    }
	  else if (k==3)
	    {
	      CutsChain.Add(StorageDirPrefix[k] + SamplesT_tw[i]);
	      AnalysisChain.Add(StorageDirPrefix[k] + SamplesT_tw[i]);
	    }
	  else if (k==4)
	    {
	      CutsChain.Add(StorageDirPrefix[k] + SamplesTbar_s[i]);
	      AnalysisChain.Add(StorageDirPrefix[k] + SamplesTbar_s[i]);
	    }
	  else if (k==5)
	    {
	      CutsChain.Add(StorageDirPrefix[k] + SamplesTbar_t[i]);
	      AnalysisChain.Add(StorageDirPrefix[k] + SamplesTbar_t[i]);
	    }
	  else if (k==6)
	    {
	      CutsChain.Add(StorageDirPrefix[k] + SamplesTbar_tw[i]);
	      AnalysisChain.Add(StorageDirPrefix[k] + SamplesTbar_tw[i]);
	    }
	  else if (k==7)
	    {
	      CutsChain.Add(StorageDirPrefix[k] + SamplesTTbar[i]);
	      AnalysisChain.Add(StorageDirPrefix[k] + SamplesTTbar[i]);
	    }
	  else if (k==11)
	    {
	      CutsChain.Add(StorageDirPrefix[k] + SamplesQCD300[i]);
	      AnalysisChain.Add(StorageDirPrefix[k] + SamplesQCD300[i]);
	    }
	  else if (k==8)
	    {
	      CutsChain.Add(StorageDirPrefix[k] + SamplesDYToCC[i]);
	      AnalysisChain.Add(StorageDirPrefix[k] + SamplesDYToCC[i]);
	    }
	  else if (k==9)
	    {
	      CutsChain.Add(StorageDirPrefix[k] + SamplesDYToBB[i]);
	      AnalysisChain.Add(StorageDirPrefix[k] + SamplesDYToBB[i]);
	    }
	  else if (k==0)
	    {
	      CutsChain.Add(StorageDirPrefix[k] + SamplesZZ[i]);
	      AnalysisChain.Add(StorageDirPrefix[k] + SamplesZZ[i]);
	    }
	  else if (k==10)
	    {
	      CutsChain.Add(StorageDirPrefix[k] + SamplesWjets_VBF[i]);
	      AnalysisChain.Add(StorageDirPrefix[k] + SamplesWjets_VBF[i]);
	    }
	}
      
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
      
      EntriePerSample[k]=EntriesPerCut[0];
      cout << PassedPerCut[NumberOfProcesses] << endl;
      gPad->Close();
      if (PassedPerCut[NumberOfProcesses]!=0)
	{
	  string A1 = Form("Reconstructed_Tprime.M() >> TprimeMass%i(155,50,1600)",k);
	  string A2 = Form("TprimeMass%i",k);
	  AnalysisChain.Draw(A1.c_str());
	  FiveJetsMass[k] = (TH1F*)gDirectory->Get(A2.c_str());
	  gPad->Close();
	  SurvivalMarker[k]=true;
	}

    }

  //Making Stack

  THStack *BKGandSignal = new THStack("BKGandSignal", "BKG and signal for 5 jets; M_{5j} GeV; Events");
  TLegend* BKGandSignallegend = new TLegend(0.75,0.55,0.90,0.9);
  //TH1F *Signal; TH1F *Singletop; TH1F *TTbar;
  for (int k=ParcialTestMin; k<ParcialTestMax; k++)
    {
      if (!SurvivalMarker[k]) continue;
      FiveJetsMass[k]->Scale(Lumi*XS[k]/EntriePerSample[k]);
      //Settings for signal
      if (k==NumberOfProcesses-1)
	{
	  FiveJetsMass[k]->SetFillColor(kSpring);
	  FiveJetsMass[k]->SetFillStyle(3444);
	  FiveJetsMass[k]->SetLineWidth(3);
	  TFile f("Signal.root", "RECREATE");
	  FiveJetsMass[k]->Write();
	  //TprimeMass[0]->Write();
	  //Signal->Add(FiveJetsMass[k]);
	}
      //Settings for Single Top
      else if (k<=6 && k>=1) 
	{
	  FiveJetsMass[k]->SetFillColor(kBlack);
	  FiveJetsMass[k]->SetFillStyle(3305);
	  TH1F *SingleTopTprimeMass=(TH1F*)gDirectory->Get("TprimeMass1");
	  if (k!=1) SingleTopTprimeMass->Add(FiveJetsMass[k]);
	  if (k==6)
	    {
	      TFile f("SingleTop.root", "RECREATE");
	      SingleTopTprimeMass->Write();
	    }
	  //Singletop->Add(FiveJetsMass[k]);
	}
      //Settings for TTbar
      else if (k==7) 
	{
	  FiveJetsMass[k]->SetFillColor(kRed);
	  FiveJetsMass[k]->SetFillStyle(3345);
	  TFile f("TTJets.root", "RECREATE");
	  FiveJetsMass[k]->Write();
	  //TTbar->Add(FiveJetsMass[k]);
	}
      //Settings for QCD
      else if (k==11) {FiveJetsMass[k]->SetFillColor(kViolet); TFile f("QCD.root", "RECREATE"); FiveJetsMass[k]->Write();}
      //Settings for DiBoson
      else if (k==0) {FiveJetsMass[k]->SetFillColor(kWhite); TFile f("Diboson.root", "RECREATE"); FiveJetsMass[k]->Write();}
      //Setting for Zjets
      else if (k==8 || k==9) 
	{
	  FiveJetsMass[k]->SetFillColor(kBlue); 
	  TH1F *ZjetsTprimeMass=(TH1F*)gDirectory->Get("TprimeMass9");
	  //if (k!=8) ZjetsTprimeMass->Add(FiveJetsMass[k]);
	  TFile f("Zjets.root", "RECREATE"); 
	  ZjetsTprimeMass->Write();
	}
      //Setting for Wjets
      if (k==10) 
	{
	  FiveJetsMass[k]->SetFillColor(kYellow+1);
	  FiveJetsMass[k]->SetFillStyle(3354);
	  TFile f("Wjets.root", "RECREATE");
	  FiveJetsMass[k]->Write();
	}
      BKGandSignal->Add(FiveJetsMass[k]);
      BKGandSignallegend->AddEntry(FiveJetsMass[k]);
    }

  //BKGandSignal->Add(Singletop);
  //BKGandSignallegend->AddEntry(Singletop);
  //BKGandSignal->Add(Signal);
  //BKGandSignallegend->AddEntry(Signal);

  /*TCanvas *MyPlot = new TCanvas("MyPlot","Single t prime to top Higgs with backgrounds",600,600);
  MyPlot->Clear();
  MyPlot->cd(1);
  BKGandSignal->Draw("hist");
  //BKGandSignallegend->Draw();
  //FiveJetsMass[0]->Draw("histsame");
  gPad->SetLogy();
  //BKGandSignal->SetMinimum(FiveJetsMass[0]->GetMinimum()+1.);  
  gPad->RedrawAxis("");
  gPad->Update();
  MyPlot->Update();*/
  
}
