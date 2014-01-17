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
#include "FilesNames_DRTH_Tighten.C"

using namespace std;

// Global Parameters
// Location of root files

const TString StorageDirPrefix[NumberOfProcesses]={"file:/gridgroup/cms/jruizalv/Extracted_MC/DRTH_PTNormMass_tighten/ZZ/", "file:/gridgroup/cms/jruizalv/Extracted_MC/DRTH_PTNormMass_tighten/WZ/", "file:/gridgroup/cms/jruizalv/Extracted_MC/DRTH_PTNormMass_tighten/WW/", "file:/gridgroup/cms/jruizalv/Extracted_MC/DRTH_PTNormMass_tighten/T-s/", "file:/gridgroup/cms/jruizalv/Extracted_MC/DRTH_PTNormMass_tighten/T-t/","file:/gridgroup/cms/jruizalv/Extracted_MC/DRTH_PTNormMass_tighten/T-tw/", "file:/gridgroup/cms/jruizalv/Extracted_MC/DRTH_PTNormMass_tighten/Tbar-s/", "file:/gridgroup/cms/jruizalv/Extracted_MC/DRTH_PTNormMass_tighten/Tbar-t/", "file:/gridgroup/cms/jruizalv/Extracted_MC/DRTH_PTNormMass_tighten/Tbar-tw/", "file:/gridgroup/cms/jruizalv/Extracted_MC/DRTH_PTNormMass_tighten/TTJets/", "file:/gridgroup/cms/jruizalv/Extracted_MC/DRTH_PTNormMass_tighten/QCD_PT_170_300/", "file:/gridgroup/cms/jruizalv/Extracted_MC/DRTH_PTNormMass_tighten/QCD_PT_300_470/", "file:/gridgroup/cms/jruizalv/Extracted_MC/DRTH_PTNormMass_tighten/QCD_PT_470_600/", "file:/gridgroup/cms/jruizalv/Extracted_MC/DRTH_PTNormMass_tighten/QCD_PT_600_800/", "file:/gridgroup/cms/jruizalv/Extracted_MC/DRTH_PTNormMass_tighten/QCD_PT_800_1000/", "file:/home/cms/jruizalv/work/CMSSW_5_3_9_patch2/src/Extractors/PatExtractor/test/"};

const int NumberOfHistos=20;
const TString Histos[NumberOfHistos] = {"Cut_0", "Cut_1", "Cut_2", "Cut_3", "Cut_4", "Cut_5", "Cut_6", "Cut_7", "Cut_8", "Cut_9", "Cut_10", "Cut_11", "Cut_12", "Cut_13", "Cut_14", "Cut_15", "Cut_16", "Cut_17", "Cut_18", "Cut_19"};

//Computing weights (Everything in pb)
const float XS[NumberOfProcesses]={8.258, 33.72, 54.838, 2.82, 47, 10.7, 1.57, 25, 10.7, 234.0, 34138.15, 1759.549, 113.8791, 26.9921, 3.550036, 0.15}; 
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

float Lumi=20000;

void GeneralAnalyzer()
{
  int ParcialTestMax=NumberOfProcesses;
  int ParcialTestMin=0;
  TH1F *FiveJetsMass[NumberOfProcesses];
  TH1F *LeadingJetPT[NumberOfProcesses];
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
  int EntriePerSample[NumberOfProcesses];
  bool SurvivalMarker[NumberOfProcesses];

  for (int k=ParcialTestMin; k<ParcialTestMax; k++)
    {
      TChain CutsChain("cuts");
      TChain AnalysisChain("stp");
      TH1F *ALLCuts[NumberOfHistos];
      for (int i=0; i<NumberOfSamples[k]; i++)
	{
	  if (k==0)
	    {
	      CutsChain.Add(StorageDirPrefix[k] + SamplesZZ[i]);
	      AnalysisChain.Add(StorageDirPrefix[k] + SamplesZZ[i]);
	    }
	  else if (k==1)
	    {
	      CutsChain.Add(StorageDirPrefix[k] + SamplesWZ[i]);
	      AnalysisChain.Add(StorageDirPrefix[k] + SamplesWZ[i]);
	    }
	  else if (k==2)
	    {
	      CutsChain.Add(StorageDirPrefix[k] + SamplesWW[i]);
	      AnalysisChain.Add(StorageDirPrefix[k] + SamplesWW[i]);
	    }
	  else if (k==3)
	    {
	      CutsChain.Add(StorageDirPrefix[k] + SamplesT_s[i]);
	      AnalysisChain.Add(StorageDirPrefix[k] + SamplesT_s[i]);
	    }
	  else if (k==4)
	    {
	      CutsChain.Add(StorageDirPrefix[k] + SamplesT_t[i]);
	      AnalysisChain.Add(StorageDirPrefix[k] + SamplesT_t[i]);
	    }
	  else if (k==5)
	    {
	      CutsChain.Add(StorageDirPrefix[k] + SamplesT_tw[i]);
	      AnalysisChain.Add(StorageDirPrefix[k] + SamplesT_tw[i]);
	    }
	  else if (k==6)
	    {
	      CutsChain.Add(StorageDirPrefix[k] + SamplesTbar_s[i]);
	      AnalysisChain.Add(StorageDirPrefix[k] + SamplesTbar_s[i]);
	    }
	  else if (k==7)
	    {
	      CutsChain.Add(StorageDirPrefix[k] + SamplesTbar_t[i]);
	      AnalysisChain.Add(StorageDirPrefix[k] + SamplesTbar_t[i]);
	    }
	  else if (k==8)
	    {
	      CutsChain.Add(StorageDirPrefix[k] + SamplesTbar_tw[i]);
	      AnalysisChain.Add(StorageDirPrefix[k] + SamplesTbar_tw[i]);
	    }
	  else if (k==9)
	    {
	      CutsChain.Add(StorageDirPrefix[k] + SamplesTTbar[i]);
	      AnalysisChain.Add(StorageDirPrefix[k] + SamplesTTbar[i]);
	    }
	  else if (k==10)
	    {
	      CutsChain.Add(StorageDirPrefix[k] + SamplesQCD170[i]);
	      AnalysisChain.Add(StorageDirPrefix[k] + SamplesQCD170[i]);
	    }
	  else if (k==11)
	    {
	      CutsChain.Add(StorageDirPrefix[k] + SamplesQCD300[i]);
	      AnalysisChain.Add(StorageDirPrefix[k] + SamplesQCD300[i]);
	    }
	  else if (k==12)
	    {
	      CutsChain.Add(StorageDirPrefix[k] + SamplesQCD470[i]);
	      AnalysisChain.Add(StorageDirPrefix[k] + SamplesQCD470[i]);
	    }
	  else if (k==13)
	    {
	      CutsChain.Add(StorageDirPrefix[k] + SamplesQCD600[i]);
	      AnalysisChain.Add(StorageDirPrefix[k] + SamplesQCD600[i]);
	    }
	  else if (k==14)
	    {
	      CutsChain.Add(StorageDirPrefix[k] + SamplesQCD800[i]);
	      AnalysisChain.Add(StorageDirPrefix[k] + SamplesQCD800[i]);
	    }
	  else if (k==NumberOfProcesses-1)
	    {
	      CutsChain.Add(StorageDirPrefix[k] + SamplesSignal[i]);
	      AnalysisChain.Add(StorageDirPrefix[k] + SamplesSignal[i]);
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
      cout << PassedPerCut[NumberOfHistos-1] << endl;
      gPad->Close();
      if (PassedPerCut[NumberOfHistos-1]!=0)
	{
	  string A1 = Form("Reconstructed_Tprime.M() >> TprimeMass%i(60,400,1600)",k);
	  string A2 = Form("TprimeMass%i",k);
	  AnalysisChain.Draw(A1.c_str());
	  FiveJetsMass[k] = (TH1F*)gDirectory->Get(A2.c_str());
	  gPad->Close();
	  SurvivalMarker[k]=true;
	}
      /////////////////////
      //Saving Histograms//
      /////////////////////
      string LJPT1 = Form("jet1_pt >> jet1_pt%i(100,10,500)",k);
      string LJPT2 = Form("jet1_pt%i",k);
      AnalysisChain.Draw(LJPT1.c_str());
      LeadingJetPT[k] = (TH1F*)gDirectory->Get(LJPT2.c_str());
      gPad->Close();
      string THT1 = Form("THT >> THT%i(100,600,1400)",k);
      string THT2 = Form("THT%i",k);
      AnalysisChain.Draw(THT1.c_str());
      THT[k] = (TH1F*)gDirectory->Get(THT2.c_str());
      gPad->Close();

    }

  //Making Stack

  THStack *BKGandSignal = new THStack("BKGandSignal", "BKG and signal for 5 jets; M_{5j} GeV; Events");
  TLegend* BKGandSignallegend = new TLegend(0.75,0.55,0.90,0.9);
  //TprimeMass histo recollection
  TH1F *SingleTopTprimeMass=(TH1F*)gDirectory->Get("TprimeMass3");
  TH1F *QCDTprimeMass=(TH1F*)gDirectory->Get("TprimeMass10");  
  TH1F *DibosonTprimeMass=(TH1F*)gDirectory->Get("TprimeMass0");
  //LeadingJetPT histo recollection
  TH1F *SingleTopTprimeMassLJPT=(TH1F*)gDirectory->Get("jet1_pt3");
  TH1F *QCDTprimeMassLJPT=(TH1F*)gDirectory->Get("jet1_pt10");  
  TH1F *DibosonTprimeMassLJPT=(TH1F*)gDirectory->Get("jet1_pt0");
  //THT histo recollection
  TH1F *SingleTopTprimeMassTHT=(TH1F*)gDirectory->Get("THT3");
  TH1F *QCDTprimeMassTHT=(TH1F*)gDirectory->Get("THT10");  
  TH1F *DibosonTprimeMassTHT=(TH1F*)gDirectory->Get("THT0");
  for (int k=ParcialTestMin; k<ParcialTestMax; k++)
    {
      if (!SurvivalMarker[k]) continue;
      FiveJetsMass[k]->Scale(Lumi*XS[k]/EntriePerSample[k]);
      //Settings for signal
      if (k==NumberOfProcesses-1)
	{
	  //TprimeMass
	  FiveJetsMass[k]->SetFillColor(kSpring);
	  FiveJetsMass[k]->SetFillStyle(3444);
	  FiveJetsMass[k]->SetLineWidth(3);
	  TFile f("Signal.root", "RECREATE");
	  FiveJetsMass[k]->Write();
	  //LJPT
	  LeadingJetPT[k]->SetFillColor(kSpring);
	  LeadingJetPT[k]->SetFillStyle(3444);
	  LeadingJetPT[k]->SetLineWidth(3);
	  LeadingJetPT[k]->Write();
	  //THT
	  THT[k]->SetFillColor(kSpring);
	  THT[k]->SetFillStyle(3444);
	  THT[k]->SetLineWidth(3);
	  THT[k]->Write();
	}
      //Settings for Single Top
      else if (k<=8 && k>=3) 
	{
	  FiveJetsMass[k]->SetFillColor(kBlack);
	  FiveJetsMass[k]->SetFillStyle(3305);
	  LeadingJetPT[k]->SetFillColor(kBlack);
	  LeadingJetPT[k]->SetFillStyle(3305);
	  THT[k]->SetFillColor(kBlack);
	  THT[k]->SetFillStyle(3305);
	  if (k!=3) 
	    {
	      SingleTopTprimeMass->Add(FiveJetsMass[k]); 
	      SingleTopTprimeMassLJPT->Add(LeadingJetPT[k]);
	      SingleTopTprimeMassTHT->Add(THT[k]);
	    }
	  if (k==8)
	    {
	      TFile f("SingleTop.root", "RECREATE");
	      SingleTopTprimeMass->Write();
	      SingleTopTprimeMassLJPT->Write();
	      SingleTopTprimeMassTHT->Write();
	    }
	  //Singletop->Add(FiveJetsMass[k]);
	}
      //Settings for TTbar
      else if (k==9) 
	{
	  FiveJetsMass[k]->SetFillColor(kRed);
	  FiveJetsMass[k]->SetFillStyle(3345);
	  LeadingJetPT[k]->SetFillColor(kRed);
	  LeadingJetPT[k]->SetFillStyle(3345);
	  THT[k]->SetFillColor(kRed);
	  THT[k]->SetFillStyle(3345);
	  TFile f("TTJets.root", "RECREATE");
	  FiveJetsMass[k]->Write();
	  LeadingJetPT[k]->Write();
	  THT[k]->Write();
	}
      //Settings for QCD
      else if (k>=10 && k<=14) 
	{
	  FiveJetsMass[k]->SetFillColor(kViolet); 
	  LeadingJetPT[k]->SetFillColor(kViolet); 
	  THT[k]->SetFillColor(kViolet); 
	  if (k!=10) 
	    {
	      QCDTprimeMass->Add(FiveJetsMass[k]);
	      QCDTprimeMassLJPT->Add(LeadingJetPT[k]);
	      QCDTprimeMassTHT->Add(THT[k]);
	    }
	  if (k==14)
	    {
	      TFile f("QCD.root", "RECREATE"); 
	      QCDTprimeMass->Write();
	      QCDTprimeMassLJPT->Write();
	      QCDTprimeMassTHT->Write();
	    }
	}
      //Settings for DiBoson
      else if (k<=2 && k>=0) 
	{
	  FiveJetsMass[k]->SetFillColor(kWhite); 
	  LeadingJetPT[k]->SetFillColor(kWhite); 
	  THT[k]->SetFillColor(kWhite); 
	  if (k!=0) 
	    {
	      DibosonTprimeMass->Add(FiveJetsMass[k]);
	      DibosonTprimeMassLJPT->Add(LeadingJetPT[k]);
	      DibosonTprimeMassTHT->Add(THT[k]);
	    }
	  if (k==2)
	    {
	      TFile f("Diboson.root", "RECREATE"); 
	      DibosonTprimeMass->Write();
	      DibosonTprimeMassLJPT->Write();
	      DibosonTprimeMassTHT->Write();
	    }
	}
      BKGandSignal->Add(FiveJetsMass[k]);
      BKGandSignallegend->AddEntry(FiveJetsMass[k]);
    }
  
}
