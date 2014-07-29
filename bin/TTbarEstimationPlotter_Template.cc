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
const int NOS=3;
//const TString Samples[NOS] = {"Diboson.root", "SingleTop.root", "TTJets.root", "Zjets.root", "Wjets.root", "QCD.root", "Signal.root"};
//const TString Samples[NOS] = {"Diboson.root", "SingleTop.root", "TTJets.root", "DY.root", "QCD.root", "Signal.root"};
//const TString Samples[NOS] = {"Diboson.root", "SingleTop.root", "TTJets.root", "QCD.root", "Signal.root"};
const TString Samples[NOS] = {"TTJets.root", "QCD.root", "Signal.root"};

// Breit-Wigner Gaussian fitting function definition
double GaussBreitW(Double_t *xx, Double_t *par)
//----------------------------------------------------------------------
{
  double M = xx[0];
  double M0 = par[1], M02=M0*M0;
  double G0 = 2.49*M0/175.16, G02=G0*G0;
  //double G0 = 2.49*M0/par[3], G02=G0*G0;
  double sigma = par[2], sigma2=sigma*sigma;
  
  double sum = 0;
  for(int i=0 ; i<=15 ; i++)
    {
      double m = i*sigma/5;
      double u = (M+m)*(M+m)-M02;
      double v = (M-m)*(M-m)-M02;
      double BW = 1/(u*u+G02*M02) + 1/(v*v+G02*M02);
      if( i==0 ) BW /= 2;
      sum += BW*exp(-m*m/sigma2/2);
    }
  sum *= G02*M02/31/sqrt(2*acos(-1.))/sigma;
  return par[0]*sum;
}

double Gaussian(Double_t *x, Double_t *par)
{
  return (par[0]*exp(-0.5*((x[0]-par[1])/par[2])*((x[0]-par[1])/par[2])));
}

double bkg(Double_t *x, Double_t *par)
{
  return par[0]+par[1]*x[0];
}

double fit_function(Double_t *x, Double_t *par)
{
  //return GaussBreitW(x,&par[2]);
  //return bkg(x,par)+GaussBreitW(x,&par[2]);
  return bkg(x,par)+Gaussian(x,&par[2]);
}

int theNumber(TFile * theFile){
  //TH1F*histogram;
  for(int k=0;k<=19;k++){
    string name="TprimeMass";
    char kchar[20];sprintf(kchar,"%d",k);
    if((TH1F*)(theFile->Get((name+kchar).c_str()))) return k;
  }
  return -1;
}

void TTbarEstimationPlotter()
{
  int ReBin=2;

  TFile *CurrentFile[NOS];
  TH1F *TprimeHistos[NOS];
  TH1F *FiveJetsMassBE[NOS];
  TH1F *FiveJetsMassLC[NOS];
  TH1F *TopMassFromHiggs[NOS];
  TH1F *FiveJetsMassLCoverBE[NOS];
  
  TH1F *FJMLCoverBE; TH1F *FJMBE; TH1F *TFH;
  TH1D *TFH_TTbar;

  //QCD estimation
  TH1F *HMLC[NOS]; TH1F *HMLCFull; TH1F * HMLCQCD;
  TH1F *FJLCQCDBE[NOS]; TH1F *FJLCQCDBEstim;
  TH1F *FJLCQCDLC[NOS];
  TH1F *FJLCQCDLCoverBE;

  for (int i=0; i<NOS; i++)
    {
      CurrentFile[i] = new TFile(StorageDirPrefix + Samples[i], "READ");
      if ( CurrentFile[i]->IsOpen() ) printf( Samples[i] + " File opened successfully\n");
      int NumberSuffix=theNumber(CurrentFile[i]);
      cout << NumberSuffix << endl;
      //
      string A=Form("TprimeMass%i",NumberSuffix);
      TprimeHistos[i]->SetDefaultSumw2(); TprimeHistos[i]= (TH1F*)gDirectory->Get(A.c_str()); TprimeHistos[i]->Rebin(ReBin);
      //
      A=Form("TprimeMassBkgE%i",NumberSuffix);
      FiveJetsMassBE[i]->SetDefaultSumw2(); FiveJetsMassBE[i]= (TH1F*)gDirectory->Get(A.c_str()); FiveJetsMassBE[i]->Rebin(ReBin);
      //cout << "FiveJetsMassBE Nbins: " << FiveJetsMassBE[i]->GetNbinsX() << endl;
      if (i==0) {FJMBE->SetDefaultSumw2(); FJMBE= (TH1F*)FiveJetsMassBE[i]->Clone("FJMBE");}
      else {FJMBE->Add(FiveJetsMassBE[i]);}
      //
      A=Form("TprimeMassLC%i",NumberSuffix);
      FiveJetsMassLC[i]->SetDefaultSumw2(); FiveJetsMassLC[i]= (TH1F*)gDirectory->Get(A.c_str()); FiveJetsMassLC[i]->Rebin(ReBin);
      //cout << "FiveJetsMass Nbins: " << FiveJetsMassLC[i]->GetNbinsX() << endl;
      if (i==0) {FJMLCoverBE->SetDefaultSumw2(); FJMLCoverBE= (TH1F*)FiveJetsMassLC[i]->Clone("FJMLCoverBE");}
      //else {FJMLCoverBE->Add(FiveJetsMassBE[i]);}
      //
      //
      A=Form("TopFromHiggsChi2Mass%i",NumberSuffix);
      TopMassFromHiggs[i]->SetDefaultSumw2(); TopMassFromHiggs[i]= (TH1F*)gDirectory->Get(A.c_str()); //TopMassFromHiggs[i]->Rebin(ReBin);
      //cout << "TopMassFromHiggs Nbins: " << TopMassFromHiggs[i]->GetNbinsX() << endl;
      if (i==0) {TFH->SetDefaultSumw2(); TFH= (TH1F*)TopMassFromHiggs[i]->Clone("TFH"); TFH_TTbar->SetDefaultSumw2(); TFH_TTbar= (TH1D*)TopMassFromHiggs[i]->Clone("TFH_TTbar");}
      else {TFH->Add(TopMassFromHiggs[i]);}
      //
      A=Form("TprimeMassLC%i",NumberSuffix);
      FiveJetsMassLCoverBE[i]->SetDefaultSumw2(); FiveJetsMassLCoverBE[i]= (TH1F*)gDirectory->Get(A.c_str()); FiveJetsMassLCoverBE[i]->Rebin(ReBin);
      //CurrentFile[i]->Close();
    }

  //cout << TFH->GetXaxis()->FindBin(220) << " " << TFH->GetBinContent(TFH->GetXaxis()->FindBin(220)) << endl;
  //TFH->SetBinContent(TFH->GetXaxis()->FindBin(220),(TFH->GetXaxis()->FindBin(210)+TFH->GetXaxis()->FindBin(230))/2);

  /*//Calculation of Estimator
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
  cout << "Number of signal events: " << sigFR << ", Number of BKG events: " << bkgFR[0]+bkgFR[1]+bkgFR[2] << ", S/sqrt(S+B)=" << sigFR/sqrt(sigFR+bkgFR[0]+bkgFR[1]+bkgFR[2]) << endl;*/

  THStack *TprimeHistosStack = new THStack("TprimeHistosStack", "BKG and signal for 5 jets; M_{5j} GeV; Events");
  THStack *FiveJetsMassBEStack = new THStack("FiveJetsMassBEStack", "BKG and signal for 5 jets; M_{5j} GeV; Events");
  THStack *FiveJetsMassLCStack = new THStack("FiveJetsMassLCStack", "BKG and signal for 5 jets; M_{5j} GeV; Events");
  THStack *TopFromHiggsChi2MassStack = new THStack("TopFromHiggsChi2MassStack", "BKG and signal for top from Higgs; M_{t}^{H} GeV; Events");
  
  cout << "1 Marker" << endl;
  for (int i=2; i>-1; i+=-1)
    {
      //cout << i << endl;
      TprimeHistosStack->Add(TprimeHistos[i]);
      FiveJetsMassBEStack->Add(FiveJetsMassBE[i]);
      FiveJetsMassLCStack->Add(FiveJetsMassLC[i]);
      TopFromHiggsChi2MassStack->Add(TopMassFromHiggs[i]);
      /*if (i!=0)
	{
	  FJMLCoverBE->Add(FiveJetsMassLC[i]);
	  cout << "a" << endl;
	  FJMBE->Add(FiveJetsMassBE[i]);
	  cout << "a" << endl;
	  TFH->Add(TopMassFromHiggs[i]);
	  cout << "a" << endl;
	  }*/
    }

  /*TCanvas *TFSCanvas1 = new TCanvas("TFSCanvas1","TFSCanvas1",600,800);
  TFSCanvas1->Clear(); TFSCanvas1->cd(1);
  TopFromHiggsChi2MassStack->Draw();
  TH1F *TFHM= (TH1F*)TopFromHiggsChi2MassStack->GetHistogram();
  TFSCanvas1->Close();*/

  TH1F *TFHC=(TH1F*)TFH->Clone("TFHClone");
  TFHC->GetXaxis()->SetRangeUser(140,230);
  TF1 *fit = new TF1("fitting",fit_function,120,800,5);
  TF1 *background = new TF1("bkgr",bkg,120,800,2);
  TF1 *peak = new TF1("Gauss",Gaussian,140,230,3);
  background->SetLineColor(1);
  peak->SetLineColor(4);
  fit->SetLineColor(2);
  Double_t par[5];
  fit->SetParameter(3,TFHC->GetMean());
  //fit->SetParameter(4,173.3);
  cout << "Mean value of TFH: " << TFHC->GetMean() << endl;
  fit->SetParameter(4,TFHC->GetRMS());
  cout << "RMS value of TFH: " << TFHC->GetRMS() << endl;
  TFH->Fit("fitting","emr");
  //TFH->Fit("gaus","wm");
  fit->GetParameters(par);
  background->SetParameters(par);
  peak->SetParameters(&par[2]);
  //background->Draw("same");
  //peak->Draw("same");

  float lowerlim=par[3]-3*par[4];
  float upperlim=par[3]+3*par[4];
  double ErrorFullIntegral=0;
  double HistoInt=TFH->IntegralAndError(TFH->GetXaxis()->FindBin(lowerlim),TFH->GetXaxis()->FindBin(2000.),ErrorFullIntegral);
  cout << "Integral of the full fit function all samples: " << fit->Integral(lowerlim,upperlim) << ", with error: " << fit->IntegralError(lowerlim,upperlim) << endl;
  cout << "Integral of Full function - flat background: " << fit->Integral(lowerlim,upperlim)-background->Integral(lowerlim,upperlim) << ", with error: " << 2*fit->IntegralError(lowerlim,upperlim) << endl;
  cout << "Integral of the peak of the fit all samples: " << peak->Integral(lowerlim,upperlim) << ", with error: " << 2*fit->IntegralError(lowerlim,upperlim) << endl;
  double chi2=fit->GetChisquare(); cout << "Chi2 of the full fitting: " << chi2 << endl;

  //Fitting only ttbar
  TF1 *peak_ttbar = new TF1("Gauss",Gaussian,140,230,3);
  Double_t par_ttbar[3];
  TH1F *TFH_TTbarC=(TH1F*)TFH_TTbar->Clone("TFH_TTbarC");
  TFH_TTbarC->GetXaxis()->SetRangeUser(140,230);
  peak_ttbar->SetParameter(2,TFH_TTbarC->GetMean());
  peak_ttbar->SetParameter(3,TFH_TTbarC->GetRMS());
  TFH_TTbarC->Fit("Gauss","emr");
  peak_ttbar->GetParameters(par_ttbar);
  //peak_ttbar->SetParameters(par_ttbar);
  par_ttbar[2]=abs(par_ttbar[2]);
  float lowerlimit=par_ttbar[1]-3*par_ttbar[2];
  float upperlimit=par_ttbar[1]+3*par_ttbar[2];
  double ErrorFullHistoIntegral=0;
  double HistoIntegral=TFH_TTbar->IntegralAndError(TFH_TTbar->GetXaxis()->FindBin(lowerlim),TFH_TTbar->GetXaxis()->FindBin(2000.),ErrorFullHistoIntegral);
  double ErrorFullHistoIntegralWin=0;
  double HistoIntegralWin=TFH_TTbar->IntegralAndError(TFH_TTbar->GetXaxis()->FindBin(lowerlimit),TFH_TTbar->GetXaxis()->FindBin(upperlimit),ErrorFullHistoIntegralWin);
  cout << "Integration limits ttbar only: " << lowerlimit << " " << upperlimit << endl;
  cout << "Full integral of the histo ttbar only: " << HistoIntegral << ", with error: " << ErrorFullHistoIntegral << endl; 
  cout << "Peak integral of the histo ttbar only: " << HistoIntegralWin << ", with error: " << ErrorFullHistoIntegralWin << endl; 
  cout << "Integral of the peak of the fit ttbar only: " << peak_ttbar->Integral(lowerlimit,upperlimit) << ", with error: " << peak_ttbar->IntegralError(lowerlimit,upperlimit) << endl;
  double chi2_peakttbar=peak_ttbar->GetChisquare(); cout << "Chi2 of the full fitting: " << chi2_peakttbar << endl;

  //double peak_histointe_ratio=((peak_ttbar->Integral(lowerlimit,upperlimit)/10)/HistoIntegral); //Calculating ratio form ttbar MC
  double peak_histointe_ratio=(HistoIntegralWin/HistoIntegral); //Calculating ratio form ttbar MC
  cout << "Ratio R=" << peak_histointe_ratio << endl;

  //cout << "2 Marker" << endl;

  double Np=(fit->Integral(lowerlim,upperlim)/10); //Calculating Number of events on the peak of the fit
  cout << "Weight factor for ttbar estimation: " << (Np*((1-peak_histointe_ratio)/peak_histointe_ratio))/FJMBE->Integral() << endl;
  FJMBE->Scale((Np*((1-peak_histointe_ratio)/peak_histointe_ratio))/FJMBE->Integral());

  //TH1F *FakeDY; FakeDY->SetFillColor(kBlue);
  //TH1F *FakeST; FakeST->SetFillColor(kBlack); FakeST->SetFillStyle(3305);
  //TH1F *FakeDB; FakeDB->SetFillColor(kWhite);

  //PostScript Plotting
  //cout << "MARKER" << endl;
  TLegend* BKGandSignallegend = new TLegend(0.75,0.65,0.90,0.9);
  if (NOS==5) 
    {
      BKGandSignallegend->AddEntry(TprimeHistos[3], "QCD", "f");
      //BKGandSignallegend->AddEntry(FakeDY, "Zjets", "f");      
      BKGandSignallegend->AddEntry(TprimeHistos[2], "TTbar", "f");
      BKGandSignallegend->AddEntry(TprimeHistos[1], "SingleT", "f");
      BKGandSignallegend->AddEntry(TprimeHistos[0], "Diboson", "f");
      BKGandSignallegend->AddEntry(TprimeHistos[4], "Signal", "f");
    }
  else if (NOS==3)
    {
      BKGandSignallegend->AddEntry(TprimeHistos[1], "QCD", "f");
      //BKGandSignallegend->AddEntry(FakeDY, "Zjets", "f");      
      BKGandSignallegend->AddEntry(TprimeHistos[0], "TTbar", "f");
      //BKGandSignallegend->AddEntry(FakeST, "SingleT", "f");
      //BKGandSignallegend->AddEntry(FakeDB, "Diboson", "f");
      BKGandSignallegend->AddEntry(TprimeHistos[2], "Signal", "f");
    }
  else 
    {
      BKGandSignallegend->AddEntry(TprimeHistos[4], "QCD", "f");
      BKGandSignallegend->AddEntry(TprimeHistos[3], "Zjets", "f");
      BKGandSignallegend->AddEntry(TprimeHistos[2], "TTbar", "f");
      BKGandSignallegend->AddEntry(TprimeHistos[1], "SingleT", "f");
      BKGandSignallegend->AddEntry(TprimeHistos[0], "Diboson", "f");
      BKGandSignallegend->AddEntry(TprimeHistos[5], "Signal", "f");
    }
  //cout << "MARKER" << endl;
  char FileName[100];
  sprintf(FileName,"TTbar_estimation_SUFFIX.pdf");
  //TPostScript *ps = new TPostScript(FileName,111);
  TPDF *ps = new TPDF(FileName,111);
  TCanvas *MyPlot = new TCanvas("MyPlot","Single t prime to top Higgs with backgrounds",600,800);
  TPaveText *p1;
  p1 = new TPaveText(3,0.5,8,3.5);
  char dateandtime[50];
  sprintf(dateandtime,"date: %s, time: %s",__DATE__,__TIME__);
  p1->AddText(dateandtime);

  //cout << "3 Marker" << endl;
  //////////////
  //First Page//
  //////////////
  MyPlot->Clear();
  MyPlot->cd(1);
  ps->NewPage();
  gStyle->SetOptStat(0);//Remove the Stat Box
  //TprimeHistos[NOS-1]->Draw("hist");
  //TprimeHistos[NOS-1]->GetYaxis()->SetRangeUser(0.1,BKGandSignal->GetMaximum());
  TprimeHistosStack->Draw("hist");
  //BKGandSignal->GetYaxis()->SetRangeUser(0.1,BKGandSignal->GetMaximum());
  //TprimeHistos[NOS-1]->Draw("histsame");
  //gPad->SetLogy();
  //TprimeHistos[NOS-1]->Draw("histsame"); //Preserve double drawing of signal
  TprimeHistosStack->SetMinimum(0.1);
  BKGandSignallegend->Draw();
  gPad->RedrawAxis("");
  gPad->Update();
  MyPlot->Update();
  cout << "First done" << endl;
  ///////////////
  //Second Page//
  ///////////////
  MyPlot->Clear();
  MyPlot->cd(1);
  ps->NewPage();
  FiveJetsMassBEStack->Draw("hist");
  //gPad->SetLogy();
  FiveJetsMassBEStack->SetMinimum(0.1);
  BKGandSignallegend->Draw();
  gPad->RedrawAxis("");
  gPad->Update();
  MyPlot->Update();
  cout << "Second done" << endl;
  //
  MyPlot->Clear();
  MyPlot->cd(1);
  ps->NewPage();
  FiveJetsMassLCStack->Draw("hist");
  //gPad->SetLogy();
  FiveJetsMassLCStack->SetMinimum(0.1);
  BKGandSignallegend->Draw();
  gPad->RedrawAxis("");
  gPad->Update();
  MyPlot->Update();
  cout << "Third done" << endl;
  //
  MyPlot->Clear();
  MyPlot->cd(1);
  ps->NewPage();
  TopFromHiggsChi2MassStack->Draw("hist");
  background->Draw("same");
  peak->Draw("same");
  fit->Draw("same");
  //gPad->SetLogy();
  TopFromHiggsChi2MassStack->SetMinimum(0.1);
  BKGandSignallegend->Draw();
  gPad->RedrawAxis("");
  gPad->Update();
  MyPlot->Update();
  cout << "Fourth done" << endl;
  //
  MyPlot->Clear();
  MyPlot->cd(1);
  ps->NewPage();
  TFH->Draw("hist");
  //gPad->SetLogy();
  TFH->SetMinimum(0.1);
  BKGandSignallegend->Draw();
  gPad->RedrawAxis("");
  gPad->Update();
  MyPlot->Update();
  cout << "Fifth done" << endl;
  //
  MyPlot->Clear();
  MyPlot->cd(1);
  ps->NewPage();
  TFH_TTbar->Draw("hist");
  peak_ttbar->Draw("same");
  gPad->Update();
  MyPlot->Update();
  cout << "Sixth done" << endl;
  //  
  MyPlot->Clear();
  MyPlot->cd(1); 
  ps->NewPage();
  FJMBE->Draw("hist");
  gPad->SetLogy();
  FJMBE->SetMinimum(0.1);
  //BKGandSignallegend->Draw();
  gPad->RedrawAxis("");
  gPad->Update();
  MyPlot->Update();
  cout << "Seventh done" << endl;
  //
  TPad *pad1 = new TPad("pad1","pad1",0,0.3,1,1);
  MyPlot->Clear();
  MyPlot->cd(1);
  ps->NewPage();
  pad1->SetBottomMargin(0); 
  pad1->Draw();
  pad1->cd(); 
  gStyle->SetOptStat(0);
  FJMLCoverBE->DrawCopy("e0 hist");
  FJMLCoverBE->SetMinimum(0.5);
  FJMLCoverBE->SetMaximum(1.5);
  //FiveJetsMassLCStack->Draw("hist same");
  //FiveJetsMassLCStack->SetMinimum(1);
  FJMBE->Draw("E1 same");
  gPad->SetLogy();
  gPad->Update();
  BKGandSignallegend->Draw();
  MyPlot->cd();
  TPad *pad2 = new TPad("pad2","pad2",0,0,1,0.3);
  pad2->SetTopMargin(0);
  pad2->Draw();
  pad2->cd();
  FJMLCoverBE->Sumw2();
  FJMLCoverBE->SetStats(0);
  FJMLCoverBE->SetTitle("");
  FJMLCoverBE->Divide(FJMBE);
  //FJMLCoverBE->Divide(FJMLCoverBE,FJMBE,1,1,"w");
  FJMLCoverBE->SetMarkerStyle(21);
  FJMLCoverBE->Draw("ep");
  FJMLCoverBE->SetMaximum(1.5);
  FJMLCoverBE->SetMinimum(0.5);
  TLine *L1 = new TLine(450,1.3,1600,1.3);
  TLine *L2 = new TLine(450,1.0,1600,1.0);
  TLine *L3 = new TLine(450,0.7,1600,0.7);
  L1->Draw(); L2->Draw(); L3->Draw();
  MyPlot->cd();
  MyPlot->Update();

  /*MyPlot->Clear();
  MyPlot->cd(1);
  ps->NewPage();
  BKGL4JPT->Draw("hist");
  Leading4JetPT[NOS-1]->Draw("histsame");
  gPad->SetLogy();
  BKGL4JPT->SetMinimum(0.1);
  BKGandSignallegend->Draw();
  gPad->RedrawAxis("");
  gPad->Update();
  MyPlot->Update();*/
  //
  ps->Close();

  exit(0);
}
