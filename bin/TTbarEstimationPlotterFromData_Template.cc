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

#include "TFitResultPtr.h"
#include "TMatrixDSym.h"
#include "TMath.h"
#include "TF1.h"

using namespace std;

const TString StorageDirPrefix="file:/sps/cms/ruizalva/CMSSW_5_3_9_patch3/src/Extractors/PatExtractor/bin/WithTrigger/MCSUFFIX/";
const TString DataStorageDirPrefix="file:/sps/cms/ruizalva/CMSSW_5_3_9_patch3/src/Extractors/PatExtractor/bin/WithTrigger/DATASUFFIX/";
const int NOS=3;
//const TString Samples[NOS] = {"Diboson.root", "SingleTop.root", "TTJets.root", "Zjets.root", "Wjets.root", "QCD.root", "Signal.root"};
//const TString Samples[NOS] = {"Diboson.root", "SingleTop.root", "TTJets.root", "DY.root", "QCD.root", "Signal.root"};
//const TString Samples[NOS] = {"Diboson.root", "SingleTop.root", "TTJets.root", "QCD.root", "Signal.root"};
const TString Samples[NOS] = {"TTJets.root", "QCD.root", "Signal.root"};
const TString DataSamples = "Data.root";

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

double Gaussian2(Double_t *x, Double_t *par)
{
  return (par[0]*exp(-0.5*((x[0]-par[1])/par[2])*((x[0]-par[1])/par[2])));
}

double bkg(Double_t *x, Double_t *par)
{
  return par[0]+par[1]*x[0];
}

double bkg2(Double_t *x, Double_t *par)
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

void TTbarEstimationPlotterFromData()
{
  int ReBin=1;

  TFile *CurrentFile[NOS];
  TH1F *TprimeHistos[NOS];
  TH1F *FiveJetsMassBE[NOS];
  TH1F *FiveJetsMassLC[NOS];
  TH1F *TopMassFromHiggs[NOS];
  TH1F *FiveJetsMassLCoverBE[NOS];
  
  TH1F *FJMLCoverBE; TH1F *FJMBE; TH1F *TFH;
  TH1D *TFH_TTbar;

  TH1F *TopMass_LC;

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
      //
      A=Form("TMass_LC_HMW%i",NumberSuffix);
      if (i==0) {TopMass_LC->SetDefaultSumw2(); TopMass_LC= (TH1F*)gDirectory->Get(A.c_str()); TopMass_LC->Rebin(ReBin);}      
      //CurrentFile[i]->Close();      
    }

  //////////////////////////////////////////////////////
  //Reading Data Files//////////////////////////////////
  //////////////////////////////////////////////////////
  TFile *DataFile = new TFile(DataStorageDirPrefix + DataSamples, "READ");
  if ( DataFile->IsOpen() ) printf( DataSamples + " File opened successfully\n");
  string A=Form("TprimeMassBkgE");
  TH1F *FJM_ttbar_Estim= (TH1F*)gDirectory->Get(A.c_str()); FJM_ttbar_Estim->Rebin(ReBin);
  A=Form("TopFromHiggsChi2Mass");
  TH1F *TopMassFromHiggsData= (TH1F*)gDirectory->Get(A.c_str());
  A=Form("TMass_LC");
  TH1F *TopMass_LC_Data= (TH1F*)gDirectory->Get(A.c_str());
  A=Form("TMass_BE");
  TH1F *TopMass_BE_Data= (TH1F*)gDirectory->Get(A.c_str());

  /////////////////////////////////////////////////////  

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
    }

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
  cout << "Integral of the full fit function all samples: " << fit->Integral(lowerlim,upperlim)/10 << ", with error: " << fit->IntegralError(lowerlim,upperlim)/10 << endl;
  cout << "Integral of Full function - flat background: " << (fit->Integral(lowerlim,upperlim)-background->Integral(lowerlim,upperlim))/10 << ", with error: " << 2*fit->IntegralError(lowerlim,upperlim)/10 << endl;
  cout << "Integral of the peak of the fit all samples: " << peak->Integral(lowerlim,upperlim)/10 << ", with error: " << 2*fit->IntegralError(lowerlim,upperlim)/10 << endl;
  cout << "Integral of the background of the fit all samples: " << background->Integral(lowerlim,upperlim)/10 << ", with error: " << fit->IntegralError(lowerlim,upperlim)/10 << endl;
  
  double chi2=fit->GetChisquare(); double ndf_fit=fit->GetNDF(); cout << "Chi2 of the full fitting: " << chi2 << ", with number of degrees of freedom: " << ndf_fit << endl;

  //Fitting only ttbar
  TF1 *peak_ttbar = new TF1("Gauss",Gaussian,140,230,3);
  Double_t par_ttbar[3];
  TH1F *TFH_TTbarC=(TH1F*)TFH_TTbar->Clone("TFH_TTbarC");
  TFH_TTbarC->GetXaxis()->SetRangeUser(140,230);
  peak_ttbar->SetParameter(2,TFH_TTbarC->GetMean());
  cout << "Mean value of TFH_TTbarC: " << TFH_TTbarC->GetMean() << endl;
  peak_ttbar->SetParameter(3,TFH_TTbarC->GetRMS());
  cout << "RMS value of TFH_TTbarC: " << TFH_TTbarC->GetRMS() << endl;
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
  cout << "Integral of the peak of the fit ttbar only: " << peak_ttbar->Integral(lowerlimit,upperlimit)/10 << ", with error: " << peak_ttbar->IntegralError(lowerlimit,upperlimit)/10 << endl;
  double chi2_peakttbar=peak_ttbar->GetChisquare(); double ndf_peakttbar=peak_ttbar->GetNDF(); cout << "Chi2 of the ttbar peak fitting: " << chi2_peakttbar << ", with number of degrees of freedom: " << ndf_peakttbar << endl;

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

  //////////////////////////
  //Fit on Data!////////////
  //////////////////////////  
  TH1F *TFHCD=(TH1F*)TopMassFromHiggsData->Clone("TFHCloneD");
  TH1F *TFHCDT=(TH1F*)TopMassFromHiggsData->Clone("TFHCloneDT");
  TFHCD->GetXaxis()->SetRangeUser(140,230);
  TFHCDT->GetXaxis()->SetRangeUser(230,1000);
  //Fit on peak + tail
  TF1 *datafit = new TF1("data_fitting",fit_function,100,1000,5);
  TF1 *databackgroundF = new TF1("databkgrF",bkg,100,1000,2);
  TF1 *datapeakF = new TF1("dataGaussF",Gaussian,140,230,3);
  databackgroundF->SetLineColor(1);
  datapeakF->SetLineColor(1);
  datafit->SetLineColor(1);
  Double_t dataparF[5];
  Double_t datapeakparF[3];
  Double_t databkgparF[2];  
  datafit->SetParameter(3,TFHCD->GetMean());
  cout << "Mean value of TFHCD: " << TFHCD->GetMean() << endl;
  datafit->SetParameter(4,TFHCD->GetRMS());
  cout << "RMS value of TFHCD: " << TFHCD->GetRMS() << endl;
  ////TFHCD->Fit("data_fitting","emr");
  TopMassFromHiggsData->Fit("data_fitting","emr");
  datafit->GetParameters(dataparF);
  databackgroundF->SetParameters(dataparF);
  datapeakF->SetParameters(&dataparF[2]);
  float datalowerlimitF=dataparF[3]-3*dataparF[4];
  float dataupperlimitF=dataparF[3]+3*dataparF[4];
  double data_bkg_integralF=databackgroundF->Integral(datalowerlimitF,dataupperlimitF)/10;
  double data_peak_integralF=(datapeakF->Integral(datalowerlimitF,dataupperlimitF)/10);
  //////////////////////
  //Separate fit on peak and tail
  TF1 *databackground = new TF1("databkgr",bkg2,100,1000,2);
  TF1 *datapeak = new TF1("dataGauss",Gaussian2,140,230,3);
  datapeak->SetLineColor(4);
  databackground->SetLineColor(2);
  Double_t datapar[5];
  Double_t datapeakpar[3];
  Double_t databkgpar[2];
  datapeak->SetParameter(2,TFHCD->GetMean());
  datapeak->SetParameter(3,TFHCD->GetRMS());
  TFHCD->Fit("dataGauss","emr");
  TFHCDT->Fit("databkgr","emr");
  datapeak->GetParameters(datapeakpar);
  databackground->GetParameters(databkgpar);
  datapar[0]=databkgpar[0]; datapar[1]=databkgpar[1]; datapar[2]=datapeakpar[0]; datapar[3]=datapeakpar[1]; datapar[4]=datapeakpar[2];
  //datafit->SetParameters(datapar);
  float datalowerlimit=datapar[3]-3*datapar[4];
  float dataupperlimit=datapar[3]+3*datapar[4];
  double data_bkg_integral=databackground->Integral(datalowerlimit,dataupperlimit)/10;
  double data_peak_integral=(datapeak->Integral(datalowerlimit,dataupperlimit)/10);
  //
  cout << "################################USING SAME METHOD THAN FOR MC################################" << endl;
  cout << "Integral from the peak of fit on data: " << data_peak_integralF << ", with limits: (" <<  datalowerlimitF << "," <<  dataupperlimitF << ")" << endl;
  cout << "Integral from the bkg fit on data: " << data_bkg_integralF << ", with limits: (" <<  datalowerlimitF << "," <<  dataupperlimitF << ")" << endl;
  cout << "Integral from the (peak-bkg) of fit on data: " << data_peak_integralF-data_bkg_integralF << ", with limits: (" <<  datalowerlimitF << "," <<  dataupperlimitF << ")" << endl;
  cout << "Integral of the full fit function all samples: " << datafit->Integral(datalowerlimitF,dataupperlimitF)/10 << ", with error: " << datafit->IntegralError(datalowerlimitF,dataupperlimitF)/10 << endl;
  cout << "Integral of Full function - flat background: " << (datafit->Integral(datalowerlimitF,dataupperlimitF)-databackgroundF->Integral(datalowerlimitF,dataupperlimitF))/10 << ", with error: " << 2*datafit->IntegralError(datalowerlimitF,dataupperlimitF)/10 << endl;
  cout << "############################USING FLAT AND GAUSSIAN INTERSECTIONS############################" << endl;
  double intersections[2] = {148.456,233.101};
  double data_bkg_integral_inter=databackground->Integral(intersections[0],intersections[1])/10;
  double data_peak_integral_inter=(datapeak->Integral(intersections[0],intersections[1])/10);
  
  cout << "Value of gaussian and flat function at intersections: " << datapeak->Eval(intersections[0]) << " " << databackground->Eval(intersections[0]) << " " << datapeak->Eval(intersections[1]) << " " << databackground->Eval(intersections[1]) << endl;
  cout << "Integral from the peak of fit on data: " << data_peak_integral_inter << ", with error: " << datapeak->IntegralError(intersections[0],intersections[1])/10 << ", with limits: (" << intersections[0] << "," << intersections[1] << ")" << endl;
  cout << "Integral from the bkg fit on data: " << data_bkg_integral_inter << ", with error: " << databackground->IntegralError(intersections[0],intersections[1])/10 << ", with limits: (" << intersections[0] << "," << intersections[1] << ")" << endl;
  cout << "Integral from the (peak-bkg) of fit on data: " << data_peak_integral_inter-data_bkg_integral_inter  << ", with error: " << (datapeak->IntegralError(intersections[0],intersections[1])+databackground->IntegralError(intersections[0],intersections[1]))/10 << ", with limits: (" << intersections[0] << "," << intersections[1] << ")" << endl;
  cout << "#############################################################################################" << endl;
  //  
  cout << "Chi2 of peak fitting on data: " << datapeak->GetChisquare() << ", with number of degrees of freedom: " << datapeak->GetNDF() << endl;
  cout << "Chi2 of tail fitting on data: " << databackground->GetChisquare()  << ", with number of degrees of freedom: " << databackground->GetNDF() << endl;

  double data_pure_peak=data_peak_integral_inter-data_bkg_integral_inter;
  /////////////////////////////////////////////////////////////////////////////////////////
  TH1F *FJM_Estim_DATAoverMC=(TH1F*)FJM_ttbar_Estim->Clone("FJMEstim_DATAoverMC"); FJM_Estim_DATAoverMC->Rebin(ReBin);
  peak_histointe_ratio=(data_pure_peak/FJM_Estim_DATAoverMC->Integral())*(1./((FiveJetsMassLC[0]->Integral()/FJM_Estim_DATAoverMC->Integral())+(data_pure_peak/FJM_Estim_DATAoverMC->Integral()))); //0.41388; //Correction of normalization of estimation!!!!!!!!!!!!!!!!
  cout << "Corrected R= " << peak_histointe_ratio << ", with correction factor: " << (HistoIntegralWin/HistoIntegral)/peak_histointe_ratio << endl;
  /////////////////////////////////////////////////////////////////////////////////////////
  //FJM_ttbar_Estim->Scale((data_pure_peak*((1-peak_histointe_ratio)/peak_histointe_ratio))/FJM_ttbar_Estim->Integral());
  cout << "Overall factor between data and MC before rescale: " << FJM_Estim_DATAoverMC->Integral()/FiveJetsMassLC[0]->Integral() << endl;
  cout << data_pure_peak << " " << peak_histointe_ratio << " " << FJM_Estim_DATAoverMC->Integral() << " " << data_pure_peak*((1.-peak_histointe_ratio)/peak_histointe_ratio) << ", with rescale factor: " << (data_pure_peak*((1.-peak_histointe_ratio)/peak_histointe_ratio))/FJM_Estim_DATAoverMC->Integral() << endl;
  FJM_Estim_DATAoverMC->Scale((data_pure_peak*((1-peak_histointe_ratio)/peak_histointe_ratio))/FJM_Estim_DATAoverMC->Integral());
  cout << "Overall factor between data and MC after rescale: " << FJM_Estim_DATAoverMC->Integral()/FiveJetsMassLC[0]->Integral() << endl;
  
  ////////////////////QCD ESTIMATION///////////////////////

  TH1F *TopMass_TTbarSubs=(TH1F*)TopMass_LC_Data->Clone("TopMass_TTbarSubs"); TopMass_TTbarSubs->Rebin(ReBin);
  cout << "Overall factor between data and MC in first top mass before rescale: " << TopMass_BE_Data->Integral()/TopMass_LC->Integral() << ", with rescale factor: " << (data_pure_peak*((1.-peak_histointe_ratio)/peak_histointe_ratio))/TopMass_BE_Data->Integral() << endl;  
  TopMass_BE_Data->Scale((data_pure_peak*((1-peak_histointe_ratio)/peak_histointe_ratio))/TopMass_BE_Data->Integral());
  cout << "Overall factor between data and MC in first top mass after rescale: " << TopMass_BE_Data->Integral()/TopMass_LC->Integral() << endl;
  
  TopMass_TTbarSubs->Add(TopMass_BE_Data,-1);

  TH1F *FirstTopMass=(TH1F*)TopMass_TTbarSubs->Clone("FirstTopMassSubs");
  FirstTopMass->GetXaxis()->SetRangeUser(230,2000);  
  TF1 *databackgroundtop1 = new TF1("databkgrtop1",bkg2,230,2000,2);
  databackgroundtop1->SetLineColor(6);
  Double_t databkgtop1par[2];
  FirstTopMass->Fit("databkgrtop1","emr");
  databackgroundtop1->GetParameters(databkgtop1par);

  cout << "QCD normalization inside the first top mass peak: " <<  databackgroundtop1->Integral(140.0,230.0) << endl;

  //////////////////////////

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
  sprintf(FileName,"Data_TTbar_estimation_MCSUFFIX_DATASUFFIX.pdf");
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
  /*MyPlot->Clear();
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
  cout << "Seventh done" << endl;*/
  //
  //PLOTING TOP MASS FROM HIGGS
  //
  MyPlot->Clear();
  MyPlot->cd(1);
  ps->NewPage();
  TopMassFromHiggsData->Draw("hist");
  databackground->Draw("same");
  datapeak->Draw("same");
  datafit->Draw("same");
  gPad->Update();
  MyPlot->Update();
  //
  //
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
  //gPad->SetLogy();
  gPad->Update();
  //BKGandSignallegend->Draw();
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
  //
  //PLOTING TTBAR ESTIMATION FROM DATA COMPARED TO MC
  //
  TPad *pad3 = new TPad("pad1","pad1",0,0.3,1,1);
  MyPlot->Clear();
  MyPlot->cd(1);
  ps->NewPage();
  pad3->SetBottomMargin(0); 
  pad3->Draw();
  pad3->cd(); 
  gStyle->SetOptStat(0);
  FJM_Estim_DATAoverMC->DrawCopy("e0 hist");
  FJM_Estim_DATAoverMC->SetMinimum(0.5);
  FJM_Estim_DATAoverMC->SetMaximum(1.5);
  //FiveJetsMassLCStack->Draw("hist same");
  //FiveJetsMassLCStack->SetMinimum(1);
  FiveJetsMassLC[0]->Draw("E1 same");
  gPad->SetLogy();
  gPad->Update();
  //BKGandSignallegend->Draw();
  MyPlot->cd();
  TPad *pad4 = new TPad("pad2","pad2",0,0,1,0.3);
  pad4->SetTopMargin(0);
  pad4->Draw();
  pad4->cd();
  FJM_Estim_DATAoverMC->Sumw2();
  FJM_Estim_DATAoverMC->SetStats(0);
  FJM_Estim_DATAoverMC->SetTitle(";M_{5j} GeV;Data/MC");
  FJM_Estim_DATAoverMC->Divide(FiveJetsMassLC[0]);
  //FJMLCoverBE->Divide(FJMLCoverBE,FJMBE,1,1,"w");
  FJM_Estim_DATAoverMC->SetMarkerStyle(21);
  FJM_Estim_DATAoverMC->Draw("ep");
  FJM_Estim_DATAoverMC->SetMaximum(1.5);
  FJM_Estim_DATAoverMC->SetMinimum(0.5);
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
