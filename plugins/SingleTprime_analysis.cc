#include "SingleTprime_analysis.h"

#include <iostream>
#include <vector>
#include <limits>
#include <algorithm>

#include "TH2.h"
#include "TString.h"
#include "TLorentzVector.h"
#include "TVector3.h"
#include <TMatrixDSym.h>
#include <TMatrixDSymEigen.h>

#include "Extractors/PatExtractor/interface/AnalysisSettings.h"
#include "Extractors/PatExtractor/interface/MCExtractor.h"
#include "Extractors/PatExtractor/interface/HLTExtractor.h"
#include "Extractors/PatExtractor/interface/MuonExtractor.h"
#include "Extractors/PatExtractor/interface/ElectronExtractor.h"
#include "Extractors/PatExtractor/interface/METExtractor.h"
#include "Extractors/PatExtractor/interface/VertexExtractor.h"
#include "Extractors/PatExtractor/interface/KinFit.h"
#include "Extractors/PatExtractor/interface/EventExtractor.h"
#include "Extractors/PatExtractor/interface/PatExtractor.h"

#include "CondFormats/JetMETObjects/interface/JetCorrectionUncertainty.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectorParameters.h"
#include "JetMETCorrections/Objects/interface/JetCorrectionsRecord.h"
#include <FWCore/Framework/interface/EventSetup.h>
#include "FWCore/ParameterSet/interface/ParameterSet.h"

//Constants 
/*int NumberOfGoodJets=5; //////////////////////////////////////////
int NumberOfBadJets=6;                                          //
float LeadingJetPt=150;                                         //
float THTcut=630;                                               //
int MinB_tags=3;
float MaxChi2=50;
float DeltaRHiggsJets=1.2;                                      //
float DeltaRWJets=3.0;                                          //
float HiggsPt=200;                                              //
float TopPt=200;                                                //
float MinDeltaRWH=2.7;                                          //
float MaxDeltaRWH=3.5;                                          //
float DeltaPhiHiggsJets=1.2; //            SET OF CUTS          //
float DeltaPhiTopJetW=1.2;                                      //
int JetMultiplicity=8;                                          //
float DeltaPhiWjets=1.3;                                        //
float MinHiggsMass=110;                                         //
float MaxHiggsMass=140;                                         //
float RelHT=0.65;                                               //
float Aplanarity=0.06;                                          //
float MaxDeltaRTH=3.3;                                          //
float MinDeltaRTH=2.8;                                          // 
float RelMassMaxCut=0.5;                                        //
float RelMassMinCut=0.3;                                        //
float MotherPTNormalizedMassCut=10.;                            //
float PTNormalizedMassCut=0.7;                                  //
int NumberOfTopsCut=2;                                          //
float NumberofLooseBtag=2;                                      //
float DeltaPhiLeadingJets=2.8;                                  //
float HMoverTM=0.6;              /////////////////////////////////
*/
float HiggsMass=125.0;
float HiggsMassWindow=2000.0;
float WMass=80.3;
float WMassWindow=2000.0;
float TopMass=172.5;
float TopMassWindow=2000.0;
float TriggerDiJet1_2=90;
float TriggerDiJet3_4=70;
float TriggerDiJet5_6=30;
float TopHadrWidth=17.35;
float WHadrWidth=10.12;
float HBBWidth=12.4;
float TopMassChi2=175.16;
float WMassChi2=84.06;
float HMassChi2=125.0;

using namespace std;

namespace patextractor {

  SingleTprime_analysis::SingleTprime_analysis(const edm::ParameterSet& cmsswSettings): Plugin(cmsswSettings)/*,
    jetEnergyResolutionScaleFactors_  (cmsswSettings.getParameter<std::vector<double> >("jetEnergyResolutionScaleFactors")),
    jetEnergyResolutionEtaBinning_    (cmsswSettings.getParameter<std::vector<double> >("jetEnergyResolutionEtaBinning"))*/
{
  std::cout << "Entering SingleTprime analysis" << std::endl;

  // Lorentz vectors

  ReconstructedHiggs = new TLorentzVector(0., 0., 0., 0.);
  ReconstructedW = new TLorentzVector(0., 0., 0., 0.);
  ReconstructedTop = new TLorentzVector(0., 0., 0., 0.);
  ReconstructedTprime = new TLorentzVector(0., 0., 0., 0.);
  TrueHiggs = new TLorentzVector(0., 0., 0., 0.);
  TrueW = new TLorentzVector(0., 0., 0., 0.);
  TrueTop = new TLorentzVector(0., 0., 0., 0.);
  TrueTprime = new TLorentzVector(0., 0., 0., 0.);
  TrueTprimeAcompainingJet = new TLorentzVector(0., 0., 0., 0.);
  FirstHiggsJet = new TLorentzVector(0., 0., 0., 0.);
  SecondHiggsJet = new TLorentzVector(0., 0., 0., 0.);
  FirstWJet = new TLorentzVector(0., 0., 0., 0.);
  SecondWJet = new TLorentzVector(0., 0., 0., 0.);
  TopJet = new TLorentzVector(0., 0., 0., 0.);
  FirstTrueHiggsJet = new TLorentzVector(0., 0., 0., 0.);
  SecondTrueHiggsJet = new TLorentzVector(0., 0., 0., 0.);
  FirstTrueWJet = new TLorentzVector(0., 0., 0., 0.);
  SecondTrueWJet = new TLorentzVector(0., 0., 0., 0.);
  TopTrueJet = new TLorentzVector(0., 0., 0., 0.);
  TprimeFromMatching = new TLorentzVector(0., 0., 0., 0.);

  /// Tree definition containing your analysis results

  m_tree_stp = new TTree("stp","Single Analysis info");  
  m_tree_cuts = new TTree("cuts","Cuts efficiencies");
  
  /// Branches definition

  m_tree_stp->Branch("Reconstructed_Higgs", &ReconstructedHiggs);
  m_tree_stp->Branch("Reconstructed_W", &ReconstructedW);
  m_tree_stp->Branch("Reconstructed_Top", &ReconstructedTop);
  m_tree_stp->Branch("Reconstructed_Tprime", &ReconstructedTprime);  
  m_tree_stp->Branch("True_Higgs", &TrueHiggs);
  m_tree_stp->Branch("True_W", &TrueW);
  m_tree_stp->Branch("True_Top", &TrueTop);
  m_tree_stp->Branch("True_Tprime", &TrueTprime);
  m_tree_stp->Branch("True_Tprime_Acompaining_Jet", &TrueTprimeAcompainingJet);
  m_tree_stp->Branch("First_Higgs_Jet", &FirstHiggsJet);
  m_tree_stp->Branch("Second_Higgs_Jet", &SecondHiggsJet);
  m_tree_stp->Branch("First_W_Jet", &FirstWJet);
  m_tree_stp->Branch("Second_W_Jet", &SecondWJet);
  m_tree_stp->Branch("Top_Jet", &TopJet);
  m_tree_stp->Branch("First_True_Higgs_Jet", &FirstTrueHiggsJet);
  m_tree_stp->Branch("Second_True_Higgs_Jet", &SecondTrueHiggsJet);
  m_tree_stp->Branch("First_True_W_Jet", &FirstTrueWJet);
  m_tree_stp->Branch("Second_True_W_Jet", &SecondTrueWJet);
  m_tree_stp->Branch("True_Top_Jet", &TopTrueJet);
  m_tree_stp->Branch("Tprime_From_Matching", &TprimeFromMatching);
  m_tree_stp->Branch("evt",         &m_evt    ,"evt/I");      // Simple evt number or event ID
  m_tree_stp->Branch("PU",         &m_nPU    ,"m_nPU/I");
  m_tree_stp->Branch("THT", &m_THT,"THT/F");
  m_tree_stp->Branch("jet1_pt",  &m_jet1pt   ,"jet1_pt/F");
  m_tree_stp->Branch("jet2_pt",  &m_jet2pt   ,"jet2_pt/F");
  m_tree_stp->Branch("jet3_pt",  &m_jet3pt   ,"jet3_pt/F");
  m_tree_stp->Branch("jet4_pt",  &m_jet4pt   ,"jet4_pt/F");
  m_tree_stp->Branch("jet5_pt",  &m_jet5pt   ,"jet5_pt/F");
  m_tree_stp->Branch("jet6_pt",  &m_jet6pt   ,"jet6_pt/F");
  m_tree_stp->Branch("Number_CSVLbtagged_jets",  &m_NBtaggedJets_CSVL   ,"m_NBtaggedJets_CSVL/I");
  m_tree_stp->Branch("Number_CSVMbtagged_jets",  &m_NBtaggedJets_CSVM   ,"m_NBtaggedJets_CSVM/I");
  m_tree_stp->Branch("Number_CSVTbtagged_jets",  &m_NBtaggedJets_CSVT   ,"m_NBtaggedJets_CSVT/I");
  m_tree_stp->Branch("DeltaR_of_Higgs_Jets",  &m_DRHiggsJets   ,"DRHiggsJets/F");
  m_tree_stp->Branch("DeltaR_of_W_Jets",  &m_DRWJets   ,"DRWJets/F");
  m_tree_stp->Branch("DeltaR_of_Top_Higgs",  &m_DRTopHiggs   ,"DRTopHiggs/F");
  m_tree_stp->Branch("DeltaR_of_W_Higgs",  &m_DRWHiggs   ,"DRWHiggs/F");
  m_tree_stp->Branch("Relative_THT",  &m_RelTHT   ,"RelTHT/F");
  m_tree_stp->Branch("Correct_Tprime",  &CorrectTprime   ,"CorrectTprime/I");
  m_tree_stp->Branch("Correct_Higgs",  &CorrectH   ,"CorrectH/I");
  m_tree_stp->Branch("Correct_W",  &CorrectW   ,"CorrectW/I");
  m_tree_stp->Branch("Correct_Top",  &CorrectTop   ,"CorrectTop/I");
  m_tree_stp->Branch("Correct_Higgs_Jet",  &CorrectHiggsJet   ,"CorrectHiggsJet/I");
  m_tree_stp->Branch("Correct_W_Jet",  &CorrectWJet   ,"CorrectWJet/I");
  m_tree_stp->Branch("Correct_Top_Jet",  &CorrectTopJet   ,"CorrectTopJet/I");
  m_tree_stp->Branch("DeltaR_TrueHiggs_RecoHiggs",  &m_DRTrueHiggsRecoHiggs   ,"DRTrueHiggsRecoHiggs/F");
  m_tree_stp->Branch("DeltaR_TrueW_RecoW",  &m_DRTrueWRecoW   ,"DRTrueWRecoW/F");
  m_tree_stp->Branch("DeltaR_TrueTop_RecoTop",  &m_DRTrueTopRecoTop   ,"DRTrueTopRecoTop/F");
  m_tree_stp->Branch("DeltaR_TrueTprime_RecoTprime",  &m_DRTrueTprimeRecoTprime   ,"DRTrueTprimeRecoTprime/F");
  m_tree_stp->Branch("DeltaR_TrueFirstHiggsJet_RecoJet",  &m_DRTrueFirstHiggsJetRecoJet   ,"DRTrueFirstHiggsJetRecoJet/F");
  m_tree_stp->Branch("DeltaR_TrueSecondHiggsJet_RecoJet",  &m_DRTrueSecondHiggsJetRecoJet   ,"DRTrueSecondHiggsJetRecoJet/F");
  m_tree_stp->Branch("DeltaR_TrueFirstWJet_RecoJet",  &m_DRTrueFirstWJetRecoJet   ,"DRTrueFirstWJetRecoJet/F");
  m_tree_stp->Branch("DeltaR_TrueSecondWJet_RecoJet",  &m_DRTrueSecondWJetRecoJet   ,"DRTrueSecondWJetRecoJet/F");
  m_tree_stp->Branch("DeltaR_TrueTopJet_RecoJet",  &m_DRTrueTopJetRecoJet   ,"DRTrueTopJetRecoJet/F");
  m_tree_stp->Branch("DeltaR_TrueWJets",  &m_DRTrueWJets   ,"DRTrueWJets/F");
  m_tree_stp->Branch("DeltaR_MatchedWJets",  &m_DRMatchedWJets   ,"DRMatchedWJets/F");
  m_tree_stp->Branch("DeltaPhi_TrueWJets",  &m_DPhiTrueWJets   ,"DRTrueWJets/F");
  m_tree_stp->Branch("DeltaPhi_MatchedWJets",  &m_DPhiMatchedWJets   ,"DRMatchedTopJets/F");
  m_tree_stp->Branch("Number_of_Matched_Higgs_Jets",  &NumbMatchedHiggsJets   ,"NumbMatchedHiggsJets/I");
  m_tree_stp->Branch("Number_of_Matched_W_Jets",  &NumbMatchedWJets   ,"NumbMatchedWJets/I");
  m_tree_stp->Branch("Number_of_Matched_Top_Jets",  &NumbMatchedTopJets   ,"NumbMatchedTopJets/I");
  m_tree_stp->Branch("Sphericity",  &m_Sphericity   ,"Sphericity/F");
  m_tree_stp->Branch("Aplanarity",  &m_Aplanarity   ,"Aplanarity/F");
  m_tree_stp->Branch("Mqq_VBF_DeltaETAcut",  &m_VBF_M   ,"Mqq/F");

  m_tree_stp->Branch("trigger_passed", &m_trigger_passed, "trigger_passed/I");

  m_tree_stp->Branch("Relative_Mass",  &RelMass   ,"RelMass/F");
  m_tree_stp->Branch("PT_Normalized_Mass",  &PTNormalizedMass   ,"PTNormalizedMass/F");
  m_tree_stp->Branch("Mother_PT_Normalized_Mass",  &MotherPTNormalizedMass   ,"MotherPTNormalizedMass/F");
  m_tree_stp->Branch("Number_of_Tops",  &NumberOfTops   ,"NumberOfTops/I");
  m_tree_stp->Branch("Number_of_Loose_and_non_med_b_tags",  &LooseNoMedBtags   ,"LooseNoMedBtags/I");
  m_tree_stp->Branch("Number_of_Loose_and_non_med_b_tags",  &LooseNoMedBtags   ,"LooseNoMedBtags/I");
  m_tree_stp->Branch("DeltaPTwoLeadingJets",  &m_DP2LeadingJets   ,"DP2LeadingJets/F");
  m_tree_stp->Branch("ChiSquaredSorting",  &chi2   ,"chi2/F");

  // Weights and errors from differents scale factoHJ.DeltaR(WJ)<MinDeltaRWH || HJ.DeltaR(WJ)>MaxDeltaRWHrs
  m_tree_stp->Branch("weight", &m_weight, "weight/F");
  m_tree_stp->Branch("weight_error_low", &m_weight_error_low, "weight_error_low/F");
  m_tree_stp->Branch("weight_error_high", &m_weight_error_high, "weight_error_high/F");

  //Proving non-boosted regime
  m_tree_stp->Branch("Njets_for_TpPT_0to100", &Njets_TpPT_0to100, "N_jets0-100/F");
  m_tree_stp->Branch("Njets_for_TpPT_100to200", &Njets_TpPT_100to200, "N_jets100-200/F");
  m_tree_stp->Branch("Njets_for_TpPT_200to300", &Njets_TpPT_200to300, "N_jets200-300/F");
  m_tree_stp->Branch("Njets_for_TpPT_300to400", &Njets_TpPT_300to400, "N_jets300-400/F");
  m_tree_stp->Branch("Njets_for_TpPT_400to500", &Njets_TpPT_400to500, "N_jets400-500/F");
  m_tree_stp->Branch("Njets_for_TpPT_500to600", &Njets_TpPT_500to600, "N_jets500-600/F");
  m_tree_stp->Branch("Njets_for_TpPT_600to700", &Njets_TpPT_600to700, "N_jets600-700/F");
  m_tree_stp->Branch("Njets_for_TpPT_700to800", &Njets_TpPT_700to800, "N_jets700-800/F");

  //Quark content
  m_tree_stp->Branch("U Quark Content",  &UQuarkContent   ,"UQuarkContent/I");
  m_tree_stp->Branch("D Quark Content",  &DQuarkContent   ,"DQuarkContent/I");
  m_tree_stp->Branch("S Quark Content",  &SQuarkContent   ,"SQuarkContent/I");
  m_tree_stp->Branch("C Quark Content",  &CQuarkContent   ,"CQuarkContent/I");
  m_tree_stp->Branch("B Quark Content",  &BQuarkContent   ,"BQuarkContent/I");

  //TTbar System PT
  m_tree_stp->Branch("PT_Difference_ttbar",  &PtTTbarDifference   ,"ttbarPTdif/F");

  //Cuts
  m_tree_cuts->Branch("Trigger_cut", &m_triggercut, "Trigger_passed/I");
  m_tree_cuts->Branch("Cut_0", &m_Cut0, "Cut0_passed/I");
  m_tree_cuts->Branch("Cut_1", &m_Cut1, "Cut1_passed/I");
  m_tree_cuts->Branch("Cut_2", &m_Cut2, "Cut2_passed/I");
  m_tree_cuts->Branch("Cut_3", &m_Cut3, "Cut3_passed/I");
  m_tree_cuts->Branch("Cut_chi2", &m_CutChi2, "CutChi2_passed/I");
  m_tree_cuts->Branch("Cut_4", &m_Cut4, "Cut4_passed/I");
  m_tree_cuts->Branch("Cut_5", &m_Cut5, "Cut5_passed/I");
  m_tree_cuts->Branch("Cut_6", &m_Cut6, "Cut6_passed/I");
  m_tree_cuts->Branch("Cut_7", &m_Cut7, "Cut7_passed/I");
  m_tree_cuts->Branch("Cut_8", &m_Cut8, "Cut8_passed/I");
  m_tree_cuts->Branch("Cut_9", &m_Cut9, "Cut9_passed/I");
  m_tree_cuts->Branch("Cut_10", &m_Cut10, "Cut10_passed/I");
  m_tree_cuts->Branch("Cut_11", &m_Cut11, "Cut11_passed/I");
  m_tree_cuts->Branch("Cut_12", &m_Cut12, "Cut12_passed/I");
  m_tree_cuts->Branch("Cut_13", &m_Cut13, "Cut13_passed/I");
  m_tree_cuts->Branch("Cut_14", &m_Cut14, "Cut14_passed/I");
  m_tree_cuts->Branch("Cut_15", &m_Cut15, "Cut15_passed/I");
  m_tree_cuts->Branch("Cut_16", &m_Cut16, "Cut16_passed/I");
  m_tree_cuts->Branch("Cut_17", &m_Cut17, "Cut17_passed/I");
  m_tree_cuts->Branch("Cut_18", &m_Cut18, "Cut18_passed/I");
  m_tree_cuts->Branch("Cut_19", &m_Cut19, "Cut19_passed/I");
  m_tree_cuts->Branch("Cut_20", &m_Cut20, "Cut20_passed/I");
  m_tree_cuts->Branch("Cut_21", &m_Cut21, "Cut21_passed/I");

  // Initialize the analysis parameters using the ParameterSet iConfig
  //int an_option = iConfig.getUntrackedParameter<int>("an_option", 0);
  m_jet_Ptcut = cmsswSettings.getParameter<edm::ParameterSet>("jets").getParameter<double>("pt_min"); //30;                               // Default val
  m_jet_EtaMaxcut = cmsswSettings.getParameter<edm::ParameterSet>("jets").getParameter<double>("eta_max"); //4.5;
  m_jet_EtaAccepcut = cmsswSettings.getParameter<edm::ParameterSet>("jets").getParameter<double>("eta_accept"); //2.5;
  m_jet_OverlapAccep = cmsswSettings.getParameter<edm::ParameterSet>("jets").getParameter<double>("eta_overlap"); //2.5;
  m_JET_btag_CSVL = cmsswSettings.getParameter<edm::ParameterSet>("jets").getParameter<double>("btag_CSVL"); //0.244;
  m_JET_btag_CSVM = cmsswSettings.getParameter<edm::ParameterSet>("jets").getParameter<double>("btag_CSVM");
  m_JET_btag_CSVT = cmsswSettings.getParameter<edm::ParameterSet>("jets").getParameter<double>("btag_CSVT");
  evt_num = 0;
  m_DRMatching=0.3;
  m_DPtMatching=10.0;

  Cut0=cmsswSettings.getParameter<edm::ParameterSet>("cuts").getParameter<bool>("cut0");
  Cut1=cmsswSettings.getParameter<edm::ParameterSet>("cuts").getParameter<bool>("cut1");
  Cut2=cmsswSettings.getParameter<edm::ParameterSet>("cuts").getParameter<bool>("cut2");
  Cut3=cmsswSettings.getParameter<edm::ParameterSet>("cuts").getParameter<bool>("cut3");
  CutChi2=cmsswSettings.getParameter<edm::ParameterSet>("cuts").getParameter<bool>("cutChi2");
  Cut4=cmsswSettings.getParameter<edm::ParameterSet>("cuts").getParameter<bool>("cut4");
  Cut5=cmsswSettings.getParameter<edm::ParameterSet>("cuts").getParameter<bool>("cut5");
  Cut6=cmsswSettings.getParameter<edm::ParameterSet>("cuts").getParameter<bool>("cut6");
  Cut7=cmsswSettings.getParameter<edm::ParameterSet>("cuts").getParameter<bool>("cut7");
  Cut8=cmsswSettings.getParameter<edm::ParameterSet>("cuts").getParameter<bool>("cut8");
  Cut9=cmsswSettings.getParameter<edm::ParameterSet>("cuts").getParameter<bool>("cut9");
  Cut10=cmsswSettings.getParameter<edm::ParameterSet>("cuts").getParameter<bool>("cut10");
  Cut11=cmsswSettings.getParameter<edm::ParameterSet>("cuts").getParameter<bool>("cut11");
  Cut12=cmsswSettings.getParameter<edm::ParameterSet>("cuts").getParameter<bool>("cut12");
  Cut13=cmsswSettings.getParameter<edm::ParameterSet>("cuts").getParameter<bool>("cut13");
  Cut14=cmsswSettings.getParameter<edm::ParameterSet>("cuts").getParameter<bool>("cut14");
  Cut15=cmsswSettings.getParameter<edm::ParameterSet>("cuts").getParameter<bool>("cut15");
  Cut16=cmsswSettings.getParameter<edm::ParameterSet>("cuts").getParameter<bool>("cut16");
  Cut17=cmsswSettings.getParameter<edm::ParameterSet>("cuts").getParameter<bool>("cut17");
  Cut18=cmsswSettings.getParameter<edm::ParameterSet>("cuts").getParameter<bool>("cut18");
  Cut19=cmsswSettings.getParameter<edm::ParameterSet>("cuts").getParameter<bool>("cut19");
  Cut20=cmsswSettings.getParameter<edm::ParameterSet>("cuts").getParameter<bool>("cut20");
  Cut21=cmsswSettings.getParameter<edm::ParameterSet>("cuts").getParameter<bool>("cut21");

  NumberOfGoodJets=cmsswSettings.getParameter<edm::ParameterSet>("selection").getParameter<double>("NumberOfGoodJets");
  NumberOfBadJets=cmsswSettings.getParameter<edm::ParameterSet>("selection").getParameter<double>("NumberOfBadJets");
  LeadingJetPt=cmsswSettings.getParameter<edm::ParameterSet>("selection").getParameter<double>("LeadingJetPt");
  THTcut=cmsswSettings.getParameter<edm::ParameterSet>("selection").getParameter<double>("THTcut");
  MinB_tags=cmsswSettings.getParameter<edm::ParameterSet>("selection").getParameter<double>("MinB_tags");
  MinLooseB_tags=cmsswSettings.getParameter<edm::ParameterSet>("selection").getParameter<double>("MinLooseB_tags");
  MaxChi2=cmsswSettings.getParameter<edm::ParameterSet>("selection").getParameter<double>("MaxChi2");
  DeltaRHiggsJets=cmsswSettings.getParameter<edm::ParameterSet>("selection").getParameter<double>("DeltaRHiggsJets");
  DeltaRWJets=cmsswSettings.getParameter<edm::ParameterSet>("selection").getParameter<double>("DeltaRWJets");
  HiggsPt=cmsswSettings.getParameter<edm::ParameterSet>("selection").getParameter<double>("HiggsPt");
  TopPt=cmsswSettings.getParameter<edm::ParameterSet>("selection").getParameter<double>("TopPt");
  MinDeltaRWH=cmsswSettings.getParameter<edm::ParameterSet>("selection").getParameter<double>("MinDeltaRWH");
  MaxDeltaRWH=cmsswSettings.getParameter<edm::ParameterSet>("selection").getParameter<double>("MaxDeltaRWH");
  DeltaPhiHiggsJets=cmsswSettings.getParameter<edm::ParameterSet>("selection").getParameter<double>("DeltaPhiHiggsJets");
  DeltaPhiTopJetW=cmsswSettings.getParameter<edm::ParameterSet>("selection").getParameter<double>("DeltaPhiTopJetW");
  JetMultiplicity=cmsswSettings.getParameter<edm::ParameterSet>("selection").getParameter<double>("JetMultiplicity");
  DeltaPhiWjets=cmsswSettings.getParameter<edm::ParameterSet>("selection").getParameter<double>("DeltaPhiWjets");
  MinHiggsMass=cmsswSettings.getParameter<edm::ParameterSet>("selection").getParameter<double>("MinHiggsMass");
  MaxHiggsMass=cmsswSettings.getParameter<edm::ParameterSet>("selection").getParameter<double>("MaxHiggsMass");
  RelHT=cmsswSettings.getParameter<edm::ParameterSet>("selection").getParameter<double>("RelHT");
  Aplanarity=cmsswSettings.getParameter<edm::ParameterSet>("selection").getParameter<double>("Aplanarity");
  MaxDeltaRTH=cmsswSettings.getParameter<edm::ParameterSet>("selection").getParameter<double>("MaxDeltaRTH");
  MinDeltaRTH=cmsswSettings.getParameter<edm::ParameterSet>("selection").getParameter<double>("MinDeltaRTH");
  RelMassMaxCut=cmsswSettings.getParameter<edm::ParameterSet>("selection").getParameter<double>("RelMassMaxCut");
  RelMassMinCut=cmsswSettings.getParameter<edm::ParameterSet>("selection").getParameter<double>("RelMassMinCut");
  MotherPTNormalizedMassCut=cmsswSettings.getParameter<edm::ParameterSet>("selection").getParameter<double>("MotherPTNormalizedMassCut");
  PTNormalizedMassCut=cmsswSettings.getParameter<edm::ParameterSet>("selection").getParameter<double>("PTNormalizedMassCut");
  NumberOfTopsCut=cmsswSettings.getParameter<edm::ParameterSet>("selection").getParameter<double>("NumberOfTopsCut");
  NumberofLooseBtag=cmsswSettings.getParameter<edm::ParameterSet>("selection").getParameter<double>("NumberofLooseBtag");
  DeltaPhiLeadingJets=cmsswSettings.getParameter<edm::ParameterSet>("selection").getParameter<double>("DeltaPhiLeadingJets");
  HMoverTM=cmsswSettings.getParameter<edm::ParameterSet>("selection").getParameter<double>("HMoverTM");

  DoMCMatching=cmsswSettings.getParameter<bool>("DoMatching");

  DoChi2Sorting=cmsswSettings.getParameter<bool>("DoChi2");

  m_weight = 1.;
}

  SingleTprime_analysis::~SingleTprime_analysis(){}

  int SingleTprime_analysis::SingleTprime_Sel() //Main function for the analysis
{
  if (!DoMCMatching)
    {
      if (!m_trigger_passed) return 0;
      m_triggercut=1;
      
      //Correction of trigger//
      TLorentzVector *jetTrigger[6];
      
      if (m_jetMet->getSize()>=6)
	{
	  jetTrigger[0] = m_jetMet->getP4(0);
	  jetTrigger[1] = m_jetMet->getP4(1);
	  jetTrigger[2] = m_jetMet->getP4(2);
	  jetTrigger[3] = m_jetMet->getP4(3);
	  jetTrigger[4] = m_jetMet->getP4(4);
	  jetTrigger[5] = m_jetMet->getP4(5);
	}
      
      if (m_jetMet->getSize()<6 || jetTrigger[0]->Pt()<TriggerDiJet1_2 ||jetTrigger[1]->Pt()<TriggerDiJet1_2 || jetTrigger[2]->Pt()<TriggerDiJet3_4 || jetTrigger[3]->Pt()<TriggerDiJet3_4 || jetTrigger[4]->Pt()<TriggerDiJet5_6 || jetTrigger[5]->Pt()<TriggerDiJet5_6) return 0;  
      /////////////////////////
    }

  int n_jets = m_jetMet->getSize();
  //cout << "Number of jets " << n_jets << endl;
  bool JetsInAcceptance[n_jets]; //Mas for keeping track of jets inside acceptance
  bool jetIsBTagged[n_jets]; //Mask for keeping track of b-tagged jets
  bool jetIsLooseBTagged[n_jets]; //Mask for keeping track of loose b-tagged jets
  int CountingGoodJets=0;
  int CountingBadJets=0;
  float TotalHT=0;
  TLorentzVector PentaJet={0,0,0,0}; int FiveJetCounter=0;
  TLorentzVector LeadingJet={0,0,0,0};
  int BtagCounter=0;
  int LooseBtagCounter=0;

  ScaleFactor jetSF[n_jets];

  TLorentzVector AllJets[n_jets];

  //Loop over jets

  /*int n_MC = m_MC->getSize();

  if (!n_MC) return 0;

  for (int i = 0; i < n_MC ; ++i)
    {
      if (abs(m_MC->getType(i))==0) continue;
      if (abs(m_MC->getType(i)) <= 5 || abs(m_MC->getType(i))==21) 
	{
	  cout << "HERE1" << endl;
	  //if (m_MC->getStatus(i)!=1) continue;
	  TLorentzVector *jetMC= m_MC->getP4(i); 
	  if (jetMC->Pt()==0 || fabs(jetMC->Eta())<=0.01) continue;
	  cout << "HERE2" << endl;
	  cout << jetMC->Pt() << endl;
	  //jeti->SetPxPyPzE(m_MC->getPx(i),m_MC->getPy(i),m_MC->getPz(i),m_MC->getE(i));
	  if (isJetForwSel(jetMC)) CountingBadJets++;
	  if (isJetAccepSel(jetMC)) CountingGoodJets++;
	}
    }
	
  if (CountingGoodJets>=NumberOfGoodJets && CountingBadJets>=NumberOfBadJets) cout << "Good Event" << endl;
  else return 0; */
	  
  for (int i=0;i<n_jets;++i)
    {
      TLorentzVector *jeti = m_jetMet->getP4(i);
      
      //cout << "Pt of jet " << i << " is " << jeti->Pt() << endl;

      TotalHT+=fabs(jeti->Pt());
      AllJets[i].SetPxPyPzE(jeti->Px(),jeti->Py(),jeti->Pz(),jeti->E());
      if (m_isMC) jetSF[i] = m_jetMet->getScaleFactor(i);
      //Finding B-tagged jets with loose requirement based on CSV algorithm                              
      /*if ((m_jetMet->getJetBTagProb_CSV(i)) > m_JET_btag_CSVL)
        {
          ++m_NBtaggedJets_CSVL;
          jetIsBTagged[i] = true;
	  ++BtagCounter;
        }
      else
        {
          jetIsBTagged[i] = false;
	  if ((m_jetMet->getJetBTagProb_CSV(i)) > 0.244)
		{
		  jetIsLooseBTagged[i] = true;
		  ++LooseBtagCounter;
		}
	      else
		{
		  jetIsLooseBTagged[i] = false;
		}
		}*/
      
      //Third Loose B-tag
      if ((m_jetMet->getJetBTagProb_CSV(i)) < m_JET_btag_CSVL && (m_jetMet->getJetBTagProb_CSV(i)) > 0.244) ++LooseNoMedBtags;

      //if (!isJetSel(jeti)) continue; // apply the pt cut
      if (isJetForwSel(jeti)/* && !JetsInAcceptance[i]*/) CountingBadJets++;
      jetIsBTagged[i] = false; jetIsLooseBTagged[i] = false; JetsInAcceptance[i]=false;
      if (isJetAccepSel(jeti)) 
	{
	  JetsInAcceptance[i]=true; CountingGoodJets++; if (FiveJetCounter==0) {LeadingJet=AllJets[i];}; if (FiveJetCounter<5) {PentaJet+=AllJets[i]; FiveJetCounter++;}
	  if ((m_jetMet->getJetBTagProb_CSV(i)) > m_JET_btag_CSVT)
	    {
	      ++m_NBtaggedJets_CSVT;
	    }
	  if ((m_jetMet->getJetBTagProb_CSV(i)) > m_JET_btag_CSVM)
	    {
	      ++m_NBtaggedJets_CSVM;
	      jetIsBTagged[i] = true;
	      ++BtagCounter;
	      jetIsLooseBTagged[i] = true;
	      ++LooseBtagCounter;
	    }
	  else
	    {
	      jetIsBTagged[i] = false;
	      if ((m_jetMet->getJetBTagProb_CSV(i)) > 0.244)
		{
		  ++m_NBtaggedJets_CSVL;
		  jetIsLooseBTagged[i] = true;
		  ++LooseBtagCounter;
		}
	      else
		{
		  jetIsLooseBTagged[i] = false;
		}
	    }	  
	}
      //if (isJetAccepSel(jeti)) cout << "The pt of InJet " << i << " is " << jeti->Pt() << endl;
      else {JetsInAcceptance[i]=false;}
      //if (isJetForwSel(jeti) && !JetsInAcceptance[i]) cout << "The pt of OutJet " << i << " is " << jeti->Pt() << endl;
    }

  //To study VBF cut
  if (true)
    {
      float DeltaEta=0;
      int VBFcounter=0;
      for (int i=0;i<n_jets;++i)
	{
	  TLorentzVector jeti; jeti.SetPxPyPzE(AllJets[i].Px(),AllJets[i].Py(),AllJets[i].Pz(),AllJets[i].E());
	  for (int j=i+1;j<n_jets;++j)
	    {
	      TLorentzVector jetj; jetj.SetPxPyPzE(AllJets[j].Px(),AllJets[j].Py(),AllJets[j].Pz(),AllJets[j].E());
	      TLorentzVector DijetVBF=jeti+jetj;
	      if (fabs(jeti.Eta()-jetj.Eta())>2.2)
		{
		  ++VBFcounter; DeltaEta=fabs(jeti.Eta()-jetj.Eta()); m_VBF_M=DijetVBF.M();
		  if (VBFcounter!=0 && DeltaEta<fabs(jeti.Eta()-jetj.Eta())) {DeltaEta=fabs(jeti.Eta()-jetj.Eta()); m_VBF_M=DijetVBF.M();}
		}
	    }
	}
    }

  // First check the number of jets in and out the acceptance
  //cout << CountingGoodJets << " " << CountingBadJets << endl;

  /////////
  //Cut 0//
  /////////

  if(Cut0) {
  if (CountingGoodJets>=NumberOfGoodJets && CountingBadJets>=NumberOfBadJets) cout << "Good Event" << endl;
  else return 0; }
  m_Cut0=1;
  
  /////////
  //Cut 1//
  /////////

  if(Cut1) {if (LeadingJet.Pt()<LeadingJetPt) return 0;}
  m_Cut1=1;

  //cout << "The HT of the event is " << TotalHT << endl;
  m_THT = TotalHT;
  if (n_jets>0) m_jet1pt = AllJets[0].Pt();
  if (n_jets>1) m_jet2pt = AllJets[1].Pt();
  if (n_jets>2) m_jet3pt = AllJets[2].Pt();
  if (n_jets>3) m_jet4pt = AllJets[3].Pt();
  if (n_jets>4) m_jet5pt = AllJets[4].Pt();
  if (n_jets>5) m_jet6pt = AllJets[5].Pt();

  m_DP2LeadingJets=fabs(AllJets[1].Phi()-AllJets[0].Phi());

  /////////
  //Cut 2//
  /////////

  if (Cut2) {if (TotalHT<THTcut) return 0;}
  m_Cut2=1;

  /////////
  //Cut 3//
  /////////

  cout << "Exclusive loose b-tags " << LooseBtagCounter-BtagCounter << endl;
  if (MinLooseB_tags!=0)
    {
      if (Cut3) {if (BtagCounter<MinB_tags || (LooseBtagCounter-BtagCounter)<MinLooseB_tags) return 0;}
    }
  else
    {
      if (Cut3) {if (BtagCounter<MinB_tags) return 0;}
    }

  m_Cut3=1;


  //////////////////////////
  //Chi2 Sorting algorithm//
  //////////////////////////

  TLorentzVector FHJ;				
  TLorentzVector SHJ; 
  TLorentzVector FWJ; 
  TLorentzVector SWJ; 
  int IndexHiggsJets[2]={0,0};
  int IndexWJets[2]={0,0};
  int IndexTopJet=0;  

  int B_tag_counter=0;
  int Bsindexes[(BtagCounter-1)*BtagCounter/2][2];
  float Hchi2[(BtagCounter-1)*BtagCounter/2];
  float Tchi2=0;
  float Wchi2=0;
  float CurrentChi2=-1;
  int IndexesFromchi2[5]={-1,-1,-1,-1,-1}; //Higgs, W, Top

  if (DoChi2Sorting)
    {
      for (int i=0;i<n_jets;++i)
	{
	  if (!JetsInAcceptance[i]) continue;
	  for (int j=i+1;j<n_jets;++j)
	    {
	      if (!JetsInAcceptance[j]) continue; 
	      if (!jetIsBTagged[i]) continue; 
	      if (!jetIsBTagged[j]) continue; //if (B_tag_counter==2) break; B_tag_counter++;
	      TLorentzVector jeti; jeti.SetPxPyPzE(AllJets[i].Px(),AllJets[i].Py(),AllJets[i].Pz(),AllJets[i].E());
	      TLorentzVector jetj; jetj.SetPxPyPzE(AllJets[j].Px(),AllJets[j].Py(),AllJets[j].Pz(),AllJets[j].E());
	      TLorentzVector DiBjet = jeti+jetj;
	      Hchi2[B_tag_counter]=0;
	      Hchi2[B_tag_counter]=(DiBjet.M()-HMassChi2)*(DiBjet.M()-HMassChi2)/((HBBWidth)*(HBBWidth));
	      Bsindexes[B_tag_counter][0]=0; Bsindexes[B_tag_counter][1]=0;
	      Bsindexes[B_tag_counter][0]=i; Bsindexes[B_tag_counter][1]=j; //cout << "Entries B_tags " << i << " " << j << endl;
	      B_tag_counter++;
	    }
	}
      
      for (int l=0;l<(BtagCounter-1)*BtagCounter/2;++l)
	{
	  for (int i=0;i<n_jets;++i)
	    {
	      if (!JetsInAcceptance[i]) continue;
	      for (int j=i+1;j<n_jets;++j)
		{
		  if (!JetsInAcceptance[j]) continue;
		  for (int k=0;k<n_jets;++k)
		    {
		      if (!JetsInAcceptance[k] || k==i || k==j) continue;
		      //cout << "FIRST CHECK " << i << j << k << endl;
		      if (i==Bsindexes[l][0] || i==Bsindexes[l][1] || j==Bsindexes[l][0] || j==Bsindexes[l][1] || k==Bsindexes[l][0] || k==Bsindexes[l][1]) continue;
		      //cout << "INDEXES!!!!!!!!! " << Bsindexes[l][0] << Bsindexes[l][1] << i << j << k << endl;		  
		      TLorentzVector jeti; jeti.SetPxPyPzE(AllJets[i].Px(),AllJets[i].Py(),AllJets[i].Pz(),AllJets[i].E());
		      TLorentzVector jetj; jetj.SetPxPyPzE(AllJets[j].Px(),AllJets[j].Py(),AllJets[j].Pz(),AllJets[j].E());
		      TLorentzVector Dijet = jeti+jetj;
		      Wchi2=(Dijet.M()-WMassChi2)*(Dijet.M()-WMassChi2)/((WHadrWidth)*(WHadrWidth));
		      TLorentzVector jetk; jetk.SetPxPyPzE(AllJets[k].Px(),AllJets[k].Py(),AllJets[k].Pz(),AllJets[k].E());
		      TLorentzVector Trijet = Dijet+jetk;
		      Tchi2=(Trijet.M()-TopMassChi2)*(Trijet.M()-TopMassChi2)/((TopHadrWidth)*(TopHadrWidth));
		      //cout << "Temporal chi2 " << Hchi2[l]+Wchi2+Tchi2 << endl;
		      if (CurrentChi2==-1) {CurrentChi2=Hchi2[l]+Wchi2+Tchi2; IndexesFromchi2[0]=Bsindexes[l][0]; IndexesFromchi2[1]=Bsindexes[l][1]; IndexesFromchi2[2]=i; IndexesFromchi2[3]=j; IndexesFromchi2[4]=k;}
		      else
			{
			  if (CurrentChi2>(Hchi2[l]+Wchi2+Tchi2)) {CurrentChi2=Hchi2[l]+Wchi2+Tchi2; IndexesFromchi2[0]=Bsindexes[l][0]; IndexesFromchi2[1]=Bsindexes[l][1]; IndexesFromchi2[2]=i; IndexesFromchi2[3]=j; IndexesFromchi2[4]=k;}
			}
		    }
		}
	    }
	  //cout << "Current Chi2 " << CurrentChi2 << endl;
	  //cout << "Current Indexes " << IndexesFromchi2[0] << " " << IndexesFromchi2[1] << " " << IndexesFromchi2[2] << " " << IndexesFromchi2[3] << " " << IndexesFromchi2[4] << endl;
	}
      
      cout << "Minimal chi2 for this events: " << CurrentChi2 << endl;

      //Reconstructing particles from the algorithm result
      FirstHiggsJet->SetPxPyPzE(AllJets[IndexesFromchi2[0]].Px(),AllJets[IndexesFromchi2[0]].Py(),AllJets[IndexesFromchi2[0]].Pz(),AllJets[IndexesFromchi2[0]].E()); SecondHiggsJet->SetPxPyPzE(AllJets[IndexesFromchi2[1]].Px(),AllJets[IndexesFromchi2[1]].Py(),AllJets[IndexesFromchi2[1]].Pz(),AllJets[IndexesFromchi2[1]].E());
      ReconstructedHiggs->SetPxPyPzE(AllJets[IndexesFromchi2[0]].Px()+AllJets[IndexesFromchi2[1]].Px(),AllJets[IndexesFromchi2[0]].Py()+AllJets[IndexesFromchi2[1]].Py(),AllJets[IndexesFromchi2[0]].Pz()+AllJets[IndexesFromchi2[1]].Pz(),AllJets[IndexesFromchi2[0]].E()+AllJets[IndexesFromchi2[1]].E());
      FHJ.SetPxPyPzE(FirstHiggsJet->Px(), FirstHiggsJet->Py(), FirstHiggsJet->Pz(), FirstHiggsJet->E());
      SHJ.SetPxPyPzE(SecondHiggsJet->Px(), SecondHiggsJet->Py(), SecondHiggsJet->Pz(), SecondHiggsJet->E());
      m_DRHiggsJets=FHJ.DeltaR(SHJ);
      IndexHiggsJets[0]=IndexesFromchi2[0]; IndexHiggsJets[1]=IndexesFromchi2[1]; 
      //
      FirstWJet->SetPxPyPzE(AllJets[IndexesFromchi2[2]].Px(),AllJets[IndexesFromchi2[2]].Py(),AllJets[IndexesFromchi2[2]].Pz(),AllJets[IndexesFromchi2[2]].E()); SecondWJet->SetPxPyPzE(AllJets[IndexesFromchi2[3]].Px(),AllJets[IndexesFromchi2[3]].Py(),AllJets[IndexesFromchi2[3]].Pz(),AllJets[IndexesFromchi2[3]].E());
      ReconstructedW->SetPxPyPzE(AllJets[IndexesFromchi2[2]].Px()+AllJets[IndexesFromchi2[3]].Px(),AllJets[IndexesFromchi2[2]].Py()+AllJets[IndexesFromchi2[3]].Py(),AllJets[IndexesFromchi2[2]].Pz()+AllJets[IndexesFromchi2[3]].Pz(),AllJets[IndexesFromchi2[2]].E()+AllJets[IndexesFromchi2[3]].E());
      FWJ.SetPxPyPzE(FirstWJet->Px(), FirstWJet->Py(), FirstWJet->Pz(), FirstWJet->E());
      SWJ.SetPxPyPzE(SecondWJet->Px(), SecondWJet->Py(), SecondWJet->Pz(), SecondWJet->E());
      m_DRWJets=FWJ.DeltaR(SWJ);
      IndexWJets[0]=IndexesFromchi2[2]; IndexWJets[1]=IndexesFromchi2[3]; 
      //
      ReconstructedTop->SetPxPyPzE(AllJets[IndexesFromchi2[2]].Px()+AllJets[IndexesFromchi2[3]].Px()+AllJets[IndexesFromchi2[4]].Px(),AllJets[IndexesFromchi2[2]].Py()+AllJets[IndexesFromchi2[3]].Py()+AllJets[IndexesFromchi2[4]].Py(),AllJets[IndexesFromchi2[2]].Pz()+AllJets[IndexesFromchi2[3]].Pz()+AllJets[IndexesFromchi2[4]].Pz(),AllJets[IndexesFromchi2[2]].E()+AllJets[IndexesFromchi2[3]].E()+AllJets[IndexesFromchi2[4]].E());
      TopJet->SetPxPyPzE(AllJets[IndexesFromchi2[4]].Px(),AllJets[IndexesFromchi2[4]].Py(),AllJets[IndexesFromchi2[4]].Pz(),AllJets[IndexesFromchi2[4]].E());
      IndexTopJet=IndexesFromchi2[4];

      if (CutChi2) {if (CurrentChi2>MaxChi2) return 0;}
      m_CutChi2=1;

      if (Cut4) {if (FHJ.DeltaR(SHJ)>DeltaRHiggsJets) return 0;}
      m_Cut4=1;

      if (Cut5) {if (FWJ.DeltaR(SWJ)>DeltaRWJets) return 0;}
      m_Cut5=1;

    }
  
  chi2=CurrentChi2;
  //////////////////////////

  ///////////////////////////////////////////////
  //Reconstructing the Higgs from b-tagged jets//
  ///////////////////////////////////////////////
  //cout << "Entering Higgs reconstruction" << endl;
  bool EventWithHiggs=false;
  int HiggsCounter=0;
  float DiffWithHiggMass=0;
  if (!DoChi2Sorting)
    {
  for (int i=0;i<n_jets;++i)
    {
      if (MinLooseB_tags!=0) {if (!jetIsLooseBTagged[i] || !JetsInAcceptance[i]) continue;}
      else {if (!jetIsBTagged[i] || !JetsInAcceptance[i]) continue;}
      TLorentzVector jeti; jeti.SetPxPyPzE(AllJets[i].Px(),AllJets[i].Py(),AllJets[i].Pz(),AllJets[i].E());
      for (int j=i+1;j<n_jets;++j)
	{
	  //if (!jetIsBTagged[j] || !JetsInAcceptance[j]) continue;
	  if (MinLooseB_tags!=0) {if (!jetIsLooseBTagged[j] || !JetsInAcceptance[j]) continue;}
	  else {if (!jetIsBTagged[j] || !JetsInAcceptance[j]) continue;}
	  TLorentzVector jetj; jetj.SetPxPyPzE(AllJets[j].Px(),AllJets[j].Py(),AllJets[j].Pz(),AllJets[j].E());
	  TLorentzVector DiBjet = jeti+jetj;
	  //DiBjet.SetPxPyPzE(jeti.Px()+jetj.Px(),jeti.Py()+jetj.Py(),jeti.Pz()+jetj.Pz(),jeti.E()+jetj.E());
	  if (fabs(DiBjet.M()-HiggsMass)>HiggsMassWindow) continue;
	  if (jeti.DeltaR(jetj)>DeltaRHiggsJets) continue;
	  //cout << "Mass of jet i: " << jeti.M() << " and j: " << jetj.M() << endl;
	  //cout << "Mass of the b-tagged jets couple " << i << j << " is " << DiBjet.M() << endl;
	  EventWithHiggs=true;
	  if (HiggsCounter==0)
	    {
	      ++HiggsCounter; DiffWithHiggMass=fabs(DiBjet.M()-HiggsMass);
	      IndexHiggsJets[0]=i; IndexHiggsJets[1]=j; ReconstructedHiggs->SetPxPyPzE(jeti.Px()+jetj.Px(),jeti.Py()+jetj.Py(),jeti.Pz()+jetj.Pz(),jeti.E()+jetj.E());
	      FirstHiggsJet->SetPxPyPzE(jeti.Px(),jeti.Py(),jeti.Pz(),jeti.E()); SecondHiggsJet->SetPxPyPzE(jetj.Px(),jetj.Py(),jetj.Pz(),jetj.E());
	    }
	  else
	    {
	      ++HiggsCounter;
	      if (DiffWithHiggMass>fabs(DiBjet.M()-HiggsMass)) 
		{
		  IndexHiggsJets[0]=i; IndexHiggsJets[1]=j; 
		  DiffWithHiggMass=fabs(DiBjet.M()-HiggsMass); ReconstructedHiggs->SetPxPyPzE(jeti.Px()+jetj.Px(),jeti.Py()+jetj.Py(),jeti.Pz()+jetj.Pz(),jeti.E()+jetj.E()); 
		  FirstHiggsJet->SetPxPyPzE(jeti.Px(),jeti.Py(),jeti.Pz(),jeti.E()); SecondHiggsJet->SetPxPyPzE(jetj.Px(),jetj.Py(),jetj.Pz(),jetj.E());
		}
	    }
	}
    }

  //cout << "Mass of selected Higgs is " << ReconstructedHiggs->M() << " from couple " << IndexHiggsJets[0] << IndexHiggsJets[1] << endl;
  //TLorentzVector BJetCouple;  BJetCouple.SetPxPyPzE(FirstHiggsJet->Px()+SecondHiggsJet->Px(),FirstHiggsJet->Py()+SecondHiggsJet->Py(),FirstHiggsJet->Pz()+SecondHiggsJet->Pz(),FirstHiggsJet->E()+SecondHiggsJet->E());
  //cout << "Mass of selected couple of b jets is " << BJetCouple.M() << endl;

  if (Cut4) {if (!EventWithHiggs) return 0;}
  m_Cut4=1;
  FHJ.SetPxPyPzE(FirstHiggsJet->Px(), FirstHiggsJet->Py(), FirstHiggsJet->Pz(), FirstHiggsJet->E());
  SHJ.SetPxPyPzE(SecondHiggsJet->Px(), SecondHiggsJet->Py(), SecondHiggsJet->Pz(), SecondHiggsJet->E());
  m_DRHiggsJets=FHJ.DeltaR(SHJ);
    }

  //////////////////////////////////////////////////////////////////
  //Reconstruction of W from full jets collection minus Higgs jets//
  //////////////////////////////////////////////////////////////////
  //cout << "Entering W reconstruction" << endl;
  bool EventWithW=false;
  int WCounter=0;
  float DiffWithWMass=0;  
  if (!DoChi2Sorting)
    {
  for (int i=0;i<n_jets;++i)
    {
      if (IndexHiggsJets[0]==i || IndexHiggsJets[1]==i || !JetsInAcceptance[i]) continue;
      TLorentzVector jeti; jeti.SetPxPyPzE(AllJets[i].Px(),AllJets[i].Py(),AllJets[i].Pz(),AllJets[i].E());
      for (int j=i+1;j<n_jets;++j)
	{
	  if (IndexHiggsJets[0]==j || IndexHiggsJets[1]==j || !JetsInAcceptance[j]) continue;
	  TLorentzVector jetj; jetj.SetPxPyPzE(AllJets[j].Px(),AllJets[j].Py(),AllJets[j].Pz(),AllJets[j].E());
	  TLorentzVector Dijet = jeti+jetj;
	  if (fabs(Dijet.M()-WMass)>WMassWindow) continue;
	  if (jeti.DeltaR(jetj)>DeltaRWJets) continue;
	  //cout << "Mass of jet i: " << jeti.M() << " and j: " << jetj.M() << endl;
	  //cout << "Mass of the b-tagged jets couple " << i << j << " is " << DiBjet.M() << endl;
	  EventWithW=true;
	  if (WCounter==0)
	    {
	      ++WCounter; DiffWithWMass=fabs(Dijet.M()-WMass);
	      IndexWJets[0]=i; IndexWJets[1]=j; ReconstructedW->SetPxPyPzE(jeti.Px()+jetj.Px(),jeti.Py()+jetj.Py(),jeti.Pz()+jetj.Pz(),jeti.E()+jetj.E());
	      FirstWJet->SetPxPyPzE(jeti.Px(),jeti.Py(),jeti.Pz(),jeti.E()); SecondWJet->SetPxPyPzE(jetj.Px(),jetj.Py(),jetj.Pz(),jetj.E());
	    }
	  else
	    {
	      ++WCounter;
	      if (DiffWithWMass>fabs(Dijet.M()-WMass)) 
		{
		  IndexWJets[0]=i; IndexWJets[1]=j; 
		  DiffWithWMass=fabs(Dijet.M()-WMass); ReconstructedW->SetPxPyPzE(jeti.Px()+jetj.Px(),jeti.Py()+jetj.Py(),jeti.Pz()+jetj.Pz(),jeti.E()+jetj.E()); 
		  FirstWJet->SetPxPyPzE(jeti.Px(),jeti.Py(),jeti.Pz(),jeti.E()); SecondWJet->SetPxPyPzE(jetj.Px(),jetj.Py(),jetj.Pz(),jetj.E());
		}
	    }
	}
    }

  if (Cut5) {if (!EventWithW) return 0;}
  m_Cut5=1;
  //cout << "Higgs Jets are: " << IndexHiggsJets[0] << IndexHiggsJets[1] << " W Jets are: " << IndexWJets[0] << IndexWJets[1] << endl;
  FWJ.SetPxPyPzE(FirstWJet->Px(), FirstWJet->Py(), FirstWJet->Pz(), FirstWJet->E());
  SWJ.SetPxPyPzE(SecondWJet->Px(), SecondWJet->Py(), SecondWJet->Pz(), SecondWJet->E());
  m_DRWJets=FWJ.DeltaR(SWJ);  
    }

  /////////////////////////
  //Reconstruction of Top//
  /////////////////////////
  //cout << "Entering Top reconstruction" << endl;
  bool EventWithTop=false;
  int TopCounter=0;
  int RestrictedTopCounter=0;
  float DiffWithTopMass=0; 
  if (!DoChi2Sorting)
    {
  for (int i=0;i<n_jets;++i)
    {
      if (IndexHiggsJets[0]==i || IndexHiggsJets[1]==i || IndexWJets[0]==i || IndexWJets[1]==i || !JetsInAcceptance[i]) continue;
      TLorentzVector jeti; jeti.SetPxPyPzE(AllJets[i].Px(),AllJets[i].Py(),AllJets[i].Pz(),AllJets[i].E());
      TLorentzVector wjet; wjet.SetPxPyPzE(ReconstructedW->Px(),ReconstructedW->Py(),ReconstructedW->Pz(),ReconstructedW->E());
      TLorentzVector Trijet = jeti+wjet;
      if (fabs(Trijet.M()-TopMass)>TopMassWindow) continue;
      //cout << "Mass of jet i: " << jeti.M() << " and j: " << jetj.M() << endl;
      //cout << "Mass of the b-tagged jets couple " << i << j << " is " << DiBjet.M() << endl;
      if (fabs(Trijet.M()-TopMass)>30) RestrictedTopCounter++;
      EventWithTop=true;
      if (TopCounter==0)
	{
	  ++TopCounter; DiffWithTopMass=fabs(Trijet.M()-TopMass);
	  IndexTopJet=i; ReconstructedTop->SetPxPyPzE(jeti.Px()+wjet.Px(),jeti.Py()+wjet.Py(),jeti.Pz()+wjet.Pz(),jeti.E()+wjet.E());
	  TopJet->SetPxPyPzE(jeti.Px(),jeti.Py(),jeti.Pz(),jeti.E());
	}
      else
	{
	  ++TopCounter;
	  if (DiffWithTopMass>fabs(Trijet.M()-TopMass)) 
	    {
		  IndexTopJet=i;
		  DiffWithTopMass=fabs(Trijet.M()-TopMass); ReconstructedTop->SetPxPyPzE(jeti.Px()+wjet.Px(),jeti.Py()+wjet.Py(),jeti.Pz()+wjet.Pz(),jeti.E()+wjet.E()); 
		  TopJet->SetPxPyPzE(jeti.Px(),jeti.Py(),jeti.Pz(),jeti.E());
	    }
	}
    }

  if (!EventWithTop) return 0;
    }
  //cout << "Higgs Jets are: " << IndexHiggsJets[0] << IndexHiggsJets[1] << " W Jets are: " << IndexWJets[0] << IndexWJets[1] << " Top jet is: " << IndexTopJet << endl;
  TLorentzVector HJ; HJ.SetPxPyPzE(ReconstructedHiggs->Px(), ReconstructedHiggs->Py(), ReconstructedHiggs->Pz(), ReconstructedHiggs->E());
  TLorentzVector TJ; TJ.SetPxPyPzE(ReconstructedTop->Px(), ReconstructedTop->Py(), ReconstructedTop->Pz(), ReconstructedTop->E());
  TLorentzVector WJ; WJ.SetPxPyPzE(ReconstructedW->Px(), ReconstructedW->Py(), ReconstructedW->Pz(), ReconstructedW->E());
  m_DRTopHiggs=HJ.DeltaR(TJ);
  m_DRWHiggs=HJ.DeltaR(WJ);
  m_RelTHT=(ReconstructedHiggs->Pt()+ReconstructedTop->Pt())/TotalHT; //Relative total hadronic energy

  ReconstructedTprime->SetPxPyPzE(ReconstructedHiggs->Px()+ReconstructedTop->Px(),ReconstructedHiggs->Py()+ReconstructedTop->Py(),ReconstructedHiggs->Pz()+ReconstructedTop->Pz(),ReconstructedHiggs->E()+ReconstructedTop->E());

  RelMass=(HJ.M()+TJ.M())/ReconstructedTprime->M();
  PTNormalizedMass=(HJ.M()+TJ.M())/(HJ.Pt()+TJ.Pt());
  MotherPTNormalizedMass=(HJ.M()+TJ.M())/ReconstructedTprime->Pt();
  NumberOfTops=RestrictedTopCounter;

  /////////
  //Cut 6//
  /////////
  
  if (Cut6) {if (HJ.Pt()<HiggsPt || TJ.Pt()<TopPt) return 0;}
  m_Cut6=1;

  /////////
  //Cut 7//
  /////////

  if (Cut7) {if (HJ.DeltaR(WJ)<MinDeltaRWH || HJ.DeltaR(WJ)>MaxDeltaRWH) return 0;}
  m_Cut7=1;

  /////////
  //Cut 8//
  /////////

  if (Cut8) {if (fabs(FHJ.Phi()-SHJ.Phi())>DeltaPhiHiggsJets || fabs(TopJet->Phi()-WJ.Phi())>DeltaPhiTopJetW) return 0;}
  m_Cut8=1;

  /////////
  //Cut 9//
  /////////

  if (Cut9) {if (CountingBadJets>JetMultiplicity) return 0;}
  m_Cut9=1;

  //////////
  //Cut 10//
  //////////

  if (Cut10) {if (fabs(FHJ.Phi()-SHJ.Phi())>DeltaPhiHiggsJets || fabs(FWJ.Phi()-SWJ.Phi())>DeltaPhiWjets) return 0;}
  m_Cut10=1;

  //////////
  //Cut 11//
  //////////

  if (Cut11) {if (HJ.M()>MaxHiggsMass || HJ.M()<MinHiggsMass) return 0;}
  m_Cut11=1;

  //////////
  //Cut 12//
  //////////

  if (Cut12) {if (m_RelTHT<RelHT) return 0;}
  m_Cut12=1;

  /////////////////////////////////////////////
  //Extracting info for Sphericity/Aplanarity//
  /////////////////////////////////////////////
  TLorentzVector ForDecayRestFrame;
  for (int i=0;i<n_jets;++i)
    {
      TLorentzVector *jeti = m_jetMet->getP4(i);
      TLorentzVector temp;
      if (i==0) ForDecayRestFrame.SetPxPyPzE(jeti->Px(),jeti->Py(),jeti->Pz(),jeti->E());
      else {temp=ForDecayRestFrame; ForDecayRestFrame.SetPxPyPzE(jeti->Px()+temp.Px(),jeti->Py()+temp.Py(),jeti->Pz()+temp.Pz(),jeti->E()+temp.E());}
    }

  double Sxx=0;
  double Sxy=0;
  double Sxz=0;
  double Syy=0;
  double Syz=0;
  double Szz=0;
  double Sxxnum=0;
  double Sxynum=0;
  double Sxznum=0;
  double Syynum=0;
  double Syznum=0;
  double Szznum=0;
  double Sden=0;
    
  for (int i=0;i<n_jets;++i)
    {
      TLorentzVector *jeti = m_jetMet->getP4(i);
      Sxxnum+=jeti->Px()*jeti->Px();
      Sxynum+=jeti->Px()*jeti->Py();
      Sxznum+=jeti->Px()*jeti->Pz();
      Syynum+=jeti->Py()*jeti->Py();
      Syznum+=jeti->Py()*jeti->Pz();
      Szznum+=jeti->Pz()*jeti->Pz();
      Sden+=jeti->P()*jeti->P();
    }
  
  Sxx=Sxxnum/Sden;
  Sxy=Sxynum/Sden;
  Sxz=Sxznum/Sden;
  Syy=Syynum/Sden;
  Syz=Syznum/Sden;
  Szz=Szznum/Sden;
  const int N = 3;
  double s[N*N] = {
    Sxx, Sxy, Sxz,
    Sxy, Syy, Syz,
    Sxz, Syz, Szz
  };
  TMatrixDSym S(N, s);
  TMatrixDSymEigen SEigen(S);
  TVectorD eigenval = SEigen.GetEigenValues();
  
  m_Sphericity=(eigenval[1]+eigenval[2])*3./2.;
  m_Aplanarity=eigenval[2]*3./2.;

  //////////
  //Cut 13//
  //////////

  if (Cut13) {if (m_Aplanarity>Aplanarity) return 0;}
  m_Cut13=1;

  //////////
  //Cut 14//
  //////////

  if (Cut14) {if (HJ.DeltaR(TJ)<MinDeltaRTH || HJ.DeltaR(TJ)>MaxDeltaRTH) return 0;}
  m_Cut14=1;

  //////////
  //Cut 15//
  //////////

  if (Cut15) {if (RelMass>RelMassMaxCut || RelMass<RelMassMinCut) return 0;}
  m_Cut15=1;

  //////////
  //Cut 16//
  //////////

  if (Cut16) {if (PTNormalizedMass>PTNormalizedMassCut) return 0;}
  m_Cut16=1;

  //////////
  //Cut 17//
  //////////

  if (Cut17) {if (MotherPTNormalizedMass>MotherPTNormalizedMassCut) return 0;}
  m_Cut17=1;

  //////////
  //Cut 18//
  //////////

  if (Cut18) {if (NumberOfTops>=NumberOfTopsCut) return 0;}
  m_Cut18=1;

  //////////
  //Cut 19//
  //////////

  if (Cut19) {if (LooseNoMedBtags>NumberofLooseBtag) return 0;}
  m_Cut19=1;

  //////////
  //Cut 20//
  //////////

  if (Cut20) {if (m_DP2LeadingJets>DeltaPhiLeadingJets) return 0;}
  m_Cut20=1;
  
  //////////
  //Cut 21//
  //////////
  
  if (Cut21) {if (HJ.M()/TJ.M()<HMoverTM) return 0;}
  m_Cut21=1;

  ///////////////////
  //TTbar system PT//
  ///////////////////
  int NewTopsCounter=0;
  int NumberOfBs=0; if (MinLooseB_tags!=0) {NumberOfBs=BtagCounter;} else {NumberOfBs=LooseBtagCounter;}
  TLorentzVector NewTops[NumberOfBs];
  for (int i=0;i<n_jets;++i)
    {
      if (MinLooseB_tags!=0) {if (!jetIsLooseBTagged[i] || !JetsInAcceptance[i]) continue;}
      else {if (!jetIsBTagged[i] || !JetsInAcceptance[i]) continue;}
      TLorentzVector jeti; jeti.SetPxPyPzE(AllJets[i].Px(),AllJets[i].Py(),AllJets[i].Pz(),AllJets[i].E());
      for (int j=i+1;j<n_jets;++j)
	{
	  if (!JetsInAcceptance[j]) continue;
	  TLorentzVector jetj; jetj.SetPxPyPzE(AllJets[j].Px(),AllJets[j].Py(),AllJets[j].Pz(),AllJets[j].E());
	  for (int k=j+1;k<n_jets;++k)
	    {
	      if (!JetsInAcceptance[k]) continue;
	      TLorentzVector jetk; jetk.SetPxPyPzE(AllJets[k].Px(),AllJets[k].Py(),AllJets[k].Pz(),AllJets[k].E());
	      TLorentzVector Trijet = jeti+jetj+jetk;
	      if (fabs(Trijet.M()-TopMass)>20) continue;
	      if (NewTopsCounter<=NumberOfBs) {NewTops[NewTopsCounter].SetPxPyPzE(Trijet.Px(),Trijet.Py(),Trijet.Pz(),Trijet.E()); NewTopsCounter++;}
	      else break;

	    }
	}
    }

  cout << NewTopsCounter << " " << NumberOfBs << endl;

  float MinimalPtTTbar=0;
  bool MarkerToDoCut=false;

  if (NewTopsCounter>=2)
    {
      for (int i=0;i<NewTopsCounter;++i)
	{
	  if (i==0) MinimalPtTTbar=abs(NewTops[i].Pt()-NewTops[i+1].Pt());
	  else if (i<NewTopsCounter-1)
	    {
	      if (MinimalPtTTbar>abs(NewTops[i].Pt()-NewTops[i+1].Pt())) MinimalPtTTbar=abs(NewTops[i].Pt()-NewTops[i+1].Pt());
	      cout << i << endl;
	    }
	}   
      MarkerToDoCut=true;
    }

  if (MarkerToDoCut) PtTTbarDifference=MinimalPtTTbar;

  ///////////////////////////
  //Comparing with MC truth//
  ///////////////////////////

  if (DoMCMatching) //Switch to true to obtain comparison of our method with MC matching
    {
      int MatchedHiggsJets=0;
      int MatchedFirstHiggsJets=0;
      int MatchedSecondHiggsJets=0;
      int MatchedWJets=0;
      int MatchedFirstWJets=0;
      int MatchedSecondWJets=0;
      int MatchedTopJets=0;
      int MatchedJetsIndexes[5]={-1,-1,-1,-1,-1}; //First Higgs Jet, Second Higgs Jet, First W Jet, Second W Jet, Top Jet
      
      //int n_jets = m_jetMet->getSize();
      for (int i=0;i<n_jets;++i)
	{
	  int jetMatch = m_jetMet->getJetMCIndex(i);
	  if (jetMatch!=-10) cout << "Jet " << i << " has index of MC particle matched " << jetMatch << " which has pdg " << m_MC->getType(jetMatch) << " mother pdg index " << m_MC->getType(patIndexToExtractorIndex(m_MC->getMom1Index(jetMatch))) << " and grandmother pdg index " << m_MC->getType(patIndexToExtractorIndex(m_MC->getMom1Index(patIndexToExtractorIndex(m_MC->getMom1Index(jetMatch))))) << endl;
	  if (jetMatch!=-10 && fabs(m_MC->getType(jetMatch))==5 && fabs(m_MC->getType(patIndexToExtractorIndex(m_MC->getMom1Index(jetMatch))))==25) ++MatchedHiggsJets;
	  if (jetMatch!=-10 && m_MC->getType(jetMatch)==5 && fabs(m_MC->getType(patIndexToExtractorIndex(m_MC->getMom1Index(jetMatch))))==25) {++MatchedFirstHiggsJets; MatchedJetsIndexes[0]=i;}
	  if (jetMatch!=-10 && m_MC->getType(jetMatch)==-5 && fabs(m_MC->getType(patIndexToExtractorIndex(m_MC->getMom1Index(jetMatch))))==25) {++MatchedSecondHiggsJets; MatchedJetsIndexes[1]=i;}
	  if (jetMatch!=-10 && fabs(m_MC->getType(jetMatch))<=5 && fabs(m_MC->getType(jetMatch))!=0 && fabs(m_MC->getType(patIndexToExtractorIndex(m_MC->getMom1Index(jetMatch))))==24) ++MatchedWJets;
	  if (jetMatch!=-10 && fabs(m_MC->getType(jetMatch))<=5 && m_MC->getType(jetMatch)>0 && fabs(m_MC->getType(patIndexToExtractorIndex(m_MC->getMom1Index(jetMatch))))==24) {++MatchedFirstWJets; MatchedJetsIndexes[2]=i;}
	  if (jetMatch!=-10 && fabs(m_MC->getType(jetMatch))<=5 && m_MC->getType(jetMatch)<0 && fabs(m_MC->getType(patIndexToExtractorIndex(m_MC->getMom1Index(jetMatch))))==24) {++MatchedSecondWJets; MatchedJetsIndexes[3]=i;}
	  if (jetMatch!=-10 && fabs(m_MC->getType(jetMatch))==5 && fabs(m_MC->getType(patIndexToExtractorIndex(m_MC->getMom1Index(jetMatch))))==6) {++MatchedTopJets; MatchedJetsIndexes[4]=i;}
	  
	}
      NumbMatchedHiggsJets=MatchedHiggsJets;
      NumbMatchedWJets=MatchedWJets;
      NumbMatchedTopJets=MatchedTopJets;
      
      //Taking into account only events where the MC indformation was obtained correctly
      if (FirstTrueHiggsJet->Pt()==0 || SecondTrueHiggsJet->Pt()==0 || FirstTrueWJet->Pt()==0 || SecondTrueWJet->Pt()==0 || TopTrueJet->Pt()==0 || TrueW->Pt()==0 || TrueTprimeAcompainingJet->Pt()==0) return 0;
      /*
      //MODIFICATION!!!!!!//

      int GoodMatchedJets=0;
      bool AllDiffIndex=true;
      
      for (int j=0; j<5; j++) 
	{
	  for (int k=j+1; k<5; k++) {if (MatchedJetsIndexes[j]==MatchedJetsIndexes[k]) AllDiffIndex=false;}
	}
      
      if (AllDiffIndex)
	{
	  for (int j=0; j<5; j++)
	    {
	      if (IndexesFromchi2[0]==MatchedJetsIndexes[j] || IndexesFromchi2[1]==MatchedJetsIndexes[j] || IndexesFromchi2[2]==MatchedJetsIndexes[j] || IndexesFromchi2[3]==MatchedJetsIndexes[j] || IndexesFromchi2[4]==MatchedJetsIndexes[j]) ++GoodMatchedJets;
	    }
	}
      else cout << "Shared index between decay products of tprime " << MatchedJetsIndexes[0] << MatchedJetsIndexes[1] << MatchedJetsIndexes[2] << MatchedJetsIndexes[3] << MatchedJetsIndexes[4] << endl;
      
      if (GoodMatchedJets==5) CorrectTprime=1;
      
      int HiggsGoodMatches=0;
      int WGoodMatches=0;
      for (int j=0; j<2; j++)
	{
	  if (IndexesFromchi2[0]==MatchedJetsIndexes[j] || IndexesFromchi2[1]==MatchedJetsIndexes[j]) HiggsGoodMatches++;
	  if (IndexesFromchi2[2]==MatchedJetsIndexes[j+2] || IndexesFromchi2[3]==MatchedJetsIndexes[j+2]) WGoodMatches++;
	}
      
      if (HiggsGoodMatches==2) CorrectH=1;
      if (WGoodMatches==2) CorrectW=1;
      if (IndexesFromchi2[4]==MatchedJetsIndexes[4]) CorrectTop=1;

      /////////////////////
      */
      //Disambiguation with Pt
      TLorentzVector TrHFJ; TrHFJ.SetPxPyPzE(FirstTrueHiggsJet->Px(), FirstTrueHiggsJet->Py(), FirstTrueHiggsJet->Pz(), FirstTrueHiggsJet->E());
      TLorentzVector TrHSJ; TrHSJ.SetPxPyPzE(SecondTrueHiggsJet->Px(), SecondTrueHiggsJet->Py(), SecondTrueHiggsJet->Pz(), SecondTrueHiggsJet->E());
      TLorentzVector TWFJ; TWFJ.SetPxPyPzE(FirstTrueWJet->Px(), FirstTrueWJet->Py(), FirstTrueWJet->Pz(), FirstTrueWJet->E());
      TLorentzVector TWSJ; TWSJ.SetPxPyPzE(SecondTrueWJet->Px(), SecondTrueWJet->Py(), SecondTrueWJet->Pz(), SecondTrueWJet->E());
      TLorentzVector TTJ; TTJ.SetPxPyPzE(TopTrueJet->Px(), TopTrueJet->Py(), TopTrueJet->Pz(), TopTrueJet->E());
      
      m_DRTrueWJets=TWFJ.DeltaR(TWSJ);
      m_DPhiTrueWJets=fabs(TWFJ.Phi()-TWSJ.Phi());
      
      if (MatchedFirstHiggsJets>2)
	{
	  float DeltaPtMatchedAndTrueHiggsJets=0;
	  int CounterDisam=0;
	  for (int i=0;i<n_jets;++i)
	    {
	      int jetMatch = m_jetMet->getJetMCIndex(i);
	      TLorentzVector jeti; jeti.SetPxPyPzE(AllJets[i].Px(),AllJets[i].Py(),AllJets[i].Pz(),AllJets[i].E());
	      if (CounterDisam==0 && jetMatch!=-10 && m_MC->getType(jetMatch)==5 && fabs(m_MC->getType(patIndexToExtractorIndex(m_MC->getMom1Index(jetMatch))))==25) {++CounterDisam; DeltaPtMatchedAndTrueHiggsJets=fabs(jeti.Pt()-TrHFJ.Pt()); MatchedJetsIndexes[0]=i;}
	      if (CounterDisam!=0 && jetMatch!=-10 && m_MC->getType(jetMatch)==5 && fabs(m_MC->getType(patIndexToExtractorIndex(m_MC->getMom1Index(jetMatch))))==25 && DeltaPtMatchedAndTrueHiggsJets>fabs(jeti.Pt()-TrHFJ.Pt())) {DeltaPtMatchedAndTrueHiggsJets=fabs(jeti.Pt()-TrHFJ.Pt()); MatchedJetsIndexes[0]=i;}
	    }
	}
      
      if (MatchedSecondHiggsJets>2)
	{
	  float DeltaPtMatchedAndTrueHiggsJets=0;
	  int CounterDisam=0;
	  for (int i=0;i<n_jets;++i)
	    {
	      int jetMatch = m_jetMet->getJetMCIndex(i);
	      TLorentzVector jeti; jeti.SetPxPyPzE(AllJets[i].Px(),AllJets[i].Py(),AllJets[i].Pz(),AllJets[i].E());
	      if (CounterDisam==0 && jetMatch!=-10 && m_MC->getType(jetMatch)==-5 && fabs(m_MC->getType(patIndexToExtractorIndex(m_MC->getMom1Index(jetMatch))))==25) {++CounterDisam; DeltaPtMatchedAndTrueHiggsJets=fabs(jeti.Pt()-TrHSJ.Pt()); MatchedJetsIndexes[1]=i;}
	      if (CounterDisam!=0 && jetMatch!=-10 && m_MC->getType(jetMatch)==-5 && fabs(m_MC->getType(patIndexToExtractorIndex(m_MC->getMom1Index(jetMatch))))==25 && DeltaPtMatchedAndTrueHiggsJets>fabs(jeti.Pt()-TrHSJ.Pt())) {DeltaPtMatchedAndTrueHiggsJets=fabs(jeti.Pt()-TrHSJ.Pt()); MatchedJetsIndexes[1]=i;}
	    }
	}
      
      if (MatchedFirstWJets>2)
	{
	  float DeltaPtMatchedAndTrueWJets=0;
	  int CounterDisam=0;
	  for (int i=0;i<n_jets;++i)
	    {
	      int jetMatch = m_jetMet->getJetMCIndex(i);
	      TLorentzVector jeti; jeti.SetPxPyPzE(AllJets[i].Px(),AllJets[i].Py(),AllJets[i].Pz(),AllJets[i].E());
	      if (CounterDisam==0 && jetMatch!=-10 && m_MC->getType(jetMatch)<=5 && m_MC->getType(jetMatch)>0 && fabs(m_MC->getType(patIndexToExtractorIndex(m_MC->getMom1Index(jetMatch))))==25) {++CounterDisam; DeltaPtMatchedAndTrueWJets=fabs(jeti.Pt()-TWFJ.Pt()); MatchedJetsIndexes[2]=i;}
	      if (CounterDisam!=0 && jetMatch!=-10 && m_MC->getType(jetMatch)<=5 && m_MC->getType(jetMatch)>0 && fabs(m_MC->getType(patIndexToExtractorIndex(m_MC->getMom1Index(jetMatch))))==25 && DeltaPtMatchedAndTrueWJets>fabs(jeti.Pt()-TWFJ.Pt())) {DeltaPtMatchedAndTrueWJets=fabs(jeti.Pt()-TWFJ.Pt()); MatchedJetsIndexes[2]=i;}
	    }
	}
      
      if (MatchedSecondWJets>2)
	{
	  float DeltaPtMatchedAndTrueWJets=0;
	  int CounterDisam=0;
	  for (int i=0;i<n_jets;++i)
	    {
	      int jetMatch = m_jetMet->getJetMCIndex(i);
	      TLorentzVector jeti; jeti.SetPxPyPzE(AllJets[i].Px(),AllJets[i].Py(),AllJets[i].Pz(),AllJets[i].E());
	      if (CounterDisam==0 && jetMatch!=-10 && m_MC->getType(jetMatch)>=-5 && m_MC->getType(jetMatch)<0 && fabs(m_MC->getType(patIndexToExtractorIndex(m_MC->getMom1Index(jetMatch))))==25) {++CounterDisam; DeltaPtMatchedAndTrueWJets=fabs(jeti.Pt()-TWSJ.Pt()); MatchedJetsIndexes[3]=i;}
	      if (CounterDisam!=0 && jetMatch!=-10 && m_MC->getType(jetMatch)>=-5 && m_MC->getType(jetMatch)<0 && fabs(m_MC->getType(patIndexToExtractorIndex(m_MC->getMom1Index(jetMatch))))==25 && DeltaPtMatchedAndTrueWJets>fabs(jeti.Pt()-TWSJ.Pt())) {DeltaPtMatchedAndTrueWJets=fabs(jeti.Pt()-TWSJ.Pt()); MatchedJetsIndexes[3]=i;}
	    }
	}
      
      if (MatchedTopJets>2)
	{
	  float DeltaPtMatchedAndTrueTopJets=0;
	  int CounterDisam=0;
	  for (int i=0;i<n_jets;++i)
	    {
	      int jetMatch = m_jetMet->getJetMCIndex(i);
	      TLorentzVector jeti; jeti.SetPxPyPzE(AllJets[i].Px(),AllJets[i].Py(),AllJets[i].Pz(),AllJets[i].E());
	      if (CounterDisam==0 && jetMatch!=-10 && fabs(m_MC->getType(jetMatch))==5 && fabs(m_MC->getType(patIndexToExtractorIndex(m_MC->getMom1Index(jetMatch))))==6) {++CounterDisam; DeltaPtMatchedAndTrueTopJets=fabs(jeti.Pt()-TTJ.Pt()); MatchedJetsIndexes[4]=i;}
	      if (CounterDisam!=0 && jetMatch!=-10 && fabs(m_MC->getType(jetMatch)==5) && fabs(m_MC->getType(patIndexToExtractorIndex(m_MC->getMom1Index(jetMatch))))==6 && DeltaPtMatchedAndTrueTopJets>fabs(jeti.Pt()-TTJ.Pt())) {DeltaPtMatchedAndTrueTopJets=fabs(jeti.Pt()-TTJ.Pt()); MatchedJetsIndexes[4]=i;}
	    }
	}
      
      if (MatchedJetsIndexes[0]!=-1 && MatchedJetsIndexes[1]!=-1 && MatchedJetsIndexes[2]!=-1 && MatchedJetsIndexes[3]!=-1 && MatchedJetsIndexes[4]!=-1)
	{
	  m_DRMatchedWJets=AllJets[MatchedJetsIndexes[2]].DeltaR(AllJets[MatchedJetsIndexes[3]]);
	  m_DPhiMatchedWJets=fabs(AllJets[MatchedJetsIndexes[2]].Phi()-AllJets[MatchedJetsIndexes[3]].Phi());
	  
	  TprimeFromMatching->SetPxPyPzE(AllJets[MatchedJetsIndexes[0]].Px()+AllJets[MatchedJetsIndexes[1]].Px()+AllJets[MatchedJetsIndexes[2]].Px()+AllJets[MatchedJetsIndexes[3]].Px()+AllJets[MatchedJetsIndexes[4]].Px(),AllJets[MatchedJetsIndexes[0]].Py()+AllJets[MatchedJetsIndexes[1]].Py()+AllJets[MatchedJetsIndexes[2]].Py()+AllJets[MatchedJetsIndexes[3]].Py()+AllJets[MatchedJetsIndexes[4]].Py(),AllJets[MatchedJetsIndexes[0]].Pz()+AllJets[MatchedJetsIndexes[1]].Pz()+AllJets[MatchedJetsIndexes[2]].Pz()+AllJets[MatchedJetsIndexes[3]].Pz()+AllJets[MatchedJetsIndexes[4]].Pz(),AllJets[MatchedJetsIndexes[0]].E()+AllJets[MatchedJetsIndexes[1]].E()+AllJets[MatchedJetsIndexes[2]].E()+AllJets[MatchedJetsIndexes[3]].E()+AllJets[MatchedJetsIndexes[4]].E());

	  int GoodMatchedJets=0;
	  bool AllDiffIndex=true;
	  
	  for (int j=0; j<5; j++) 
	    {
	      for (int k=j+1; k<5; k++) {if (MatchedJetsIndexes[j]==MatchedJetsIndexes[k]) AllDiffIndex=false;}
	    }
	  
	  if (AllDiffIndex)
	    {
	      for (int j=0; j<5; j++)
		{
		  if (IndexHiggsJets[0]==MatchedJetsIndexes[j] || IndexHiggsJets[1]==MatchedJetsIndexes[j] || IndexWJets[0]==MatchedJetsIndexes[j] || IndexWJets[1]==MatchedJetsIndexes[j] || IndexTopJet==MatchedJetsIndexes[j]) ++GoodMatchedJets;
		}
	    }
	  else cout << "Shared index between decay products of tprime after disambiguation " << MatchedJetsIndexes[0] << MatchedJetsIndexes[1] << MatchedJetsIndexes[2] << MatchedJetsIndexes[3] << MatchedJetsIndexes[4] << endl;
	  
	  if (GoodMatchedJets==5) CorrectTprime=1;
	  
	  int HiggsGoodMatches=0;
	  int WGoodMatches=0;
	  for (int j=0; j<2; j++)
	    {
	      if (IndexHiggsJets[0]==MatchedJetsIndexes[j] || IndexHiggsJets[1]==MatchedJetsIndexes[j]) HiggsGoodMatches++;
	      if (IndexWJets[0]==MatchedJetsIndexes[j+2] || IndexWJets[1]==MatchedJetsIndexes[j+2]) WGoodMatches++;
	    }
	  
	  if (HiggsGoodMatches==2) CorrectH=1;
	  if (WGoodMatches==2) CorrectW=1;
	  if (IndexTopJet==MatchedJetsIndexes[4]) CorrectTop=1;
	  
	}
      else return 0;
    }

  // Finally fill the analysis tree 
  //SingleTprime_analysis::fillTree();
  
  if (m_isMC) {
    // Get scale factors
    double sf = 1;
    double sf_error = 0;

    if (m_NBtaggedJets_CSVL >= 2) {

      ScaleFactor sf_jet1 = {1, 0, 0};
      ScaleFactor sf_jet2 = {1, 0, 0};
      bool gotJet1 = false;

      for (int i = 0; i < n_jets; i++)
      {
        if (jetIsBTagged[i]) {
          if (! gotJet1) {
            gotJet1 = true;
            sf_jet1 = jetSF[i];
          } else {
            sf_jet2 = jetSF[i];
            break;
          }
        }
      }

      sf = sf_jet1.getValue() * sf_jet2.getValue();
      sf_error = sf_jet2.getValue() * sf_jet2.getValue() * sf_jet1.getErrorHigh() * sf_jet1.getErrorHigh() +
        sf_jet1.getValue() * sf_jet1.getValue() * sf_jet2.getErrorHigh() * sf_jet2.getErrorHigh();
    }

    m_weight *= sf;

    double squared_sf_error = sf_error * sf_error;
    m_weight_error_low += squared_sf_error;
    m_weight_error_high += squared_sf_error;
  }

  //if (TrueTprime->Pt()==0) return 0;
  //TLorentzVector TrTp; TrTp.SetPxPyPzE(TrueTprime->Px(), TrueTprime->Py(), TrueTprime->Pz(), TrueTprime->E());
  //if (TrTp.Pt()<100) Njets_TpPT_0to100=n_jets;
  //else if (TrTp.Pt()>=100 && TrTp.Pt()<200) Njets_TpPT_100to200=n_jets;
  //else if (TrTp.Pt()>=200 && TrTp.Pt()<300) Njets_TpPT_200to300=n_jets;
  //else if (TrTp.Pt()>=300 && TrTp.Pt()<400) Njets_TpPT_300to400=n_jets;
  //else if (TrTp.Pt()>=400 && TrTp.Pt()<500) Njets_TpPT_400to500=n_jets;
  //else if (TrTp.Pt()>=500 && TrTp.Pt()<600) Njets_TpPT_500to600=n_jets;
  //else if (TrTp.Pt()>=600 && TrTp.Pt()<700) Njets_TpPT_600to700=n_jets;
  //else if (TrTp.Pt()>=700 && TrTp.Pt()<800) Njets_TpPT_700to800=n_jets;

  return 1;
}

bool SingleTprime_analysis::isJetSel(TLorentzVector *jet) //Function to select Jets of interest
{  
  double jetpt = jet->Pt();
  double jeteta = jet->Eta();

  if (jetpt<m_jet_Ptcut || fabs(jeteta)>m_jet_EtaMaxcut) return false;
  
  return true;
}

bool SingleTprime_analysis::isJetAccepSel(TLorentzVector *jet) //Function to select jets inside the acceptance
{  
  double jetpt = jet->Pt();
  double jeteta = jet->Eta();

  if (jetpt>m_jet_Ptcut && fabs(jeteta)<m_jet_EtaAccepcut) return true;

  return false;
}

bool SingleTprime_analysis::isJetForwSel(TLorentzVector *jet) // Function to select jets outside acceptance
{  
  double jetpt = jet->Pt();
  double jeteta = jet->Eta();

  if (jetpt>(m_jet_Ptcut-10) && fabs(jeteta)>=(m_jet_EtaAccepcut-m_jet_OverlapAccep) && fabs(jeteta)<m_jet_EtaMaxcut) return true;

  return false;
}

#define   MTT_TRIGGER_NOT_FOUND   1000

void SingleTprime_analysis::analyze(const edm::Event& event, const edm::EventSetup& iSetup, PatExtractor& extractor) {
  analyze(iSetup, extractor);
}

void SingleTprime_analysis::analyze(const edm::EventSetup& iSetup, PatExtractor& extractor)
{
  reset(); //Setting all variables to zero

  //Pointer to Extractor objects
  m_jetMet   = std::static_pointer_cast<JetMETExtractor>(extractor.getExtractor("JetMET"));
  m_event    = std::static_pointer_cast<EventExtractor>(extractor.getExtractor("event"));

  std::shared_ptr<HLTExtractor> HLT = std::static_pointer_cast<HLTExtractor>(extractor.getExtractor("HLT"));
  m_trigger_passed = HLT->isTriggerFired();

  if (m_isMC) 
    {
      m_MC = std::static_pointer_cast<MCExtractor>(extractor.getExtractor("MC"));
      MCidentification();
    }

  m_nPU = m_event->nPU();
  m_evt = m_event->n_events();
  evt_num++; cout << evt_num << endl;

  //Performing analysis
  //cout << "Going to selection" << endl;
  int ToFillTree = SingleTprime_Sel();

  if (ToFillTree==1) fillTree();
  m_tree_cuts->Fill();
}

  //MC Identification

#define ID_U (1)
#define ID_D (2)
#define ID_S (3)
#define ID_C (4)
#define ID_B (5)
#define ID_T (6)
#define ID_H (25)
#define ID_W (24)
#define ID_Tp (6000006)

int SingleTprime_analysis::patIndexToExtractorIndex(int patIndex) const {

  for (int i = 0; i < m_MC->getSize() ; i++) {
    if (m_MC->getPatIndex(i) == patIndex)
      return i;
  }

  return -1;
}

void SingleTprime_analysis::MCidentification()
{

  int n_MC = m_MC->getSize();

  if (!n_MC) return;

  for (int i = 0; i < n_MC ; ++i)
    {

      int motherIndex = patIndexToExtractorIndex(m_MC->getMom1Index(i));
      int grandMotherIndex = -1;
      if (motherIndex != -1) grandMotherIndex = patIndexToExtractorIndex(m_MC->getMom1Index(motherIndex));

      if (abs(m_MC->getType(i)) == ID_H) 
	{
	  TrueHiggs->SetPxPyPzE(m_MC->getPx(i),m_MC->getPy(i),m_MC->getPz(i),m_MC->getE(i));
	}
            
      if (abs(m_MC->getType(i)) == ID_T) 
	{
	  TrueTop->SetPxPyPzE(m_MC->getPx(i),m_MC->getPy(i),m_MC->getPz(i),m_MC->getE(i));
	}

      if (abs(m_MC->getType(i)) == ID_Tp) 
	{
	  TrueTprime->SetPxPyPzE(m_MC->getPx(i),m_MC->getPy(i),m_MC->getPz(i),m_MC->getE(i));
	  for (int j = 0; j < n_MC ; ++j)
	    {
	      
	      int motherIndex1 = patIndexToExtractorIndex(m_MC->getMom1Index(j));
	      int grandMotherIndex1 = -1;
	      if (motherIndex1 != -1) grandMotherIndex1 = patIndexToExtractorIndex(m_MC->getMom1Index(motherIndex1));
	     
	      if (abs(m_MC->getType(j)) < ID_B && abs(m_MC->getType(j))>0 && abs(m_MC->getType(motherIndex)) < ID_B && m_MC->getType(motherIndex)==m_MC->getType(motherIndex1)) 
		{
		  TrueTprimeAcompainingJet->SetPxPyPzE(m_MC->getPx(j),m_MC->getPy(j),m_MC->getPz(j),m_MC->getE(j));
		}
	      
	    }
	}

      if (motherIndex == -1) continue;
      TLorentzVector CurrentMCParticle; CurrentMCParticle.SetPxPyPzE(m_MC->getPx(i),m_MC->getPy(i),m_MC->getPz(i),m_MC->getE(i));
      if (CurrentMCParticle.Pt()==0) continue;

      if (abs(m_MC->getType(i)) == ID_W && abs(m_MC->getType(motherIndex)) == ID_T) 
	{
	  TrueW->SetPxPyPzE(m_MC->getPx(i),m_MC->getPy(i),m_MC->getPz(i),m_MC->getE(i));
	}
      
      if (m_MC->getType(i) == ID_B && abs(m_MC->getType(motherIndex)) == ID_H) 
	{
	  FirstTrueHiggsJet->SetPxPyPzE(m_MC->getPx(i),m_MC->getPy(i),m_MC->getPz(i),m_MC->getE(i));
	}

      if (m_MC->getType(i) == -1*ID_B && abs(m_MC->getType(motherIndex)) == ID_H) 
	{
	  SecondTrueHiggsJet->SetPxPyPzE(m_MC->getPx(i),m_MC->getPy(i),m_MC->getPz(i),m_MC->getE(i));
	}

      if (m_MC->getType(i) <= ID_B && m_MC->getType(i) > 0 && abs(m_MC->getType(motherIndex)) == ID_W) 
	{
	  FirstTrueWJet->SetPxPyPzE(m_MC->getPx(i),m_MC->getPy(i),m_MC->getPz(i),m_MC->getE(i));
	}

      if (m_MC->getType(i) >= -1*ID_B && m_MC->getType(i) < 0 && abs(m_MC->getType(motherIndex)) == ID_W) 
	{
	  SecondTrueWJet->SetPxPyPzE(m_MC->getPx(i),m_MC->getPy(i),m_MC->getPz(i),m_MC->getE(i));
	}

      if (abs(m_MC->getType(i)) == ID_B && abs(m_MC->getType(motherIndex)) == ID_T) 
	{
	  TopTrueJet->SetPxPyPzE(m_MC->getPx(i),m_MC->getPy(i),m_MC->getPz(i),m_MC->getE(i));
	}

      if (false) 
	{
	  std::cout << "Type: " << m_MC->getType(i) << std::endl;
	  std::cout << "Mother type: " << m_MC->getType(motherIndex) << std::endl;
	  if (grandMotherIndex != -1) std::cout << "Grandmother type: " << m_MC->getType(grandMotherIndex) << std::endl;
	}
      
      if (abs(m_MC->getType(i)) == ID_U)
	{
	  UQuarkContent+=1;
	}
      else if (abs(m_MC->getType(i)) == ID_D)
	{
	  DQuarkContent+=1;
	}
      else if (abs(m_MC->getType(i)) == ID_C)
	{
	  SQuarkContent+=1;
	}
      else if (abs(m_MC->getType(i)) == ID_S)
	{
	  CQuarkContent+=1;
	}
      else if (abs(m_MC->getType(i)) == ID_B)
	{
	  BQuarkContent+=1;
	}

      
    }

}

void SingleTprime_analysis::reset()
{
  m_evt = 0;
  m_nPU = 0;

  m_THT = 0.;
  m_jet1pt = 0.;
  m_jet2pt = 0.;
  m_jet3pt = 0.;
  m_jet4pt = 0.;
  m_jet5pt = 0.;
  m_jet6pt = 0.;
  m_NBtaggedJets_CSVL=0;
  m_NBtaggedJets_CSVM=0;
  m_NBtaggedJets_CSVT=0;
  m_DRHiggsJets=0.;
  m_DRWJets=0.;
  m_DRTopHiggs=0.;
  m_DRWHiggs=0.;
  m_RelTHT=0.;
  m_Sphericity=0.;
  m_Aplanarity=0.;
  m_DRTrueHiggsRecoHiggs=0.;
  m_DRTrueWRecoW=0.;
  m_DRTrueTopRecoTop=0.;
  m_DRTrueTprimeRecoTprime=0.;
  m_DRTrueFirstHiggsJetRecoJet=0;
  m_DRTrueSecondHiggsJetRecoJet=0;
  m_DRTrueFirstWJetRecoJet=0;
  m_DRTrueSecondWJetRecoJet=0; 
  m_DRTrueTopJetRecoJet=0;

  RelMass=0.;
  PTNormalizedMass=0.;
  MotherPTNormalizedMass=0.;
  NumberOfTops=0;
  LooseNoMedBtags=0;
  m_DP2LeadingJets=0;
  chi2=0;

  UQuarkContent=0;
  DQuarkContent=0;
  SQuarkContent=0;
  CQuarkContent=0;
  BQuarkContent=0;

  PtTTbarDifference=-10;

  Njets_TpPT_0to100=0;
  Njets_TpPT_100to200=0;
  Njets_TpPT_200to300=0;
  Njets_TpPT_300to400=0;
  Njets_TpPT_400to500=0;
  Njets_TpPT_500to600=0;
  Njets_TpPT_600to700=0;
  Njets_TpPT_700to800=0;
  
  m_DRTrueWJets=0;   
  m_DRMatchedWJets=0;
  m_DPhiTrueWJets=0;
  m_DPhiMatchedWJets=0;

  m_VBF_M=0;

  CorrectTprime = 0;
  CorrectH = 0;
  CorrectW = 0;
  CorrectTop = 0;
  CorrectHiggsJet = 0;
  CorrectWJet = 0;
  CorrectTopJet = 0;
  NumbMatchedHiggsJets = 0;
  NumbMatchedWJets = 0;
  NumbMatchedTopJets = 0;

  m_triggercut= 0;
  m_Cut0= 0;
  m_Cut1= 0;
  m_Cut2= 0;
  m_Cut3= 0;
  m_CutChi2= 0;
  m_Cut4= 0;
  m_Cut5= 0;
  m_Cut6= 0;
  m_Cut7= 0;
  m_Cut8= 0;
  m_Cut9= 0;
  m_Cut10= 0;
  m_Cut11= 0;
  m_Cut12= 0;
  m_Cut13= 0;
  m_Cut14= 0;
  m_Cut15= 0;
  m_Cut16= 0;
  m_Cut17= 0;
  m_Cut18= 0;
  m_Cut19= 0;
  m_Cut20= 0;
  m_Cut21= 0;

  m_trigger_passed = false;

  m_weight = 1.;
  m_weight_error_low = 0.;
  m_weight_error_high = 0.;

  ReconstructedHiggs->SetPxPyPzE(0., 0., 0., 0.);
  ReconstructedW->SetPxPyPzE(0., 0., 0., 0.);
  ReconstructedTop->SetPxPyPzE(0., 0., 0., 0.);
  ReconstructedTprime->SetPxPyPzE(0., 0., 0., 0.);
  TrueHiggs->SetPxPyPzE(0., 0., 0., 0.);
  TrueW->SetPxPyPzE(0., 0., 0., 0.);
  TrueTop->SetPxPyPzE(0., 0., 0., 0.);
  TrueTprime->SetPxPyPzE(0., 0., 0., 0.);
  TrueTprimeAcompainingJet->SetPxPyPzE(0., 0., 0., 0.);
  FirstTrueHiggsJet->SetPxPyPzE(0., 0., 0., 0.);
  SecondTrueHiggsJet->SetPxPyPzE(0., 0., 0., 0.);
  FirstTrueWJet->SetPxPyPzE(0., 0., 0., 0.);
  SecondTrueWJet->SetPxPyPzE(0., 0., 0., 0.);
  TopTrueJet->SetPxPyPzE(0., 0., 0., 0.);
  TprimeFromMatching->SetPxPyPzE(0., 0., 0., 0.);
}

// Fill the root tree containing analysis results

void SingleTprime_analysis::fillTree()
{
  m_weight_error_low = sqrt(m_weight_error_low);
  m_weight_error_high = sqrt(m_weight_error_high);
  m_tree_stp->Fill(); 
}

}

// Register the plugin inside the factory
DEFINE_EDM_PLUGIN(PatExtractorPluginFactory, patextractor::SingleTprime_analysis, "SingleTprime_analysis");
