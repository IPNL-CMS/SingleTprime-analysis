#include <Extractors/PatExtractor/interface/ExtractorPlugin.h>

#include <vector>
#include <boost/regex.hpp>

#include <TTree.h>
#include <TLorentzVector.h>
#include <TRef.h>

class AnalysisSettings;
class EventExtractor;
class MuonExtractor;
class ElectronExtractor;
class METExtractor;
class VertexExtractor;
class KinFit;
class HLTExtractor;
class PatExtractor;
class MCExtractor;
class JetMETExtractor;

namespace edm {
  class EventSetup;
}

namespace patextractor {

  //class SingleTprime_analysis: patextractor::Plugin {
  class SingleTprime_analysis: public Plugin {
  public:
    SingleTprime_analysis(const edm::ParameterSet& iConfig);
    ~SingleTprime_analysis();
    
    virtual void analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup, PatExtractor& extractor);
    virtual void analyze(const edm::EventSetup& iSetup, PatExtractor& extractor);

    //Functions
    bool isJetSel(TLorentzVector *jet);
    bool isJetAccepSel(TLorentzVector *jet);
    bool isJetForwSel(TLorentzVector *jet);
    int  SingleTprime_Sel();
    int  SingleTprime_Sel2M1L();
    int patIndexToExtractorIndex(int patIndex) const;

    void MCidentification();
    void reset();
    void fillTree();
    void fillTree2M1L();

  private:
    //int j;

    TTree* m_tree_stp;
    TTree* m_tree_cuts;
    TTree* m_tree_stp2M1L;
    //TTree* m_tree_cuts2M1L;    

    int m_evt;
    int m_nPU; //PileUp
    int evt_num;
    int n_vtx;
    int n_TrueInteractions;

    int n_jets;    
    float m_THT;
    float m_jet1pt;
    float m_jet2pt;
    float m_jet3pt;
    float m_jet4pt;
    float m_jet5pt;
    float m_jet6pt;
    float m_DRHiggsJets;
    float m_DRWJets;
    float m_DRTopHiggs;
    float m_DRWHiggs;
    float m_RelTHT;
    float m_Sphericity;
    float m_Aplanarity;
    float m_DRTrueHiggsRecoHiggs;
    float m_DRTrueWRecoW;
    float m_DRTrueTopRecoTop;
    float m_DRTrueTprimeRecoTprime;
    float m_DRTrueFirstHiggsJetRecoJet;
    float m_DRTrueSecondHiggsJetRecoJet;
    float m_DRTrueFirstWJetRecoJet;
    float m_DRTrueSecondWJetRecoJet; 
    float m_DRTrueTopJetRecoJet;   
    float m_DP2LeadingJets;

    float DRWjetsFromHiggsChi2;
    float DRWbFromHiggsChi2;

    float RelMass;
    float PTNormalizedMass;
    float MotherPTNormalizedMass;
    int NumberOfTops;
    int NumberOfHiggses;
    int LooseNoMedBtags;
    float chi2;

    int UQuarkContent;
    int DQuarkContent;
    int SQuarkContent;
    int CQuarkContent;
    int BQuarkContent;
    int RealLB_tag_ToB;
    int FakeLB_tag_ToLight;
    int FakeLB_tag_ToC;
    int RealMB_tag_ToB;
    int FakeMB_tag_ToLight;
    int FakeMB_tag_ToC;
    int RealTB_tag_ToB;
    int FakeTB_tag_ToLight;
    int FakeTB_tag_ToC;

    float PtTTbarDifference;

    //Proving non-boosted regime//
    float Njets_TpPT_0to100;
    float Njets_TpPT_100to200;
    float Njets_TpPT_200to300;
    float Njets_TpPT_300to400;
    float Njets_TpPT_400to500;
    float Njets_TpPT_500to600;
    float Njets_TpPT_600to700;
    float Njets_TpPT_700to800;
    //////////////////////////////

    float m_DRTrueWJets;
    float m_DRMatchedWJets;
    float m_DPhiTrueWJets;   
    float m_DPhiMatchedWJets;

    float CSV3MMinimum;
    float CSV3MSum;
    float CSV3MAverage;

    float DeltaEtaTpMaxEtaJet;

    //float WFromHiggs;
    //float TopFromHiggs;
    int LeadingAfterRecoJetMother;

    //bool isMC;
    float m_jet_Ptcut;
    float m_jet_EtaMaxcut;
    float m_jet_EtaAccepcut;
    float m_jet_OverlapAccep;
    float m_jet_MultInAcceptance;
    float m_jet_MultOutAcceptance;
    float m_DRMatching;
    float m_DPtMatching;

    //Booleans to activate each cut
    bool Cut0; //Event preselection
    bool Cut1; //Leading jet PT
    bool Cut2; //HT
    bool Cut3; //B-tags
    bool CutChi2; //Chi2
    bool Cut4; //DR(jj^{H})
    bool Cut5; //DR(jj^{W})
    bool Cut6; //H pt and Top pt
    bool Cut7; //DR(WH)
    bool Cut8; //DPH and DPT
    bool Cut9; //Jetmultiplicity
    bool Cut10; //DPH DPW
    bool Cut11; //MH
    bool Cut12; //RelHT
    bool Cut13; //Aplanarity
    bool Cut14; //DR(TH)
    bool Cut15; //Mh related to Mtop
    bool Cut16; //Normalized Pt Mass
    bool Cut17; //Normalized Mother Pt Mass
    bool Cut18; //One Top
    bool Cut19; //Third loose b-tag
    bool Cut20; //Delta Phi between two leading jets
    bool Cut21; //HMass over TMass
    bool Cut22; //WMass from Higgs with Chi2

    int NumberOfGoodJets; //////////////////////////////////////////
    int NumberOfBadJets;                                          //
    float LeadingJetPt;                                           //
    float THTcut;                                                 //
    int MinB_tags;
    int MinLooseB_tags;
    int MinTightB_tags;
    float MaxChi2;
    float DeltaRHiggsJets;                                        //
    float DeltaRWJets;                                            //
    float HiggsPt;                                                //
    float TopPt;                                                  //
    float MinDeltaRWH;                                            //
    float MaxDeltaRWH;                                            //
    float DeltaPhiHiggsJets; //              SET OF CUTS          //
    float DeltaPhiTopJetW;                                        //
    int JetMultiplicity;                                          //
    float DeltaPhiWjets;                                          //
    float MinHiggsMass;                                           //
    float MaxHiggsMass;                                           //
    float RelHT;                                                  //
    float Aplanarity;                                             //
    float MaxDeltaRTH;                                            //
    float MinDeltaRTH;                                            // 
    float RelMassMaxCut;                                          //
    float RelMassMinCut;                                          //
    float MotherPTNormalizedMassCut;                              //
    float PTNormalizedMassCut;                                    //
    int NumberOfTopsCut;                                          //
    float NumberofLooseBtag;                                      //
    float DeltaPhiLeadingJets;                                    //
    float HMoverTM;                /////////////////////////////////
    float Wchi2Min;
    float Wchi2Max;
 
    bool DoMCMatching;
    bool DoChi2Sorting;

    //Cuts
    int m_triggercut;
    int m_Cut0;
    int m_Cut1;
    int m_Cut2;
    int m_Cut3;
    int m_CutChi2;
    int m_Cut4;
    int m_Cut5;
    int m_Cut6;
    int m_Cut7;
    int m_Cut8;
    int m_Cut9;
    int m_Cut10;
    int m_Cut11;
    int m_Cut12;
    int m_Cut13;
    int m_Cut14;
    int m_Cut15;
    int m_Cut16;
    int m_Cut17;
    int m_Cut18;
    int m_Cut19;
    int m_Cut20;
    int m_Cut21;
    int m_Cut22;
 
    //Reconstructed particles
    TLorentzVector* ReconstructedHiggs;
    TLorentzVector* ReconstructedW;
    TLorentzVector* ReconstructedTop;
    TLorentzVector* ReconstructedTprime;
    TLorentzVector* TrueHiggs;
    TLorentzVector* TrueW;
    TLorentzVector* TrueTop;
    TLorentzVector* TrueTprime;
    TLorentzVector* TrueTprimeAcompainingJet;
    TLorentzVector* FirstHiggsJet;
    TLorentzVector* SecondHiggsJet;
    TLorentzVector* FirstWJet;
    TLorentzVector* SecondWJet;
    TLorentzVector* TopJet;
    TLorentzVector* FirstTrueHiggsJet;
    TLorentzVector* SecondTrueHiggsJet;
    TLorentzVector* FirstTrueWJet;
    TLorentzVector* SecondTrueWJet;
    TLorentzVector* TopTrueJet;
    TLorentzVector* TprimeFromMatching;
    TLorentzVector* WFromHiggs;
    TLorentzVector* TopFromHiggs;
    TLorentzVector* WFromHiggsChi2;
    TLorentzVector* TopFromHiggsChi2;
    TLorentzVector* ReconstructedHiggs2M1L;
    TLorentzVector* ReconstructedW2M1L;
    TLorentzVector* ReconstructedTop2M1L;
    TLorentzVector* ReconstructedTprime2M1L;
    TLorentzVector* FirstHiggsJet2M1L;
    TLorentzVector* SecondHiggsJet2M1L;
    TLorentzVector* FirstWJet2M1L;
    TLorentzVector* SecondWJet2M1L;
    TLorentzVector* TopJet2M1L;
    TLorentzVector* WFromHiggsChi22M1L;
    TLorentzVector* TopFromHiggsChi22M1L;

    //Matching with MC truth
    int CorrectTprime;
    int CorrectH;
    int CorrectW;
    int CorrectTop;
    int CorrectHiggsJet;
    int CorrectWJet;
    int CorrectTopJet;
    int NumbMatchedHiggsJets;
    int NumbMatchedWJets;
    int NumbMatchedTopJets;
    int HiggsTTbarReco;

    //Linking extractors
    std::shared_ptr<EventExtractor> m_event;
    std::shared_ptr<MCExtractor>    m_MC;
    std::shared_ptr<VertexExtractor> m_vertex;
    std::shared_ptr<JetMETExtractor> m_jetMet;
    float m_JET_btag_CSVL;
    float m_JET_btag_CSVM;
    float m_JET_btag_CSVT;

    //Counter for b-tagged jets
    int m_NBtaggedJets_CSVL;
    int m_NBtaggedJets_CSVM;
    int m_NBtaggedJets_CSVT;

    /// scale factors for jet energy resolution
    std::vector<double> jetEnergyResolutionScaleFactors_;
    std::vector<double> jetEnergyResolutionEtaBinning_;

    bool m_trigger_passed;	

    float m_weight;	
    float m_weight_error_low;
    float m_weight_error_high;

    float m_VBF_M;

    float PU_weight;

  };
  
}
