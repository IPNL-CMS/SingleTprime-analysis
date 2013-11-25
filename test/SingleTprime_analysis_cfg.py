####################################
# Variable to run only the analysis#
####################################
OnlyAnalysis = True
####################################

import FWCore.ParameterSet.Config as cms

def readFile(file):
  return cms.untracked.string(open(file).read())

def createExtractorProcess(isMC, isSemiMu, useShiftCorrectedMET, globalTag):
  process = cms.Process("PATextractor2")

  process.load('Configuration/StandardSequences/Services_cff')
  process.load('Configuration/StandardSequences/GeometryIdeal_cff')
  process.load('Configuration/StandardSequences/MagneticField_38T_cff')
  process.load('Configuration/StandardSequences/EndOfProcess_cff')
  process.load('Configuration/StandardSequences/FrontierConditions_GlobalTag_cff')
  process.load("FWCore.MessageLogger.MessageLogger_cfi")
  process.load("Extractors.PatExtractor.PAT_extractor_cff")
  
  #Global tag and data type choice
  process.GlobalTag.globaltag = '%s::All' % globalTag
  process.PATextraction.isMC  = isMC
  process.PATextraction.doMC  = isMC

  #Input PAT file to extract
  if not OnlyAnalysis:
    process.source = cms.Source("PoolSource",fileNames = cms.untracked.vstring(),duplicateCheckMode = cms.untracked.string( 'noDuplicateCheck' ))
    process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(1000))
  else:
    process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(1))
    process.PATextraction.n_events=100000
    process.source = cms.Source("EmptySource")  
    process.PATextraction.fillTree   = False
    process.PATextraction.inputRootFile=cms.string('Signal_extracted_mc.root')

  #Output extracted file name
  if isMC:
    if not OnlyAnalysis: process.PATextraction.extractedRootFile = cms.string('extracted_mc.root')
    else: process.PATextraction.extractedRootFile = cms.string('analyzed.root')
  else:
    process.PATextraction.extractedRootFile = cms.string('extracted.root')

  #########################################
  #
  #  Loading my analysis
  #
  #########################################

  process.PATextraction.plugins = cms.PSet( # <1>
    SingleTprime_analysis = cms.PSet(
        jets = cms.PSet(
          pt_min = cms.double(30),
          eta_max = cms.double(4.5),
          eta_accept = cms.double(2.5),
          eta_overlap = cms.double(2.5),
          btag_CSVL = cms.double(0.679)
          ),
	cuts = cms.PSet(
	  cut0 = cms.bool(True), #Event preselection
	  cut1 = cms.bool(True), #Leading jet PT
	  cut2 = cms.bool(True), #HT
          cut3 = cms.bool(True), #2btags
          cut4 = cms.bool(True), #Higgs jets DR<2.5
          cut5 = cms.bool(True), #W jets DR<3
          cut6 = cms.bool(True), #H pt and Top pt
          cut7 = cms.bool(True), #DR(WH)
          cut8 = cms.bool(True), #DPH and DPT
          cut9 = cms.bool(True), #Jet multiplicity
          cut10 = cms.bool(True), #DPH and DPW
          cut11 = cms.bool(True), #MH
          cut12 = cms.bool(True), #RelHT
          cut13 = cms.bool(True), #Aplanarity
          cut14 = cms.bool(False), #DR(TH)
          cut15 = cms.bool(False)  #RelMass
          ),
        DoMatching = cms.bool(False)
      )
    )
  
  #########################################
  #
  # PAT extractor main options statements
  #
  #########################################
  # Jets correction : needs a valid global tags, or an external DB where JEC are stored
  process.PATextraction.jet_PF.redoJetCorrection = True
  
  if isMC:
    process.PATextraction.jet_PF.jetCorrectorLabel = "ak5PFchsL1FastL2L3"
  else:
    process.PATextraction.jet_PF.jetCorrectorLabel = "ak5PFchsL1FastL2L3Residual"

  process.PATextraction.jet_PF.doJER = True # Disable automatically on data

  process.PATextraction.doJet      = True
  process.PATextraction.doMET      = True
  process.PATextraction.doHLT      = True
  process.PATextraction.doElectron = True
  process.PATextraction.doMuon     = True
  process.PATextraction.doVertex   = True

  #if not isMC: process.PATextraction.triggersXML = readFile("triggers_jets.xml")

  process.PATextraction.triggersXML = readFile("triggers_jets.xml")
  
  # JER systematics:
  # Use -1 for 1-sigma down, 0 for nominal correction, and 1 for 1-sigma up
  process.PATextraction.jet_PF.jerSign = 0

  from Extractor_MTT_ScaleFactors import loadMuonScaleFactor, loadBTagScaleFactors, loadElectronScaleFactor
  loadBTagScaleFactors(process)

  # Scale factors
  process.PATextraction.muon_scale_factors = loadMuonScaleFactor("Muon_ID_iso_Efficiencies_Run_2012ABCD_53X.pkl")
  process.PATextraction.electron_scale_factors = loadElectronScaleFactor("Electron_scale_factors.json")

  #########################################
  #
  # Launch the job
  #
  #########################################

  process.p = cms.Path(process.PATextraction)
  process.MessageLogger.cerr.FwkReport.reportEvery = 1000

  return process

if __name__ == "__main__":

  process = createExtractorProcess(True, False, useShiftCorrectedMET = False, globalTag = "START53_V27")
  
  # To build Transient Tracks
  process.load("TrackingTools.TransientTrack.TransientTrackBuilder_cfi")
  
  if not OnlyAnalysis:
    FileList=[]
    BASE="file:/home/cms/jruizalv/work/CMSSW_5_3_9_patch2/src/Extractors/SignalSamples/734GeV/Signal-PAT_tuple_"
    for i in xrange(1000):
      if i!=175: FileList.append(BASE+str(i)+".root")
    process.source.fileNames = cms.untracked.vstring( FileList
    #    'file:/home/cms/jruizalv/work/CMSSW_5_3_9_patch2/src/Extractors/Small_Sample_patTuple.root'
    #'file:/home/cms/jruizalv/work/CMSSW_5_3_9_patch2/src/Extractors/PatExtractor/test/extracted_mc.root'      
        )
    process.maxEvents.input = -1
