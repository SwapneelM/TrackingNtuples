import FWCore.ParameterSet.Config as cms

process = cms.Process("RECO")

process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.MessageLogger = cms.Service(
    "MessageLogger",
    destinations   = cms.untracked.vstring(
       'detailedInfo'
         ,'critical'
    ),
    ,critical      = cms.untracked.PSet(
                   , threshold = cms.untracked.string('ERROR') 
    ),
    detailedInfo   = cms.untracked.PSet(
                  threshold = cms.untracked.string('DEBUG')
    ),
    debugModules  = cms.untracked.vstring("*")
)

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
#        'file:/data/ml/smehta/cmssw-hep/CMSSW_10_2_4_Patatrack/src/10824.5_TTbar_13+TTbar_13TeV_TuneCUETP8M1_2018_GenSimFull+DigiFull_2018+RecoFull_pixelTrackingOnly_2018+HARVESTFull_pixelTrackingOnly_2018/step2.root',
        'file:/data/ml/smehta/cmssw-hep/CMSSW_10_2_4_Patatrack/src/10824.5_TTbar_13+TTbar_13TeV_TuneCUETP8M1_2018_GenSimFull+DigiFull_2018+RecoFull_pixelTrackingOnly_2018+HARVESTFull_pixelTrackingOnly_2018/step3.root'
    )
)

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(10)
)

process.EventContentAnalyzer = cms.EDAnalyzer('EventContentAnalyzer',
  indentation = cms.untracked.string('++'),
  verbose = cms.untracked.bool(False),
  verboseIndentation = cms.untracked.string('  '),
  verboseForModuleLabels = cms.untracked.vstring(),
  getData = cms.untracked.bool(False),
  getDataForModuleLabels = cms.untracked.vstring(),
  listContent = cms.untracked.bool(True)
)

process.path = cms.Path(process.EventContentAnalyzer)

