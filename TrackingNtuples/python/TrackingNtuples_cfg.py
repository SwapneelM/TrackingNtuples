import FWCore.ParameterSet.Config as cms

process = cms.Process("RECO")

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
        'file:step2.root',
        'file:step3.root'
    )
)

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(1)
)

process.printContent = cms.EDAnalyzer('EventContentAnalyzer',
  indentation = cms.untracked.string('++'),
  verbose = cms.untracked.bool(False),
  verboseIndentation = cms.untracked.string('  '),
  verboseForModuleLabels = cms.untracked.vstring(),
  getData = cms.untracked.bool(False),
  getDataForModuleLabels = cms.untracked.vstring(),
  listContent = cms.untracked.bool(True)
)

process.path = cms.Path(process.printContent)
