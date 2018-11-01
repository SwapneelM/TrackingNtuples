import FWCore.ParameterSet.Config as cms

maxLuminosityBlocks = cms.untracked.PSet(
            input = cms.untracked.int32(-1)
)

options = cms.untracked.PSet(
	SkipEvent	= cms.untracked.vstring('ProductNotFound')
)

maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(-1)
)

TFileService=cms.Service(
    'TFileService',
    fileName=cms.string('outfile.root')
    )

ntuples = cms.EDAnalyzer(
  'MyTrackingNtuples',
  pixelTracks = cms.InputTag('pixelTracks'),
  matchedRecHits = cms.InputTag("siStripMatchedRecHits","matchedRecHit"),
  rphiRecHits = cms.InputTag("siStripMatchedRecHits","rphiRecHit"),
  stereoRecHits = cms.InputTag("siStripMatchedRecHits","stereoRecHit"),
  siPixelRecHits = cms.InputTag("siPixelRecHits")
)

path = cms.Path(process.ntuples)
