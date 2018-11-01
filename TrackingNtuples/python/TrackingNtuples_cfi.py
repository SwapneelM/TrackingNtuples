import FWCore.ParameterSet.Config as cms

process = cms.Process("TrackingNTuples")

process.maxLuminosityBlocks = cms.untracked.PSet(
            input = cms.untracked.int32(-1)
)

process.options = cms.untracked.PSet(
	SkipEvent	= cms.untracked.vstring('ProductNotFound')
)

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(-1)
)

process.TFileService=cms.Service(
    'TFileService',
    fileName=cms.string('outfile.root')
    )

process.ntuples = cms.EDAnalyzer(
  'MyTrackingNtuples',
  pixelTracks = cms.InputTag('pixelTracks'),
  matchedRecHits = cms.InputTag("siStripMatchedRecHits","matchedRecHit"),
  rphiRecHits = cms.InputTag("siStripMatchedRecHits","rphiRecHit"),
  stereoRecHits = cms.InputTag("siStripMatchedRecHits","stereoRecHit"),
  siPixelRecHits = cms.InputTag("siPixelRecHits")
)

process.path = cms.Path(process.ntuples)
