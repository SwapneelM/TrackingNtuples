import FWCore.ParameterSet.Config as cms

process = cms.Process("TrackingNTuples")

'''
process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.MessageLogger.cerr.threshold = 'INFO'
process.MessageLogger.cerr.INFO = cms.untracked.PSet(
    limit = cms.untracked.int32(-1)
)
process.MessageLogger = cms.Service(
    "MessageLogger",
    destinations   = cms.untracked.vstring(
       'detailedInfo',
       'critical'
    ),
    critical      = cms.untracked.PSet(
                  threshold = cms.untracked.string('ERROR')
    ),
    detailedInfo   = cms.untracked.PSet(
                  threshold = cms.untracked.string('INFO')
    )
)
process.MessageLogger.categories.append('Tracks', 'TrackExtra')
'''

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
