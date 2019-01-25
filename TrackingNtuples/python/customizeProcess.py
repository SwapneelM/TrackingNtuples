import FWCore.ParameterSet.Config as cms
# from SimTracker.TrackHistory.TrackClassifier_cff import *

# This is needed to find the "quickTrackAssociatorByHits"
from SimTracker.TrackAssociatorProducers.quickTrackAssociatorByHits_cfi import *
from SimGeneral.TrackingAnalysis.simHitTPAssociation_cfi import *
from SimTracker.TrackAssociation.trackingParticleRecoTrackAsssociation_cfi import *

def customize(process, outfile='outfile.root'):
    '''Portable function to run the custom tracking ntuples in a RECOStep CMSSW config'''

    process.TFileService = cms.Service(
        'TFileService',
        fileName=cms.string(outfile)
    )

    process.ntuples = cms.EDAnalyzer(
        'MyTrackingNtuples',
        pixelTracks = cms.InputTag('pixelTracks'),
        matchedRecHits = cms.InputTag("siStripMatchedRecHits","matchedRecHit"),
        rphiRecHits = cms.InputTag("siStripMatchedRecHits","rphiRecHit"),
        stereoRecHits = cms.InputTag("siStripMatchedRecHits","stereoRecHit"),
        siPixelRecHits = cms.InputTag("siPixelRecHits"),
        # associator = cms.InputTag("trackingParticleRecoTrackAsssociation")
        simHitTPMap = cms.InputTag("simHitTPAssocProducer"),
        associator = cms.InputTag("swapneelAssociator"),
        clusterTPMap = cms.InputTag("tpClusterProducer")
    )
    process.swapneelAssociator = process.trackingParticleRecoTrackAsssociation.clone(
                                    label_tr = cms.InputTag("pixelTracks")
                                )

    process.reconstruction_pixelTrackingOnly *= process.reconstruction
    process.reconstruction_pixelTrackingOnly *= process.simHitTPAssocProducer
    process.reconstruction_pixelTrackingOnly *= process.tpClusterProducer
    process.reconstruction_pixelTrackingOnly *= process.quickTrackAssociatorByHits
    process.reconstruction_pixelTrackingOnly *= process.swapneelAssociator
    process.reconstruction_pixelTrackingOnly *= process.ntuples

