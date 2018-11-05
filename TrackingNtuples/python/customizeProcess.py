import FWCore.ParameterSet.Config as cms
from SimTracker.TrackHistory.TrackClassifier_cff import *
# from SimTracker.TrackAssociation.trackingParticleRecoTrackAsssociation_cfi import *

def customize(process, outfile='outfile.root'):
    '''Portable function to run the custom tracking ntuples in a RECOStep CMSSW config'''

    process.load("SimTracker.TrackAssociation.trackingParticleRecoTrackAsssociation_cfi")

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
        associator = cms.InputTag('trackingParticleRecoTrackAsssociation')
    )

    process.trackingParticleRecoTrackAsssociation.label_tr = cms.InputTag("pixelTracks")

    # Overriding the variables in track association
    # This was obtained from the 'process.load' above
    #process.trackreconstruction = process.trackingParticleRecoTrackAsssociation.clone(
    #            associator = 'QuickTrackAssociatorByHits'
    #            )


    process.reconstruction_pixelTrackingOnly *= process.reconstruction
    # process.reconstruction_pixelTrackingOnly *= process.trackreconstruction
    process.reconstruction_pixelTrackingOnly *= process.ntuples
