import FWCore.ParameterSet.Config as cms

def customize(process, outfile='outfile.root'):
    '''Protable funtion to run our tracking ntuples in a RECOStep CMSSW config'''
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
        siPixelRecHits = cms.InputTag("siPixelRecHits")
    )

    # Overriding the variables in track association
    process.trackingParticleRecoTrackAsssociation.label_tr = cms.InputTag("pixelTracks")

    process.reconstruction_pixelTrackingOnly *= process.reconstruction
    process.reconstruction_pixelTrackingOnly *= process.ntuples
