// -*- C++ -*-
//
// Class:      TrackingNtuples
//
/**\class TrackingNtuples TrackingNtuples.cc TrackingNtuples/TrackingNtuples/plugins/TrackingNtuples.cc

Description: [one line class summary]

         Implementation:
             [Notes on implementation]
*/
//
// Original Author:  Swapneel Sundeep Mehta
//         Created:  Mon, 24 Sep 2018 09:17:20 GMT
//
//


// system include files
#include <memory>
#include <utility>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"

// Including TrackReco/TrackBase and related files to retrieve the covariance matrix
// and the parameter set and store it in corresponding branches of the 
// output ROOT tree
#include "Rtypes.h"
#include "DataFormats/TrackReco/interface/TrackBase.h"
#include "DataFormats/TrackReco/interface/fillCovariance.h"
#include <algorithm>

// If the analyzer does not use TFileService, please remove
// the template argument to the base class so the class inherits
// from  edm::one::EDAnalyzer<>
// This will improve performance in multithreaded jobs.
#include "FWCore/Utilities/interface/InputTag.h"
#include "FWCore/ServiceRegistry/interface/Service.h"

#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "TMatrixD.h"
#include "TROOT.h"
#include "TTree.h"
#include "TFile.h"

// Adding Message Logging Capabilities
#include "FWCore/MessageLogger/interface/MessageLogger.h"

// Adding DataFormats for tracking Rechits
#include "DataFormats/TrackerRecHit2D/interface/SiPixelRecHitCollection.h"
#include "DataFormats/TrackerRecHit2D/interface/SiStripRecHit2D.h"
#include "DataFormats/TrackerRecHit2D/interface/SiStripRecHit2DCollection.h"

#include "DataFormats/TrackerRecHit2D/interface/SiStripMatchedRecHit2D.h"
#include "DataFormats/TrackerRecHit2D/interface/SiStripMatchedRecHit2DCollection.h"

#include "DataFormats/TrackerRecHit2D/interface/SiPixelRecHitCollection.h"
#include "DataFormats/TrajectorySeed/interface/TrajectorySeedCollection.h"
#include "DataFormats/TrackCandidate/interface/TrackCandidateCollection.h"

#include "TrackingTools/TransientTrackingRecHit/interface/TransientTrackingRecHit.h"
#include "DataFormats/TrackingRecHit/interface/TrackingRecHit.h"

// Include the library file required for localpoint
#include "DataFormats/GeometrySurface/interface/GloballyPositioned.h"

// #include "DataFormats/Common/interface/DetSet.h"

// Include the Dataformat for track association
#include "DataFormats/RecoCandidate/interface/TrackAssociation.h"

// Include data format for simhits
#include "SimDataFormats/TrackingHit/interface/PSimHit.h"
#include "SimDataFormats/TrackingHit/interface/PSimHitContainer.h"

#include "DataFormats/Common/interface/DetSetVector.h"
#include "DataFormats/SiPixelDigi/interface/PixelDigi.h"
#include "SimDataFormats/TrackerDigiSimLink/interface/PixelDigiSimLink.h"
#include "SimDataFormats/TrackerDigiSimLink/interface/StripDigiSimLink.h"

// Tracking Particle Clustering
#include "SimTracker/TrackHistory/interface/TrackClassifier.h"
#include "SimTracker/TrackerHitAssociation/interface/ClusterTPAssociation.h"
#include "SimTracker/TrackerHitAssociation/interface/TrackerHitAssociator.h"

#include "SimDataFormats/Associations/interface/TrackToTrackingParticleAssociator.h"

// Possibly needed for the SiPixelCluster definition
#include "RecoPixelVertexing/PixelLowPtUtilities/interface/ClusterData.h"

// Required for the Tracking Parrticle Association
#include "DataFormats/TrackerRecHit2D/interface/OmniClusterRef.h"

// TODO: Check if this is necessary
// #include "DataFormats/TrackerRecHit2D/interface/Phase2TrackerRecHit1D.h"
// #include "DataFormats/TrackerRecHit2D/interface/TkCloner.h"

// For the Geometry - is this required? Probably at a later point.
/*
#include "Geometry/TrackerGeometryBuilder/interface/TrackerGeometry.h"
#include "Geometry/Records/interface/TrackerDigiGeometryRecord.h"
#include "Geometry/TrackerNumberingBuilder/interface/GeometricDet.h"
#include "Geometry/TrackerGeometryBuilder/interface/PixelGeomDetType.h"
#include "Geometry/TrackerGeometryBuilder/interface/StripGeomDetType.h"
#include "DataFormats/DetId/interface/DetId.h"
#include "DataFormats/GeometryVector/interface/LocalPoint.h"
#include "DataFormats/GeometryVector/interface/GlobalPoint.h"
*/

//
// class declaration
//

class MyTrackingNtuples : public edm::one::EDAnalyzer<edm::one::SharedResources>  {
   public:
      explicit MyTrackingNtuples(const edm::ParameterSet&);
      ~MyTrackingNtuples();
   
   private:
      virtual void beginJob() override;
      virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
      virtual void endJob() override;
//      static bool compare(const std::pair<OmniClusterRef, TrackingParticleRef>&, const std::pair<OmniClusterRef, TrackingParticleRef>&);

      // ----------member data ---------------------------
      int nevent_ = 0;
      int nlumi_ = 0;
      int nrun_ = 0;

      //used to select what tracks to read from configuration file
      edm::EDGetTokenT<edm::View<reco::Track> > tracksToken_; 
      
      //used to select extra track information to read 
      edm::EDGetTokenT<reco::TrackExtraCollection> trackExtraToken_;
      
      edm::Service<TFileService> fs_;
      TTree *tree_;
      
      // Track properties
      std::vector<double> eta_;
      std::vector<double> phi_;
      std::vector<double> qoverp_;
      std::vector<double> dxy_;
      std::vector<double> dsz_;
      std::vector<double> dz_;
      
      std::vector<double> eta_Error_;
      std::vector<double> phi_Error_;
      std::vector<double> qoverp_Error_;
      std::vector<double> dxy_Error_;
      std::vector<double> dsz_Error_;
      std::vector<double> dz_Error_;
      

      // Store the return values of parameter and 
      // covariance matrix functions
      reco::TrackBase::CovarianceMatrix covariance_mat_;
      reco::TrackBase::ParameterVector track_parameters_;
      
      // Handle the Reshaping of the Covariance Matrix
      std::vector<double> reshaped_cov_mat_;
      std::vector< std::vector<double> > covariance_array_;
          
      // Temporary Variables
      std::vector< std::vector<double> > tmpMatrix;
      std::vector<double> tmpVector1;
      std::vector<double> tmpVector2;
      std::vector<double> tmpVector3;
      std::vector<double> tmpVector4;
      std::vector<double> tmpVector5;

      // Declare default buffer size for SiStripRecHit Tracking
      int bufsize = 32000;

      // -------------------- Testing Rechit Retrieval --------------------------
      //edm::EDGetTokenT<bool> doPixel_;
      //edm::EDGetTokenT<bool> doStrip_;
      
      // doPixel_( conf.getParameter<bool>("associatePixel") );
      // doStrip_( conf.getParameter<bool>("associateStrip") );
      
      edm::EDGetTokenT< SiStripMatchedRecHit2DCollection > matchedRecHitToken_;
      edm::EDGetTokenT< SiStripRecHit2DCollection > rphiRecHitToken_;
      edm::EDGetTokenT< SiStripRecHit2DCollection > stereoRecHitToken_;
      edm::EDGetTokenT< SiPixelRecHitCollection > siPixelRecHitsToken_;
      // edm::EDGetTokenT< edm::AssociationMap< edm::OneToManyWithQualityGeneric< edm::View< reco::Track >, TrackingParticleCollection, double > > > association_;

      edm::EDGetTokenT< reco::RecoToSimCollection > association_;
      // edm::EDGetTokenT<edm::DetSetVector<Phase2TrackerRecHit1D> > siPhase2RecHitsToken_;

      // Tracking Particle Association Code
      edm::EDGetTokenT< ClusterTPAssociation > clusterTPMapToken_;
      // TrackClassifier classifier_;

      std::vector<float> stereo_x_;
      std::vector<float> stereo_y_;
      std::vector<float> stereo_z_;

      std::vector<float> rphi_x_;
      std::vector<float> rphi_y_;
      std::vector<float> rphi_z_;

      /*
      std::vector<float> pixel_x_;
      std::vector<float> pixel_y_;
      std::vector<float> pixel_z_;
      */

      std::vector<float> stereo_r;
      std::vector<float> stereo_phi;
      std::vector<float> stereo_eta;

      // TODO: Try using a struct to save time and effort
      struct customEventData {  
          std::vector<double> eta_;
          std::vector<double> phi_;
          std::vector<double> qoverp_;
          std::vector<double> dxy_;
          std::vector<double> dsz_;
          
          reco::TrackBase::CovarianceMatrix covariance_mat_;
          reco::TrackBase::ParameterVector track_parameters_;
        };
};

//
// constants, enums and typedefs
//

//
// static data member definitions
//


//
// constructors and destructor
//
MyTrackingNtuples::MyTrackingNtuples(const edm::ParameterSet& iConfig)
 :
//  src_(iConfig.getParameter<edm::InputTag>( "src" )),
  tracksToken_(consumes< edm::View<reco::Track> >(iConfig.getParameter<edm::InputTag>("pixelTracks"))),
  trackExtraToken_(consumes<reco::TrackExtraCollection>(iConfig.getParameter<edm::InputTag>("pixelTracks"))),
//  rphiRecHits_(consumes<std::vector<SiStripRecHit2D>&>(iConfig.getParameter<edm::InputTag>("rphiRecHits"))),
//  stereoRecHits_(consumes<std::vector<SiStripRecHit2D>&>(iConfig.getParameter<edm::InputTag>("stereoRecHits"))),
//  matchedRecHits_(consumes<std::vector<SiStripRecHit2D>&>(iConfig.getParameter<edm::InputTag>("matchedRecHits"))),
  matchedRecHitToken_(consumes<SiStripMatchedRecHit2DCollection>(iConfig.getParameter<edm::InputTag>("matchedRecHits"))),
  rphiRecHitToken_(consumes<SiStripRecHit2DCollection>(iConfig.getParameter<edm::InputTag>("rphiRecHits"))),   
//  siPhase2RecHitsToken_(consumes<edm::DetSetVector<Phase2TrackerRecHit1D> >(iConfig.getParameter<edm::InputTag>("siPhase2RecHits"))),
  stereoRecHitToken_(consumes<SiStripRecHit2DCollection>(iConfig.getParameter<edm::InputTag>("stereoRecHits"))),
  siPixelRecHitsToken_(consumes<SiPixelRecHitCollection>(iConfig.getParameter<edm::InputTag>("siPixelRecHits"))),
  association_(consumes< reco::RecoToSimCollection >(iConfig.getParameter<edm::InputTag>("associator"))),
  clusterTPMapToken_(consumes< ClusterTPAssociation >(iConfig.getParameter<edm::InputTag>("clusterTPMap")))
//  classifier_(iConfig, consumesCollector())
{
    gROOT->Reset();
    usesResource("TFileService");
    
    //now do what ever initialization is needed
    tree_ = fs_->make<TTree>("tree","tree");
    tree_->Branch("nevent", &nevent_, "nevent/I");
    tree_->Branch("nlumi", &nlumi_, "nlumi/I");
    tree_->Branch("nrun", &nrun_, "nrun/I");
    
    tree_->Branch("trackEta", &eta_);
    tree_->Branch("trackPhi", &phi_);
    tree_->Branch("qoverp", &eta_);
    tree_->Branch("dxy", &dxy_);
    tree_->Branch("dsz", &dsz_);

    tree_->Branch("trackEtaError", &eta_Error_);
    tree_->Branch("trackPhiError", &phi_Error_);
    tree_->Branch("qoverpError", &eta_Error_);
    tree_->Branch("dxyError", &dxy_Error_);
    tree_->Branch("dszError", &dsz_Error_);
    
    tree_->Branch("trackParameters", &track_parameters_);
    tree_->Branch("covarianceArray", &covariance_array_);

    //tree_->Branch("NStereoHits",&NStereoHits_,"NStereoHits/I");
    tree_->Branch("StereoHitX",&stereo_x_);
    tree_->Branch("StereoHitY",&stereo_y_);
    tree_->Branch("StereoHitZ",&stereo_z_);

    tree_->Branch("MonoHitX",&rphi_x_);
    tree_->Branch("MonoHitY",&rphi_y_);
    tree_->Branch("MonoHitZ",&rphi_z_);

    // tree_->Branch("PixelHitX",&pixel_x_);
    // tree_->Branch("PixelHitY",&pixel_y_);
    // tree_->Branch("PixelHitZ",&pixel_z_);
    /* 
    tree_->Branch("StereoHitSigX",&StereoHitSigX_,"StereoHitSigX[1000]/F");
    tree_->Branch("StereoHitSigY",&StereoHitSigY_,"StereoHitSigY[1000]/F");
    tree_->Branch("StereoHitCorr",&StereoHitCorr_,"StereoHitCorr[1000]/F");
 
    tree_->Branch("StereoHitSignal",&StereoHitSignal_,"StereoHitSignal[1000]/F");
    tree_->Branch("StereoHitNoise",&StereoHitNoise_,"StereoHitNoise[1000]/F");
    tree_->Branch("StereoHitWidth",&StereoHitWidth_,"StereoHitWidth[1000]/I");
    */

    // siRphiHitCollection_ = iConfig.getParameter<std::string>("siRphiHitCollection");
    // siStereoHitCollection_ = iConfig.getParameter<std::string>("siStereoHitCollection");
    // siMatchedHitCollection_ = iConfig.getParameter<std::string>("siMatchedHitCollection");

}


MyTrackingNtuples::~MyTrackingNtuples()
{
    // do anything here that needs to be done at destruction time
    // (e.g. close files, deallocate resources etc.)
}


//
// member functions
//

// Custom comparison function for Tracking Particle
// using P = std::pair<OmniClusterRef, TrackingParticleRef>;

/*
static bool MyTrackingNtuples::compare(const std::pair<OmniClusterRef, TrackingParticleRef>& i, const std::pair<OmniClusterRef, TrackingParticleRef>& j) {
    return i.second.index() > j.second.index();
}
*/

// ------------ method called for each event  ------------

void MyTrackingNtuples::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
    using namespace edm;
    using namespace std;
    using namespace reco;

    nevent_ = iEvent.id().event();
    nlumi_ = iEvent.id().luminosityBlock();
    nrun_ = iEvent.id().run();
    
    std::cout << "Event Number: " << nevent_ << std::endl;
    std::cout << "Luminosity Block: " << nlumi_ << std::endl;
    std::cout << "Run Number: " << nrun_ << std::endl;
    
    // Set a counter to check the number of hits
    int numtracks_ = 0;
    
    eta_.clear();
    eta_Error_.clear();
    phi_.clear();
    phi_Error_.clear();
    qoverp_.clear();
    qoverp_Error_.clear();
    dxy_.clear();
    dxy_Error_.clear();
    dsz_.clear();
    dsz_Error_.clear();

    // TODO: Is this clearing the internal vectors 
    // as well as the container?
    covariance_array_.clear();

    // Get the information from the pixeltrack branches
    Handle< edm::View<reco::Track> > tracks_;
    iEvent.getByToken(tracksToken_, tracks_);
    
    edm::Handle<reco::RecoToSimCollection> association;
    iEvent.getByToken(association_, association);    

    // Retrieving Cluster to Tracking Particle Association
    edm::Handle<ClusterTPAssociation> pCluster2TPListH;
    iEvent.getByToken(clusterTPMapToken_, pCluster2TPListH);
    const ClusterTPAssociation& clusterToTPMap_ = *pCluster2TPListH;

    for(size_t track_idx=0; track_idx<tracks_->size(); ++track_idx) {
        
        edm::RefToBase<reco::Track> trk_(tracks_, track_idx);
        
          //find the tracking particle based on the track
          auto gen_match = association->find(trk_);
          if(gen_match != association->end()) {
            auto tracking_particle = gen_match->val.front().first;
            std::cout << "TP PDG ID: " << tracking_particle->pdgId() << std::endl;
            std::cout << "TP Index: " << tracking_particle.index() << std::endl;
            std::cout << "TP Number of Hits: " << tracking_particle->numberOfHits() << std::endl;
            
            
            /*
            // Incorrectly using equal_range - don't need to do this
            // Instead just pass the clusterref to equal_range and pick the
            // tracking particle in the first element as the reference for the rechit
            auto clusterRange_ = clusterToTPMap.equal_range(// add OmniClusterRef here); 
            for(auto clusterIterator_ = clusterRange_.first; clusterIterator_ != clusterRange_.second;
                    clusterIterator_++) {
                customTPMap
            }
            */

            // Custom methods for Track Association

            /*auto clusterTPmap = clusterToTPMap.map();
            std::sort(clusterTPmap.begin(), clusterTPmap.end(), compare);
            auto clusterRange = std::equal_range(clusterTPmap.begin(), clusterTPmap.end(),std::make_pair(OmniClusterRef(), tracking_particle), compare);      
           
            int TrkTruthnHitAll=-1;
            int TrkTruthnHitPixel=-1;
            int TrkTruthnHitStrip=-1;
       
            if( clusterRange.first != clusterRange.second ) {
                for( auto ip=clusterRange.first; ip != clusterRange.second; ++ip ) {
                    const OmniClusterRef& cluster = ip->first;
                    if (cluster.isPixel() && cluster.isValid()){ TrkTruthnHitPixel+=1;}
                    if (cluster.isStrip() && cluster.isValid()){ TrkTruthnHitStrip+=1;}
            }
             TrkTruthnHitAll = TrkTruthnHitPixel + TrkTruthnHitStrip;
             std::cout << "All hits: " << TrkTruthnHitAll << std::endl;
            */

		  } else { //no matching
              std::cout << "Match not found" << std::endl;
            //TODO
          }
        
        // Add all the Track parameters and corresponding errors to the vectors
        // to be put into the TTree
        eta_.push_back(trk_->eta());
        eta_Error_.push_back(trk_->etaError());
        phi_.push_back(trk_->phi());
        phi_Error_.push_back(trk_->phi());
        qoverp_.push_back(trk_->qoverp());
        dsz_.push_back(trk_->dsz());
        dsz_Error_.push_back(trk_->dszError());
        dxy_.push_back(trk_->dxy());
        dxy_Error_.push_back(trk_->dxyError());
        dz_.push_back(trk_->dz());
        dz_Error_.push_back(trk_->dzError());
        // std::cout << "Jet Data: " << trk_->eta() << trk_->phi() << trk_->qoverp() << trk_->dxy() << trk_->dsz() << std::endl;
        
        covariance_mat_ = trk_->covariance();
        track_parameters_ = trk_->parameters();

        // TODO #1: Figure out how to store covariance
        // matrix in the TTree - Ntuples?
        // Edit: DONE -> std::vector< double >

        std::cout << "Track Parameters: ";
        // Print the collected parameters from the parameter set
        for (int i_ = 0; i_ < track_parameters_.kSize; i_++) {
            std::cout << track_parameters_.At(i_);
        }
        std::cout << std::endl;

        // Reset counter variable
        // i_ = 0;

        // Print the covariance matrix (assume 5 x 5), eliminate the upper triangular half
        // and reshape it to store it in a fixed-dimension vector of doubles
        for (int i_ = 0; i_ < 5; i_++) {
          for (int j_ = 0; j_ <= i_; j_++) {
            // std::cout << "i: " << i_ << " j: " << j_ << std::endl;
            // std::cout << covariance_mat_[i_][j_] << " | ";
            reshaped_cov_mat_.push_back(covariance_mat_[i_][j_]);
          }
          // std::cout << std::endl;
        }

        // std::cout << "Reshaped Covariance Matrix: " << reshaped_cov_mat_.size() << std::endl;
        for (int i_ = 0; (unsigned)i_ < reshaped_cov_mat_.size(); i_++) {
            std::cout << reshaped_cov_mat_.at(i_) << " ";
        }
        std::cout << std::endl;

        // Push the reshaped vector into a container of multiple vectors
        // which will contain all the track covariances for an event
        covariance_array_.push_back(reshaped_cov_mat_);
        reshaped_cov_mat_.clear();
        
        // std::cout << "Track Covariance and Parameter Set found" << std::endl;

        numtracks_ ++;
        LogInfo("Tracks") << "Found " << numtracks_ << " tracks" << std::endl;;       
        // do something with track parameters, e.g, plot the charge.
    }
    
    // Reset number of tracks to zero
    numtracks_ = 0;
    
    // Get the extra information from the pixeltrack branches
    // We haven't really done anything with this though
    Handle<reco::TrackExtraCollection> trackExtra_;
    iEvent.getByToken(trackExtraToken_, trackExtra_);
            
    /*  for(reco::TrackExtraCollection::const_iterator itTrackExtra_ = trackExtra_->begin();
        itTrackExtra_ != trackExtra_->end();
        ++itTrackExtra_) {
        numtracks_ ++;
        LogInfo("TrackExtra") << "Found extra info on " << numtracks_ << " tracks" << std::endl;
    }*/

    // std::cout << "Run in TrackExtra section " << numtracks_ << "times." << std::endl;

    // TODO #2: Print extra information about the tracks

// ----------------------- Get the RecHits from the Data -------------------

    // int pixelcounter_ = 0;
    // int stripcounter_ = 0;
    
    // Find the collections of rechits and pixelhits
    // TODO: We don't need pixel hits since we already have the tracks
    // Remove the last part corresponding to pixel hits
    edm::Handle<SiPixelRecHitCollection> pixelrechitColl_;
    edm::Handle<SiStripRecHit2DCollection> rphirechitColl_;
    edm::Handle<SiStripRecHit2DCollection> stereorechitColl_;
    edm::Handle<SiStripMatchedRecHit2DCollection> matchedrechitColl_;

    iEvent.getByToken(siPixelRecHitsToken_, pixelrechitColl_);
    iEvent.getByToken(rphiRecHitToken_, rphirechitColl_);
    iEvent.getByToken(stereoRecHitToken_, stereorechitColl_);
    iEvent.getByToken(matchedRecHitToken_, matchedrechitColl_);
    
    // TODO: Write a function that initializes the value of all of the variables 
    // initNtuple();

    // Print size of rphirechits
    std::cout << "RPhiRecHitColl Data Size: " << (rphirechitColl_.product())->dataSize() << std::endl;
    std::cout << "StereoRecHitColl Data Size: " << (stereorechitColl_.product())->dataSize() << std::endl;
    std::cout << "PixelRecHitColl Data Size: " << (pixelrechitColl_.product())->dataSize() << std::endl;

    // Iterate over rphirechits and check if the begin/end methods work
    /*for(std::vector<SiStripRecHit2D>::const_iterator hiter = rphirechits_.begin();
        hiter != rphirechits_.end();
        ++hiter) {
      std::cout << ".";
    }
    std::cout << std::endl;*/

    // Different approach to iterating over the stereorechits
    if((stereorechitColl_.product())->dataSize() > 0) {
        SiStripRecHit2DCollection::const_iterator stereorecHitIdIterator = (stereorechitColl_.product())->begin();
        SiStripRecHit2DCollection::const_iterator stereorecHitIdIteratorEnd = (stereorechitColl_.product())->end();
      
        for(SiStripRecHit2DCollection::const_iterator stereo_detunit_iterator_ = stereorecHitIdIterator;
            stereo_detunit_iterator_ != stereorecHitIdIteratorEnd; stereo_detunit_iterator_++) {
        
            SiStripRecHit2DCollection::DetSet stereo_detset_ = *stereo_detunit_iterator_;
            SiStripRecHit2DCollection::DetSet::const_iterator stereorechitRangeIteratorBegin = stereo_detset_.begin();
            SiStripRecHit2DCollection::DetSet::const_iterator stereorechitRangeIteratorEnd = stereo_detset_.end();
            SiStripRecHit2DCollection::DetSet::const_iterator stereo_iterRecHit_;
            
            for (stereo_iterRecHit_ = stereorechitRangeIteratorBegin; 
                  stereo_iterRecHit_ != stereorechitRangeIteratorEnd; ++stereo_iterRecHit_) {

                auto stereo_cluster_ = stereo_iterRecHit_->firstClusterRef();
                // OmniClusterRef::ClusterStripRef stereo_cluster_ = stereo_iterRecHit_->stripCluster();
                // auto clusterRange_ = clusterToTPMap_.equal_range(stereo_cluster_);
                auto clusterTPMapIter_ = clusterToTPMap_.equal_range(stereo_cluster_);
                for (auto iter__ = clusterTPMapIter_.first; iter__ != clusterTPMapIter_.second;
                        iter__++){
                    std::cout << "Iterating through cluster mapping" << std::endl;
                }
                // Obtain the local position in terms of coordinates and store it in the vector
                LocalPoint stereo_lp = stereo_iterRecHit_->localPosition();
                stereo_x_.push_back(stereo_lp.x());
                stereo_y_.push_back(stereo_lp.y());
                stereo_z_.push_back(stereo_lp.z());        
            }
        }
    }


    // Different approach to iterating over the rphirechits/monorechits
    if((rphirechitColl_.product())->dataSize() > 0) {
        SiStripRecHit2DCollection::const_iterator rphirecHitIdIterator = (rphirechitColl_.product())->begin();
        SiStripRecHit2DCollection::const_iterator rphirecHitIdIteratorEnd = (rphirechitColl_.product())->end();

        for(SiStripRecHit2DCollection::const_iterator rphi_detunit_iterator_ = rphirecHitIdIterator;
            rphi_detunit_iterator_ != rphirecHitIdIteratorEnd; rphi_detunit_iterator_++) {
        
            SiStripRecHit2DCollection::DetSet rphi_detset = *rphi_detunit_iterator_;

            SiStripRecHit2DCollection::DetSet::const_iterator rechitRangeIteratorBegin = rphi_detset.begin();
            SiStripRecHit2DCollection::DetSet::const_iterator rechitRangeIteratorEnd   = rphi_detset.end();
            SiStripRecHit2DCollection::DetSet::const_iterator rphi_iterRecHit_;
        
            for (rphi_iterRecHit_ = rechitRangeIteratorBegin; 
                rphi_iterRecHit_ != rechitRangeIteratorEnd; ++rphi_iterRecHit_) {
          
                // Test get rechit cluster
                // SiPixelRecHit::ClusterRef const& clust = iterRecHit_.cluster();
		        // std::cout << "Rechit Stripcluster: " << recHit.stripCluster() << std::endl;
		        const SiStripCluster& rphi_cluster_ = rphi_iterRecHit_->stripCluster();
                
	            LocalPoint rphi_lp = rphi_iterRecHit_->localPosition();
                rphi_x_.push_back(rphi_lp.x());
                rphi_y_.push_back(rphi_lp.y());
                rphi_z_.push_back(rphi_lp.z());        
            }
        }
    }
      
    /*
    const SiPixelRecHitCollection *hits = pixelrechitColl_.product();

    for(SiPixelRecHitCollection::DataContainer::const_iterator hit = hits->data().begin(),
        end = hits->data().end(); hit != end; ++hit) {

      // Is this a different (x, y) location ?

      // std::vector< SiPixelCluster::Pixel > pixels(hit->cluster()->pixels());
      // bool pixelOnEdge = false;
      // for(std::vector< SiPixelCluster::Pixel >::const_iterator pixel = pixels.begin();
        // pixel != pixels.end(); ++pixel) {
        // int pixelX = pixel->x;
        // int pixelY = pixel->y;
        
        // Here lay code to check if it is an edge pixel
        // }
	  	SiPixelRecHit::ClusterRef const& clust = hit->cluster();

	    LocalPoint lp = hit->localPosition();
        pixel_x.push_back(lp.x());
        pixel_y.push_back(lp.y());
        pixel_z.push_back(lp.z());
    }
    */
  // ------------------------------ Fill and Print the Tree -------------------------

  tree_->Fill();
  
  std::cout << "Tree filled" << std::endl;
  // tree_->Print();
}


// ------------ method called once each job just before starting event loop  ------------
void
MyTrackingNtuples::beginJob()
{
  // Defined a custom check for file pointer to store root ntuple
  if( !fs_ ){
        throw edm::Exception( edm::errors::Configuration,
                "TFile Service is not registered in cfg file" );
    }
}

// ------------ method called once each job just after ending the event loop  ------------
void
MyTrackingNtuples::endJob()
{
    // fs_->Write();
    
    // The Covariance Matrix and Track Parameters are cleared by gROOT->Reset()
    tmpMatrix.swap(covariance_array_);

    /*
    tmpVector1.swap(eta_);
    tmpVector2.swap(phi_);
    tmpVector3.swap(qoverp_);
    tmpVector4.swap(dxy_);
    tmpVector5.swap(dsz_);
    
    TVectorD tmpVector1;
    tmpVector1.swap(eta_);
    TVectorD tmpVector2;
    tmpVector2.swap(phi_);
    TVectorD tmpVector3;
    tmpVector3.swap(qoverp_);
    */
    
    //std::vector<covariance_mat_>.swap(covariance_mat_);
    //std::vector<track_parameters_>.swap(track_parameters_);

    std::cout << ">> Ending job." << std::endl;
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------

//define this as a plug-in
DEFINE_FWK_MODULE(MyTrackingNtuples);
