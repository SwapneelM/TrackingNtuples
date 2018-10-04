// -*- C++ -*-
//
// Package:    TrackingNtuples/TrackingNtuples
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

// If the analyzer does not use TFileService, please remove
// the template argument to the base class so the class inherits
// from  edm::one::EDAnalyzer<>
// This will improve performance in multithreaded jobs.
#include "FWCore/Utilities/interface/InputTag.h"
#include "FWCore/ServiceRegistry/interface/Service.h"

#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "TTree.h"
#include "TFile.h"
//
// Adding Message Logging Capabilities
//
#include "FWCore/MessageLogger/interface/MessageLogger.h"

//
// class declaration
//

class MyTrackingNtuples : public edm::one::EDAnalyzer<edm::one::SharedResources>  {
   public:
      explicit MyTrackingNtuples(const edm::ParameterSet&){
          usesResource("TFileService");
      };
      ~MyTrackingNtuples();
   
   private:
      virtual void beginJob() override;
      virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
      virtual void endJob() override;

      // ----------member data ---------------------------
      int nevent_ = 0;
      int nlumi_ = 0;
      int nrun_ = 0;

      edm::EDGetTokenT<reco::TrackCollection> tracksToken_;  //used to select what tracks to read from configuration file
      edm::EDGetTokenT<reco::TrackExtraCollection> trackExtraToken_;  //used to select extra track information to read 
      
      edm::Service<TFileService> fs_;
      TTree *tree_;
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
  tree_(nullptr),
  tracksToken_(consumes<reco::TrackCollection>(iConfig.getParameter<edm::InputTag>("pixelTracks"))),
  trackExtraToken_(consumes<reco::TrackExtraCollection>(iConfig.getParameter<edm::InputTag>("pixelTracks")))
{
   //now do what ever initialization is needed
  tree_ = fs->make<TTree>("tree","tree");
  tree_->Branch("nevent",&nevent_,"nevent/I");
  tree_->Branch("nlumi",&nlumi_,"nlumi/I");
  tree_->Branch("nrun",&nrun_,"nrun/I");

}


MyTrackingNtuples::~MyTrackingNtuples()
{

   // do anything here that needs to be done at destruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called for each event  ------------
void
MyTrackingNtuples::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
    using namespace edm;
    using namespace std;
    using namespace reco;

    nevent_=iEvent.id().event();
    nlumi_=iEvent.id().luminosityBlock();
    nrun_=iEvent.id().run();

    // Set a counter to check the number of hits
    int numhits_ = 0;
    
    Handle<reco::TrackCollection> tracks_;
    iEvent.getByToken(tracksToken_, tracks_);

    for(reco::TrackCollection::const_iterator itTrack_ = tracks_->begin();
        itTrack_ != tracks_->end(); 
        ++itTrack_) {
        numhits_ ++;
        LogDebug("TrackFinder") << "Found tracks: " << numhits_ << " hits\n";       // do something with track parameters, e.g, plot the charge.
    }

    // Reset number of hits to zero
    numhits_ = 0;
    
    Handle<reco::TrackExtraCollection> trackExtra_;
    iEvent.getByToken(trackExtraToken_, trackExtra_);

    for(reco::TrackExtraCollection::const_iterator itTrackExtra_ = trackExtra_->begin();
        itTrackExtra_ != trackExtra_->end();
        ++itTrackExtra_) {
        numhits_ ++;
        LogDebug("TrackExtraFinder") << "Found extra info on " << numhits_ << " tracks\n";
    }

    tree_->Fill();

    /* #ifdef THIS_IS_AN_EVENT_EXAMPLE
       Handle<ExampleData> pIn;
       iEvent.getByLabel("example",pIn);
    #endif

    #ifdef THIS_IS_AN_EVENTSETUP_EXAMPLE
       ESHandle<SetupData> pSetup;
       iSetup.get<SetupRecord>().get(pSetup);
    #endif */
}


// ------------ method called once each job just before starting event loop  ------------
void
MyTrackingNtuples::beginJob()
{
  if( !fs_ ){
        throw edm::Exception( edm::errors::Configuration,
                "TFile Service is not registered in cfg file" );
    }
}

// ------------ method called once each job just after ending the event loop  ------------
void
MyTrackingNtuples::endJob()
{
    tree_->Write();
    std::cout << ">> Ending job." << std::endl;
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------

//define this as a plug-in
DEFINE_FWK_MODULE(MyTrackingNtuples);
