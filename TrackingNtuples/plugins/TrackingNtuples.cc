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
      explicit MyTrackingNtuples(const edm::ParameterSet&);
      ~MyTrackingNtuples();

   private:
      virtual void beginJob() override;
      virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
      virtual void endJob() override;

      // ----------member data ---------------------------
      int nevent, nlumi, nrun;

      edm::EDGetTokenT<reco::TrackCollection> tracksToken_;  //used to select what tracks to read from configuration file
      edm::EDGetTokenT<reco::TrackExtraCollection> trackExtraToken_;  //used to select extra track information to read 
      
      edm::Service<TFileService> fs;
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
  tree(nullptr), tracksToken_(consumes<reco::TrackCollection>(iConfig.getParameter<edm::InputTag>("pixelTracks"))),
  trackExtraToken_(consumes<reco::TrackExtraCollection>(iConfig.getParameter<edm::InputTag>("pixelTracks")))

{
   //now do what ever initialization is needed
  tree = fs->make<TTree>("tree","tree");
  tree->Branch("nevent",&nevent,"nevent/I");
  tree->Branch("nlumi",&nlumi,"nlumi/I");
  tree->Branch("nrun",&nrun,"nrun/I");

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

    nevent=iEvent.id().event();
    nlumi=iEvent.id().luminosityBlock();
    nrun=iEvent.id().run();

    // Set a counter to check the number of hits
    int numhits = 0;
    
    Handle<reco::TrackCollection> tracks;
    iEvent.getByToken(tracksToken_, tracks);
    for(reco::TrackCollection::const_iterator itTrack = tracks->begin();
        itTrack != tracks->end();
        ++itTrack) {
        numhits ++;
        LogDebug("TrackFinder") << "Found tracks: " << numhits << " hits\n";       // do something with track parameters, e.g, plot the charge.
    }

    // Reset number of hits to zero
    numhits = 0;
    
    Handle<reco::TrackExtraCollection> trackExtra;
    iEvent.getByToken(trackExtraToken_, trackExtra);
    for(reco::TrackExtraCollection::const_iterator itTrack = trackExtra->begin();
        itTrack != trackExtra->end();
        ++itTrack) {
        numhits ++;
        LogDebug("TrackFinder") << "Found extra info on " << numhits << " tracks\n";
    }

    tree=fs->make<TTree>("tree","tree");
    tree->Branch("nevent",&nevent,"nevent/I");
    tree->Branch("nlumi",&nlumi,"nlumi/I");
    tree->Branch("nrun",&nrun,"nrun/I");

    tree->Fill();

    TFile outputFile ("outputFile.root","RECREATE");
    tree.Write()
    outputFile.ls();
    outputFile.Close();
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
  if( !fs ){
        throw edm::Exception( edm::errors::Configuration,
                "TFile Service is not registered in cfg file" );
    }
}

// ------------ method called once each job just after ending the event loop  ------------
void
MyTrackingNtuples::endJob()
{
    cout << ">> Ending job." << endl;
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------

//define this as a plug-in
DEFINE_FWK_MODULE(MyTrackingNtuples);
