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

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
//
// class declaration
//

// If the analyzer does not use TFileService, please remove
// the template argument to the base class so the class inherits
// from  edm::one::EDAnalyzer<>
// This will improve performance in multithreaded jobs.

//
// Custom: Adding Message Logging Capabilities
//
#include "FWCore/MessageLogger/interface/MessageLogger.h"

class TrackingNtuples : public edm::one::EDAnalyzer<edm::one::SharedResources>  {
   public:
      explicit TrackingNtuples(const edm::ParameterSet&);
      ~TrackingNtuples();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


   private:
      virtual void beginJob() override;
      virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
      virtual void endJob() override;

      // ----------member data ---------------------------
      edm::EDGetTokenT<reco::TrackCollection> tracksToken_;  //used to select what tracks to read from configuration file
      edm::EDGetTokenT<reco::TrackExtraCollection> trackExtraToken_;  //used to select extra track information to read 
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
TrackingNtuples::TrackingNtuples(const edm::ParameterSet& iConfig)
 :
  tracksToken_(consumes<reco::TrackCollection>(iConfig.getParameter<edm::InputTag>("pixelTracks"))),
  trackExtraToken_(consumes<reco::TrackExtraCollection>(iConfig.getParameter<edm::InputTag>("pixelTracks")))

{
   //now do what ever initialization is needed

}


TrackingNtuples::~TrackingNtuples()
{

   // do anything here that needs to be done at destruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called for each event  ------------
void
TrackingNtuples::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   using namespace edm;
    int numhits = 0;
    Handle<reco::TrackCollection> tracks;
    iEvent.getByToken(tracksToken_, tracks);
    for(reco::TrackCollection::const_iterator itTrack = tracks->begin();
        itTrack != tracks->end();
        ++itTrack) {
        numhits ++;
        LogDebug("TrackFinder") << "Found tracks: " << numhits << " hits\n";       // do something with track parameters, e.g, plot the charge.
      // int charge = itTrack->charge();
    }
    numhits = 0;
    Handle<reco::TrackExtraCollection> trackExtra;
    iEvent.getByToken(trackExtraToken_, trackExtra);
    for(reco::TrackExtraCollection::const_iterator itTrack = trackExtra->begin();
        itTrack != trackExtra->end();
        ++itTrack) {
        numhits ++;
        LogDebug("TrackFinder") << "Found extra info on " << numhits << " tracks\n";
    }

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
TrackingNtuples::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void
TrackingNtuples::endJob()
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
TrackingNtuples::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.addUntracked<edm::InputTag>("pixelTracks");
  descriptions.addDefault(desc);

  //Specify that only 'tracks' is allowed
  //To use, remove the default given above and uncomment below
  //ParameterSetDescription desc;
  //desc.addUntracked<edm::InputTag>("tracks","ctfWithMaterialTracks");
  //descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(TrackingNtuples);
