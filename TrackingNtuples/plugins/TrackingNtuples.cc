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
#include "TROOT.h"
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
      int nevent_ = 0;
      int nlumi_ = 0;
      int nrun_ = 0;

      edm::EDGetTokenT<reco::TrackCollection> tracksToken_;  //used to select what tracks to read from configuration file
      edm::EDGetTokenT<reco::TrackExtraCollection> trackExtraToken_;  //used to select extra track information to read 
      
      edm::Service<TFileService> fs_;
      TTree *tree_;

      std::vector<double> jet_eta_;
      std::vector<double> jet_phi_;
      std::vector<double> qoverp_;
      reco::TrackBase::CovarianceMatrix covariance_mat_;
      reco::TrackBase::ParameterVector track_parameters_;

      struct customEventData {  
        std::vector<double> jet_eta_;
        std::vector<double> jet_phi_;
        std::vector<double> qoverp_;
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
  tracksToken_(consumes<reco::TrackCollection>(iConfig.getParameter<edm::InputTag>("pixelTracks"))),
  trackExtraToken_(consumes<reco::TrackExtraCollection>(iConfig.getParameter<edm::InputTag>("pixelTracks")))
{
    gROOT->Reset();
    usesResource("TFileService");
    
    //now do what ever initialization is needed
    tree_ = fs_->make<TTree>("tree","tree");
    tree_->Branch("nevent", &nevent_, "nevent/I");
    tree_->Branch("nlumi", &nlumi_, "nlumi/I");
    tree_->Branch("nrun", &nrun_, "nrun/I");
    
    tree_->Branch("jeteta", &jet_eta_);
    tree_->Branch("jetphi", &jet_phi_);
    tree_->Branch("qoverp", &jet_eta_);
    //tree_->Branch("covarianceMatrix", &covariance_mat_);
    tree_->Branch("trackparameters", &track_parameters_);
}


MyTrackingNtuples::~MyTrackingNtuples()
{
    // do anything here that needs to be done at destruction time
    // (e.g. close files, deallocate resources etc.)
    tree_->Write();
    tree_->ResetBranchAddresses();

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

    std::cout << "Event Number: " << nevent_ << std::endl;
    std::cout << "Luminosity Block: " << nlumi_ << std::endl;
    std::cout << "Run Number: " << nrun_ << std::endl;
    
    // Set a counter to check the number of hits
    int numtracks_ = 0;
    
    jet_eta_.clear();
    jet_phi_.clear();
    qoverp_.clear();

    // Get the information from the pixeltrack branches
    Handle<reco::TrackCollection> tracks_;
    iEvent.getByToken(tracksToken_, tracks_);

    for(reco::TrackCollection::const_iterator itTrack_ = tracks_->begin();
        itTrack_ != tracks_->end(); 
        ++itTrack_) {

        reco::Track trk_ = *itTrack_;
        
        jet_eta_.push_back(trk_.eta());
        jet_phi_.push_back(trk_.phi());
        qoverp_.push_back(trk_.qoverp());
        std::cout << 'Jet Data: ' << trk_.eta() << trk_.phi() << trk_.qoverp_() << endl;
        // covariance_mat_ = trk_.covariance();
        track_parameters_ = trk_.parameters();

        // TODO #1: Print the covariance matrix and 
        // figure out how to store it in the TTree

        // Print the collected parameters from the parameter set
        for (int i_=0; i_ < track_parameters_.kSize; i_++) {
            std::cout << 'Track Param: ' << track_parameters_.At(i_) << std::endl;
        }
        std::cout << "Track Covariance and Parameter Set found ?" << std::endl;
        //std::cout << "Track Phi" << trk_.phi();

        numtracks_ ++;
        LogInfo("Tracks") << "Found " << numtracks_ << " tracks\n";       // do something with track parameters, e.g, plot the charge.
    }
    
    // Reset number of hits to zero
    numtracks_ = 0;
    
    // Get the extra information from the pixeltrack branches
    Handle<reco::TrackExtraCollection> trackExtra_;
    iEvent.getByToken(trackExtraToken_, trackExtra_);

    for(reco::TrackExtraCollection::const_iterator itTrackExtra_ = trackExtra_->begin();
        itTrackExtra_ != trackExtra_->end();
        ++itTrackExtra_) {
        numtracks_ ++;
        LogInfo("TrackExtra") << "Found extra info on " << numtracks_ << " tracks\n";
    }

    std::cout << "Run in TrackExtra section " << numtracks_ << "times." << std::endl;

    // TODO #2: Print extra information about the tracks
    
    tree_->Fill();
    tree_->Print();
    
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
    
    // The Covariance Matrix and Track Parameters are cleared by gROOT->Reset()
    // covariance_mat_.clear();
    //track_parameters_.clear();

    std::vector<double> tmpVector1;
    tmpVector1.swap(jet_eta_);
    std::vector<double> tmpVector2;
    tmpVector2.swap(jet_phi_);
    std::vector<double> tmpVector3;
    tmpVector3.swap(qoverp_);
    //std::vector<covariance_mat_>.swap(covariance_mat_);
    //std::vector<track_parameters_>.swap(track_parameters_);

    std::cout << ">> Ending job." << std::endl;
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------

//define this as a plug-in
DEFINE_FWK_MODULE(MyTrackingNtuples);
