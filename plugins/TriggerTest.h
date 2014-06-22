// -*- C++ -*-
//
// Package:    hugues/TriggerTest
// Class:      TriggerTest
// 
/**\class TriggerTest TriggerTest.cc hugues/TriggerTest/plugins/TriggerTest.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Hugues Louis Brun
//         Created:  Mon, 09 Jun 2014 15:02:11 GMT
//
//


// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"


#include "SimDataFormats/GeneratorProducts/interface/HepMCProduct.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"
#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"
#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"
#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"

#include "HLTrigger/HLTcore/interface/HLTConfigProvider.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/HLTReco/interface/TriggerEventWithRefs.h"


#include "DataFormats/HLTReco/interface/TriggerEvent.h"
#include "DataFormats/HLTReco/interface/TriggerObject.h"

#include "DataFormats/RecoCandidate/interface/RecoEcalCandidate.h"
#include "DataFormats/RecoCandidate/interface/RecoEcalCandidateFwd.h"

#include "DataFormats/EgammaCandidates/interface/ElectronIsolationAssociation.h"
#include "DataFormats/RecoCandidate/interface/RecoChargedCandidateFwd.h"
#include "DataFormats/RecoCandidate/interface/RecoChargedCandidate.h"
#include "DataFormats/Candidate/interface/CandMatchMap.h"
#include "DataFormats/Candidate/interface/CandidateFwd.h"
#include "DataFormats/Common/interface/ValueMap.h"
#include "DataFormats/EgammaCandidates/interface/ElectronIsolationAssociation.h"

#include "DataFormats/EgammaCandidates/interface/Electron.h"
#include "DataFormats/EgammaCandidates/interface/ElectronFwd.h"

#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectronFwd.h"
#include "DataFormats/EgammaReco/interface/ElectronSeedFwd.h"
#include "DataFormats/EgammaReco/interface/ElectronSeed.h"

#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonFwd.h"


// root stuff !
#include "TH1D.h"
#include <map>
#include "TFile.h"
#include <math.h>
#include <sstream>
#include <string>
#include <stdlib.h>
#include <string.h>
#include "TLorentzVector.h"
#include "TVector3.h"
#include "TString.h"
#include "TTree.h"


//
// class declaration
//

class TriggerTest : public edm::EDAnalyzer {
   public:
      explicit TriggerTest(const edm::ParameterSet&);
      ~TriggerTest();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


   private:
      virtual void beginJob() override;
      virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
      virtual void endJob() override;
      virtual void beginEvent();
      virtual void endEvent();
      virtual bool hasWasMother(const reco::GenParticle);
      virtual bool hasTauasMother(const reco::GenParticle);

    
    
    std::vector<int> triggerBits_;
    
      //virtual void beginRun(edm::Run const&, edm::EventSetup const&) override;
      //virtual void endRun(edm::Run const&, edm::EventSetup const&) override;
      //virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;
      //virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;

      // ----------member data ---------------------------
    HLTConfigProvider hltConfig;
    
    //edm::InputTag electronsCollection_;
    //edm::InputTag muonsCollection_;

    edm::InputTag triggerResultsTag_;
    edm::InputTag triggerSummaryLabel_;
    
    std::vector<std::string> pathsToSave_;
    std::vector<std::string> filterToMatch_;
    std::vector<std::string> mapsValues_;
    
    std::string HLTprocess_;
    std::string outputFile_; // output file
    
    // ---------- output ROOT file
    TFile*  rootFile_;
    
    // ---------- tree declaration
    TTree *mytree_;
    
    // -----------tree variables
    //Events
    int T_Event_RunNumber;
    int T_Event_EventNumber;
    int T_Event_LuminosityBlock;
    
    std::vector<int> *T_Event_pathsFired;
    
    
    std::vector<float> *T_Trig_Eta;
    std::vector<float> *T_Trig_Pt;
    std::vector<float> *T_Trig_Phi;
    std::vector<int> *T_Trig_Leg;
    std::vector<float> *T_Trig_Value;
    std::vector<float> *T_Trig_Value2;
    
    std::vector<float> *T_Gen_Eta;
    std::vector<float> *T_Gen_Phi;
    std::vector<float> *T_Gen_Pt;
    std::vector<int> *T_Gen_pdgID;
    std::vector<int> *T_Gen_MotherID;
    std::vector<int> *T_Gen_FromW;
    std::vector<int> *T_Gen_FromTau;
    
    
   /* std::vector<float> *T_Elec_Eta;
    std::vector<float> *T_Elec_Phi;
    std::vector<float> *T_Elec_Pt;

    
    std::vector<float> *T_Muon_Eta;
    std::vector<float> *T_Muon_Phi;
    std::vector<float> *T_Muon_Pt;*/
};


typedef std::vector< edm::Handle<trigger::TriggerFilterObjectWithRefs> > TrigFiltVect;
typedef edm::AssociationMap<edm::OneToValue<std::vector<reco::RecoEcalCandidate>,float,unsigned int> > recoEcalCandidateMap;



