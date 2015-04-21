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
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "FWCore/Framework/interface/ConsumesCollector.h"

#include "SimDataFormats/GeneratorProducts/interface/HepMCProduct.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"
#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"
#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"
#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"

#include "HLTrigger/HLTcore/interface/HLTConfigProvider.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/HLTReco/interface/TriggerEventWithRefs.h"
#include "DataFormats/HLTReco/interface/TriggerRefsCollections.h"

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

#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"

#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonFwd.h"
#include "DataFormats/MuonReco/interface/MuonSelectors.h"

#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"
#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"

#include "DataFormats/RecoCandidate/interface/RecoChargedCandidate.h"
#include "DataFormats/RecoCandidate/interface/RecoChargedCandidateFwd.h"

#include "DataFormats/RecoCandidate/interface/IsoDeposit.h"
#include "DataFormats/RecoCandidate/interface/IsoDepositFwd.h"

#include "DataFormats/RecoCandidate/interface/RecoChargedCandidateIsolation.h"

// extract isolation deposit
#include "PhysicsTools/IsolationAlgos/interface/IsoDepositExtractor.h"
#include "PhysicsTools/IsolationAlgos/interface/IsoDepositExtractor.h"
#include "PhysicsTools/IsolationAlgos/interface/IsoDepositExtractorFactory.h"

#include "DataFormats/TrackReco/interface/Track.h"
#include "SimDataFormats/TrackingAnalysis/interface/TrackingParticle.h"



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
typedef std::vector<edm::InputTag> vtag;
typedef edm::AssociationMap<edm::OneToValue<std::vector<reco::Muon>, float > > MuonIsolationMap;
typedef std::vector< edm::Handle< double > >   rhoHandles;


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
    reco::isodeposit::IsoDepositExtractor * trkExtractor;
   // reco::isodeposit::IsoDepositExtractor * caloExtractor;
    
    //edm::InputTag electronsCollection_;
    //edm::InputTag muonsCollection_;
    
    bool doMuon_;

    vtag muonProducers_;

    edm::InputTag primaryVertexInputTag_;
    edm::InputTag triggerResultsTag_;
    edm::InputTag triggerSummaryLabel_;
    
    std::vector<std::string> pathsToSave_;
    std::vector<std::string> filterToMatch_;
    std::vector<std::string> mapsValues_;
    edm::InputTag rhoTag_;
    
    edm::InputTag muonECALpfIsoTag_;
    edm::InputTag muonHCALpfIsoTag_;
    
    std::vector<edm::InputTag> rhoInputTags_;

    
    std::string HLTprocess_;
    std::string outputFile_; // output file
    
    edm::ParameterSet trkExtractorPSet_;
   // edm::ParameterSet caloExtractorPSet_;
    
    // ---------- output ROOT file
    TFile*  rootFile_;
    
    // ---------- tree declaration
    TTree *mytree_;
    
    // -----------tree variables
    //Events
    int T_Event_RunNumber;
    int T_Event_EventNumber;
    int T_Event_LuminosityBlock;
    
    int T_Event_nPU;
    float T_Event_nTruePU;
    int T_Event_nPUm;
    int T_Event_nPUp;
    float T_Event_AveNTruePU;

    std::vector<float> * T_Event_Rho;

    
    int antiMuEnrichementVeto;
    
    std::vector<int> *T_Event_pathsFired;
    
    
    std::vector<float> *T_Trig_Eta;
    std::vector<float> *T_Trig_Pt;
    std::vector<float> *T_Trig_Phi;
    std::vector<int> *T_Trig_Leg;
    std::vector<float> *T_Trig_TkIsoVeto;
    std::vector<float> *T_Trig_trackerIso;
    std::vector<float> *T_Trig_detBasedHCAL;
    std::vector<float> *T_Trig_detBasedECAL;
    std::vector<float> *T_Trig_detBasedCALO;
    std::vector<float> *T_Trig_PFecal;
    std::vector<float> *T_Trig_PFhcal;
    std::vector<float> *T_Trig_PFecalUnseeded;
    std::vector<float> *T_Trig_PFhcalUnseeded;
    
    std::vector<float> *T_Trig_rho;
    
    std::vector<float> *T_Gen_Eta;
    std::vector<float> *T_Gen_Phi;
    std::vector<float> *T_Gen_Pt;
    std::vector<float> *T_Gen_Px;
    std::vector<float> *T_Gen_Py;
    std::vector<float> *T_Gen_Pz;
    std::vector<float> *T_Gen_Energy;
    std::vector<int> *T_Gen_pdgID;
    std::vector<int> *T_Gen_MotherID;
    std::vector<int> *T_Gen_FromW;
    std::vector<int> *T_Gen_FromTau;
    
    
    
    //trigger leg
    std::vector<int> *T_Muon_TriggerLeg;
    
    std::vector<float>*T_Muon_Eta;
    std::vector<float>*T_Muon_Phi;
    std::vector<float>*T_Muon_Energy;
    std::vector<float>*T_Muon_Et;
    std::vector<float>*T_Muon_Pt;
    std::vector<float>*T_Muon_Px;
    std::vector<float>*T_Muon_Py;
    std::vector<float>*T_Muon_Pz;
    std::vector<float>*T_Muon_Mass;
    
    
    std::vector<bool> *T_Muon_IsGlobalMuon;
    std::vector<bool> *T_Muon_IsTrackerMuon;
    std::vector<bool> *T_Muon_IsPFMuon;
    std::vector<bool> *T_Muon_IsCaloMuon;
    std::vector<bool> *T_Muon_IsStandAloneMuon;
    std::vector<bool> *T_Muon_IsMuon;
    std::vector<bool> *T_Muon_IsGlobalMuon_PromptTight;
    std::vector<bool> *T_Muon_IsTrackerMuonArbitrated;
    std::vector<int>  *T_Muon_numberOfChambers;
    std::vector<int>  *T_Muon_numberOfChambersRPC;
    std::vector<int>  *T_Muon_numberOfMatches;
    std::vector<int>  *T_Muon_numberOfMatchedStations;
    std::vector<int>  *T_Muon_charge;
    
    
    std::vector<bool> *T_Muon_TMLastStationTight;
    std::vector<float> *T_Muon_globalTrackChi2;
    std::vector<int>  *T_Muon_validMuonHits;
    std::vector<float> *T_Muon_trkKink;
    std::vector<int>  *T_Muon_trkNbOfTrackerLayers;
    std::vector<int>  *T_Muon_trkNbOfValidTrackeHits;
    std::vector<int>  *T_Muon_trkValidPixelHits;
    std::vector<float> *T_Muon_trkError;
    std::vector<float> *T_Muon_dB;
    std::vector<float> *T_Muon_dzPV;
    std::vector<float> *T_Muon_dBstop;
    std::vector<float> *T_Muon_dzstop;
    
    // PF isolation
    std::vector<float> *T_Muon_chargedHadronIsoR04;
    std::vector<float> *T_Muon_neutralHadronIsoR04;
    std::vector<float> *T_Muon_photonIsoR04;
    std::vector<float> *T_Muon_chargedHadronIsoPUR04;
    
    std::vector<float> *T_Muon_chargedHadronIsoR03;
    std::vector<float> *T_Muon_neutralHadronIsoR03;
    std::vector<float> *T_Muon_photonIsoR03;
    std::vector<float> *T_Muon_chargedHadronIsoPUR03;
    
    std::vector<float> *T_Muon_isoR03_emEt;
    std::vector<float> *T_Muon_isoR03_hadEt;
    std::vector<float> *T_Muon_isoR03_hoEt;
    std::vector<float> *T_Muon_isoR03_sumPt;
    std::vector<int> *T_Muon_isoR03_nTracks;
    std::vector<int> *T_Muon_isoR03_nJets;
    std::vector<float> *T_Muon_isoRingsMVA;
    
    std::vector<float> *T_Muon_ecalPFiso;
    std::vector<float> *T_Muon_hcalPFiso;
    
  
    
    
   /* std::vector<float> *T_Elec_Eta;
    std::vector<float> *T_Elec_Phi;
    std::vector<float> *T_Elec_Pt;

    
    std::vector<float> *T_Muon_Eta;
    std::vector<float> *T_Muon_Phi;
    std::vector<float> *T_Muon_Pt;*/
};


typedef std::vector< edm::Handle<trigger::TriggerFilterObjectWithRefs> > TrigFiltVect;
typedef edm::AssociationMap<edm::OneToValue<std::vector<reco::RecoEcalCandidate>,float,unsigned int> > recoEcalCandidateMap;



