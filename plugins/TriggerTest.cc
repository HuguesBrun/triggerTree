#include "TriggerTest.h"





TriggerTest::TriggerTest(const edm::ParameterSet& iConfig)

{
    //now do what ever initialization is needed
    //   electronsCollection_      = iConfig.getParameter<edm::InputTag>("electronsCollection");
    //  muonsCollection_      = iConfig.getParameter<edm::InputTag>("muonsCollection");
    triggerResultsTag_ = iConfig.getParameter<edm::InputTag>("triggerResultTag");
    triggerSummaryLabel_= iConfig.getParameter<edm::InputTag>("triggerSummaryTag");
    pathsToSave_ = iConfig.getParameter<std::vector<std::string> >("pathsToSave");
    filterToMatch_ = iConfig.getParameter<std::vector<std::string> >("filterToMatch");
    mapsValues_ = iConfig.getParameter<std::vector<std::string> >("mapsValues");
    HLTprocess_   = iConfig.getParameter<std::string>("HLTprocess");
    rhoTag_ = iConfig.getParameter<edm::InputTag>("rhoTag");
    outputFile_   = iConfig.getParameter<std::string>("outputFile");
    trkExtractorPSet_ = iConfig.getParameter<edm::ParameterSet>("TrkExtractorPSet");
    // caloExtractorPSet_ = iConfig.getParameter<edm::ParameterSet>("CaloExtractorPSet");
    rootFile_ = TFile::Open(outputFile_.c_str(),"RECREATE");
    
    
    std::string trkExtractorName = trkExtractorPSet_.getParameter<std::string>("ComponentName");
    trkExtractor = IsoDepositExtractorFactory::get()->create( trkExtractorName, trkExtractorPSet_, consumesCollector());
    
}


TriggerTest::~TriggerTest()
{
    
    // do anything here that needs to be done at desctruction time
    // (e.g. close files, deallocate resources etc.)
    
}


//
// member functions
//

// ------------ method called for each event  ------------
void
TriggerTest::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
    using namespace edm;
    using namespace std;
    
    beginEvent(); //create the vectors
    
    /* edm::Handle<reco::GsfElectronCollection> electronsCollection;
     iEvent.getByLabel(electronsCollection_ , electronsCollection);
     
     
     edm::Handle < std::vector <reco::Muon> > recoMuons;
     iEvent.getByLabel(muonsCollection_, recoMuons);*/
    
    edm::Handle<GenEventInfoProduct> genEvent;
    iEvent.getByLabel("generator", genEvent);
    
    
    edm::Handle<reco::GenParticleCollection> genParticles;
    iEvent.getByLabel( "genParticles", genParticles );
    
    T_Event_RunNumber = iEvent.id().run();
    T_Event_EventNumber = iEvent.id().event();
    T_Event_LuminosityBlock = iEvent.id().luminosityBlock();
    
    float truePu=0.;
    Handle<std::vector< PileupSummaryInfo > > puInfo;
    try {
        iEvent.getByLabel("addPileupInfo",puInfo);
        std::vector<PileupSummaryInfo>::const_iterator PVI;
        //The in-time crossing is getBunchCrossing = 0; negative ones are early, positive ones are late.
        for(PVI = puInfo->begin(); PVI != puInfo->end(); ++PVI) {
            
            //    std::cout << " Pileup Information: bunchXing, nvtx: " << PVI->getBunchCrossing() << " " << PVI->getPU_NumInteractions() << std::endl;
            if(PVI->getBunchCrossing()==0){
                T_Event_nPU =PVI->getPU_NumInteractions();
                T_Event_nTruePU=PVI->getTrueNumInteractions();
                
            }
            
            else if(PVI->getBunchCrossing()==-1){
                T_Event_nPUm=PVI->getPU_NumInteractions();
            }
            else if(PVI->getBunchCrossing()==1){
                T_Event_nPUp=PVI->getPU_NumInteractions();
            }
            truePu += PVI->getTrueNumInteractions();
        }
    } catch (...) {}
    T_Event_AveNTruePU=truePu;
    
    
    bool changedConfig = false;
    if (!hltConfig.init(iEvent.getRun(), iSetup, HLTprocess_.c_str(), changedConfig)) {
        cout << "Initialization of HLTConfigProvider failed!!" << endl;
        return;
    }
    if (changedConfig){
        unsigned int nbPaths = pathsToSave_.size();
        
        std::cout << "the curent menu is " << hltConfig.tableName() << std::endl;
        for (unsigned int itePath=0 ; itePath<nbPaths ; itePath++){
            for (size_t j = 0; j < hltConfig.triggerNames().size(); j++) {
                if (TString(hltConfig.triggerNames()[j]).Contains(pathsToSave_.at(itePath))){
                    triggerBits_.push_back(j);
                    cout << "found the path " << pathsToSave_.at(itePath) << endl;
                }
                
            }
        }
        if (triggerBits_.size() <nbPaths) cout << "an HLT paths is not found ! ! " << endl;
        
    }
    
    edm::Handle<edm::TriggerResults> triggerResults;
    iEvent.getByLabel(triggerResultsTag_, triggerResults);
    
    for (unsigned int itePath = 0 ; itePath < triggerBits_.size() ; itePath++){
        if (triggerResults->accept(triggerBits_.at(itePath))) {
            edm::LogVerbatim("TriggerTest") << "passing path " << itePath << endl;
            T_Event_pathsFired->push_back(1);
        }
        else T_Event_pathsFired->push_back(0);
    }
    
    
    //std::string caloExtractorName = caloExtractorPSet_.getParameter<std::string>("ComponentName");
    //caloExtractor = IsoDepositExtractorFactory::get()->create( caloExtractorName, caloExtractorPSet_, consumesCollector());
    
    //if (mapsValues_.size() != filterToMatch_.size()) cout << "warning MAP and filters are not matching 1 to 1 !!!" << endl;
    
    // get rho
    edm::Handle<double> rhoHandle;
    iEvent.getByLabel(rhoTag_, rhoHandle);
    double rho = -1;
    
    edm::Handle<trigger::TriggerEvent> triggerSummary;
    iEvent.getByLabel(triggerSummaryLabel_, triggerSummary);
    trigger::TriggerObjectCollection allTriggerObjects = triggerSummary->getObjects();
    trigger::TriggerObjectCollection legObjects;
    std::vector<unsigned int> legRefs;
    // find the ref of the legs
    for (size_t iteFilter=0; iteFilter<filterToMatch_.size(); iteFilter++) {
        edm::InputTag filterTag = edm::InputTag(filterToMatch_.at(iteFilter), "", HLTprocess_.c_str());
        size_t filterIndex = (*triggerSummary).filterIndex(filterTag);
        if (filterIndex < (*triggerSummary).sizeFilters()) { //check if the trigger object is present
            //save the trigger objects corresponding to muon leg
            //cout << "found the filter " << filterToMatch_.at(iteFilter) << endl;
            const trigger::Keys &keys = (*triggerSummary).filterKeys(filterIndex);
            for (size_t j = 0; j < keys.size(); j++) {
                trigger::TriggerObject foundObject = (allTriggerObjects)[keys[j]];
                legObjects.push_back(foundObject);
                legRefs.push_back(iteFilter);
            }
        }
    }
    
    
    
    /* for (unsigned k = 0 ; k < legObjects.size() ; k++){
     trigger::TriggerObject theObject = legObjects.at(k);
     T_Trig_Pt->push_back(theObject.pt());
     T_Trig_Eta->push_back(theObject.eta());
     T_Trig_Phi->push_back(theObject.phi());
     T_Trig_Leg->push_back(legRefs.at(k));
     }*/
    
    TrigFiltVect hltFilters( filterToMatch_.size() );
    for (unsigned int i = 0 ; i <  filterToMatch_.size() ; i++){
        // cout << "i=" << i << endl;
        iEvent.getByLabel(edm::InputTag(filterToMatch_.at(i),"",HLTprocess_.c_str()), hltFilters[i]);
    }
    
    for (unsigned int i = 0 ; i <  filterToMatch_.size() ; i++){
        if (!(hltFilters[i].isValid())) continue;
        TString nameOfFilter = filterToMatch_.at(i);
        
        
        if (nameOfFilter.Contains("hltL3fL1sMu16L1f0L2f16QL3Filtered24Q")){
            std::vector<reco::RecoChargedCandidateRef> prevMuonRefs;
            reco::RecoChargedCandidateRef ref;
            hltFilters[i]->getObjects(trigger::TriggerMuon, prevMuonRefs);
            edm::Handle<reco::IsoDepositMap> IsoCaloMap;
            iEvent.getByLabel(edm::InputTag("hltL3CaloMuonCorrectedIsolations","",HLTprocess_.c_str()),IsoCaloMap);
            edm::Handle<reco::IsoDepositMap> IsoCaloMapNoECAL;
            iEvent.getByLabel(edm::InputTag("hltL3CaloMuonCorrectedIsolationsNoECAL","",HLTprocess_.c_str()),IsoCaloMapNoECAL);
            edm::Handle<reco::IsoDepositMap> IsoCaloMapNoHCAL;
            iEvent.getByLabel(edm::InputTag("hltL3CaloMuonCorrectedIsolationsNoHCAL","",HLTprocess_.c_str()),IsoCaloMapNoHCAL);
            
            edm::Handle<reco::RecoChargedCandidateIsolationMap> hcalIsoMap;
            iEvent.getByLabel(edm::InputTag("hltMuonHcalPFClusterIso","",HLTprocess_.c_str()),hcalIsoMap);
            
            edm::Handle<reco::RecoChargedCandidateIsolationMap> ecalIsoMap;
            iEvent.getByLabel(edm::InputTag("hltMuonEcalPFClusterIso","",HLTprocess_.c_str()),ecalIsoMap);
            
            edm::Handle<reco::RecoChargedCandidateIsolationMap> hcalIsoMapUnSeeded;
            iEvent.getByLabel(edm::InputTag("hltMuonHcalPFClusterIsoUnseeded","",HLTprocess_.c_str()),hcalIsoMapUnSeeded);
            
            edm::Handle<reco::RecoChargedCandidateIsolationMap> ecalIsoMapUnSeeded;
            iEvent.getByLabel(edm::InputTag("hltMuonEcalPFClusterIsoUnseeded","",HLTprocess_.c_str()),ecalIsoMapUnSeeded);
            
            edm::Handle<reco::IsoDepositMap> IsoTkMap;
            iEvent.getByLabel(edm::InputTag("hltL3MuonCombRelIsolationsIterTrkRegIter02","trkIsoDeposits",HLTprocess_.c_str()),IsoTkMap);
            
            unsigned int nMuons = prevMuonRefs.size();
            
            reco::IsoDeposit::Vetos trkVetos(nMuons);
            std::vector<reco::IsoDeposit> trkDeps(nMuons);
            
            //     reco::IsoDeposit::Vetos caloVetos(nMuons);
            //    std::vector<reco::IsoDeposit> caloDeps(nMuons);
            
            for (size_t j = 0 ; j < prevMuonRefs.size() ; j++){
                ref = prevMuonRefs[j];
                reco::TrackRef mu = ref->track();
                trkDeps[j] = trkExtractor->deposit(iEvent, iSetup, *mu);
                trkVetos[j] = trkDeps[j].veto();
                
                //     caloDeps[j] = caloExtractor->deposit(iEvent, iSetup, *mu);
                //     caloVetos[j] = caloDeps[j].veto();
                
            }
            
            for (size_t j = 0 ; j < prevMuonRefs.size() ; j++){
                ref = prevMuonRefs[j];
                reco::TrackRef mu = ref->track();
                //    reco::IsoDeposit trackDep = trkExtractor->deposit(iEvent, iSetup, *mu);
                //    reco::IsoDeposit::Veto trkVeto = trackDep.veto();
                float trkIso = trkDeps[j].depositWithin(0.3, trkVetos, -1.);
                T_Trig_TkIsoVeto->push_back(trkIso);
                //    float caloIsoSum = caloDeps[j].depositWithin(0.3, caloVetos);
                //    cout << "caloIso=" << caloIsoSum << endl;
                T_Trig_Eta->push_back(ref->eta());
                T_Trig_Pt->push_back(ref->pt());
                T_Trig_Phi->push_back(ref->phi());
                T_Trig_Leg->push_back(i);
                if (rhoHandle.isValid()) rho = *(rhoHandle.product());
                T_Trig_rho->push_back(rho);
                if (IsoTkMap.isValid()) {
                    reco::IsoDeposit theTkIsolation = (*IsoTkMap)[ref];
                    T_Trig_trackerIso->push_back(theTkIsolation.depositWithin(0.3));
                }
                else T_Trig_trackerIso->push_back(-1);
                if (IsoCaloMap.isValid()){
                    reco::IsoDeposit theCaloIsolation = (*IsoCaloMap)[ref];
                    T_Trig_detBasedCALO->push_back(theCaloIsolation.depositWithin(0.3));
                }
                else T_Trig_detBasedCALO->push_back(-1);
                if (IsoCaloMapNoHCAL.isValid()){
                    reco::IsoDeposit theCaloIsolationNoHCAL = (*IsoCaloMapNoHCAL)[ref];
                    T_Trig_detBasedECAL->push_back(theCaloIsolationNoHCAL.depositWithin(0.3));
                }
                else T_Trig_detBasedECAL->push_back(-1);
                if (IsoCaloMapNoECAL.isValid()){
                    reco::IsoDeposit theCaloIsolationNoECAL = (*IsoCaloMapNoECAL)[ref];
                    T_Trig_detBasedHCAL->push_back(theCaloIsolationNoECAL.depositWithin(0.3));
                }
                else T_Trig_detBasedHCAL->push_back(-1);
                
                if (hcalIsoMap.isValid()){
                    reco::RecoChargedCandidateIsolationMap::const_iterator valHcalIso = (*hcalIsoMap).find(ref);
                    T_Trig_PFhcal->push_back(valHcalIso->val);
                }
                else T_Trig_PFhcal->push_back(-1.);
                if (ecalIsoMap.isValid()){
                    reco::RecoChargedCandidateIsolationMap::const_iterator valEcalIso = (*ecalIsoMap).find(ref);
                    T_Trig_PFecal->push_back(valEcalIso->val);
                }
                else T_Trig_PFecal->push_back(-1.);
                
                if (hcalIsoMapUnSeeded.isValid()){
                    reco::RecoChargedCandidateIsolationMap::const_iterator valHcalIsoUnSeeded = (*hcalIsoMapUnSeeded).find(ref);
                    T_Trig_PFhcalUnseeded->push_back(valHcalIsoUnSeeded->val);
                }
                else T_Trig_PFhcalUnseeded->push_back(-1.);
                if (ecalIsoMapUnSeeded.isValid()){
                    reco::RecoChargedCandidateIsolationMap::const_iterator valEcalIsoUnSeeded = (*ecalIsoMapUnSeeded).find(ref);
                    T_Trig_PFecalUnseeded->push_back(valEcalIsoUnSeeded->val);
                }
                else T_Trig_PFecalUnseeded->push_back(-1.);
                
            }
        }
        else if (nameOfFilter.Contains("hltMuonHcalIsoFilter")){
            std::vector<reco::RecoChargedCandidateRef> prevMuonRefs;
            reco::RecoChargedCandidateRef ref;
            hltFilters[i]->getObjects(trigger::TriggerMuon, prevMuonRefs);
            cout << "nb of candidates =" << prevMuonRefs.size() << endl;
            for (size_t j = 0 ; j < prevMuonRefs.size() ; j++){
                ref = prevMuonRefs[j];
            }
            
        }
        else {
            std::vector<reco::RecoChargedCandidateRef> prevMuonRefs;
            reco::RecoChargedCandidateRef ref;
            hltFilters[i]->getObjects(trigger::TriggerMuon, prevMuonRefs);
            for (size_t j = 0 ; j < prevMuonRefs.size() ; j++){
                ref = prevMuonRefs[j];
                T_Trig_Eta->push_back(ref->eta());
                T_Trig_Pt->push_back(ref->pt());
                T_Trig_Phi->push_back(ref->phi());
                T_Trig_Leg->push_back(i);
                T_Trig_trackerIso->push_back(-1);
                T_Trig_detBasedCALO->push_back(-1);
                T_Trig_rho->push_back(-1);
                T_Trig_detBasedECAL->push_back(-1);
                T_Trig_detBasedHCAL->push_back(-1);
                T_Trig_PFhcal->push_back(-1);
                T_Trig_PFecal->push_back(-1);
                T_Trig_PFhcalUnseeded->push_back(-1);
                T_Trig_PFecalUnseeded->push_back(-1);
            }
        }
        
    }
    
    
    
    
    
    
    int theNbOfGen = genParticles->size();
    for (int i=0 ; i < theNbOfGen; i++){
        const reco::GenParticle & genMuon = (*genParticles)[i];
        if ((fabs(genMuon.pdgId())!=13)&&(fabs(genMuon.pdgId())!=11)) continue;
        if (!(genMuon.status()==1)) continue;
        
        const reco::Candidate  *MomPart =genMuon.mother();
        
        T_Gen_Eta->push_back(genMuon.eta());
        T_Gen_Phi->push_back(genMuon.phi());
        T_Gen_Pt->push_back(genMuon.pt());
        
        T_Gen_pdgID->push_back(genMuon.pdgId());
        T_Gen_MotherID->push_back(MomPart->pdgId());
        T_Gen_FromW->push_back(hasWasMother(genMuon));
        T_Gen_FromTau->push_back(hasTauasMother(genMuon));
        
        
    }
    
    
    
    mytree_->Fill();
    
    endEvent();
    
}


// ------------ method called once each job just before starting event loop  ------------
void
TriggerTest::beginJob()
{
    mytree_ = new TTree("eventsTree","");
    mytree_->Branch("T_Event_RunNumber", &T_Event_RunNumber, "T_Event_RunNumber/I");
    mytree_->Branch("T_Event_EventNumber", &T_Event_EventNumber, "T_Event_EventNumber/I");
    mytree_->Branch("T_Event_LuminosityBlock", &T_Event_LuminosityBlock, "T_Event_LuminosityBlock/I");
    
    mytree_->Branch("T_Event_nPU", &T_Event_nPU, "T_Event_nPU/I");
    mytree_->Branch("T_Event_nTruePU", &T_Event_nTruePU, "T_Event_nTruePU/F");
    mytree_->Branch("T_Event_nPUm", &T_Event_nPUm, "T_Event_nPUm/I");
    mytree_->Branch("T_Event_nPUp", &T_Event_nPUp, "T_Event_nPUp/I");
    mytree_->Branch("T_Event_AveNTruePU", &T_Event_AveNTruePU, "T_Event_AveNTruePU/F");
    
    
    mytree_->Branch("T_Event_pathsFired", "std::vector<int>", &T_Event_pathsFired);
    
    mytree_->Branch("T_Trig_Eta", "std::vector<float>", &T_Trig_Eta);
    mytree_->Branch("T_Trig_Pt", "std::vector<float>", &T_Trig_Pt);
    mytree_->Branch("T_Trig_Phi", "std::vector<float>", &T_Trig_Phi);
    mytree_->Branch("T_Trig_Leg", "std::vector<int>", &T_Trig_Leg);
    mytree_->Branch("T_Trig_TkIsoVeto", "std::vector<float>", &T_Trig_TkIsoVeto);
    mytree_->Branch("T_Trig_trackerIso", "std::vector<float>", &T_Trig_trackerIso);
    mytree_->Branch("T_Trig_detBasedHCAL", "std::vector<float>", &T_Trig_detBasedHCAL);
    mytree_->Branch("T_Trig_detBasedECAL", "std::vector<float>", &T_Trig_detBasedECAL);
    mytree_->Branch("T_Trig_detBasedCALO", "std::vector<float>", &T_Trig_detBasedCALO);
    // mytree_->Branch("T_Trig_PFecal", "std::vector<float>", &T_Trig_PFecal);
    //mytree_->Branch("T_Trig_PFhcal", "std::vector<float>", &T_Trig_PFhcal);
    mytree_->Branch("T_Trig_PFecalUnseeded", "std::vector<float>", &T_Trig_PFecalUnseeded);
    mytree_->Branch("T_Trig_PFhcalUnseeded", "std::vector<float>", &T_Trig_PFhcalUnseeded);
    mytree_->Branch("T_Trig_rho", "std::vector<float>", &T_Trig_rho);
    
    
    mytree_->Branch("T_Gen_Eta", "std::vector<float>", &T_Gen_Eta);
    mytree_->Branch("T_Gen_Phi", "std::vector<float>", &T_Gen_Phi);
    mytree_->Branch("T_Gen_Pt", "std::vector<float>", &T_Gen_Pt);
    mytree_->Branch("T_Gen_pdgID", "std::vector<int>", &T_Gen_pdgID);
    mytree_->Branch("T_Gen_MotherID", "std::vector<int>", &T_Gen_MotherID);
    mytree_->Branch("T_Gen_FromW", "std::vector<int>", &T_Gen_FromW);
    mytree_->Branch("T_Gen_FromTau", "std::vector<int>", &T_Gen_FromTau);
    
    /* mytree_->Branch("T_Elec_Eta", "std::vector<float>", &T_Elec_Eta);
     mytree_->Branch("T_Elec_Phi", "std::vector<float>", &T_Elec_Phi);
     mytree_->Branch("T_Elec_Pt", "std::vector<float>", &T_Elec_Pt);
     
     mytree_->Branch("T_Muon_Eta", "std::vector<float>", &T_Muon_Eta);
     mytree_->Branch("T_Muon_Phi", "std::vector<float>", &T_Muon_Phi);
     mytree_->Branch("T_Muon_Pt", "std::vector<float>", &T_Muon_Pt);*/
    
    
}

// ------------ method called once each job just after ending the event loop  ------------
void
TriggerTest::endJob()
{
    rootFile_->Write();
    rootFile_->Close();
}

void
TriggerTest::beginEvent()
{
    T_Event_pathsFired = new std::vector<int>;
    
    T_Trig_Eta = new std::vector<float>;
    T_Trig_Pt = new std::vector<float>;
    T_Trig_Phi = new std::vector<float>;
    T_Trig_Leg = new std::vector<int>;
    T_Trig_TkIsoVeto = new std::vector<float>;
    T_Trig_trackerIso = new std::vector<float>;
    T_Trig_detBasedHCAL = new std::vector<float>;
    T_Trig_detBasedECAL = new std::vector<float>;
    T_Trig_detBasedCALO = new std::vector<float>;
    T_Trig_PFecal = new std::vector<float>;
    T_Trig_PFhcal = new std::vector<float>;
    T_Trig_PFecalUnseeded = new std::vector<float>;
    T_Trig_PFhcalUnseeded = new std::vector<float>;
    T_Trig_rho = new std::vector<float>;
    
    T_Gen_Eta = new std::vector<float>;
    T_Gen_Phi = new std::vector<float>;
    T_Gen_Pt = new std::vector<float>;
    T_Gen_pdgID = new std::vector<int>;
    T_Gen_MotherID = new std::vector<int>;
    T_Gen_FromW = new std::vector<int>;
    T_Gen_FromTau = new std::vector<int>;
    
    /* T_Elec_Eta = new std::vector<float>;
     T_Elec_Phi = new std::vector<float>;
     T_Elec_Pt = new std::vector<float>;
     
     T_Muon_Eta = new std::vector<float>;
     T_Muon_Phi = new std::vector<float>;
     T_Muon_Pt = new std::vector<float>;*/
    
}

void
TriggerTest::endEvent()
{
    delete T_Event_pathsFired;
    
    delete T_Trig_Eta;
    delete T_Trig_Pt;
    delete T_Trig_Phi;
    delete T_Trig_Leg;
    delete T_Trig_TkIsoVeto;
    delete T_Trig_trackerIso;
    delete T_Trig_detBasedHCAL;
    delete T_Trig_detBasedECAL;
    delete T_Trig_detBasedCALO;
    delete T_Trig_PFecal;
    delete T_Trig_PFhcal;
    delete T_Trig_PFecalUnseeded;
    delete T_Trig_PFhcalUnseeded;
    delete T_Trig_rho;
    
    delete T_Gen_Eta;
    delete T_Gen_Phi;
    delete T_Gen_Pt;
    delete T_Gen_pdgID;
    delete T_Gen_MotherID;
    delete T_Gen_FromW;
    delete T_Gen_FromTau;
    
    /* delete T_Muon_Eta;
     delete T_Muon_Phi;
     delete T_Muon_Pt;
     
     
     delete T_Elec_Eta;
     delete T_Elec_Phi;
     delete T_Elec_Pt;*/
    
}

// ------------ method called when starting to processes a run  ------------
/*
 void
 TriggerTest::beginRun(edm::Run const&, edm::EventSetup const&)
 {
 }
 */

// ------------ method called when ending the processing of a run  ------------
/*
 void
 TriggerTest::endRun(edm::Run const&, edm::EventSetup const&)
 {
 }
 */

// ------------ method called when starting to processes a luminosity block  ------------
/*
 void
 TriggerTest::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
 {
 }
 */

// ------------ method called when ending the processing of a luminosity block  ------------
/*
 void
 TriggerTest::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
 {
 }
 */

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
TriggerTest::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
    
    
}

bool
TriggerTest::hasWasMother(const reco::GenParticle  p)
{
    bool foundW = false;
    if (p.numberOfMothers()==0) return foundW;
    const reco::Candidate  *part = (p.mother());
    // loop on the mother particles to check if is has a W has mother
    while ((part->numberOfMothers()>0)) {
        const reco::Candidate  *MomPart =part->mother();
        if (fabs(MomPart->pdgId())==24){
            foundW = true;
            break;
        }
        part = MomPart;
    }
    return foundW;
}

bool
TriggerTest::hasTauasMother(const reco::GenParticle  p)
{
    bool foundW = false;
    if (p.numberOfMothers()==0) return foundW;
    const reco::Candidate  *part = (p.mother());
    // loop on the mother particles to check if is has a W has mother
    while ((part->numberOfMothers()>0)) {
        const reco::Candidate  *MomPart =part->mother();
        if (fabs(MomPart->pdgId())==15){
            foundW = true;
            break;
        }
        part = MomPart;
    }
    return foundW;
}

//define this as a plug-in
DEFINE_FWK_MODULE(TriggerTest);
