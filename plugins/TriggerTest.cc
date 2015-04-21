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
    primaryVertexInputTag_    = iConfig.getParameter<edm::InputTag>("primaryVertexInputTag");
    doMuon_                   = iConfig.getParameter<bool>("doMuon");
    muonProducers_			= iConfig.getParameter<vtag>("muonProducer");
    muonECALpfIsoTag_       = iConfig.getParameter<edm::InputTag>("muonECALpfIsoTag");
    muonHCALpfIsoTag_       = iConfig.getParameter<edm::InputTag>("muonHCALpfIsoTag");
    rhoInputTags_           = iConfig.getParameter<std::vector<edm::InputTag> >("rhoTags");
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
     */
     
    
    edm::Handle<GenEventInfoProduct> genEvent;
    iEvent.getByLabel("generator", genEvent);
    
    
    edm::Handle<reco::GenParticleCollection> genParticles;
    iEvent.getByLabel( "genParticles", genParticles );
    
    // load the vertices collection
    edm::Handle<reco::VertexCollection> vtx_h;
    iEvent.getByLabel(primaryVertexInputTag_, vtx_h);
    
    //muon collection :
    edm::Handle < std::vector <reco::Muon> > recoMuons;
    edm::InputTag muonProducer = muonProducers_.at(0);
    iEvent.getByLabel(muonProducer, recoMuons);
    
    // muon PF clustering isolation

    edm::Handle< edm::ValueMap<float> > muonECALIsoMapH;
    iEvent.getByLabel(muonECALpfIsoTag_,muonECALIsoMapH);
    const edm::ValueMap<float> muonECALIsoMap  = *(muonECALIsoMapH);
    
    edm::Handle< edm::ValueMap<float> > muonHCALIsoMapH;
    iEvent.getByLabel(muonHCALpfIsoTag_,muonHCALIsoMapH);
    const edm::ValueMap<float> muonHCALIsoMap  = *(muonHCALIsoMapH);
    

    
    T_Event_RunNumber = iEvent.id().run();
    T_Event_EventNumber = iEvent.id().event();
    T_Event_LuminosityBlock = iEvent.id().luminosityBlock();
    
    
    edm::InputTag trackingParticlesTag = edm::InputTag("mix","MergedTrackTruth");
    
    edm::Handle<TrackingParticleCollection>  TPCollectionH ;
    TrackingParticleCollection tPC;
    iEvent.getByLabel(trackingParticlesTag,TPCollectionH);
    if (TPCollectionH.isValid()) tPC   = *(TPCollectionH.product());
    else cout << "not found tracking particles collection" << endl;
    
    
    reco::Vertex dummy;
    const reco::Vertex *pv = &dummy;
    if (vtx_h->size() != 0) {
        pv = &*vtx_h->begin();
    } else { // create a dummy PV
        reco::Vertex::Error e;
        e(0, 0) = 0.0015 * 0.0015;
        e(1, 1) = 0.0015 * 0.0015;
        e(2, 2) = 15. * 15.;
        reco::Vertex::Point p(0, 0, 0);
        dummy = reco::Vertex(p, e, 0, 0, 0);
    }
    
    
   // cout << "event=" << T_Event_EventNumber << endl;
    
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
                //if (TString(hltConfig.triggerNames()[j]).Contains(pathsToSave_.at(itePath))){
                if (TString(hltConfig.triggerNames()[j]) == pathsToSave_.at(itePath)){
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
            //edm::LogVerbatim("TriggerTest") << "passing path " << itePath << endl;
//	    cout  << "passing path " << itePath << endl;
            T_Event_pathsFired->push_back(1);
        }
        else T_Event_pathsFired->push_back(0);
    }
    
    
    //std::string caloExtractorName = caloExtractorPSet_.getParameter<std::string>("ComponentName");
    //caloExtractor = IsoDepositExtractorFactory::get()->create( caloExtractorName, caloExtractorPSet_, consumesCollector());
    
    //if (mapsValues_.size() != filterToMatch_.size()) cout << "warning MAP and filters are not matching 1 to 1 !!!" << endl;
    
    //fill rho offline
    rhoHandles rhos(rhoInputTags_.size());
    for (unsigned int iteRho = 0 ; iteRho < rhoInputTags_.size() ; iteRho++){
        iEvent.getByLabel(rhoInputTags_[iteRho], rhos[iteRho]);
        T_Event_Rho->push_back(*rhos[iteRho]);
    }
    
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
        //    cout << "found the filter " << filterToMatch_.at(iteFilter) << endl;
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
        
        
        if (nameOfFilter.Contains("hltL3fL1sMu16L1f0L2f10QL3Filtered20Q")){
            std::vector<reco::RecoChargedCandidateRef> prevMuonRefs;
            reco::RecoChargedCandidateRef ref;
            hltFilters[i]->getObjects(trigger::TriggerMuon, prevMuonRefs);
            edm::Handle< edm::ValueMap<float> > IsoCaloMap;
            iEvent.getByLabel(edm::InputTag("hltL3CaloMuonCorrectedIsolations","",HLTprocess_.c_str()),IsoCaloMap);
            edm::Handle< edm::ValueMap<float> > IsoCaloMapECAL;
            iEvent.getByLabel(edm::InputTag("hltL3CaloMuonCorrectedIsolationsECAL","",HLTprocess_.c_str()),IsoCaloMapECAL);
            edm::Handle< edm::ValueMap<float> > IsoCaloMapHCAL;
            iEvent.getByLabel(edm::InputTag("hltL3CaloMuonCorrectedIsolationsHCAL","",HLTprocess_.c_str()),IsoCaloMapHCAL);
            
            edm::Handle<reco::RecoChargedCandidateIsolationMap> hcalIsoMap;
            iEvent.getByLabel(edm::InputTag("hltMuonHcalPFClusterIsoForMuons","",HLTprocess_.c_str()),hcalIsoMap);
            
            edm::Handle<reco::RecoChargedCandidateIsolationMap> ecalIsoMap;
            iEvent.getByLabel(edm::InputTag("hltMuonEcalPFClusterIsoForMuons","",HLTprocess_.c_str()),ecalIsoMap);
            
            /*edm::Handle<reco::RecoChargedCandidateIsolationMap> hcalIsoMapUnSeeded;
            iEvent.getByLabel(edm::InputTag("hltMuonHcalPFClusterIsoUnseeded","",HLTprocess_.c_str()),hcalIsoMapUnSeeded);
            
            edm::Handle<reco::RecoChargedCandidateIsolationMap> ecalIsoMapUnSeeded;
            iEvent.getByLabel(edm::InputTag("hltMuonEcalPFClusterIsoUnseeded","",HLTprocess_.c_str()),ecalIsoMapUnSeeded);*/
            
            edm::Handle<reco::IsoDepositMap> IsoTkMap;
            iEvent.getByLabel(edm::InputTag("hltMuonTkRelIsolationCut0p09Map","trkIsoDeposits",HLTprocess_.c_str()),IsoTkMap);
            
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
                //    reco::IsoDeposit::Vet$i trkVeto = trackDep.veto();
                //   
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
		//float myIsoCalo = -1;
                if (IsoCaloMap.isValid()){
                    float theCaloIsolation = (*IsoCaloMap)[ref];
                    T_Trig_detBasedCALO->push_back(theCaloIsolation);
     		  //  myIsoCalo = (theCaloIsolation>0 ? theCaloIsolation : 0);
                }
                else T_Trig_detBasedCALO->push_back(-1);
   		//cout << "eta=" << ref->eta() << " phi=" << ref->phi() << " trkIso=" << trkIso << " caloIso=" << myIsoCalo  << "pt=" << ref->pt() << endl;
		//cout << "comb=" << ((trkIso+myIsoCalo)/ref->pt()) << " pass=" << ((trkIso+myIsoCalo)/ref->pt()<0.15) << endl;
                if (IsoCaloMapHCAL.isValid()){
                    float theCaloIsolationHCAL = (*IsoCaloMapHCAL)[ref];
                    T_Trig_detBasedHCAL->push_back(theCaloIsolationHCAL);
                }
                else T_Trig_detBasedHCAL->push_back(-1);
                if (IsoCaloMapECAL.isValid()){
                    float theCaloIsolationECAL = (*IsoCaloMapECAL)[ref];
                    T_Trig_detBasedECAL->push_back(theCaloIsolationECAL);
                }
                else T_Trig_detBasedECAL->push_back(-1);
                
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
                
                /*if (hcalIsoMapUnSeeded.isValid()){
                    reco::RecoChargedCandidateIsolationMap::const_iterator valHcalIsoUnSeeded = (*hcalIsoMapUnSeeded).find(ref);
                    T_Trig_PFhcalUnseeded->push_back(valHcalIsoUnSeeded->val);
                }
                else T_Trig_PFhcalUnseeded->push_back(-1.);
                if (ecalIsoMapUnSeeded.isValid()){
                    reco::RecoChargedCandidateIsolationMap::const_iterator valEcalIsoUnSeeded = (*ecalIsoMapUnSeeded).find(ref);
                    T_Trig_PFecalUnseeded->push_back(valEcalIsoUnSeeded->val);
                }
                else T_Trig_PFecalUnseeded->push_back(-1.);*/
                
            }
        }
        else if (nameOfFilter.Contains("hltMuonHcalIsoFilter")){
            std::vector<reco::RecoChargedCandidateRef> prevMuonRefs;
            reco::RecoChargedCandidateRef ref;
            hltFilters[i]->getObjects(trigger::TriggerMuon, prevMuonRefs);
            //cout << "nb of candidates =" << prevMuonRefs.size() << endl;
            for (size_t j = 0 ; j < prevMuonRefs.size() ; j++){
                ref = prevMuonRefs[j];
            }
            
        }
	else if (nameOfFilter.Contains("hltL3crIsoL1sMu16L1f0L2f10QL3f20QL3trkIsoFiltered0p09")){
	    std::vector<reco::RecoChargedCandidateRef> prevMuonRefs;
            reco::RecoChargedCandidateRef ref;
            hltFilters[i]->getObjects(trigger::TriggerMuon, prevMuonRefs);
            //cout << "pass the isolation=" << endl;
            for (size_t j = 0 ; j < prevMuonRefs.size() ; j++){
                ref = prevMuonRefs[j];
		//cout << "eta=" << ref->eta() << " phi=" << ref->phi() << endl;
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
        
        T_Gen_Px->push_back(genMuon.px());
        T_Gen_Py->push_back(genMuon.py());
        T_Gen_Pz->push_back(genMuon.pz());
        T_Gen_Energy->push_back(genMuon.energy());

        
        
        T_Gen_pdgID->push_back(genMuon.pdgId());
        T_Gen_MotherID->push_back(MomPart->pdgId());
        T_Gen_FromW->push_back(hasWasMother(genMuon));
        T_Gen_FromTau->push_back(hasTauasMother(genMuon));
        
        
    }
    
    
    

   
   //starting loop on tracking particles
    
    bool needToVeto = false;
   for (TrackingParticleCollection::size_type i=0; i<tPC.size(); i++) {
   TrackingParticleRef trpart(TPCollectionH, i);
   
      /* cout << "tp eta=" << trpart->eta() << endl;
       cout << "tp phi=" << trpart->phi() << endl;
       
       cout << "tp pT=" << trpart->pt() << "pdgID=" << trpart->pdgId() << " vx=" << trpart->vx() << " vy=" << trpart->vy() << endl;*/
       if (fabs(trpart->pdgId())!=13) continue;
       float rhoVtx = sqrt(pow(trpart->vx(),2)+pow(trpart->vy(),2));
       float zVtx = fabs(trpart->vz());
       if ((trpart->pt()>5)&&(fabs(trpart->eta())<2.5)&&(rhoVtx<2000)&&(zVtx<4000)) needToVeto = true;
   }
    
    antiMuEnrichementVeto = needToVeto;
    
    if (doMuon_){
        int nbMuons = recoMuons->size();
     //   cout << "il y a " << nbMuons << " muons " << endl;
        //loop on the muons in the event
        for (int k = 0 ; k < nbMuons ; k++){
            
            const reco::Muon* muon = &((*recoMuons)[k]);
            //  cout << "le muon" << k << " : eta=" << muon->eta() << " phi=" << muon->phi() << endl;

            T_Muon_Eta->push_back(muon->eta());
            T_Muon_Phi->push_back(muon->phi());
            T_Muon_IsGlobalMuon->push_back(muon->isGlobalMuon());
            T_Muon_IsPFMuon->push_back(muon->isPFMuon());
            T_Muon_IsTrackerMuon->push_back(muon->isTrackerMuon());
            T_Muon_IsCaloMuon->push_back(muon->isCaloMuon());
            T_Muon_IsStandAloneMuon->push_back(muon->isStandAloneMuon());
            T_Muon_IsMuon->push_back(muon->isMuon());
            T_Muon_Energy->push_back(muon->energy());
            T_Muon_Et->push_back(muon->et());
            T_Muon_Pt->push_back(muon->pt());
            T_Muon_Px->push_back(muon->px());
            T_Muon_Py->push_back(muon->py());
            T_Muon_Pz->push_back(muon->pz());
            T_Muon_Mass->push_back(muon->mass());
            T_Muon_charge->push_back(muon->charge());
            
            T_Muon_numberOfChambers->push_back(muon->numberOfChambers());
            T_Muon_numberOfChambersRPC->push_back(muon->numberOfChambersNoRPC());
            T_Muon_numberOfMatches->push_back(muon->numberOfMatches());
            T_Muon_numberOfMatchedStations->push_back(muon->numberOfMatchedStations());
            bool isMatchTheStation = muon::isGoodMuon(*muon, muon::TMOneStationTight);
            bool isGlobalMuonPT = muon::isGoodMuon(*muon, muon::GlobalMuonPromptTight);
            bool isGlobalMuonArbitrated = muon::isGoodMuon(*muon, muon::TrackerMuonArbitrated);
            T_Muon_TMLastStationTight->push_back(isMatchTheStation);
            T_Muon_IsGlobalMuon_PromptTight->push_back(isGlobalMuonPT);
            T_Muon_IsTrackerMuonArbitrated->push_back(isGlobalMuonArbitrated);
            
            if (muon->globalTrack().isNull()) T_Muon_globalTrackChi2->push_back(-1); else T_Muon_globalTrackChi2->push_back(muon->globalTrack()->normalizedChi2());
            if (muon->globalTrack().isNull()) T_Muon_validMuonHits->push_back(-1); else T_Muon_validMuonHits->push_back(muon->globalTrack()->hitPattern().numberOfValidMuonHits());
            T_Muon_trkKink->push_back(muon->combinedQuality().trkKink);
            if (muon->muonBestTrack().isNull()) {
                T_Muon_trkNbOfTrackerLayers->push_back(-1);
                T_Muon_trkError->push_back(-1);
                T_Muon_dB->push_back(-1);
                T_Muon_dBstop->push_back(-1);
                T_Muon_dzPV->push_back(-1);
                T_Muon_dzstop->push_back(-1);
                T_Muon_trkValidPixelHits->push_back(-1);
                T_Muon_trkNbOfValidTrackeHits->push_back(-1);
            }
            else {
                T_Muon_trkNbOfTrackerLayers->push_back(muon->muonBestTrack()->hitPattern().trackerLayersWithMeasurement());
                T_Muon_trkError->push_back(muon->muonBestTrack()->ptError());
                T_Muon_trkValidPixelHits->push_back(muon->muonBestTrack()->hitPattern().numberOfValidPixelHits());
                T_Muon_dB->push_back(fabs(muon->muonBestTrack()->dxy(pv->position())));
                T_Muon_dBstop->push_back(fabs(muon->muonBestTrack()->dxy(pv->position())));
                T_Muon_dzPV->push_back(fabs(muon->muonBestTrack()->dz(pv->position())));
                T_Muon_dzstop->push_back(fabs(muon->muonBestTrack()->dz(pv->position())));
                T_Muon_trkNbOfValidTrackeHits->push_back(muon->muonBestTrack()->hitPattern().numberOfValidTrackerHits());
            }
            T_Muon_isoR03_emEt->push_back(muon->isolationR03().emEt);
            T_Muon_isoR03_hadEt->push_back(muon->isolationR03().hadEt);
            T_Muon_isoR03_hoEt->push_back(muon->isolationR03().hoEt);
            T_Muon_isoR03_sumPt->push_back(muon->isolationR03().sumPt);
            T_Muon_isoR03_nTracks->push_back(muon->isolationR03().nTracks);
            T_Muon_isoR03_nJets->push_back(muon->isolationR03().nJets);
            T_Muon_chargedHadronIsoR04->push_back(muon->pfIsolationR04().sumChargedHadronPt);
            T_Muon_neutralHadronIsoR04->push_back(muon->pfIsolationR04().sumNeutralHadronEt);
            T_Muon_photonIsoR04->push_back(muon->pfIsolationR04().sumPhotonEt);
            T_Muon_chargedHadronIsoPUR04->push_back(muon->pfIsolationR04().sumPUPt);
            T_Muon_chargedHadronIsoR03->push_back(muon->pfIsolationR03().sumChargedHadronPt);
            T_Muon_neutralHadronIsoR03->push_back(muon->pfIsolationR03().sumNeutralHadronEt);
            T_Muon_photonIsoR03->push_back(muon->pfIsolationR03().sumPhotonEt);
            T_Muon_chargedHadronIsoPUR03->push_back(muon->pfIsolationR03().sumPUPt);
            

            //PF clustering isolation
            reco::MuonRef thisMuonRef = reco::MuonRef(recoMuons, k);
            T_Muon_ecalPFiso->push_back(muonECALIsoMap[thisMuonRef]);
            T_Muon_hcalPFiso->push_back(muonHCALIsoMap[thisMuonRef]);
            
        }
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
    
    mytree_->Branch("T_Event_Rho", "std::vector<float>", &T_Event_Rho);

    
    mytree_->Branch("antiMuEnrichementVeto", &antiMuEnrichementVeto, "antiMuEnrichementVeto/I");

    mytree_->Branch("T_Event_pathsFired", "std::vector<int>", &T_Event_pathsFired);
    
    mytree_->Branch("T_Trig_Eta", "std::vector<float>", &T_Trig_Eta);
    mytree_->Branch("T_Trig_Pt", "std::vector<float>", &T_Trig_Pt);
    mytree_->Branch("T_Trig_Phi", "std::vector<float>", &T_Trig_Phi);
    mytree_->Branch("T_Trig_Leg", "std::vector<int>", &T_Trig_Leg);
    mytree_->Branch("T_Trig_TkIsoVeto", "std::vector<float>", &T_Trig_TkIsoVeto);
    mytree_->Branch("T_Trig_trackerIso", "std::vector<float>", &T_Trig_trackerIso);
  //  mytree_->Branch("T_Trig_detBasedHCAL", "std::vector<float>", &T_Trig_detBasedHCAL);
   // mytree_->Branch("T_Trig_detBasedECAL", "std::vector<float>", &T_Trig_detBasedECAL);
  //  mytree_->Branch("T_Trig_detBasedCALO", "std::vector<float>", &T_Trig_detBasedCALO);
    mytree_->Branch("T_Trig_PFecal", "std::vector<float>", &T_Trig_PFecal);
    mytree_->Branch("T_Trig_PFhcal", "std::vector<float>", &T_Trig_PFhcal);
    //mytree_->Branch("T_Trig_PFecalUnseeded", "std::vector<float>", &T_Trig_PFecalUnseeded);
    //mytree_->Branch("T_Trig_PFhcalUnseeded", "std::vector<float>", &T_Trig_PFhcalUnseeded);
    mytree_->Branch("T_Trig_rho", "std::vector<float>", &T_Trig_rho);
    
    
    mytree_->Branch("T_Gen_Eta", "std::vector<float>", &T_Gen_Eta);
    mytree_->Branch("T_Gen_Phi", "std::vector<float>", &T_Gen_Phi);
    mytree_->Branch("T_Gen_Pt", "std::vector<float>", &T_Gen_Pt);
    mytree_->Branch("T_Gen_Px", "std::vector<float>", &T_Gen_Px);
    mytree_->Branch("T_Gen_Py", "std::vector<float>", &T_Gen_Py);
    mytree_->Branch("T_Gen_Pz", "std::vector<float>", &T_Gen_Pz);
    mytree_->Branch("T_Gen_Energy", "std::vector<float>", &T_Gen_Energy);
    
    mytree_->Branch("T_Gen_pdgID", "std::vector<int>", &T_Gen_pdgID);
    mytree_->Branch("T_Gen_MotherID", "std::vector<int>", &T_Gen_MotherID);
    mytree_->Branch("T_Gen_FromW", "std::vector<int>", &T_Gen_FromW);
    mytree_->Branch("T_Gen_FromTau", "std::vector<int>", &T_Gen_FromTau);
    
    if (doMuon_){
        
        mytree_->Branch("T_Muon_Eta", "std::vector<float>", &T_Muon_Eta);
        mytree_->Branch("T_Muon_Phi", "std::vector<float>", &T_Muon_Phi);
        mytree_->Branch("T_Muon_Energy", "std::vector<float>", &T_Muon_Energy);
        mytree_->Branch("T_Muon_Et", "std::vector<float>", &T_Muon_Et);
        mytree_->Branch("T_Muon_Pt", "std::vector<float>", &T_Muon_Pt);
        mytree_->Branch("T_Muon_Px", "std::vector<float>", &T_Muon_Px);
        mytree_->Branch("T_Muon_Py", "std::vector<float>", &T_Muon_Py);
        mytree_->Branch("T_Muon_Pz", "std::vector<float>", &T_Muon_Pz);
        mytree_->Branch("T_Muon_Mass", "std::vector<float>", &T_Muon_Mass);
        mytree_->Branch("T_Muon_IsGlobalMuon", "std::vector<bool>", &T_Muon_IsGlobalMuon);
        mytree_->Branch("T_Muon_IsTrackerMuon", "std::vector<bool>", &T_Muon_IsTrackerMuon);
        mytree_->Branch("T_Muon_IsPFMuon", "std::vector<bool>", &T_Muon_IsPFMuon);
        mytree_->Branch("T_Muon_IsCaloMuon", "std::vector<bool>", &T_Muon_IsCaloMuon);
        mytree_->Branch("T_Muon_IsStandAloneMuon", "std::vector<bool>", &T_Muon_IsStandAloneMuon);
        mytree_->Branch("T_Muon_IsMuon", "std::vector<bool>", &T_Muon_IsMuon);
        mytree_->Branch("T_Muon_IsGlobalMuon_PromptTight", "std::vector<bool>", &T_Muon_IsGlobalMuon_PromptTight);
        mytree_->Branch("T_Muon_IsTrackerMuonArbitrated", "std::vector<bool>", &T_Muon_IsTrackerMuonArbitrated);
        mytree_->Branch("T_Muon_numberOfChambers", "std::vector<int>", &T_Muon_numberOfChambers);
        mytree_->Branch("T_Muon_numberOfChambersRPC", "std::vector<int>", &T_Muon_numberOfChambersRPC);
        mytree_->Branch("T_Muon_numberOfMatches", "std::vector<int>", &T_Muon_numberOfMatches);
        mytree_->Branch("T_Muon_numberOfMatchedStations", "std::vector<int>", &T_Muon_numberOfMatchedStations);
        mytree_->Branch("T_Muon_charge", "std::vector<int>", &T_Muon_charge);
        mytree_->Branch("T_Muon_TMLastStationTight", "std::vector<bool>", &T_Muon_TMLastStationTight);
        mytree_->Branch("T_Muon_globalTrackChi2", "std::vector<float>", &T_Muon_globalTrackChi2);
        mytree_->Branch("T_Muon_validMuonHits", "std::vector<int>", &T_Muon_validMuonHits);
        mytree_->Branch("T_Muon_trkKink", "std::vector<float>", &T_Muon_trkKink);
        mytree_->Branch("T_Muon_trkNbOfTrackerLayers", "std::vector<int>", &T_Muon_trkNbOfTrackerLayers);
        mytree_->Branch("T_Muon_trkNbOfValidTrackeHits", "std::vector<int>", &T_Muon_trkNbOfValidTrackeHits);
        mytree_->Branch("T_Muon_trkValidPixelHits", "std::vector<int>", &T_Muon_trkValidPixelHits);
        mytree_->Branch("T_Muon_trkError", "std::vector<float>", &T_Muon_trkError);
        mytree_->Branch("T_Muon_dB", "std::vector<float>", &T_Muon_dB);
        mytree_->Branch("T_Muon_dzPV", "std::vector<float>", &T_Muon_dzPV);
        mytree_->Branch("T_Muon_dBstop", "std::vector<float>", &T_Muon_dBstop);
        mytree_->Branch("T_Muon_dzstop", "std::vector<float>", &T_Muon_dzstop);
        mytree_->Branch("T_Muon_chargedHadronIsoR04", "std::vector<float>", &T_Muon_chargedHadronIsoR04);
        mytree_->Branch("T_Muon_neutralHadronIsoR04", "std::vector<float>", &T_Muon_neutralHadronIsoR04);
        mytree_->Branch("T_Muon_photonIsoR04", "std::vector<float>", &T_Muon_photonIsoR04);
        mytree_->Branch("T_Muon_chargedHadronIsoPUR04", "std::vector<float>", &T_Muon_chargedHadronIsoPUR04);
        mytree_->Branch("T_Muon_chargedHadronIsoR03", "std::vector<float>", &T_Muon_chargedHadronIsoR03);
        mytree_->Branch("T_Muon_neutralHadronIsoR03", "std::vector<float>", &T_Muon_neutralHadronIsoR03);
        mytree_->Branch("T_Muon_photonIsoR03", "std::vector<float>", &T_Muon_photonIsoR03);
        mytree_->Branch("T_Muon_chargedHadronIsoPUR03", "std::vector<float>", &T_Muon_chargedHadronIsoPUR03);
        mytree_->Branch("T_Muon_isoR03_emEt", "std::vector<float>", &T_Muon_isoR03_emEt);
        mytree_->Branch("T_Muon_isoR03_hadEt", "std::vector<float>", &T_Muon_isoR03_hadEt);
        mytree_->Branch("T_Muon_isoR03_hoEt", "std::vector<float>", &T_Muon_isoR03_hoEt);
        mytree_->Branch("T_Muon_isoR03_sumPt", "std::vector<float>", &T_Muon_isoR03_sumPt);
        mytree_->Branch("T_Muon_isoR03_nTracks", "std::vector<int>", &T_Muon_isoR03_nTracks);
        mytree_->Branch("T_Muon_isoR03_nJets", "std::vector<int>", &T_Muon_isoR03_nJets);
        mytree_->Branch("T_Muon_isoRingsMVA", "std::vector<float>", &T_Muon_isoRingsMVA);
        
        
        mytree_->Branch("T_Muon_ecalPFiso", "std::vector<float>", &T_Muon_ecalPFiso);
        mytree_->Branch("T_Muon_hcalPFiso", "std::vector<float>", &T_Muon_hcalPFiso);
    }
    
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
    
    T_Event_Rho = new std::vector<float>;

    
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
    T_Gen_Px = new std::vector<float>;
    T_Gen_Py = new std::vector<float>;
    T_Gen_Pz = new std::vector<float>;
    T_Gen_Energy = new std::vector<float>;
    T_Gen_pdgID = new std::vector<int>;
    T_Gen_MotherID = new std::vector<int>;
    T_Gen_FromW = new std::vector<int>;
    T_Gen_FromTau = new std::vector<int>;
    
    T_Muon_Eta = new std::vector<float>;
    T_Muon_Phi = new std::vector<float>;
    T_Muon_Energy = new std::vector<float>;
    T_Muon_Et = new std::vector<float>;
    T_Muon_Pt = new std::vector<float>;
    T_Muon_Px = new std::vector<float>;
    T_Muon_Py = new std::vector<float>;
    T_Muon_Pz = new std::vector<float>;
    T_Muon_Mass = new std::vector<float>;
    T_Muon_IsGlobalMuon = new std::vector<bool>;
    T_Muon_IsTrackerMuon = new std::vector<bool>;
    T_Muon_IsPFMuon = new std::vector<bool>;
    T_Muon_IsCaloMuon = new std::vector<bool>;
    T_Muon_IsStandAloneMuon = new std::vector<bool>;
    T_Muon_IsMuon = new std::vector<bool>;
    T_Muon_IsGlobalMuon_PromptTight = new std::vector<bool>;
    T_Muon_IsTrackerMuonArbitrated = new std::vector<bool>;
    T_Muon_numberOfChambers = new std::vector<int>;
    T_Muon_numberOfChambersRPC = new std::vector<int>;
    T_Muon_numberOfMatches = new std::vector<int>;
    T_Muon_numberOfMatchedStations = new std::vector<int>;
    T_Muon_charge = new std::vector<int>;
    T_Muon_TMLastStationTight = new std::vector<bool>;
    T_Muon_globalTrackChi2 = new std::vector<float>;
    T_Muon_validMuonHits = new std::vector<int>;
    T_Muon_trkKink = new std::vector<float>;
    T_Muon_trkNbOfTrackerLayers = new std::vector<int>;
    T_Muon_trkNbOfValidTrackeHits = new std::vector<int>;
    T_Muon_trkValidPixelHits = new std::vector<int>;
    T_Muon_trkError = new std::vector<float>;
    T_Muon_dB = new std::vector<float>;
    T_Muon_dzPV = new std::vector<float>;
    T_Muon_dBstop = new std::vector<float>;
    T_Muon_dzstop = new std::vector<float>;
    T_Muon_chargedHadronIsoR04 = new std::vector<float>;
    T_Muon_neutralHadronIsoR04 = new std::vector<float>;
    T_Muon_photonIsoR04 = new std::vector<float>;
    T_Muon_chargedHadronIsoPUR04 = new std::vector<float>;
    T_Muon_chargedHadronIsoR03 = new std::vector<float>;
    T_Muon_neutralHadronIsoR03 = new std::vector<float>;
    T_Muon_photonIsoR03 = new std::vector<float>;
    T_Muon_chargedHadronIsoPUR03 = new std::vector<float>;
    T_Muon_isoR03_emEt = new std::vector<float>;
    T_Muon_isoR03_hadEt = new std::vector<float>;
    T_Muon_isoR03_hoEt = new std::vector<float>;
    T_Muon_isoR03_sumPt = new std::vector<float>;
    T_Muon_isoR03_nTracks = new std::vector<int>;
    T_Muon_isoR03_nJets = new std::vector<int>;
    T_Muon_isoRingsMVA = new std::vector<float>;
    
    T_Muon_ecalPFiso = new std::vector<float>;
    T_Muon_hcalPFiso = new std::vector<float>;
    
    
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
    
    delete T_Event_Rho;

    
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
    delete T_Gen_Px;
    delete T_Gen_Py;
    delete T_Gen_Pz;
    delete T_Gen_Energy;
    
    delete T_Gen_pdgID;
    delete T_Gen_MotherID;
    delete T_Gen_FromW;
    delete T_Gen_FromTau;
    
    
    delete T_Muon_Eta;
    delete T_Muon_Phi;
    delete T_Muon_Energy;
    delete T_Muon_Et;
    delete T_Muon_Pt;
    delete T_Muon_Px;
    delete T_Muon_Py;
    delete T_Muon_Pz;
    delete T_Muon_Mass;
    delete T_Muon_IsGlobalMuon;
    delete T_Muon_IsTrackerMuon;
    delete T_Muon_IsPFMuon;
    delete T_Muon_IsCaloMuon;
    delete T_Muon_IsStandAloneMuon;
    delete T_Muon_IsMuon;
    delete T_Muon_IsGlobalMuon_PromptTight;
    delete T_Muon_IsTrackerMuonArbitrated;
    delete T_Muon_numberOfChambers;
    delete T_Muon_numberOfChambersRPC;
    delete T_Muon_numberOfMatches;
    delete T_Muon_numberOfMatchedStations;
    delete T_Muon_charge;
    delete T_Muon_TMLastStationTight;
    delete T_Muon_globalTrackChi2;
    delete T_Muon_validMuonHits;
    delete T_Muon_trkKink;
    delete T_Muon_trkNbOfTrackerLayers;
    delete T_Muon_trkNbOfValidTrackeHits;
    delete T_Muon_trkValidPixelHits;
    delete T_Muon_trkError;
    delete T_Muon_dB;
    delete T_Muon_dzPV;
    delete T_Muon_dBstop;
    delete T_Muon_dzstop;
    delete T_Muon_chargedHadronIsoR04;
    delete T_Muon_neutralHadronIsoR04;
    delete T_Muon_photonIsoR04;
    delete T_Muon_chargedHadronIsoPUR04;
    delete T_Muon_chargedHadronIsoR03;
    delete T_Muon_neutralHadronIsoR03;
    delete T_Muon_photonIsoR03;
    delete T_Muon_chargedHadronIsoPUR03;
    delete T_Muon_isoR03_emEt;
    delete T_Muon_isoR03_hadEt;
    delete T_Muon_isoR03_hoEt;
    delete T_Muon_isoR03_sumPt;
    delete T_Muon_isoR03_nTracks;
    delete T_Muon_isoR03_nJets;
    delete T_Muon_isoRingsMVA;
    
    delete T_Muon_ecalPFiso;
    delete T_Muon_hcalPFiso;
    
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
