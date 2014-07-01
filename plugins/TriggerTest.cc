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
    outputFile_   = iConfig.getParameter<std::string>("outputFile");
    rootFile_ = TFile::Open(outputFile_.c_str(),"RECREATE");

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
    
    
    //if (mapsValues_.size() != filterToMatch_.size()) cout << "warning MAP and filters are not matching 1 to 1 !!!" << endl;
    
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
    
    
    
    
    
    ///loop on the trigger objects
  /*  cout << "loop on objects " << legObjects.size() << endl;
    for (size_t t = 0 ; t < legObjects.size() ; t++){
        float valueToStore = 0;
        cout << "hello an object :), eta=" << legObjects[t].eta() << " phi=" << legObjects[t].phi() << "type =" << filterToMatch_.at(legRefs[t]) << endl;
        TString nameOfFilter = filterToMatch_.at(legRefs[t]);
        if (nameOfFilter.Contains("EtFilter")){
            edm::Handle<trigger::TriggerFilterObjectWithRefs>  localFilter;
            iEvent.getByLabel(edm::InputTag(nameOfFilter,"",HLTprocess_.c_str()), localFilter);
            cout << "on a le filter EtFilter" << localFilter.isValid() << endl;
            edm::Handle<recoEcalCandidateMap> clusterShapeMap;
            iEvent.getByLabel(edm::InputTag(mapsValues_.at(legRefs[t]),"",HLTprocess_.c_str()),clusterShapeMap);
        }
        cout << "value to store=" << valueToStore << endl;
    }*/
    TrigFiltVect hltFilters( filterToMatch_.size() );
    for (unsigned int i = 0 ; i <  filterToMatch_.size() ; i++){
        // cout << "i=" << i << endl;
        iEvent.getByLabel(edm::InputTag(filterToMatch_.at(i),"",HLTprocess_.c_str()), hltFilters[i]);
    }
   
    for (unsigned int i = 0 ; i <  filterToMatch_.size() ; i++){
        if (!(hltFilters[i].isValid())) continue;
        TString nameOfFilter = filterToMatch_.at(i);
     //   cout << "filter to match=" << nameOfFilter << endl;
        if ((nameOfFilter.Contains("EtFilter"))||(nameOfFilter.Contains("ClusterShapeFilter"))||(nameOfFilter.Contains("EcalIsoFilter"))||(nameOfFilter.Contains("HEFilter"))){
        //if ((nameOfFilter.Contains("ClusterShapeFilter"))){
            edm::LogVerbatim("TriggerTest") << "et filter !!";
            edm::LogVerbatim("TriggerTest") << nameOfFilter;
            edm::Handle<recoEcalCandidateMap> clusterShapeMap;
            iEvent.getByLabel(edm::InputTag(mapsValues_.at(i),"",HLTprocess_.c_str()),clusterShapeMap);
            edm::LogVerbatim("TriggerTest") << "is the map valid" << clusterShapeMap.isValid();
            std::vector<edm::Ref<reco::RecoEcalCandidateCollection> > candsFilterEtLead;
            edm::Ref<reco::RecoEcalCandidateCollection> ref;
            hltFilters[i]->getObjects(trigger::TriggerCluster, candsFilterEtLead);
            for (size_t j = 0 ; j < candsFilterEtLead.size() ; j++){
                ref = candsFilterEtLead[j];
                T_Trig_Eta->push_back(ref->eta());
                T_Trig_Pt->push_back(ref->pt());
                T_Trig_Phi->push_back(ref->phi());
                T_Trig_Leg->push_back(i);
                if (clusterShapeMap.isValid()){
                    recoEcalCandidateMap::const_iterator valSigEta = (*clusterShapeMap).find(ref);
                    T_Trig_Value->push_back(valSigEta->val);
                }
                else T_Trig_Value->push_back(-1);
                T_Trig_Value2->push_back(-1);
            }
        }
        if (nameOfFilter.Contains("HcalIsoFilter")){
            edm::LogVerbatim("TriggerTest") << "HCAL iso";
            edm::LogVerbatim("TriggerTest") << nameOfFilter;
            edm::Handle<reco::ElectronSeedCollection> hltElectronSeed;
            iEvent.getByLabel(edm::InputTag(mapsValues_.at(i),"",HLTprocess_.c_str()),hltElectronSeed);
            edm::LogVerbatim("TriggerTest") << "is the map valid" << hltElectronSeed.isValid();
            std::vector<edm::Ref<reco::RecoEcalCandidateCollection> > candsFilterEtLead;
            edm::Ref<reco::RecoEcalCandidateCollection> ref;
            hltFilters[i]->getObjects(trigger::TriggerCluster, candsFilterEtLead);
            for (size_t j = 0 ; j < candsFilterEtLead.size() ; j++){
                ref = candsFilterEtLead[j];
                T_Trig_Eta->push_back(ref->eta());
                T_Trig_Pt->push_back(ref->pt());
                T_Trig_Phi->push_back(ref->phi());
                T_Trig_Leg->push_back(i);
                reco::SuperClusterRef recr2 = ref->superCluster();
                if (hltElectronSeed.isValid()){
                    int nmatch = 0;
                    for(reco::ElectronSeedCollection::const_iterator it = hltElectronSeed->begin(); it != hltElectronSeed->end(); it++){
                        edm::RefToBase<reco::CaloCluster> caloCluster = it->caloCluster() ;
                        reco::SuperClusterRef scRef = caloCluster.castTo<reco::SuperClusterRef>() ;
                        // std::cout<<"BB ref, e, eta, phi"<<&(*scRef)<<" "<<scRef->energy()<<" "<<scRef->eta()<<" "<<scRef->phi()<<std::endl;
                        if(&(*recr2) ==  &(*scRef)) {
                            nmatch++;
                        }
                    }
                    T_Trig_Value->push_back(nmatch);
                }
                else T_Trig_Value->push_back(-1);
                T_Trig_Value2->push_back(-1);
            }
        }
        if ((nameOfFilter.Contains("PixelMatchFilter"))){
            //if ((nameOfFilter.Contains("ClusterShapeFilter"))){
            edm::LogVerbatim("TriggerTest") << "pixel match";
            edm::LogVerbatim("TriggerTest") << nameOfFilter;
            edm::Handle<reco::ElectronCollection> hltElectronColl;
            iEvent.getByLabel(edm::InputTag(mapsValues_.at(i),"",HLTprocess_.c_str()),hltElectronColl);
            edm::LogVerbatim("TriggerTest") << "is the map valid" << hltElectronColl.isValid();
            std::vector<edm::Ref<reco::RecoEcalCandidateCollection> > candsFilterEtLead;
            edm::Ref<reco::RecoEcalCandidateCollection> ref;
            hltFilters[i]->getObjects(trigger::TriggerCluster, candsFilterEtLead);
            for (size_t j = 0 ; j < candsFilterEtLead.size() ; j++){
                ref = candsFilterEtLead[j];
                T_Trig_Eta->push_back(ref->eta());
                T_Trig_Pt->push_back(ref->pt());
                T_Trig_Phi->push_back(ref->phi());
                T_Trig_Leg->push_back(i);
                reco::SuperClusterRef recr2 = ref->superCluster();
                if (hltElectronColl.isValid()){
                    float candValEoP = -1;
                    for(reco::ElectronCollection::const_iterator iElectron = hltElectronColl->begin(); iElectron != hltElectronColl->end(); iElectron++){
                        reco::ElectronRef electronref(reco::ElectronRef(hltElectronColl,iElectron - hltElectronColl->begin()));
                        const reco::SuperClusterRef theClus = electronref->superCluster();
                        if(&(*recr2) ==  &(*theClus)) {
                            const math::XYZVector trackMom =  electronref->track()->momentum();
                            if( trackMom.R() != 0) candValEoP = fabs((1/electronref->superCluster()->energy()) - (1/trackMom.R()));
                        }
                    }
                    T_Trig_Value->push_back(candValEoP);
                }
                else T_Trig_Value->push_back(-1);
                T_Trig_Value2->push_back(-1);
            }
        }
        if ((nameOfFilter.Contains("OneOEMinusOneOPFilter"))){
            //if ((nameOfFilter.Contains("ClusterShapeFilter"))){
            edm::LogVerbatim("TriggerTest") << "OneOOE !!";
            edm::LogVerbatim("TriggerTest") << nameOfFilter;
            edm::Handle<reco::ElectronIsolationMap> hltElectrondPhiMap;
            iEvent.getByLabel(edm::InputTag(mapsValues_.at(i),"Deta",HLTprocess_.c_str()),hltElectrondPhiMap);
            edm::Handle<reco::ElectronIsolationMap> hltElectrondEtaMap;
            iEvent.getByLabel(edm::InputTag(mapsValues_.at(i),"Dphi",HLTprocess_.c_str()),hltElectrondEtaMap);
            edm::LogVerbatim("TriggerTest") << "is the map valid" << hltElectrondPhiMap.isValid() << " " << hltElectrondEtaMap.isValid();
            std::vector<edm::Ref<reco::ElectronCollection> > candsFilterEtLead;
            edm::Ref<reco::ElectronCollection> ref;
            hltFilters[i]->getObjects(trigger::TriggerElectron, candsFilterEtLead);
            for (size_t j = 0 ; j < candsFilterEtLead.size() ; j++){
                ref = candsFilterEtLead[j];
                T_Trig_Eta->push_back(ref->eta());
                T_Trig_Pt->push_back(ref->pt());
                T_Trig_Phi->push_back(ref->phi());
                T_Trig_Leg->push_back(i);
                if (hltElectrondPhiMap.isValid()){
                    reco::ElectronIsolationMap::const_iterator valSigEta = (*hltElectrondPhiMap).find(ref);
                    T_Trig_Value->push_back(valSigEta->val);
                }
                else T_Trig_Value->push_back(-1);
                if (hltElectrondEtaMap.isValid()){
                    reco::ElectronIsolationMap::const_iterator valSigEta = (*hltElectrondEtaMap).find(ref);
                    T_Trig_Value2->push_back(valSigEta->val);
                }
                else T_Trig_Value2->push_back(-1);
            }
        }
        if ((nameOfFilter.Contains("OneOEMinusOneOPFilter"))){
            //if ((nameOfFilter.Contains("ClusterShapeFilter"))){
            edm::LogVerbatim("TriggerTest") << "OneOOE !!";
            edm::LogVerbatim("TriggerTest") << nameOfFilter;
            edm::Handle<reco::ElectronIsolationMap> hltElectrondPhiMap;
            iEvent.getByLabel(edm::InputTag(mapsValues_.at(i),"Deta",HLTprocess_.c_str()),hltElectrondPhiMap);
            edm::Handle<reco::ElectronIsolationMap> hltElectrondEtaMap;
            iEvent.getByLabel(edm::InputTag(mapsValues_.at(i),"Dphi",HLTprocess_.c_str()),hltElectrondEtaMap);
            edm::LogVerbatim("TriggerTest") << "is the map valid" << hltElectrondPhiMap.isValid() << " " << hltElectrondEtaMap.isValid();
            std::vector<edm::Ref<reco::ElectronCollection> > candsFilterEtLead;
            edm::Ref<reco::ElectronCollection> ref;
            hltFilters[i]->getObjects(trigger::TriggerElectron, candsFilterEtLead);
            for (size_t j = 0 ; j < candsFilterEtLead.size() ; j++){
                ref = candsFilterEtLead[j];
                T_Trig_Eta->push_back(ref->eta());
                T_Trig_Pt->push_back(ref->pt());
                T_Trig_Phi->push_back(ref->phi());
                T_Trig_Leg->push_back(i);
                if (hltElectrondPhiMap.isValid()){
                    reco::ElectronIsolationMap::const_iterator valSigEta = (*hltElectrondPhiMap).find(ref);
                    T_Trig_Value->push_back(valSigEta->val);
                }
                else T_Trig_Value->push_back(-1);
                if (hltElectrondEtaMap.isValid()){
                    reco::ElectronIsolationMap::const_iterator valSigEta = (*hltElectrondEtaMap).find(ref);
                    T_Trig_Value2->push_back(valSigEta->val);
                }
                else T_Trig_Value2->push_back(-1);
            }
        }
        if ((nameOfFilter.Contains("DphiFilter"))){
            //if ((nameOfFilter.Contains("ClusterShapeFilter"))){
            edm::LogVerbatim("TriggerTest") << "isolation !!";
            edm::LogVerbatim("TriggerTest") << nameOfFilter;
            edm::Handle<reco::ElectronIsolationMap> hltElectrondPhiMap;
            iEvent.getByLabel(edm::InputTag(mapsValues_.at(i),"",HLTprocess_.c_str()),hltElectrondPhiMap);
            edm::LogVerbatim("TriggerTest") << "is the map valid" << hltElectrondPhiMap.isValid();
            std::vector<edm::Ref<reco::ElectronCollection> > candsFilterEtLead;
            edm::Ref<reco::ElectronCollection> ref;
            hltFilters[i]->getObjects(trigger::TriggerElectron, candsFilterEtLead);
            for (size_t j = 0 ; j < candsFilterEtLead.size() ; j++){
                ref = candsFilterEtLead[j];
                T_Trig_Eta->push_back(ref->eta());
                T_Trig_Pt->push_back(ref->pt());
                T_Trig_Phi->push_back(ref->phi());
                T_Trig_Leg->push_back(i);
                if (hltElectrondPhiMap.isValid()){
                    reco::ElectronIsolationMap::const_iterator valSigEta = (*hltElectrondPhiMap).find(ref);
                    T_Trig_Value->push_back(valSigEta->val);
                }
                else T_Trig_Value->push_back(-1);
                T_Trig_Value2->push_back(-1);

            }
        }
        if ((nameOfFilter.Contains("DphiFilter"))){
            //if ((nameOfFilter.Contains("ClusterShapeFilter"))){
            edm::LogVerbatim("TriggerTest") << "last !!";
            edm::LogVerbatim("TriggerTest") << nameOfFilter;

            std::vector<edm::Ref<reco::ElectronCollection> > candsFilterEtLead;
            edm::Ref<reco::ElectronCollection> ref;
            hltFilters[i]->getObjects(trigger::TriggerElectron, candsFilterEtLead);
            for (size_t j = 0 ; j < candsFilterEtLead.size() ; j++){
                ref = candsFilterEtLead[j];
                T_Trig_Eta->push_back(ref->eta());
                T_Trig_Pt->push_back(ref->pt());
                T_Trig_Phi->push_back(ref->phi());
                T_Trig_Leg->push_back(i);
                T_Trig_Value->push_back(-1);
                T_Trig_Value2->push_back(-1);
            }
        }
      /*  if (nameOfFilter.Contains("L3Filtered24Q")){
            //if ((nameOfFilter.Contains("ClusterShapeFilter"))){
          //  edm::LogVerbatim("TriggerTest") << "last !!";
            std::vector<reco::RecoChargedCandidateRef> prevMuonRefs;
            reco::RecoChargedCandidateRef ref;
            hltFilters[i]->getObjects(trigger::TriggerMuon, prevMuonRefs);
            edm::Handle<reco::IsoDepositMap> IsoCaloMap;
            iEvent.getByLabel(edm::InputTag("hltL3CaloMuonCorrectedIsolations","",HLTprocess_.c_str()),IsoCaloMap);
            edm::Handle<reco::IsoDepositMap> IsoTkMap;
            iEvent.getByLabel(edm::InputTag(mapsValues_.at(i),"trkIsoDeposits",HLTprocess_.c_str()),IsoTkMap);
            for (size_t j = 0 ; j < prevMuonRefs.size() ; j++){
                ref = prevMuonRefs[j];
                T_Trig_Eta->push_back(ref->eta());
                T_Trig_Pt->push_back(ref->pt());
                T_Trig_Phi->push_back(ref->phi());
                T_Trig_Leg->push_back(i);
                if (IsoTkMap.isValid()) {
                    reco::IsoDeposit theTkIsolation = (*IsoTkMap)[ref];
                    T_Trig_Value->push_back(theTkIsolation.depositWithin(0.3));
                }
                else T_Trig_Value->push_back(-1);
                if (IsoCaloMap.isValid()){
                    reco::IsoDeposit theCaloIsolation = (*IsoCaloMap)[ref];
                    T_Trig_Value2->push_back(theCaloIsolation.depositWithin(0.3));
                }
                else T_Trig_Value2->push_back(-1);

            }
        }*/
        else if (nameOfFilter.Contains("L3PreFiltered8")||nameOfFilter.Contains("L3Filtered17")||nameOfFilter.Contains("TrkFiltered8")){
            //if ((nameOfFilter.Contains("ClusterShapeFilter"))){
        //    cout << nameOfFilter << endl;
          //  edm::LogVerbatim("TriggerTest") << "last !!";
            std::vector<reco::RecoChargedCandidateRef> prevMuonRefs;
            reco::RecoChargedCandidateRef ref;
            hltFilters[i]->getObjects(trigger::TriggerMuon, prevMuonRefs);
            edm::Handle<reco::IsoDepositMap> IsoCaloMap;
            iEvent.getByLabel(edm::InputTag("hltL3CaloMuonCorrectedIsolations","",HLTprocess_.c_str()),IsoCaloMap);
            edm::Handle<reco::IsoDepositMap> IsoTkMap;
            iEvent.getByLabel(edm::InputTag(mapsValues_.at(i),"trkIsoDeposits",HLTprocess_.c_str()),IsoTkMap);
            for (size_t j = 0 ; j < prevMuonRefs.size() ; j++){
                ref = prevMuonRefs[j];
                T_Trig_Eta->push_back(ref->eta());
                T_Trig_Pt->push_back(ref->pt());
                T_Trig_Phi->push_back(ref->phi());
                T_Trig_Leg->push_back(i);
                if (IsoTkMap.isValid()) {
                    reco::IsoDeposit theTkIsolation = (*IsoTkMap)[ref];
                    T_Trig_Value->push_back(theTkIsolation.depositWithin(0.3));
                }
                else T_Trig_Value->push_back(-1);
      /*          if (IsoCaloMap.isValid()){
                    reco::IsoDeposit theCaloIsolation = (*IsoCaloMap)[ref];
                    T_Trig_Value2->push_back(theCaloIsolation.depositWithin(0.3));
                }*/
               T_Trig_Value2->push_back(-1);

            }
        }
        else {//if ((nameOfFilter.Contains("L3crIsoRhoFiltered0p15")||(nameOfFilter.Contains("L1Filtered"))||(nameOfFilter.Contains("L2Filtered"))||(nameOfFilter.Contains("L2")))){
            //if ((nameOfFilter.Contains("ClusterShapeFilter"))){
            //  edm::LogVerbatim("TriggerTest") << "last !!";
            std::vector<reco::RecoChargedCandidateRef> prevMuonRefs;
            reco::RecoChargedCandidateRef ref;
            hltFilters[i]->getObjects(trigger::TriggerMuon, prevMuonRefs);
            for (size_t j = 0 ; j < prevMuonRefs.size() ; j++){
                ref = prevMuonRefs[j];
                T_Trig_Eta->push_back(ref->eta());
                T_Trig_Pt->push_back(ref->pt());
                T_Trig_Phi->push_back(ref->phi());
                T_Trig_Leg->push_back(i);
                T_Trig_Value->push_back(-1);
                T_Trig_Value2->push_back(-1);
            }
        }
    }
    
    
    
   
    /*int theNbOfGen = genParticles->size();
    for (int i=0 ; i < theNbOfGen; i++){
        const reco::GenParticle & genParticle = (*genParticles)[i];
        if (genParticle.status()>3) continue; 
        cout << "part found=  pt=" << genParticle.pt() << " eta=" << genParticle.eta();
	cout << " pdgID=" << genParticle.pdgId() << " status=" << genParticle.status() << endl; 
    
    
    }*/
    
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
    
    
  /*  for(reco::GsfElectronCollection::const_iterator eleIt = electronsCollection->begin(); eleIt != electronsCollection->end(); eleIt++){
        //cout << "pt=" << eleIt->pt() << endl;
        T_Elec_Eta->push_back(eleIt->eta());
        T_Elec_Phi->push_back(eleIt->phi());
        T_Elec_Pt->push_back(eleIt->pt());
        
    }
    
    int nbMuons = recoMuons->size();
    for (int k = 0 ; k < nbMuons ; k++){
        
        const reco::Muon* muon = &((*recoMuons)[k]);
        
        T_Muon_Eta->push_back(muon->eta());
        T_Muon_Phi->push_back(muon->phi());
        T_Muon_Pt->push_back(muon->pt());
    
    }*/
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
    mytree_->Branch("T_Trig_Value", "std::vector<float>", &T_Trig_Value);
    mytree_->Branch("T_Trig_Value2", "std::vector<float>", &T_Trig_Value2);
    
    
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
    T_Trig_Value = new std::vector<float>;
    T_Trig_Value2 = new std::vector<float>;
    
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
    delete T_Trig_Value;
    delete T_Trig_Value2;
    
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
