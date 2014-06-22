import FWCore.ParameterSet.Config as cms

process = cms.Process("testTrigger2")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.load("Configuration.Geometry.GeometryIdeal_cff")
process.load("Configuration/StandardSequences/MagneticField_38T_cff")
process.load("Configuration/StandardSequences/FrontierConditions_GlobalTag_cff")
process.load("Configuration.StandardSequences.Reconstruction_cff")
process.load("Configuration.StandardSequences.Services_cff")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

process.source = cms.Source("PoolSource",
                            # replace 'myfile.root' with the source file you want to use
                            fileNames = cms.untracked.vstring(
                                                              #'file:/tmp/hbrun/theMuEnrichedFile.root'
                                                              'file:/tmp/hbrun/theDYRaw.root'
                                                              )
                            )


process.GlobalTag.globaltag = 'PRE_LS171V9A::All'


process.triggerTestProducer = cms.EDAnalyzer('TriggerTest',
            triggerResultTag     = cms.InputTag("TriggerResults", "", "TEST"),
            triggerSummaryTag    = cms.InputTag("hltTriggerSummaryAOD", "", "TEST"),
            pathsToSave           =cms.vstring("HLT_IsoMu24_v", "HLT_Mu24_v","HLT_Mu17_v"),
            filterToMatch           =cms.vstring("hltL1sL1Mu12EG7",
                                                 "hltL1Mu12EG7L1MuFiltered0",
                                                 "hltL1Mu12EG7L3MuFiltered17",
                                                 "hltEG8EtFilterL1Mu12EG7",
                                                 "hltMu17Ele8CaloIdTCaloIsoVLTrkIdVLTrkIsoVLClusterShapeFilter",
                                                 "hltMu17Ele8CaloIdTCaloIsoVLTrkIdVLTrkIsoVLEcalIsoFilter",
                                                 "hltMu17Ele8CaloIdTCaloIsoVLTrkIdVLTrkIsoVLHEFilter",
                                                 "hltMu17Ele8CaloIdTCaloIsoVLTrkIdVLTrkIsoVLHcalIsoFilter",
                                                 "hltMu17Ele8CaloIdTCaloIsoVLTrkIdVLTrkIsoVLPixelMatchFilter",
                                                 "hltMu17Ele8CaloIdTCaloIsoVLTrkIdVLTrkIsoVLOneOEMinusOneOPFilter",
                                                 "hltMu17Ele8CaloIdTCaloIsoVLTrkIdVLTrkIsoVLDetaFilter",
                                                 "hltMu17Ele8CaloIdTCaloIsoVLTrkIdVLTrkIsoVLDphiFilter",
                                                 "hltMu17Ele8CaloIdTCaloIsoVLTrkIdVLTrkIsoVLTrackIsoFilter"),
            mapsValues            = cms.vstring("hltL1SeededHLTClusterShape",
                                                "hltL1SeededPhotonEcalIso",
                                                "hltL1SeededPhotonHcalForHE",
                                                "hltL1SeededPhotonHcalIso",
                                                "hltL1SeededStartUpElectronPixelSeeds",
                                                "hltPixelMatch3HitElectronsL1Seeded",
                                                "hlt3HitElectronL1SeededDetaDphi",
                                                "",
                                                "hltL1Seeded3HitElectronTrackIso",
                                                "",
                                                "hltActivityPhotonClusterShape",
                                                "hltActivityPhotonEcalIso",
                                                "hltActivityPhotonHcalForHE",
                                                "hltActivityPhotonHcalIso",
                                                "hltActivityStartUpElectronPixelSeeds",
                                                "hltPixelMatch3HitElectronsActivity",
                                                "hlt3HitElectronActivityDetaDphi",
                                                "",
                                                "hlt3HitElectronActivityTrackIso",
                                                "",
                                                ""
                                                ),
            HLTprocess            = cms.string("TEST"),
            outputFile		        = cms.string("triggerTree.root")
)


process.p = cms.Path(process.triggerTestProducer)

