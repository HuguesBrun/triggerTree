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
                                                              'file:/tmp/hbrun/5.4_ZEE+ZEEFS+HARVESTFS/step1.root'
                                                              )
                            )


process.GlobalTag.globaltag = 'PRE_LS171V9A::All'


process.triggerTestProducer = cms.EDAnalyzer('TriggerTest',
            triggerResultTag     = cms.InputTag("TriggerResults", "", "HLT"),
            triggerSummaryTag    = cms.InputTag("hltTriggerSummaryAOD", "", "HLT"),
            pathsToSave           =cms.vstring("HLT_Mu17_Mu8_v"),
            filterToMatch           =cms.vstring("hltL1sL1DoubleMu10MuOpenORDoubleMu103p5",
                                                 "hltL2pfL1DoubleMu10MuOpenOR3p5L1f0L2PreFiltered0",
                                                 "hltL2fL1DoubleMu10MuOpenOR3p5L1f0L2Filtered10",
                                                 "hltL3pfL1DoubleMu10MuOpenOR3p5L1f0L2pf0L3PreFiltered8",
                                                 "hltL3fL1DoubleMu10MuOpenOR3p5L1f0L2f10L3Filtered17",
                                                 "hltDiMuonGlb17Glb8DzFiltered0p2"),
            mapsValues            = cms.vstring(),
            HLTprocess            = cms.string("HLT"),
            outputFile		        = cms.string("triggerTree.root")
)


process.p = cms.Path(process.triggerTestProducer)

