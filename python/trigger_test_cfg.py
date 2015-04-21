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
                                             doMuon               = cms.bool(False),
                                             primaryVertexInputTag   	= cms.InputTag("offlinePrimaryVertices","","RECO"),
                                             muonProducer 	         	= cms.VInputTag(cms.InputTag("muons")),
                                             triggerResultTag     = cms.InputTag("TriggerResults", "", "HLTfancy"),
                                             triggerSummaryTag    = cms.InputTag("hltTriggerSummaryAOD", "", "HLTfancy"),
                                             pathsToSave          = cms.vstring("HLT_IsoMu20_v1",
                                                                                "HLT_IsoMu24_eta2p1_v1",
                                                                                ),
                                             filterToMatch         = cms.vstring("hltL3fL1sMu16L1f0L2f10QL3Filtered20Q"),
                                             mapsValues            = cms.vstring(""),
                                             HLTprocess            = cms.string("HLTfancy"),
                                             rhoTag                = cms.InputTag("hltFixedGridRhoFastjetAllCaloForMuons","","HLTfancy"),
                                             TrkExtractorPSet = cms.PSet(
                                                                         DR_VetoPt = cms.double( 0.025 ),
                                                                         Diff_z = cms.double( 0.2 ),
                                                                         inputTrackCollection = cms.InputTag( "hltIter2L3MuonMerged" ),
                                                                         ReferenceRadius = cms.double( 6.0 ),
                                                                         BeamSpotLabel = cms.InputTag( "hltOnlineBeamSpot" ),
                                                                         ComponentName = cms.string( "PixelTrackExtractor" ),
                                                                         DR_Max = cms.double( 0.3 ),
                                                                         Diff_r = cms.double( 0.1 ),
                                                                         PropagateTracksToRadius = cms.bool( True ),
                                                                         Chi2Prob_Min = cms.double( -1.0 ),
                                                                         DR_Veto = cms.double( 0.01 ),
                                                                         NHits_Min = cms.uint32( 0 ),
                                                                         Chi2Ndof_Max = cms.double( 1.0E64 ),
                                                                         Pt_Min = cms.double( -1.0 ),
                                                                         DepositLabel = cms.untracked.string( "PXLS" ),
                                                                         BeamlineOption = cms.string( "BeamSpotFromEvent" ),
                                                                         VetoLeadingTrack = cms.bool( True ),
                                                                         PtVeto_Min = cms.double( 2.0 )
                                                                         ),
                                             outputFile		        = cms.string("triggerTree.root")
                                             )

process.p = cms.Path(process.triggerTestProducer)

