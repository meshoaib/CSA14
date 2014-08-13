import FWCore.ParameterSet.Config as cms

demo = cms.EDAnalyzer("MiniAnalyzer",
                      triggerBits = cms.InputTag("TriggerResults","","HLT"),
                      rho         = cms.InputTag("fixedGridRhoFastjetAll"),
                      vertices    = cms.InputTag("offlineSlimmedPrimaryVertices"),
                      muons       = cms.InputTag("slimmedMuons"),
                      electrons   = cms.InputTag("slimmedElectrons"),
                      jets        = cms.InputTag("slimmedJets"),
                      mets        = cms.InputTag("slimmedMETs")
                      )

