import FWCore.ParameterSet.Config as cms

miniAnalyzer = cms.EDAnalyzer("MiniAnalyzer",
                              vertices = cms.InputTag("offlineSlimmedPrimaryVertices"),
                              muons = cms.InputTag("slimmedMuons"),
                              electrons = cms.InputTag("slimmedElectrons"),
                              jets = cms.InputTag("slimmedJets"),
                              mets = cms.InputTag("slimmedMETs")
                              )
