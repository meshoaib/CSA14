import FWCore.ParameterSet.Config as cms

process = cms.Process("MiniAnalyzer")

process.load("FWCore.MessageService.MessageLogger_cfi")

from UserCode.TopAnalysis.TTJets_MG_PU_S14_POSTLS170_cfi import source as TTJets_MG_PU_S14_POSTLS170_Source

process.source = TTJets_MG_PU_S14_POSTLS170_Source

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(100) )

process.load('UserCode.TopAnalysis.miniAnalyzer_cfi')

process.p = cms.Path(process.miniAnalyzer)
