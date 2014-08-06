import FWCore.ParameterSet.Config as cms

process = cms.Process("MiniAnalyzer")

#reduce verbosity
process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.threshold = ''
process.MessageLogger.cerr.FwkReport.reportEvery = 1000

#source of events
#sourceName="TTJets_MG_PU_S14_POSTLS170"
sourceName="TTJets_MG_PU20bx25_POSTLS170"
if sourceName=="TTJets_MG_PU_S14_POSTLS170"   : process.load('UserCode.TopAnalysis.TTJets_MG_PU_S14_POSTLS170_cfi')
if sourceName=="TTJets_MG_PU20bx25_POSTLS170" : process.load('UserCode.TopAnalysis.TTJets_MG_PU20bx25_POSTLS170_cfi')
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(5000) )

#tfileservice
process.TFileService = cms.Service("TFileService", 
                                   fileName = cms.string("%s.root"%sourceName),
                                   closeFileFast = cms.untracked.bool(True)
                                   )

#analyzer
process.load('UserCode.TopAnalysis.miniAnalyzer_cfi')

#running sequence
process.p = cms.Path(process.miniAnalyzer)
