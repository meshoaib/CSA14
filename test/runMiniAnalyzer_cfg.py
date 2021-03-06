import FWCore.ParameterSet.Config as cms

process = cms.Process("Demo")

process.load("FWCore.MessageService.MessageLogger_cfi")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(10000))


#<<<<<<< HEAD
#from UserCode.TopAnalysis.csa14.TT_PU_S14_V5_cfi import source as TT_source
#from UserCode.TopAnalysis.csa14.TT_PU_S14_V7_cfi import source as TT_source
#from UserCode.TopAnalysis.csa14.TT_PU_v2_cfi import source as TT_source
#from UserCode.TopAnalysis.csa14.TT_PU_S14_cfi import source as TT_source
#from UserCode.TopAnalysis.csa14.TTJets_MG_PU20bx25_POSTLS170_cfi import source as TT_source
#from UserCode.TopAnalysis.csa14.TT_wjets_cfi import source as TT_source
#from UserCode.TopAnalysis.csa14.TT_w1234jets_V5_v1_cfi import source as TT_source
#from UserCode.TopAnalysis.csa14.TT_w1234jets_V5_v2_cfi import source as TT_source 
#from UserCode.TopAnalysis.csa14.TT_DY_cfi import source as TT_source
#from UserCode.TopAnalysis.csa14.QCD_80_120_MuEnriched_pythia8_cfi import source as TT_source
#from UserCode.TopAnalysis.csa14.QCD_120_170_MuEnriched_pythia8_cfi import source as TT_source
#from UserCode.TopAnalysis.csa14.QCD_170_300_MuEnriched_pythia8_cfi import source as TT_source
#from UserCode.TopAnalysis.csa14.QCD_300_470_MuEnriched_pythia8_cfi import source as TT_source
#from UserCode.TopAnalysis.csa14.QCD_470_600_MuEnriched_pythia8_cfi import source as TT_source
#from UserCode.TopAnalysis.csa14.QCD_600_800_MuEnriched_pythia8_cfi import source as TT_source
#from UserCode.TopAnalysis.csa14.QCD_800_1000_MuEnriched_pythia8_cfi import source as TT_source
#from UserCode.TopAnalysis.csa14.QCD_1000_MuEnriched_pythia8_cfi import source as TT_source
#from UserCode.TopAnalysis.csa14.QCD_1000_1400_Tune4C_pythia8_cfi import source as TT_source
#from UserCode.TopAnalysis.csa14.QCD_1400_1800_Tune4C_pythia8_cfi import source as TT_source
#from UserCode.TopAnalysis.csa14.QCD_1800_Tune4C_pythia8_cfi import source as TT_source
#from UserCode.TopAnalysis.csa14.QCD_15_3000_Tune4C_flat_pythia8_cfi import source as TT_source
#from UserCode.TopAnalysis.csa14.QCD_170_300_Tune4C_pythia8_cfi import source as TT_source 
#from UserCode.TopAnalysis.csa14.QCD_300_470_Tune4C_pythia8_cfi import source as TT_source
#from UserCode.TopAnalysis.csa14.QCD_470_600_Tune4C_pythia8_cfi import source as TT_source
#from UserCode.TopAnalysis.csa14.QCD_600_800_Tune4C_pythia8_cfi import source as TT_source
#from UserCode.TopAnalysis.csa14.QCD_800_1000_Tune4C_pythia8_cfi import source as TT_source 
#from UserCode.TopAnalysis.csa14.QCD_80_120_Tune4C_pythia8_cfi import source as TT_source
#from UserCode.TopAnalysis.csa14.QCD_120_170_Tune4C_pythia8_cfi import source as TT_source
#=======
from UserCode.TopAnalysis.csa14.TT_PU_v2_cfi import source as events_source
#from UserCode.TopAnalysis.csa14.TT_PU_S14_cfi import source as events_source
#from UserCode.TopAnalysis.csa14.TTJets_MG_PU20bx25_POSTLS170_cfi import source as events_source
#from UserCode.TopAnalysis.csa14.TT_wjets_cfi import source as events_source
#from UserCode.TopAnalysis.csa14.TT_w1234jets_V5_v1_cfi import source as events_source
#from UserCode.TopAnalysis.csa14.TT_w1234jets_V5_v2_cfi import source as events_source 
#from UserCode.TopAnalysis.csa14.TT_DY_cfi import source as events_source
#from UserCode.TopAnalysis.csa14.QCD_80_120_MuEnriched_pythia8_cfi import source as events_source
#from UserCode.TopAnalysis.csa14.QCD_120_170_MuEnriched_pythia8_cfi import source as events_source
#from UserCode.TopAnalysis.csa14.QCD_170_300_MuEnriched_pythia8_cfi import source as events_source
#from UserCode.TopAnalysis.csa14.QCD_300_470_MuEnriched_pythia8_cfi import source as events_source
#from UserCode.TopAnalysis.csa14.QCD_470_600_MuEnriched_pythia8_cfi import source as events_source
#from UserCode.TopAnalysis.csa14.QCD_600_800_MuEnriched_pythia8_cfi import source as events_source
#from UserCode.TopAnalysis.csa14.QCD_800_1000_MuEnriched_pythia8_cfi import source as events_source
#from UserCode.TopAnalysis.csa14.QCD_1000_MuEnriched_pythia8_cfi import source as events_source
#from UserCode.TopAnalysis.csa14.QCD_1000_1400_Tune4C_pythia8_cfi import source as events_source
#from UserCode.TopAnalysis.csa14.QCD_1400_1800_Tune4C_pythia8_cfi import source as events_source
#from UserCode.TopAnalysis.csa14.QCD_1800_Tune4C_pythia8_cfi import source as events_source
#from UserCode.TopAnalysis.csa14.QCD_15_3000_Tune4C_flat_pythia8_cfi import source as events_source
#from UserCode.TopAnalysis.csa14.QCD_170_300_Tune4C_pythia8_cfi import source as events_source 
#from UserCode.TopAnalysis.csa14.QCD_300_470_Tune4C_pythia8_cfi import source as events_source
#from UserCode.TopAnalysis.csa14.QCD_470_600_Tune4C_pythia8_cfi import source as events_source
#from UserCode.TopAnalysis.csa14.QCD_600_800_Tune4C_pythia8_cfi import source as events_source
#from UserCode.TopAnalysis.csa14.QCD_800_1000_Tune4C_pythia8_cfi import source as events_source 
#from UserCode.TopAnalysis.csa14.QCD_80_120_Tune4C_pythia8_cfi import source as events_source
#from UserCode.TopAnalysis.csa14.QCD_120_170_Tune4C_pythia8_cfi import source as events_source
#>>>>>>> b32c9b8249b94e9257f9330212c10ab47d682abc

process.source=events_source

#reduce verbosity
process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.threshold = ''
process.MessageLogger.cerr.FwkReport.reportEvery = 1000

#tfileservice
process.TFileService = cms.Service("TFileService",
        fileName = cms.string("pf_combinediso_TT_PU_S14_V5.root")
#        fileName = cms.string("chargediso_TT_PU_S14_V7.root")
#        fileName = cms.string("chargediso_TT_PU20bx25_v2.root")
#        fileName = cms.string("chargediso_TT_S14_PU40bx50.root")
#                                   fileName = cms.string('TTJets_MG_PU20bx25_POSTLS170.root')
#        fileName = cms.string("chargediso_wjets_PU20bx25.root")
#        fileName = cms.string("w1234jets_v5_v1.root")
#        fileName = cms.string("w1234jets_v5_v2.root") 
#        fileName = cms.string("chargediso_DY_PU20bx25.root")
#        fileName = cms.string("chargediso_QCD_80_120_MuEnriched_PU20bx25.root")
#        fileName = cms.string("chargediso_QCD_120_170_MuEnriched_PU20bx25.root")
#        fileName = cms.string("chargediso_QCD_170_300_MuEnriched_PU20bx25.root")
#        fileName = cms.string("chargediso_QCD_300_470_MuEnriched_PU20bx25.root")
#        fileName = cms.string("chargediso_QCD_470_600_MuEnriched_PU20bx25.root")
#        fileName = cms.string("chargediso_QCD_600_800_MuEnriched_PU20bx25.root")
#        fileName = cms.string("chargediso_QCD_800_1000_MuEnriched_PU20bx25.root")
#        fileName = cms.string("chargediso_QCD_1000_MuEnriched_PU20bx25.root")
#        fileName = cms.string("QCD_1000_1400_Tune4C_PU20bx25.root")
#        fileName = cms.string("QCD_1400_1800_Tune4C.root")
#        fileName = cms.string("QCD_1800_Tune4C.root")
#        fileName = cms.string("QCD_15_3000_Tune4C.root")
#        fileName = cms.string("QCD_170_300_Tune4C.root") 
#        fileName = cms.string("QCD_300_470_Tune4C.root")
#        fileName = cms.string("QCD_470_600_Tune4C.root")
#        fileName = cms.string("QCD_600_800_Tune4C.root")
#        fileName = cms.string("QCD_800_1000_Tune4C.root")
#        fileName = cms.string("QCD_80_120_Tune4C.root")
#        fileName = cms.string("QCD_120_170_Tune4C.root")
)

#running sequence
#process.load('UserCode.TopAnalysis.myChargedPFJets_cfi')
process.load('UserCode.TopAnalysis.miniAnalyzer_cfi')
#process.p = cms.Path(process.myChargedPFJets*process.demo)
process.p = cms.Path(process.demo)


