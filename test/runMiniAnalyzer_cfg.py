import FWCore.ParameterSet.Config as cms

process = cms.Process("Demo")

process.load("FWCore.MessageService.MessageLogger_cfi")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(3000000))

from UserCode.TopAnalysis.csa14.TT_PU_v2_cfi import source as TT_source
#from UserCode.TopAnalysis.csa14.TT_wjets_cfi import source as TT_source
#from UserCode.TopAnalysis.csa14.TT_w1234jets_V5_v1_cfi import source as TT_source
#from UserCode.TopAnalysis.csa14.TT_w1234jets_V5_v2_cfi import source as TT_source 
#from UserCode.TopAnalysis.csa14.TT_DY_cfi import source as TT_source
#from UserCode.TopAnalysis.csa14.QCD_80_120_MuEnriched_pythia8_cfi import source as TT_source
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

process.source=TT_source

#reduce verbosity
process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.threshold = ''
process.MessageLogger.cerr.FwkReport.reportEvery = 1000

#tfileservice
process.TFileService = cms.Service("TFileService",
        fileName = cms.string("TT_PU_v2.root")
#        fileName = cms.string("wjets.root")
#        fileName = cms.string("w1234jets_v5_v1.root")
#        fileName = cms.string("w1234jets_v5_v2.root") 
#        fileName = cms.string("DY.root")
#        fileName = cms.string("QCD_80_120_MuEnriched.root")
#        fileName = cms.string("QCD_170_300_MuEnriched.root")
#        fileName = cms.string("QCD_300_470_MuEnriched.root")
#        fileName = cms.string("QCD_470_600_MuEnriched.root")
#        fileName = cms.string("QCD_600_800_MuEnriched.root")
#        fileName = cms.string("QCD_800_1000_MuEnriched.root")
#        fileName = cms.string("QCD_1000_MuEnriched.root")
#        fileName = cms.string("QCD_1000_1400_Tune4C.root")
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
process.load('UserCode.TopAnalysis.miniAnalyzer_cfi')
process.p = cms.Path(process.demo)


