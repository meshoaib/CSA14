from WMCore.Configuration import Configuration
config = Configuration()

config.section_("General")
config.General.requestName = 'CSA14_MC'
config.General.workArea = 'crab_projects'
config.General.transferOutput = False

config.section_("JobType")
config.JobType.pluginName = 'Analysis'
config.JobType.psetName = 'runMiniAnalyzer_cfg.py'

config.section_("Data")
config.Data.inputDataset = '/TTJets_MSDecaysCKM_central_Tune4C_13TeV-madgraph-tauola/Spring14miniaod-PU20bx25_POSTLS170_V5-v2/MINIAODSIM'
config.Data.dbsUrl = 'global'
config.Data.splitting = 'FileBased'
config.Data.unitsPerJob = 10
config.Data.publication = False
#config.Data.publishDbsUrl = 'phys03'
#config.Data.publishDataName = 'CRAB3_tutorial_MC_analysis_test1'

config.section_("Site")
config.Site.storageSite = 'T2_CH_CERN'
