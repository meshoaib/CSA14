### Installation

cmsrel CMSSW_7_0_6_patch1
cd CMSSW_7_0_6_patch1/src
cmsenv
git clone git@github.com:pfs/CSA14.git UserCode/TopAnalysis
scram b -j 9


### Running

cmsRun UserCode/TopAnalysis/test/runMiniAnalyzer_cfg.py

Edit cfg for different miniaod sources

Put sources in python according to py generated from DAS


