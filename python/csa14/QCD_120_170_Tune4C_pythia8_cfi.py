import FWCore.ParameterSet.Config as cms

maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )
readFiles = cms.untracked.vstring()
secFiles = cms.untracked.vstring() 
source = cms.Source ("PoolSource",fileNames = readFiles, secondaryFileNames = secFiles)
readFiles.extend( [
       '/store/mc/Spring14miniaod/QCD_Pt-120to170_Tune4C_13TeV_pythia8/MINIAODSIM/castor_PU20bx25_POSTLS170_V5-v1/00000/008C7009-9006-E411-A888-002618943842.root',
       '/store/mc/Spring14miniaod/QCD_Pt-120to170_Tune4C_13TeV_pythia8/MINIAODSIM/castor_PU20bx25_POSTLS170_V5-v1/00000/048AE5B8-7D06-E411-9906-0025905A60F8.root',
       '/store/mc/Spring14miniaod/QCD_Pt-120to170_Tune4C_13TeV_pythia8/MINIAODSIM/castor_PU20bx25_POSTLS170_V5-v1/00000/06D72FE5-8F06-E411-B55B-0025905A60E4.root',
       '/store/mc/Spring14miniaod/QCD_Pt-120to170_Tune4C_13TeV_pythia8/MINIAODSIM/castor_PU20bx25_POSTLS170_V5-v1/00000/081D87ED-7906-E411-81C3-0025905A60C6.root',
       '/store/mc/Spring14miniaod/QCD_Pt-120to170_Tune4C_13TeV_pythia8/MINIAODSIM/castor_PU20bx25_POSTLS170_V5-v1/00000/086EC0A4-7906-E411-A998-0025905A6068.root',
       '/store/mc/Spring14miniaod/QCD_Pt-120to170_Tune4C_13TeV_pythia8/MINIAODSIM/castor_PU20bx25_POSTLS170_V5-v1/00000/0A0F360C-9006-E411-AE53-00261894394F.root',
       '/store/mc/Spring14miniaod/QCD_Pt-120to170_Tune4C_13TeV_pythia8/MINIAODSIM/castor_PU20bx25_POSTLS170_V5-v1/00000/0A72D207-7C06-E411-BFA2-0025905B8596.root',
       '/store/mc/Spring14miniaod/QCD_Pt-120to170_Tune4C_13TeV_pythia8/MINIAODSIM/castor_PU20bx25_POSTLS170_V5-v1/00000/0A8335CB-7B06-E411-9E17-002618FDA28E.root',
       '/store/mc/Spring14miniaod/QCD_Pt-120to170_Tune4C_13TeV_pythia8/MINIAODSIM/castor_PU20bx25_POSTLS170_V5-v1/00000/0C8796B6-7906-E411-BC86-003048FFCC0A.root',
       '/store/mc/Spring14miniaod/QCD_Pt-120to170_Tune4C_13TeV_pythia8/MINIAODSIM/castor_PU20bx25_POSTLS170_V5-v1/00000/0E1D2EFD-8F06-E411-A2F0-0025905A60A8.root',
       '/store/mc/Spring14miniaod/QCD_Pt-120to170_Tune4C_13TeV_pythia8/MINIAODSIM/castor_PU20bx25_POSTLS170_V5-v1/00000/0EFD5770-9006-E411-AAAF-0025905A6118.root',
       '/store/mc/Spring14miniaod/QCD_Pt-120to170_Tune4C_13TeV_pythia8/MINIAODSIM/castor_PU20bx25_POSTLS170_V5-v1/00000/1262BA09-9006-E411-82CF-002618943933.root',
       '/store/mc/Spring14miniaod/QCD_Pt-120to170_Tune4C_13TeV_pythia8/MINIAODSIM/castor_PU20bx25_POSTLS170_V5-v1/00000/18BE6822-7B06-E411-BD57-002618943882.root',
       '/store/mc/Spring14miniaod/QCD_Pt-120to170_Tune4C_13TeV_pythia8/MINIAODSIM/castor_PU20bx25_POSTLS170_V5-v1/00000/1A3D7D40-7C06-E411-B6D3-003048FFD7A2.root',
       '/store/mc/Spring14miniaod/QCD_Pt-120to170_Tune4C_13TeV_pythia8/MINIAODSIM/castor_PU20bx25_POSTLS170_V5-v1/00000/1AD0F044-7B06-E411-8448-0025905938A4.root',
       '/store/mc/Spring14miniaod/QCD_Pt-120to170_Tune4C_13TeV_pythia8/MINIAODSIM/castor_PU20bx25_POSTLS170_V5-v1/00000/20B9B2DD-7C06-E411-95C7-0025905B860E.root',
       '/store/mc/Spring14miniaod/QCD_Pt-120to170_Tune4C_13TeV_pythia8/MINIAODSIM/castor_PU20bx25_POSTLS170_V5-v1/00000/242F5A21-7C06-E411-8240-0025905A608C.root',
       '/store/mc/Spring14miniaod/QCD_Pt-120to170_Tune4C_13TeV_pythia8/MINIAODSIM/castor_PU20bx25_POSTLS170_V5-v1/00000/24482EB9-7906-E411-A051-0025905A60DE.root',
       '/store/mc/Spring14miniaod/QCD_Pt-120to170_Tune4C_13TeV_pythia8/MINIAODSIM/castor_PU20bx25_POSTLS170_V5-v1/00000/24D03C07-9006-E411-A2CA-0025905B85B2.root',
       '/store/mc/Spring14miniaod/QCD_Pt-120to170_Tune4C_13TeV_pythia8/MINIAODSIM/castor_PU20bx25_POSTLS170_V5-v1/00000/264BD6B6-7906-E411-9FE7-003048FFD760.root',
       '/store/mc/Spring14miniaod/QCD_Pt-120to170_Tune4C_13TeV_pythia8/MINIAODSIM/castor_PU20bx25_POSTLS170_V5-v1/00000/2816D4AF-7906-E411-95E9-002590596468.root',
       '/store/mc/Spring14miniaod/QCD_Pt-120to170_Tune4C_13TeV_pythia8/MINIAODSIM/castor_PU20bx25_POSTLS170_V5-v1/00000/304F2058-7B06-E411-9B59-0025905AA9F0.root',
       '/store/mc/Spring14miniaod/QCD_Pt-120to170_Tune4C_13TeV_pythia8/MINIAODSIM/castor_PU20bx25_POSTLS170_V5-v1/00000/3478E5B8-7D06-E411-B20A-0025905A60F8.root',
       '/store/mc/Spring14miniaod/QCD_Pt-120to170_Tune4C_13TeV_pythia8/MINIAODSIM/castor_PU20bx25_POSTLS170_V5-v1/00000/34D650BF-7B06-E411-8851-0025905A60B6.root',
       '/store/mc/Spring14miniaod/QCD_Pt-120to170_Tune4C_13TeV_pythia8/MINIAODSIM/castor_PU20bx25_POSTLS170_V5-v1/00000/367866BC-7906-E411-AC79-0025905A611C.root',
       '/store/mc/Spring14miniaod/QCD_Pt-120to170_Tune4C_13TeV_pythia8/MINIAODSIM/castor_PU20bx25_POSTLS170_V5-v1/00000/4079628A-9006-E411-8D6C-0025905A6070.root',
       '/store/mc/Spring14miniaod/QCD_Pt-120to170_Tune4C_13TeV_pythia8/MINIAODSIM/castor_PU20bx25_POSTLS170_V5-v1/00000/40DEA41A-7C06-E411-B19E-002590593872.root',
       '/store/mc/Spring14miniaod/QCD_Pt-120to170_Tune4C_13TeV_pythia8/MINIAODSIM/castor_PU20bx25_POSTLS170_V5-v1/00000/4200D5DC-7906-E411-97DE-0025905A60F2.root',
       '/store/mc/Spring14miniaod/QCD_Pt-120to170_Tune4C_13TeV_pythia8/MINIAODSIM/castor_PU20bx25_POSTLS170_V5-v1/00000/44E05EF3-7906-E411-9AFF-0025905A605E.root',
       '/store/mc/Spring14miniaod/QCD_Pt-120to170_Tune4C_13TeV_pythia8/MINIAODSIM/castor_PU20bx25_POSTLS170_V5-v1/00000/44EB1B95-7906-E411-87E5-00261894393A.root',
       '/store/mc/Spring14miniaod/QCD_Pt-120to170_Tune4C_13TeV_pythia8/MINIAODSIM/castor_PU20bx25_POSTLS170_V5-v1/00000/4A4688CC-7906-E411-987B-0025905822B6.root',
       '/store/mc/Spring14miniaod/QCD_Pt-120to170_Tune4C_13TeV_pythia8/MINIAODSIM/castor_PU20bx25_POSTLS170_V5-v1/00000/4C56A1D8-7906-E411-82DD-0025905B860C.root',
       '/store/mc/Spring14miniaod/QCD_Pt-120to170_Tune4C_13TeV_pythia8/MINIAODSIM/castor_PU20bx25_POSTLS170_V5-v1/00000/504C0DED-7C06-E411-9BC6-002354EF3BE4.root',
       '/store/mc/Spring14miniaod/QCD_Pt-120to170_Tune4C_13TeV_pythia8/MINIAODSIM/castor_PU20bx25_POSTLS170_V5-v1/00000/56A421D3-7906-E411-949B-0025905A60B2.root',
       '/store/mc/Spring14miniaod/QCD_Pt-120to170_Tune4C_13TeV_pythia8/MINIAODSIM/castor_PU20bx25_POSTLS170_V5-v1/00000/5861E8F9-7C06-E411-85F2-0025905B855E.root',
       '/store/mc/Spring14miniaod/QCD_Pt-120to170_Tune4C_13TeV_pythia8/MINIAODSIM/castor_PU20bx25_POSTLS170_V5-v1/00000/5A5FD8B7-7906-E411-BCD9-003048FFD740.root',
       '/store/mc/Spring14miniaod/QCD_Pt-120to170_Tune4C_13TeV_pythia8/MINIAODSIM/castor_PU20bx25_POSTLS170_V5-v1/00000/5E41F93B-7B06-E411-A2BF-0025905B860C.root',
       '/store/mc/Spring14miniaod/QCD_Pt-120to170_Tune4C_13TeV_pythia8/MINIAODSIM/castor_PU20bx25_POSTLS170_V5-v1/00000/6671C33E-7B06-E411-8A9F-003048FFD760.root',
       '/store/mc/Spring14miniaod/QCD_Pt-120to170_Tune4C_13TeV_pythia8/MINIAODSIM/castor_PU20bx25_POSTLS170_V5-v1/00000/6A4657E7-7906-E411-8549-0025905A609E.root',
       '/store/mc/Spring14miniaod/QCD_Pt-120to170_Tune4C_13TeV_pythia8/MINIAODSIM/castor_PU20bx25_POSTLS170_V5-v1/00000/70D1941B-9006-E411-B1D3-0025905A608A.root',
       '/store/mc/Spring14miniaod/QCD_Pt-120to170_Tune4C_13TeV_pythia8/MINIAODSIM/castor_PU20bx25_POSTLS170_V5-v1/00000/721277C6-7906-E411-AA64-002618943981.root',
       '/store/mc/Spring14miniaod/QCD_Pt-120to170_Tune4C_13TeV_pythia8/MINIAODSIM/castor_PU20bx25_POSTLS170_V5-v1/00000/76F101EF-7906-E411-A8BE-0025905A611E.root',
       '/store/mc/Spring14miniaod/QCD_Pt-120to170_Tune4C_13TeV_pythia8/MINIAODSIM/castor_PU20bx25_POSTLS170_V5-v1/00000/78A9E7D9-8F06-E411-98EA-002354EF3BD2.root',
       '/store/mc/Spring14miniaod/QCD_Pt-120to170_Tune4C_13TeV_pythia8/MINIAODSIM/castor_PU20bx25_POSTLS170_V5-v1/00000/841A47E3-7C06-E411-AB27-002618943809.root',
       '/store/mc/Spring14miniaod/QCD_Pt-120to170_Tune4C_13TeV_pythia8/MINIAODSIM/castor_PU20bx25_POSTLS170_V5-v1/00000/84D3C09D-7B06-E411-9CCB-002618943829.root',
       '/store/mc/Spring14miniaod/QCD_Pt-120to170_Tune4C_13TeV_pythia8/MINIAODSIM/castor_PU20bx25_POSTLS170_V5-v1/00000/8A993E27-7B06-E411-92FE-0025905A606A.root',
       '/store/mc/Spring14miniaod/QCD_Pt-120to170_Tune4C_13TeV_pythia8/MINIAODSIM/castor_PU20bx25_POSTLS170_V5-v1/00000/8C206DF6-8F06-E411-90DF-0025905A607E.root',
       '/store/mc/Spring14miniaod/QCD_Pt-120to170_Tune4C_13TeV_pythia8/MINIAODSIM/castor_PU20bx25_POSTLS170_V5-v1/00000/8EBC0812-9006-E411-9C90-003048FFD740.root',
       '/store/mc/Spring14miniaod/QCD_Pt-120to170_Tune4C_13TeV_pythia8/MINIAODSIM/castor_PU20bx25_POSTLS170_V5-v1/00000/9007BACF-7906-E411-AD7A-003048FFD770.root',
       '/store/mc/Spring14miniaod/QCD_Pt-120to170_Tune4C_13TeV_pythia8/MINIAODSIM/castor_PU20bx25_POSTLS170_V5-v1/00000/96B10F1D-7C06-E411-A56C-0025905A48F0.root',
       '/store/mc/Spring14miniaod/QCD_Pt-120to170_Tune4C_13TeV_pythia8/MINIAODSIM/castor_PU20bx25_POSTLS170_V5-v1/00000/9821707C-7C06-E411-8BBC-0025905A48F2.root',
       '/store/mc/Spring14miniaod/QCD_Pt-120to170_Tune4C_13TeV_pythia8/MINIAODSIM/castor_PU20bx25_POSTLS170_V5-v1/00000/9830C6EC-7906-E411-AFE5-0025905A6084.root',
       '/store/mc/Spring14miniaod/QCD_Pt-120to170_Tune4C_13TeV_pythia8/MINIAODSIM/castor_PU20bx25_POSTLS170_V5-v1/00000/98E07AC1-7C06-E411-B598-002618FDA28E.root',
       '/store/mc/Spring14miniaod/QCD_Pt-120to170_Tune4C_13TeV_pythia8/MINIAODSIM/castor_PU20bx25_POSTLS170_V5-v1/00000/9A9C67D3-7906-E411-BCEE-003048FFCB96.root',
       '/store/mc/Spring14miniaod/QCD_Pt-120to170_Tune4C_13TeV_pythia8/MINIAODSIM/castor_PU20bx25_POSTLS170_V5-v1/00000/9C0CAB47-7D06-E411-B602-003048FFD720.root',
       '/store/mc/Spring14miniaod/QCD_Pt-120to170_Tune4C_13TeV_pythia8/MINIAODSIM/castor_PU20bx25_POSTLS170_V5-v1/00000/A0424DC6-9006-E411-83D2-0025905AA9CC.root',
       '/store/mc/Spring14miniaod/QCD_Pt-120to170_Tune4C_13TeV_pythia8/MINIAODSIM/castor_PU20bx25_POSTLS170_V5-v1/00000/A2146E7C-7C06-E411-BB5C-0025905A48F2.root',
       '/store/mc/Spring14miniaod/QCD_Pt-120to170_Tune4C_13TeV_pythia8/MINIAODSIM/castor_PU20bx25_POSTLS170_V5-v1/00000/A25908B7-7B06-E411-864D-003048FFD730.root',
       '/store/mc/Spring14miniaod/QCD_Pt-120to170_Tune4C_13TeV_pythia8/MINIAODSIM/castor_PU20bx25_POSTLS170_V5-v1/00000/A2C7FA14-9006-E411-971D-003048FFCBA4.root',
       '/store/mc/Spring14miniaod/QCD_Pt-120to170_Tune4C_13TeV_pythia8/MINIAODSIM/castor_PU20bx25_POSTLS170_V5-v1/00000/A4919CA6-9006-E411-957D-0025905A608C.root',
       '/store/mc/Spring14miniaod/QCD_Pt-120to170_Tune4C_13TeV_pythia8/MINIAODSIM/castor_PU20bx25_POSTLS170_V5-v1/00000/A60446B3-7906-E411-9345-003048FFCBA8.root',
       '/store/mc/Spring14miniaod/QCD_Pt-120to170_Tune4C_13TeV_pythia8/MINIAODSIM/castor_PU20bx25_POSTLS170_V5-v1/00000/A87C67A9-7906-E411-A8B8-0025905B8596.root',
       '/store/mc/Spring14miniaod/QCD_Pt-120to170_Tune4C_13TeV_pythia8/MINIAODSIM/castor_PU20bx25_POSTLS170_V5-v1/00000/A895B81D-7D06-E411-98EE-002590596498.root',
       '/store/mc/Spring14miniaod/QCD_Pt-120to170_Tune4C_13TeV_pythia8/MINIAODSIM/castor_PU20bx25_POSTLS170_V5-v1/00000/A8D7D617-7D06-E411-AFEA-0025905B85E8.root',
       '/store/mc/Spring14miniaod/QCD_Pt-120to170_Tune4C_13TeV_pythia8/MINIAODSIM/castor_PU20bx25_POSTLS170_V5-v1/00000/AA45F1E5-7906-E411-975C-0025905A60A6.root',
       '/store/mc/Spring14miniaod/QCD_Pt-120to170_Tune4C_13TeV_pythia8/MINIAODSIM/castor_PU20bx25_POSTLS170_V5-v1/00000/AC1446CF-7906-E411-8249-002618943976.root',
       '/store/mc/Spring14miniaod/QCD_Pt-120to170_Tune4C_13TeV_pythia8/MINIAODSIM/castor_PU20bx25_POSTLS170_V5-v1/00000/ACE77EB9-7906-E411-B0A9-00259059649C.root',
       '/store/mc/Spring14miniaod/QCD_Pt-120to170_Tune4C_13TeV_pythia8/MINIAODSIM/castor_PU20bx25_POSTLS170_V5-v1/00000/B032A8DA-7B06-E411-BB7F-002354EF3BE4.root',
       '/store/mc/Spring14miniaod/QCD_Pt-120to170_Tune4C_13TeV_pythia8/MINIAODSIM/castor_PU20bx25_POSTLS170_V5-v1/00000/B08F26AC-7906-E411-B9CC-0025905938A4.root',
       '/store/mc/Spring14miniaod/QCD_Pt-120to170_Tune4C_13TeV_pythia8/MINIAODSIM/castor_PU20bx25_POSTLS170_V5-v1/00000/B2206AD3-7B06-E411-A418-0025905A48F2.root',
       '/store/mc/Spring14miniaod/QCD_Pt-120to170_Tune4C_13TeV_pythia8/MINIAODSIM/castor_PU20bx25_POSTLS170_V5-v1/00000/B47134D1-7906-E411-8FD6-0025905A60E4.root',
       '/store/mc/Spring14miniaod/QCD_Pt-120to170_Tune4C_13TeV_pythia8/MINIAODSIM/castor_PU20bx25_POSTLS170_V5-v1/00000/B4893F48-7B06-E411-8C46-0025905A6088.root',
       '/store/mc/Spring14miniaod/QCD_Pt-120to170_Tune4C_13TeV_pythia8/MINIAODSIM/castor_PU20bx25_POSTLS170_V5-v1/00000/B8BD6EDB-7B06-E411-BCDC-0025905A609E.root',
       '/store/mc/Spring14miniaod/QCD_Pt-120to170_Tune4C_13TeV_pythia8/MINIAODSIM/castor_PU20bx25_POSTLS170_V5-v1/00000/BCABA309-7C06-E411-AF56-002618943960.root',
       '/store/mc/Spring14miniaod/QCD_Pt-120to170_Tune4C_13TeV_pythia8/MINIAODSIM/castor_PU20bx25_POSTLS170_V5-v1/00000/C41766BC-7906-E411-A302-0025905A611C.root',
       '/store/mc/Spring14miniaod/QCD_Pt-120to170_Tune4C_13TeV_pythia8/MINIAODSIM/castor_PU20bx25_POSTLS170_V5-v1/00000/C47761E1-7B06-E411-A094-003048FFCC18.root',
       '/store/mc/Spring14miniaod/QCD_Pt-120to170_Tune4C_13TeV_pythia8/MINIAODSIM/castor_PU20bx25_POSTLS170_V5-v1/00000/C8FB7CF4-7906-E411-B5AC-0025905A608A.root',
       '/store/mc/Spring14miniaod/QCD_Pt-120to170_Tune4C_13TeV_pythia8/MINIAODSIM/castor_PU20bx25_POSTLS170_V5-v1/00000/CCA84C07-9106-E411-8C1B-0025905A48FC.root',
       '/store/mc/Spring14miniaod/QCD_Pt-120to170_Tune4C_13TeV_pythia8/MINIAODSIM/castor_PU20bx25_POSTLS170_V5-v1/00000/D090B229-9006-E411-9929-003048FFCBB0.root',
       '/store/mc/Spring14miniaod/QCD_Pt-120to170_Tune4C_13TeV_pythia8/MINIAODSIM/castor_PU20bx25_POSTLS170_V5-v1/00000/D0E3F8F9-8F06-E411-9E54-002618943885.root',
       '/store/mc/Spring14miniaod/QCD_Pt-120to170_Tune4C_13TeV_pythia8/MINIAODSIM/castor_PU20bx25_POSTLS170_V5-v1/00000/D20EB0FB-7C06-E411-A94B-0025905A60DE.root',
       '/store/mc/Spring14miniaod/QCD_Pt-120to170_Tune4C_13TeV_pythia8/MINIAODSIM/castor_PU20bx25_POSTLS170_V5-v1/00000/D2831DE7-7C06-E411-B601-0025905938A4.root',
       '/store/mc/Spring14miniaod/QCD_Pt-120to170_Tune4C_13TeV_pythia8/MINIAODSIM/castor_PU20bx25_POSTLS170_V5-v1/00000/D2CDD7D6-7906-E411-AA58-0025905A6076.root',
       '/store/mc/Spring14miniaod/QCD_Pt-120to170_Tune4C_13TeV_pythia8/MINIAODSIM/castor_PU20bx25_POSTLS170_V5-v1/00000/D65B55FC-7B06-E411-B5B1-0025905938A8.root',
       '/store/mc/Spring14miniaod/QCD_Pt-120to170_Tune4C_13TeV_pythia8/MINIAODSIM/castor_PU20bx25_POSTLS170_V5-v1/00000/DC216B74-9006-E411-B418-003048FFCC2C.root',
       '/store/mc/Spring14miniaod/QCD_Pt-120to170_Tune4C_13TeV_pythia8/MINIAODSIM/castor_PU20bx25_POSTLS170_V5-v1/00000/DEDCB008-7C06-E411-8B69-002618943876.root',
       '/store/mc/Spring14miniaod/QCD_Pt-120to170_Tune4C_13TeV_pythia8/MINIAODSIM/castor_PU20bx25_POSTLS170_V5-v1/00000/E02C851A-7B06-E411-8C3D-003048FFD756.root',
       '/store/mc/Spring14miniaod/QCD_Pt-120to170_Tune4C_13TeV_pythia8/MINIAODSIM/castor_PU20bx25_POSTLS170_V5-v1/00000/E06FD1E6-7906-E411-BEFD-0025905A60F8.root',
       '/store/mc/Spring14miniaod/QCD_Pt-120to170_Tune4C_13TeV_pythia8/MINIAODSIM/castor_PU20bx25_POSTLS170_V5-v1/00000/E236D945-7C06-E411-B833-0025905A60DA.root',
       '/store/mc/Spring14miniaod/QCD_Pt-120to170_Tune4C_13TeV_pythia8/MINIAODSIM/castor_PU20bx25_POSTLS170_V5-v1/00000/E272E8C5-7906-E411-9253-00261894382D.root',
       '/store/mc/Spring14miniaod/QCD_Pt-120to170_Tune4C_13TeV_pythia8/MINIAODSIM/castor_PU20bx25_POSTLS170_V5-v1/00000/E2DBD611-9006-E411-8232-0025905A611E.root',
       '/store/mc/Spring14miniaod/QCD_Pt-120to170_Tune4C_13TeV_pythia8/MINIAODSIM/castor_PU20bx25_POSTLS170_V5-v1/00000/E476FFEB-7C06-E411-9190-0025905A6090.root',
       '/store/mc/Spring14miniaod/QCD_Pt-120to170_Tune4C_13TeV_pythia8/MINIAODSIM/castor_PU20bx25_POSTLS170_V5-v1/00000/E4B5BFC7-7906-E411-AE9D-002618943919.root',
       '/store/mc/Spring14miniaod/QCD_Pt-120to170_Tune4C_13TeV_pythia8/MINIAODSIM/castor_PU20bx25_POSTLS170_V5-v1/00000/E6D56625-7B06-E411-8ED1-003048D15DE0.root',
       '/store/mc/Spring14miniaod/QCD_Pt-120to170_Tune4C_13TeV_pythia8/MINIAODSIM/castor_PU20bx25_POSTLS170_V5-v1/00000/E8B5BFAF-7906-E411-82D2-0025905938B4.root',
       '/store/mc/Spring14miniaod/QCD_Pt-120to170_Tune4C_13TeV_pythia8/MINIAODSIM/castor_PU20bx25_POSTLS170_V5-v1/00000/EA9277FD-7B06-E411-B3F7-002590596484.root',
       '/store/mc/Spring14miniaod/QCD_Pt-120to170_Tune4C_13TeV_pythia8/MINIAODSIM/castor_PU20bx25_POSTLS170_V5-v1/00000/EA9CFD64-7A06-E411-B934-003048FFCB6A.root',
       '/store/mc/Spring14miniaod/QCD_Pt-120to170_Tune4C_13TeV_pythia8/MINIAODSIM/castor_PU20bx25_POSTLS170_V5-v1/00000/EC5F3888-7A06-E411-A5D7-0025905AA9F0.root',
       '/store/mc/Spring14miniaod/QCD_Pt-120to170_Tune4C_13TeV_pythia8/MINIAODSIM/castor_PU20bx25_POSTLS170_V5-v1/00000/F02984F9-8F06-E411-9054-0026189438DF.root',
       '/store/mc/Spring14miniaod/QCD_Pt-120to170_Tune4C_13TeV_pythia8/MINIAODSIM/castor_PU20bx25_POSTLS170_V5-v1/00000/F039B684-7C06-E411-AC5A-003048FFD740.root',
       '/store/mc/Spring14miniaod/QCD_Pt-120to170_Tune4C_13TeV_pythia8/MINIAODSIM/castor_PU20bx25_POSTLS170_V5-v1/00000/F0B92F0E-9006-E411-9447-0025905A60D2.root',
       '/store/mc/Spring14miniaod/QCD_Pt-120to170_Tune4C_13TeV_pythia8/MINIAODSIM/castor_PU20bx25_POSTLS170_V5-v1/00000/F214CD13-7C06-E411-9B70-0025905B85B2.root',
       '/store/mc/Spring14miniaod/QCD_Pt-120to170_Tune4C_13TeV_pythia8/MINIAODSIM/castor_PU20bx25_POSTLS170_V5-v1/00000/F461FFAD-7906-E411-ACDC-003048D15DE0.root',
       '/store/mc/Spring14miniaod/QCD_Pt-120to170_Tune4C_13TeV_pythia8/MINIAODSIM/castor_PU20bx25_POSTLS170_V5-v1/00000/F8A84E07-7C06-E411-ABD5-0025905A60B2.root',
       '/store/mc/Spring14miniaod/QCD_Pt-120to170_Tune4C_13TeV_pythia8/MINIAODSIM/castor_PU20bx25_POSTLS170_V5-v1/00000/FA644A44-7B06-E411-908A-0025905A612A.root' ] );


secFiles.extend( [
               ] )


