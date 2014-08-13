import FWCore.ParameterSet.Config as cms

maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )
readFiles = cms.untracked.vstring()
secFiles = cms.untracked.vstring() 
source = cms.Source ("PoolSource",fileNames = readFiles, secondaryFileNames = secFiles)
readFiles.extend( [
       '/store/mc/Spring14miniaod/QCD_Pt-1800_Tune4C_13TeV_pythia8/MINIAODSIM/castor_PU20bx25_POSTLS170_V5-v1/00000/0239FB21-7A06-E411-9499-0025905A60FE.root',
       '/store/mc/Spring14miniaod/QCD_Pt-1800_Tune4C_13TeV_pythia8/MINIAODSIM/castor_PU20bx25_POSTLS170_V5-v1/00000/0C85BCA6-7906-E411-82EC-002618943972.root',
       '/store/mc/Spring14miniaod/QCD_Pt-1800_Tune4C_13TeV_pythia8/MINIAODSIM/castor_PU20bx25_POSTLS170_V5-v1/00000/0E8AE7E1-7A06-E411-9553-0025905AA9CC.root',
       '/store/mc/Spring14miniaod/QCD_Pt-1800_Tune4C_13TeV_pythia8/MINIAODSIM/castor_PU20bx25_POSTLS170_V5-v1/00000/105A18A4-7906-E411-AB3A-00261894394F.root',
       '/store/mc/Spring14miniaod/QCD_Pt-1800_Tune4C_13TeV_pythia8/MINIAODSIM/castor_PU20bx25_POSTLS170_V5-v1/00000/12933B13-7A06-E411-BDCA-0025905A609E.root',
       '/store/mc/Spring14miniaod/QCD_Pt-1800_Tune4C_13TeV_pythia8/MINIAODSIM/castor_PU20bx25_POSTLS170_V5-v1/00000/12A6856E-7B06-E411-879D-00261894389A.root',
       '/store/mc/Spring14miniaod/QCD_Pt-1800_Tune4C_13TeV_pythia8/MINIAODSIM/castor_PU20bx25_POSTLS170_V5-v1/00000/1459AD10-7A06-E411-AC83-0025905A609E.root',
       '/store/mc/Spring14miniaod/QCD_Pt-1800_Tune4C_13TeV_pythia8/MINIAODSIM/castor_PU20bx25_POSTLS170_V5-v1/00000/1A97105E-7B06-E411-9595-0026189437EB.root',
       '/store/mc/Spring14miniaod/QCD_Pt-1800_Tune4C_13TeV_pythia8/MINIAODSIM/castor_PU20bx25_POSTLS170_V5-v1/00000/1A9C6FDD-7906-E411-BE86-00248C55CC9D.root',
       '/store/mc/Spring14miniaod/QCD_Pt-1800_Tune4C_13TeV_pythia8/MINIAODSIM/castor_PU20bx25_POSTLS170_V5-v1/00000/1C39F40D-7B06-E411-9161-0025905A608E.root',
       '/store/mc/Spring14miniaod/QCD_Pt-1800_Tune4C_13TeV_pythia8/MINIAODSIM/castor_PU20bx25_POSTLS170_V5-v1/00000/1C7F8DAB-7B06-E411-9A90-0025905A608C.root',
       '/store/mc/Spring14miniaod/QCD_Pt-1800_Tune4C_13TeV_pythia8/MINIAODSIM/castor_PU20bx25_POSTLS170_V5-v1/00000/1E56903A-7A06-E411-BA28-00259059642E.root',
       '/store/mc/Spring14miniaod/QCD_Pt-1800_Tune4C_13TeV_pythia8/MINIAODSIM/castor_PU20bx25_POSTLS170_V5-v1/00000/24DE131A-7A06-E411-AD45-0025905B855C.root',
       '/store/mc/Spring14miniaod/QCD_Pt-1800_Tune4C_13TeV_pythia8/MINIAODSIM/castor_PU20bx25_POSTLS170_V5-v1/00000/264B5603-7B06-E411-9A62-002618943972.root',
       '/store/mc/Spring14miniaod/QCD_Pt-1800_Tune4C_13TeV_pythia8/MINIAODSIM/castor_PU20bx25_POSTLS170_V5-v1/00000/2E3AF2F8-7906-E411-ABB1-0026189437EB.root',
       '/store/mc/Spring14miniaod/QCD_Pt-1800_Tune4C_13TeV_pythia8/MINIAODSIM/castor_PU20bx25_POSTLS170_V5-v1/00000/30E00F1E-7B06-E411-907E-0025905938D4.root',
       '/store/mc/Spring14miniaod/QCD_Pt-1800_Tune4C_13TeV_pythia8/MINIAODSIM/castor_PU20bx25_POSTLS170_V5-v1/00000/32493016-7A06-E411-8CC4-0025905A60E4.root',
       '/store/mc/Spring14miniaod/QCD_Pt-1800_Tune4C_13TeV_pythia8/MINIAODSIM/castor_PU20bx25_POSTLS170_V5-v1/00000/36A2FCBD-7B06-E411-A6E2-0025905A6088.root',
       '/store/mc/Spring14miniaod/QCD_Pt-1800_Tune4C_13TeV_pythia8/MINIAODSIM/castor_PU20bx25_POSTLS170_V5-v1/00000/402AEF25-7B06-E411-B376-0025905A6070.root',
       '/store/mc/Spring14miniaod/QCD_Pt-1800_Tune4C_13TeV_pythia8/MINIAODSIM/castor_PU20bx25_POSTLS170_V5-v1/00000/42C9778E-7A06-E411-A49A-0025905A60B0.root',
       '/store/mc/Spring14miniaod/QCD_Pt-1800_Tune4C_13TeV_pythia8/MINIAODSIM/castor_PU20bx25_POSTLS170_V5-v1/00000/48141D0B-7B06-E411-88BA-00261894394F.root',
       '/store/mc/Spring14miniaod/QCD_Pt-1800_Tune4C_13TeV_pythia8/MINIAODSIM/castor_PU20bx25_POSTLS170_V5-v1/00000/4AB96332-7A06-E411-AC7B-00259059642A.root',
       '/store/mc/Spring14miniaod/QCD_Pt-1800_Tune4C_13TeV_pythia8/MINIAODSIM/castor_PU20bx25_POSTLS170_V5-v1/00000/5A9D73D2-8006-E411-8882-0025905A60AA.root',
       '/store/mc/Spring14miniaod/QCD_Pt-1800_Tune4C_13TeV_pythia8/MINIAODSIM/castor_PU20bx25_POSTLS170_V5-v1/00000/64326A8F-7906-E411-B5D0-00259059391E.root',
       '/store/mc/Spring14miniaod/QCD_Pt-1800_Tune4C_13TeV_pythia8/MINIAODSIM/castor_PU20bx25_POSTLS170_V5-v1/00000/66D43249-7B06-E411-9301-0026189438B1.root',
       '/store/mc/Spring14miniaod/QCD_Pt-1800_Tune4C_13TeV_pythia8/MINIAODSIM/castor_PU20bx25_POSTLS170_V5-v1/00000/6C0BAA67-7B06-E411-A9E0-0025905B85B2.root',
       '/store/mc/Spring14miniaod/QCD_Pt-1800_Tune4C_13TeV_pythia8/MINIAODSIM/castor_PU20bx25_POSTLS170_V5-v1/00000/707105EF-7906-E411-B043-0025905A6118.root',
       '/store/mc/Spring14miniaod/QCD_Pt-1800_Tune4C_13TeV_pythia8/MINIAODSIM/castor_PU20bx25_POSTLS170_V5-v1/00000/76E70D1A-7A06-E411-B5B1-0025905A60B8.root',
       '/store/mc/Spring14miniaod/QCD_Pt-1800_Tune4C_13TeV_pythia8/MINIAODSIM/castor_PU20bx25_POSTLS170_V5-v1/00000/7CBBB870-7B06-E411-87A2-0026189438CB.root',
       '/store/mc/Spring14miniaod/QCD_Pt-1800_Tune4C_13TeV_pythia8/MINIAODSIM/castor_PU20bx25_POSTLS170_V5-v1/00000/8229F4CF-7B06-E411-BB21-0025905A60FE.root',
       '/store/mc/Spring14miniaod/QCD_Pt-1800_Tune4C_13TeV_pythia8/MINIAODSIM/castor_PU20bx25_POSTLS170_V5-v1/00000/8A791CC5-7906-E411-9DB3-0026189438B1.root',
       '/store/mc/Spring14miniaod/QCD_Pt-1800_Tune4C_13TeV_pythia8/MINIAODSIM/castor_PU20bx25_POSTLS170_V5-v1/00000/90F4DF2F-7A06-E411-8D29-0026189438CB.root',
       '/store/mc/Spring14miniaod/QCD_Pt-1800_Tune4C_13TeV_pythia8/MINIAODSIM/castor_PU20bx25_POSTLS170_V5-v1/00000/9297DD1E-7A06-E411-AC2B-00261894389A.root',
       '/store/mc/Spring14miniaod/QCD_Pt-1800_Tune4C_13TeV_pythia8/MINIAODSIM/castor_PU20bx25_POSTLS170_V5-v1/00000/9404E52F-7A06-E411-BBFA-0026189438CB.root',
       '/store/mc/Spring14miniaod/QCD_Pt-1800_Tune4C_13TeV_pythia8/MINIAODSIM/castor_PU20bx25_POSTLS170_V5-v1/00000/9467917E-7B06-E411-AFE6-00259059391E.root',
       '/store/mc/Spring14miniaod/QCD_Pt-1800_Tune4C_13TeV_pythia8/MINIAODSIM/castor_PU20bx25_POSTLS170_V5-v1/00000/9A4E9939-7B06-E411-84A3-0026189438BF.root',
       '/store/mc/Spring14miniaod/QCD_Pt-1800_Tune4C_13TeV_pythia8/MINIAODSIM/castor_PU20bx25_POSTLS170_V5-v1/00000/9CC8AB10-7A06-E411-B174-0025905A6076.root',
       '/store/mc/Spring14miniaod/QCD_Pt-1800_Tune4C_13TeV_pythia8/MINIAODSIM/castor_PU20bx25_POSTLS170_V5-v1/00000/9EEBB329-7A06-E411-B6A6-0025905964C4.root',
       '/store/mc/Spring14miniaod/QCD_Pt-1800_Tune4C_13TeV_pythia8/MINIAODSIM/castor_PU20bx25_POSTLS170_V5-v1/00000/A092AF24-7B06-E411-93DB-00248C55CC9D.root',
       '/store/mc/Spring14miniaod/QCD_Pt-1800_Tune4C_13TeV_pythia8/MINIAODSIM/castor_PU20bx25_POSTLS170_V5-v1/00000/A239DC58-7B06-E411-8588-002618943920.root',
       '/store/mc/Spring14miniaod/QCD_Pt-1800_Tune4C_13TeV_pythia8/MINIAODSIM/castor_PU20bx25_POSTLS170_V5-v1/00000/A29E7F96-7B06-E411-9E21-0025905A610A.root',
       '/store/mc/Spring14miniaod/QCD_Pt-1800_Tune4C_13TeV_pythia8/MINIAODSIM/castor_PU20bx25_POSTLS170_V5-v1/00000/AA8DCD1B-7A06-E411-9BFA-0026189438BF.root',
       '/store/mc/Spring14miniaod/QCD_Pt-1800_Tune4C_13TeV_pythia8/MINIAODSIM/castor_PU20bx25_POSTLS170_V5-v1/00000/AABB6B5C-7C06-E411-8E6B-0025905A611C.root',
       '/store/mc/Spring14miniaod/QCD_Pt-1800_Tune4C_13TeV_pythia8/MINIAODSIM/castor_PU20bx25_POSTLS170_V5-v1/00000/B6B4931D-7A06-E411-A126-002590596486.root',
       '/store/mc/Spring14miniaod/QCD_Pt-1800_Tune4C_13TeV_pythia8/MINIAODSIM/castor_PU20bx25_POSTLS170_V5-v1/00000/BE1A6E1B-7A06-E411-A023-003048FFD756.root',
       '/store/mc/Spring14miniaod/QCD_Pt-1800_Tune4C_13TeV_pythia8/MINIAODSIM/castor_PU20bx25_POSTLS170_V5-v1/00000/BEC8A539-7B06-E411-BEA3-0025905A607E.root',
       '/store/mc/Spring14miniaod/QCD_Pt-1800_Tune4C_13TeV_pythia8/MINIAODSIM/castor_PU20bx25_POSTLS170_V5-v1/00000/C0297B5B-7B06-E411-877B-002590593876.root',
       '/store/mc/Spring14miniaod/QCD_Pt-1800_Tune4C_13TeV_pythia8/MINIAODSIM/castor_PU20bx25_POSTLS170_V5-v1/00000/CE5DD24A-7A06-E411-8F34-0025905A6126.root',
       '/store/mc/Spring14miniaod/QCD_Pt-1800_Tune4C_13TeV_pythia8/MINIAODSIM/castor_PU20bx25_POSTLS170_V5-v1/00000/D052CDCA-7B06-E411-84E9-003048FFD75C.root',
       '/store/mc/Spring14miniaod/QCD_Pt-1800_Tune4C_13TeV_pythia8/MINIAODSIM/castor_PU20bx25_POSTLS170_V5-v1/00000/D09007B7-7B06-E411-8868-003048FFD730.root',
       '/store/mc/Spring14miniaod/QCD_Pt-1800_Tune4C_13TeV_pythia8/MINIAODSIM/castor_PU20bx25_POSTLS170_V5-v1/00000/D6B57F96-7B06-E411-B1EF-0025905A610A.root',
       '/store/mc/Spring14miniaod/QCD_Pt-1800_Tune4C_13TeV_pythia8/MINIAODSIM/castor_PU20bx25_POSTLS170_V5-v1/00000/D8A6665D-7B06-E411-A966-0025905A48F2.root',
       '/store/mc/Spring14miniaod/QCD_Pt-1800_Tune4C_13TeV_pythia8/MINIAODSIM/castor_PU20bx25_POSTLS170_V5-v1/00000/DA6B537C-7B06-E411-B6B2-0026189438CB.root',
       '/store/mc/Spring14miniaod/QCD_Pt-1800_Tune4C_13TeV_pythia8/MINIAODSIM/castor_PU20bx25_POSTLS170_V5-v1/00000/DE161B11-7A06-E411-8AC3-002618943920.root',
       '/store/mc/Spring14miniaod/QCD_Pt-1800_Tune4C_13TeV_pythia8/MINIAODSIM/castor_PU20bx25_POSTLS170_V5-v1/00000/E021DCB6-7906-E411-8E67-0025905A6122.root',
       '/store/mc/Spring14miniaod/QCD_Pt-1800_Tune4C_13TeV_pythia8/MINIAODSIM/castor_PU20bx25_POSTLS170_V5-v1/00000/E4861CC0-7B06-E411-A44E-0026189438CB.root',
       '/store/mc/Spring14miniaod/QCD_Pt-1800_Tune4C_13TeV_pythia8/MINIAODSIM/castor_PU20bx25_POSTLS170_V5-v1/00000/EC80ABF8-7A06-E411-AE7A-0026189438B5.root',
       '/store/mc/Spring14miniaod/QCD_Pt-1800_Tune4C_13TeV_pythia8/MINIAODSIM/castor_PU20bx25_POSTLS170_V5-v1/00000/F03AAEE2-7A06-E411-8D7B-0025905B85D0.root',
       '/store/mc/Spring14miniaod/QCD_Pt-1800_Tune4C_13TeV_pythia8/MINIAODSIM/castor_PU20bx25_POSTLS170_V5-v1/00000/F21557FC-7B06-E411-BA5E-0025905938A8.root',
       '/store/mc/Spring14miniaod/QCD_Pt-1800_Tune4C_13TeV_pythia8/MINIAODSIM/castor_PU20bx25_POSTLS170_V5-v1/00000/F68CE00D-7A06-E411-B734-00261894393A.root',
       '/store/mc/Spring14miniaod/QCD_Pt-1800_Tune4C_13TeV_pythia8/MINIAODSIM/castor_PU20bx25_POSTLS170_V5-v1/00000/F800658C-7906-E411-939B-003048FFCC1E.root',
       '/store/mc/Spring14miniaod/QCD_Pt-1800_Tune4C_13TeV_pythia8/MINIAODSIM/castor_PU20bx25_POSTLS170_V5-v1/00000/F86C959D-7B06-E411-AC51-0026189438D5.root',
       '/store/mc/Spring14miniaod/QCD_Pt-1800_Tune4C_13TeV_pythia8/MINIAODSIM/castor_PU20bx25_POSTLS170_V5-v1/00000/F8FF9AA5-7B06-E411-B62B-0025905938A4.root',
       '/store/mc/Spring14miniaod/QCD_Pt-1800_Tune4C_13TeV_pythia8/MINIAODSIM/castor_PU20bx25_POSTLS170_V5-v1/00000/FC9E516A-7B06-E411-8D39-00261894393A.root' ] );


secFiles.extend( [
               ] )

