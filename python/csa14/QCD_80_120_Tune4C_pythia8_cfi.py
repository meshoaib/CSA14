import FWCore.ParameterSet.Config as cms

maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )
readFiles = cms.untracked.vstring()
secFiles = cms.untracked.vstring() 
source = cms.Source ("PoolSource",fileNames = readFiles, secondaryFileNames = secFiles)
readFiles.extend( [
       '/store/mc/Spring14miniaod/QCD_Pt-80to120_Tune4C_13TeV_pythia8/MINIAODSIM/castor_PU20bx25_POSTLS170_V5-v1/00000/025BE6D8-7E06-E411-9580-0025905A48F2.root',
       '/store/mc/Spring14miniaod/QCD_Pt-80to120_Tune4C_13TeV_pythia8/MINIAODSIM/castor_PU20bx25_POSTLS170_V5-v1/00000/0432D51B-7E06-E411-BDE4-0025905A48F2.root',
       '/store/mc/Spring14miniaod/QCD_Pt-80to120_Tune4C_13TeV_pythia8/MINIAODSIM/castor_PU20bx25_POSTLS170_V5-v1/00000/04B29E96-7D06-E411-A5FF-0025905A60D2.root',
       '/store/mc/Spring14miniaod/QCD_Pt-80to120_Tune4C_13TeV_pythia8/MINIAODSIM/castor_PU20bx25_POSTLS170_V5-v1/00000/069A4D2A-7906-E411-BB15-002618943875.root',
       '/store/mc/Spring14miniaod/QCD_Pt-80to120_Tune4C_13TeV_pythia8/MINIAODSIM/castor_PU20bx25_POSTLS170_V5-v1/00000/06BCAF73-7906-E411-9C27-002618FDA28E.root',
       '/store/mc/Spring14miniaod/QCD_Pt-80to120_Tune4C_13TeV_pythia8/MINIAODSIM/castor_PU20bx25_POSTLS170_V5-v1/00000/0C577586-7D06-E411-9358-0025905A48F0.root',
       '/store/mc/Spring14miniaod/QCD_Pt-80to120_Tune4C_13TeV_pythia8/MINIAODSIM/castor_PU20bx25_POSTLS170_V5-v1/00000/0E8147FF-7E06-E411-B4D6-0025905A60D0.root',
       '/store/mc/Spring14miniaod/QCD_Pt-80to120_Tune4C_13TeV_pythia8/MINIAODSIM/castor_PU20bx25_POSTLS170_V5-v1/00000/102D7C4F-7E06-E411-A1E3-0025905B85AE.root',
       '/store/mc/Spring14miniaod/QCD_Pt-80to120_Tune4C_13TeV_pythia8/MINIAODSIM/castor_PU20bx25_POSTLS170_V5-v1/00000/1217217B-8006-E411-A0D2-0025905A60D2.root',
       '/store/mc/Spring14miniaod/QCD_Pt-80to120_Tune4C_13TeV_pythia8/MINIAODSIM/castor_PU20bx25_POSTLS170_V5-v1/00000/142D2895-7E06-E411-B525-0025905B855E.root',
       '/store/mc/Spring14miniaod/QCD_Pt-80to120_Tune4C_13TeV_pythia8/MINIAODSIM/castor_PU20bx25_POSTLS170_V5-v1/00000/161357E8-7D06-E411-A9EE-0025905A60B2.root',
       '/store/mc/Spring14miniaod/QCD_Pt-80to120_Tune4C_13TeV_pythia8/MINIAODSIM/castor_PU20bx25_POSTLS170_V5-v1/00000/189A858A-7906-E411-B3F1-002618943894.root',
       '/store/mc/Spring14miniaod/QCD_Pt-80to120_Tune4C_13TeV_pythia8/MINIAODSIM/castor_PU20bx25_POSTLS170_V5-v1/00000/1A3C5256-7E06-E411-9E0A-0025905964B6.root',
       '/store/mc/Spring14miniaod/QCD_Pt-80to120_Tune4C_13TeV_pythia8/MINIAODSIM/castor_PU20bx25_POSTLS170_V5-v1/00000/1EECF19E-7D06-E411-AF40-0025905A60CE.root',
       '/store/mc/Spring14miniaod/QCD_Pt-80to120_Tune4C_13TeV_pythia8/MINIAODSIM/castor_PU20bx25_POSTLS170_V5-v1/00000/1EFF4254-7D06-E411-A159-003048D15DF0.root',
       '/store/mc/Spring14miniaod/QCD_Pt-80to120_Tune4C_13TeV_pythia8/MINIAODSIM/castor_PU20bx25_POSTLS170_V5-v1/00000/20499D1B-7906-E411-8EA1-003048FFD754.root',
       '/store/mc/Spring14miniaod/QCD_Pt-80to120_Tune4C_13TeV_pythia8/MINIAODSIM/castor_PU20bx25_POSTLS170_V5-v1/00000/2619D11F-7A06-E411-AD66-0025905A60F4.root',
       '/store/mc/Spring14miniaod/QCD_Pt-80to120_Tune4C_13TeV_pythia8/MINIAODSIM/castor_PU20bx25_POSTLS170_V5-v1/00000/286BA773-7906-E411-82AD-002618943829.root',
       '/store/mc/Spring14miniaod/QCD_Pt-80to120_Tune4C_13TeV_pythia8/MINIAODSIM/castor_PU20bx25_POSTLS170_V5-v1/00000/288D6881-7E06-E411-B57E-002590593878.root',
       '/store/mc/Spring14miniaod/QCD_Pt-80to120_Tune4C_13TeV_pythia8/MINIAODSIM/castor_PU20bx25_POSTLS170_V5-v1/00000/2CE9CFBA-7D06-E411-8986-0025905A613C.root',
       '/store/mc/Spring14miniaod/QCD_Pt-80to120_Tune4C_13TeV_pythia8/MINIAODSIM/castor_PU20bx25_POSTLS170_V5-v1/00000/2E7BE851-8006-E411-874E-0025905A612C.root',
       '/store/mc/Spring14miniaod/QCD_Pt-80to120_Tune4C_13TeV_pythia8/MINIAODSIM/castor_PU20bx25_POSTLS170_V5-v1/00000/32743974-7E06-E411-A2A8-0025905A48BC.root',
       '/store/mc/Spring14miniaod/QCD_Pt-80to120_Tune4C_13TeV_pythia8/MINIAODSIM/castor_PU20bx25_POSTLS170_V5-v1/00000/361F0567-7F06-E411-93A8-0025905A6136.root',
       '/store/mc/Spring14miniaod/QCD_Pt-80to120_Tune4C_13TeV_pythia8/MINIAODSIM/castor_PU20bx25_POSTLS170_V5-v1/00000/36D5BBF1-7D06-E411-A62A-0025905A606A.root',
       '/store/mc/Spring14miniaod/QCD_Pt-80to120_Tune4C_13TeV_pythia8/MINIAODSIM/castor_PU20bx25_POSTLS170_V5-v1/00000/38930B5F-7906-E411-964B-00261894390E.root',
       '/store/mc/Spring14miniaod/QCD_Pt-80to120_Tune4C_13TeV_pythia8/MINIAODSIM/castor_PU20bx25_POSTLS170_V5-v1/00000/3A0C124D-7D06-E411-AC2D-003048FFCC0A.root',
       '/store/mc/Spring14miniaod/QCD_Pt-80to120_Tune4C_13TeV_pythia8/MINIAODSIM/castor_PU20bx25_POSTLS170_V5-v1/00000/3A72F5A7-7E06-E411-AD2B-0025905A6132.root',
       '/store/mc/Spring14miniaod/QCD_Pt-80to120_Tune4C_13TeV_pythia8/MINIAODSIM/castor_PU20bx25_POSTLS170_V5-v1/00000/3C632E3C-7D06-E411-879F-0025905A613C.root',
       '/store/mc/Spring14miniaod/QCD_Pt-80to120_Tune4C_13TeV_pythia8/MINIAODSIM/castor_PU20bx25_POSTLS170_V5-v1/00000/3C67369E-7806-E411-BEE1-0025905A608A.root',
       '/store/mc/Spring14miniaod/QCD_Pt-80to120_Tune4C_13TeV_pythia8/MINIAODSIM/castor_PU20bx25_POSTLS170_V5-v1/00000/3C9C388E-7906-E411-A56D-00261894398A.root',
       '/store/mc/Spring14miniaod/QCD_Pt-80to120_Tune4C_13TeV_pythia8/MINIAODSIM/castor_PU20bx25_POSTLS170_V5-v1/00000/3CA53DE5-7806-E411-AE6F-002618943862.root',
       '/store/mc/Spring14miniaod/QCD_Pt-80to120_Tune4C_13TeV_pythia8/MINIAODSIM/castor_PU20bx25_POSTLS170_V5-v1/00000/40A723E5-7806-E411-8CED-00248C0BE005.root',
       '/store/mc/Spring14miniaod/QCD_Pt-80to120_Tune4C_13TeV_pythia8/MINIAODSIM/castor_PU20bx25_POSTLS170_V5-v1/00000/40AC4A22-7E06-E411-98B1-0025905B85B2.root',
       '/store/mc/Spring14miniaod/QCD_Pt-80to120_Tune4C_13TeV_pythia8/MINIAODSIM/castor_PU20bx25_POSTLS170_V5-v1/00000/44107568-7806-E411-9588-00261894390E.root',
       '/store/mc/Spring14miniaod/QCD_Pt-80to120_Tune4C_13TeV_pythia8/MINIAODSIM/castor_PU20bx25_POSTLS170_V5-v1/00000/46399B8F-7D06-E411-A7A6-00259059642E.root',
       '/store/mc/Spring14miniaod/QCD_Pt-80to120_Tune4C_13TeV_pythia8/MINIAODSIM/castor_PU20bx25_POSTLS170_V5-v1/00000/46769CA6-7806-E411-98A6-00259059649C.root',
       '/store/mc/Spring14miniaod/QCD_Pt-80to120_Tune4C_13TeV_pythia8/MINIAODSIM/castor_PU20bx25_POSTLS170_V5-v1/00000/4837139F-7D06-E411-A5E1-0025905B85E8.root',
       '/store/mc/Spring14miniaod/QCD_Pt-80to120_Tune4C_13TeV_pythia8/MINIAODSIM/castor_PU20bx25_POSTLS170_V5-v1/00000/48CD2087-7E06-E411-94CB-003048FFD720.root',
       '/store/mc/Spring14miniaod/QCD_Pt-80to120_Tune4C_13TeV_pythia8/MINIAODSIM/castor_PU20bx25_POSTLS170_V5-v1/00000/4CF27586-7906-E411-8F40-0026189438E0.root',
       '/store/mc/Spring14miniaod/QCD_Pt-80to120_Tune4C_13TeV_pythia8/MINIAODSIM/castor_PU20bx25_POSTLS170_V5-v1/00000/5029BD40-7E06-E411-93EF-0025905A6132.root',
       '/store/mc/Spring14miniaod/QCD_Pt-80to120_Tune4C_13TeV_pythia8/MINIAODSIM/castor_PU20bx25_POSTLS170_V5-v1/00000/507B67C8-7C06-E411-853A-0025905A60D0.root',
       '/store/mc/Spring14miniaod/QCD_Pt-80to120_Tune4C_13TeV_pythia8/MINIAODSIM/castor_PU20bx25_POSTLS170_V5-v1/00000/54F6AD15-7E06-E411-8B7E-002590596468.root',
       '/store/mc/Spring14miniaod/QCD_Pt-80to120_Tune4C_13TeV_pythia8/MINIAODSIM/castor_PU20bx25_POSTLS170_V5-v1/00000/588A4A8B-7E06-E411-8EC6-003048FFD756.root',
       '/store/mc/Spring14miniaod/QCD_Pt-80to120_Tune4C_13TeV_pythia8/MINIAODSIM/castor_PU20bx25_POSTLS170_V5-v1/00000/58F3854B-7906-E411-970B-0026189438E1.root',
       '/store/mc/Spring14miniaod/QCD_Pt-80to120_Tune4C_13TeV_pythia8/MINIAODSIM/castor_PU20bx25_POSTLS170_V5-v1/00000/5A1BFEC2-7D06-E411-8168-003048FFCB74.root',
       '/store/mc/Spring14miniaod/QCD_Pt-80to120_Tune4C_13TeV_pythia8/MINIAODSIM/castor_PU20bx25_POSTLS170_V5-v1/00000/628FEF3C-7D06-E411-BCAA-0025905964B6.root',
       '/store/mc/Spring14miniaod/QCD_Pt-80to120_Tune4C_13TeV_pythia8/MINIAODSIM/castor_PU20bx25_POSTLS170_V5-v1/00000/64077E88-7906-E411-9292-00261894393F.root',
       '/store/mc/Spring14miniaod/QCD_Pt-80to120_Tune4C_13TeV_pythia8/MINIAODSIM/castor_PU20bx25_POSTLS170_V5-v1/00000/683BFE3F-7D06-E411-8ADE-002590593872.root',
       '/store/mc/Spring14miniaod/QCD_Pt-80to120_Tune4C_13TeV_pythia8/MINIAODSIM/castor_PU20bx25_POSTLS170_V5-v1/00000/72E23188-7806-E411-A47C-003048FFD740.root',
       '/store/mc/Spring14miniaod/QCD_Pt-80to120_Tune4C_13TeV_pythia8/MINIAODSIM/castor_PU20bx25_POSTLS170_V5-v1/00000/7451F03E-7E06-E411-A7F3-003048FF9AA6.root',
       '/store/mc/Spring14miniaod/QCD_Pt-80to120_Tune4C_13TeV_pythia8/MINIAODSIM/castor_PU20bx25_POSTLS170_V5-v1/00000/78508372-7806-E411-80E2-003048FFD744.root',
       '/store/mc/Spring14miniaod/QCD_Pt-80to120_Tune4C_13TeV_pythia8/MINIAODSIM/castor_PU20bx25_POSTLS170_V5-v1/00000/785B2948-7E06-E411-BD1A-0025905A60F4.root',
       '/store/mc/Spring14miniaod/QCD_Pt-80to120_Tune4C_13TeV_pythia8/MINIAODSIM/castor_PU20bx25_POSTLS170_V5-v1/00000/7E59653B-7D06-E411-B14D-0026189438BA.root',
       '/store/mc/Spring14miniaod/QCD_Pt-80to120_Tune4C_13TeV_pythia8/MINIAODSIM/castor_PU20bx25_POSTLS170_V5-v1/00000/8450854F-7E06-E411-83EE-003048FFD740.root',
       '/store/mc/Spring14miniaod/QCD_Pt-80to120_Tune4C_13TeV_pythia8/MINIAODSIM/castor_PU20bx25_POSTLS170_V5-v1/00000/848B6C8C-7E06-E411-86B8-00261894390B.root',
       '/store/mc/Spring14miniaod/QCD_Pt-80to120_Tune4C_13TeV_pythia8/MINIAODSIM/castor_PU20bx25_POSTLS170_V5-v1/00000/86842083-7906-E411-A0BB-0025905B855E.root',
       '/store/mc/Spring14miniaod/QCD_Pt-80to120_Tune4C_13TeV_pythia8/MINIAODSIM/castor_PU20bx25_POSTLS170_V5-v1/00000/88666A88-7D06-E411-9BC1-0025905A608E.root',
       '/store/mc/Spring14miniaod/QCD_Pt-80to120_Tune4C_13TeV_pythia8/MINIAODSIM/castor_PU20bx25_POSTLS170_V5-v1/00000/88A33D14-7D06-E411-B78C-002618943807.root',
       '/store/mc/Spring14miniaod/QCD_Pt-80to120_Tune4C_13TeV_pythia8/MINIAODSIM/castor_PU20bx25_POSTLS170_V5-v1/00000/8AFA338B-7F06-E411-80C4-0025905AA9F0.root',
       '/store/mc/Spring14miniaod/QCD_Pt-80to120_Tune4C_13TeV_pythia8/MINIAODSIM/castor_PU20bx25_POSTLS170_V5-v1/00000/8C220497-7906-E411-AA9E-0026189438C0.root',
       '/store/mc/Spring14miniaod/QCD_Pt-80to120_Tune4C_13TeV_pythia8/MINIAODSIM/castor_PU20bx25_POSTLS170_V5-v1/00000/90FD6BE1-7E06-E411-97D7-003048FFD756.root',
       '/store/mc/Spring14miniaod/QCD_Pt-80to120_Tune4C_13TeV_pythia8/MINIAODSIM/castor_PU20bx25_POSTLS170_V5-v1/00000/9845260A-7E06-E411-82BE-0025905A60BE.root',
       '/store/mc/Spring14miniaod/QCD_Pt-80to120_Tune4C_13TeV_pythia8/MINIAODSIM/castor_PU20bx25_POSTLS170_V5-v1/00000/98C0C413-7E06-E411-83A3-002590596484.root',
       '/store/mc/Spring14miniaod/QCD_Pt-80to120_Tune4C_13TeV_pythia8/MINIAODSIM/castor_PU20bx25_POSTLS170_V5-v1/00000/9AB30E84-7906-E411-8092-002618943914.root',
       '/store/mc/Spring14miniaod/QCD_Pt-80to120_Tune4C_13TeV_pythia8/MINIAODSIM/castor_PU20bx25_POSTLS170_V5-v1/00000/9AC956B0-8006-E411-A07F-00259059642A.root',
       '/store/mc/Spring14miniaod/QCD_Pt-80to120_Tune4C_13TeV_pythia8/MINIAODSIM/castor_PU20bx25_POSTLS170_V5-v1/00000/9AF6A932-7E06-E411-947C-0025905A48F0.root',
       '/store/mc/Spring14miniaod/QCD_Pt-80to120_Tune4C_13TeV_pythia8/MINIAODSIM/castor_PU20bx25_POSTLS170_V5-v1/00000/9C4350C4-7E06-E411-8546-003048FFD760.root',
       '/store/mc/Spring14miniaod/QCD_Pt-80to120_Tune4C_13TeV_pythia8/MINIAODSIM/castor_PU20bx25_POSTLS170_V5-v1/00000/9C8018C4-7E06-E411-8D8C-003048FFCBA8.root',
       '/store/mc/Spring14miniaod/QCD_Pt-80to120_Tune4C_13TeV_pythia8/MINIAODSIM/castor_PU20bx25_POSTLS170_V5-v1/00000/9CE8DA6A-8006-E411-8166-002590596484.root',
       '/store/mc/Spring14miniaod/QCD_Pt-80to120_Tune4C_13TeV_pythia8/MINIAODSIM/castor_PU20bx25_POSTLS170_V5-v1/00000/9ECA543C-7906-E411-8871-00261894392F.root',
       '/store/mc/Spring14miniaod/QCD_Pt-80to120_Tune4C_13TeV_pythia8/MINIAODSIM/castor_PU20bx25_POSTLS170_V5-v1/00000/A0AEEAD8-7E06-E411-9378-0025905A48F2.root',
       '/store/mc/Spring14miniaod/QCD_Pt-80to120_Tune4C_13TeV_pythia8/MINIAODSIM/castor_PU20bx25_POSTLS170_V5-v1/00000/A43514C1-7D06-E411-BBB3-0025905B8596.root',
       '/store/mc/Spring14miniaod/QCD_Pt-80to120_Tune4C_13TeV_pythia8/MINIAODSIM/castor_PU20bx25_POSTLS170_V5-v1/00000/AA96FDCA-7D06-E411-9623-0025905A48D0.root',
       '/store/mc/Spring14miniaod/QCD_Pt-80to120_Tune4C_13TeV_pythia8/MINIAODSIM/castor_PU20bx25_POSTLS170_V5-v1/00000/AC2FA643-8006-E411-9FF8-00261894396F.root',
       '/store/mc/Spring14miniaod/QCD_Pt-80to120_Tune4C_13TeV_pythia8/MINIAODSIM/castor_PU20bx25_POSTLS170_V5-v1/00000/B405AE7E-7D06-E411-B204-0026189438B8.root',
       '/store/mc/Spring14miniaod/QCD_Pt-80to120_Tune4C_13TeV_pythia8/MINIAODSIM/castor_PU20bx25_POSTLS170_V5-v1/00000/B4D3E70C-7E06-E411-9091-00261894389E.root',
       '/store/mc/Spring14miniaod/QCD_Pt-80to120_Tune4C_13TeV_pythia8/MINIAODSIM/castor_PU20bx25_POSTLS170_V5-v1/00000/B62B4E15-7906-E411-A315-00261894394A.root',
       '/store/mc/Spring14miniaod/QCD_Pt-80to120_Tune4C_13TeV_pythia8/MINIAODSIM/castor_PU20bx25_POSTLS170_V5-v1/00000/B6C48B39-7E06-E411-9877-0026189438AE.root',
       '/store/mc/Spring14miniaod/QCD_Pt-80to120_Tune4C_13TeV_pythia8/MINIAODSIM/castor_PU20bx25_POSTLS170_V5-v1/00000/B87E47B2-7E06-E411-8DA2-003048FFCC2C.root',
       '/store/mc/Spring14miniaod/QCD_Pt-80to120_Tune4C_13TeV_pythia8/MINIAODSIM/castor_PU20bx25_POSTLS170_V5-v1/00000/BAA95C7A-7E06-E411-A17F-0025905938A4.root',
       '/store/mc/Spring14miniaod/QCD_Pt-80to120_Tune4C_13TeV_pythia8/MINIAODSIM/castor_PU20bx25_POSTLS170_V5-v1/00000/BE0A20CB-7806-E411-9999-0025905A60DE.root',
       '/store/mc/Spring14miniaod/QCD_Pt-80to120_Tune4C_13TeV_pythia8/MINIAODSIM/castor_PU20bx25_POSTLS170_V5-v1/00000/C272D00D-7E06-E411-9B42-002618943807.root',
       '/store/mc/Spring14miniaod/QCD_Pt-80to120_Tune4C_13TeV_pythia8/MINIAODSIM/castor_PU20bx25_POSTLS170_V5-v1/00000/C28C9395-7D06-E411-8A1B-0025905A6060.root',
       '/store/mc/Spring14miniaod/QCD_Pt-80to120_Tune4C_13TeV_pythia8/MINIAODSIM/castor_PU20bx25_POSTLS170_V5-v1/00000/C48D20DB-7F06-E411-9B00-0025905AA9CC.root',
       '/store/mc/Spring14miniaod/QCD_Pt-80to120_Tune4C_13TeV_pythia8/MINIAODSIM/castor_PU20bx25_POSTLS170_V5-v1/00000/C49301F6-7D06-E411-B1A6-0025905A60DA.root',
       '/store/mc/Spring14miniaod/QCD_Pt-80to120_Tune4C_13TeV_pythia8/MINIAODSIM/castor_PU20bx25_POSTLS170_V5-v1/00000/C65FA61E-7E06-E411-9F5E-002618943809.root',
       '/store/mc/Spring14miniaod/QCD_Pt-80to120_Tune4C_13TeV_pythia8/MINIAODSIM/castor_PU20bx25_POSTLS170_V5-v1/00000/C68CE927-7D06-E411-89E3-0026189438AE.root',
       '/store/mc/Spring14miniaod/QCD_Pt-80to120_Tune4C_13TeV_pythia8/MINIAODSIM/castor_PU20bx25_POSTLS170_V5-v1/00000/CAE54E7E-7806-E411-A76D-002590596468.root',
       '/store/mc/Spring14miniaod/QCD_Pt-80to120_Tune4C_13TeV_pythia8/MINIAODSIM/castor_PU20bx25_POSTLS170_V5-v1/00000/CCAB4D94-7D06-E411-8040-0025905A612E.root',
       '/store/mc/Spring14miniaod/QCD_Pt-80to120_Tune4C_13TeV_pythia8/MINIAODSIM/castor_PU20bx25_POSTLS170_V5-v1/00000/D664714B-8006-E411-B2B6-0025905A48C0.root',
       '/store/mc/Spring14miniaod/QCD_Pt-80to120_Tune4C_13TeV_pythia8/MINIAODSIM/castor_PU20bx25_POSTLS170_V5-v1/00000/D827DE89-7906-E411-8CD2-003048FFD756.root',
       '/store/mc/Spring14miniaod/QCD_Pt-80to120_Tune4C_13TeV_pythia8/MINIAODSIM/castor_PU20bx25_POSTLS170_V5-v1/00000/DE949D1B-7E06-E411-AC25-003048FFCC1E.root',
       '/store/mc/Spring14miniaod/QCD_Pt-80to120_Tune4C_13TeV_pythia8/MINIAODSIM/castor_PU20bx25_POSTLS170_V5-v1/00000/DEA4AFAD-7E06-E411-9C23-0025905A60B6.root',
       '/store/mc/Spring14miniaod/QCD_Pt-80to120_Tune4C_13TeV_pythia8/MINIAODSIM/castor_PU20bx25_POSTLS170_V5-v1/00000/DEFEB418-7D06-E411-AAAE-00261894389E.root',
       '/store/mc/Spring14miniaod/QCD_Pt-80to120_Tune4C_13TeV_pythia8/MINIAODSIM/castor_PU20bx25_POSTLS170_V5-v1/00000/E0A2F6A2-7E06-E411-BD96-003048FFD76E.root',
       '/store/mc/Spring14miniaod/QCD_Pt-80to120_Tune4C_13TeV_pythia8/MINIAODSIM/castor_PU20bx25_POSTLS170_V5-v1/00000/E244DD2C-7906-E411-A6E3-002590593878.root',
       '/store/mc/Spring14miniaod/QCD_Pt-80to120_Tune4C_13TeV_pythia8/MINIAODSIM/castor_PU20bx25_POSTLS170_V5-v1/00000/E2823481-8006-E411-A604-0025905A6068.root',
       '/store/mc/Spring14miniaod/QCD_Pt-80to120_Tune4C_13TeV_pythia8/MINIAODSIM/castor_PU20bx25_POSTLS170_V5-v1/00000/EA4230B0-7E06-E411-B0DD-003048FFCBFC.root',
       '/store/mc/Spring14miniaod/QCD_Pt-80to120_Tune4C_13TeV_pythia8/MINIAODSIM/castor_PU20bx25_POSTLS170_V5-v1/00000/EC20313E-7906-E411-83A5-0025905B85E8.root',
       '/store/mc/Spring14miniaod/QCD_Pt-80to120_Tune4C_13TeV_pythia8/MINIAODSIM/castor_PU20bx25_POSTLS170_V5-v1/00000/F0A0A742-7F06-E411-8F59-0025905A6094.root',
       '/store/mc/Spring14miniaod/QCD_Pt-80to120_Tune4C_13TeV_pythia8/MINIAODSIM/castor_PU20bx25_POSTLS170_V5-v1/00000/F4143F8E-7906-E411-A178-0026189438B1.root',
       '/store/mc/Spring14miniaod/QCD_Pt-80to120_Tune4C_13TeV_pythia8/MINIAODSIM/castor_PU20bx25_POSTLS170_V5-v1/00000/F44A4758-8006-E411-A829-002590596486.root',
       '/store/mc/Spring14miniaod/QCD_Pt-80to120_Tune4C_13TeV_pythia8/MINIAODSIM/castor_PU20bx25_POSTLS170_V5-v1/00000/F44FB177-7E06-E411-9691-0025905B860E.root',
       '/store/mc/Spring14miniaod/QCD_Pt-80to120_Tune4C_13TeV_pythia8/MINIAODSIM/castor_PU20bx25_POSTLS170_V5-v1/00000/F83138E3-7D06-E411-8BEF-002618FDA28E.root' ] );


secFiles.extend( [
               ] )


