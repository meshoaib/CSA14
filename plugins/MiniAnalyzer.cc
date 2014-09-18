// -*- C++ -*-
//
// Package:    UserCode/MiniAnalyzer
// Class:      MiniAnalyzer
// 
/**\class MiniAnalyzer MiniAnalyzer.cc Test/MiniAnalyzer/plugins/MiniAnalyzer.cc

   Description: [one line class summary]

   Implementation:
   [Notes on implementation]
*/
//
// Original Author:  Qamar Ul Hassan
//         Created:  Sun, 13 Jul 2014 06:22:18 GMT
//
//

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "FWCore/Common/interface/TriggerNames.h"

#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/MET.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"
#include "DataFormats/PatCandidates/interface/TriggerObjectStandAlone.h"
#include "DataFormats/PatCandidates/interface/PackedTriggerPrescales.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/Common/interface/TriggerResults.h"

#include "EgammaAnalysis/ElectronTools/interface/EGammaCutBasedEleId.h"

#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "UserCode/TopAnalysis/interface/MiniEvent.h"

#include "TH1.h"
#include "TH1F.h"
#include "TH2F.h"

#include <vector>
#include <unordered_map>
#include <memory>
#include <cmath>


using namespace edm;
using namespace std;
using namespace reco;
using namespace pat; 

//
// class declaration
//

class MiniAnalyzer : public edm::EDAnalyzer {
public:
  explicit MiniAnalyzer(const edm::ParameterSet&);
  ~MiniAnalyzer();

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


private:
  virtual void beginJob() override;
  virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
  virtual void endJob() override;

  //virtual void beginRun(edm::Run const&, edm::EventSetup const&) override;
  //virtual void endRun(edm::Run const&, edm::EventSetup const&) override;
  //virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;
  //virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;

  // ----------member data ---------------------------
      
  edm::EDGetTokenT<edm::TriggerResults> triggerBits_;
  edm::EDGetTokenT<reco::VertexCollection> vtxToken_;
  edm::EDGetTokenT<double> rhoToken_;
  edm::EDGetTokenT<pat::MuonCollection> muonToken_;
  edm::EDGetTokenT<pat::ElectronCollection> electronToken_;
  edm::EDGetTokenT<pat::JetCollection> jetToken_;
  edm::EDGetTokenT<pat::METCollection> metToken_;
  edm::EDGetTokenT<pat::PackedCandidateCollection> pfToken_;

  //TH1F* fHistnew_Histo;

  std::unordered_map<std::string,TH1F*> histContainer_;
  std::unordered_map<std::string,TH2F*> histContainer2d_; 

  TTree *tree_;
  MiniEvent_t ev_;
};

//
// constants, enums and typedefs
//

//
// static data member definitions
//


//
// constructors and destructor
//
MiniAnalyzer::MiniAnalyzer(const edm::ParameterSet& iConfig) :
  triggerBits_(consumes<edm::TriggerResults>(iConfig.getParameter<edm::InputTag>("triggerBits"))),
  vtxToken_(consumes<reco::VertexCollection>(iConfig.getParameter<edm::InputTag>("vertices"))),
  rhoToken_(consumes<double>(iConfig.getParameter<edm::InputTag>("rho"))),
  muonToken_(consumes<pat::MuonCollection>(iConfig.getParameter<edm::InputTag>("muons"))),
  electronToken_(consumes<pat::ElectronCollection>(iConfig.getParameter<edm::InputTag>("electrons"))),
  jetToken_(consumes<pat::JetCollection>(iConfig.getParameter<edm::InputTag>("jets"))),
  metToken_(consumes<pat::METCollection>(iConfig.getParameter<edm::InputTag>("mets"))),
 pfToken_(consumes<pat::PackedCandidateCollection>(iConfig.getParameter<edm::InputTag>("pfCands")))
{
  //now do what ever initialization is needed

}


MiniAnalyzer::~MiniAnalyzer()
{
 
  // do anything here that needs to be done at desctruction time
  // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called for each event  ------------
void
MiniAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  //  bool isData = iEvent.isRealData();
  ev_.run     = iEvent.id().run();
  ev_.lumi    = iEvent.luminosityBlock();
  ev_.event   = iEvent.id().event();

  edm::Handle<edm::TriggerResults> triggerBits;
  iEvent.getByToken(triggerBits_, triggerBits);
  /*
    const edm::TriggerNames &names = iEvent.triggerNames(*triggerBits);
    std::cout << "\n === TRIGGER PATHS === " << std::endl;
    for (unsigned int i = 0, n = triggerBits->size(); i < n; ++i) {
    std::cout << "Trigger " << names.triggerName(i) <<
    //                ", prescale " << triggerPrescales->getPrescaleForIndex(i) <<
    ": " << (triggerBits->accept(i) ? "PASS" : "fail (or not run)")
    << std::endl;
    }
  */
  edm::Handle<reco::VertexCollection> vertices;
  iEvent.getByToken(vtxToken_, vertices);
  if (vertices->empty()) return; // skip the event if no PV found
  const reco::Vertex &primVtx = vertices->front();
  ev_.nvtx=vertices->size();
  
  edm::Handle< double > rhoH;
  iEvent.getByToken(rhoToken_,rhoH);
  float rho=*rhoH;
  ev_.rho=rho;
	    
  histContainer_["cutflow"]->Fill(0);
  histContainer_["ecutflow"]->Fill(0);
  histContainer_["mucutflow"]->Fill(0);


  //
  // MUONS
  // cf. https://twiki.cern.ch/twiki/bin/view/CMSPublic/SWGuideMuonId  
  //
  float leptonpt(0), leptonphi(0);
  edm::Handle<pat::MuonCollection> muons;
  iEvent.getByToken(muonToken_, muons);
  histContainer_["nrecomuons"]->Fill(muons->size());
  std::vector<const pat::Muon *> selectedMuons,vetoMuons;        
  for (const pat::Muon &mu : *muons) { 
    
    //kinematics
    bool passPt( mu.pt() > 26 );
    bool passVetoPt( mu.pt()>10 );
    bool passEta(fabs(mu.eta()) < 2.1 );
    
    //distance to the PV
    float dz(fabs( mu.vertex().z() - primVtx.z())); 
    bool passDB( mu.dB()<0.2 && dz<0.5 );
    
    //isolation
    float relchIso((mu.chargedHadronIso())/mu.pt());
    bool passIso( relchIso<0.05 );

    if( mu.isPFMuon() 
	&& mu.isGlobalMuon() 
	&& mu.normChi2() < 10 
	&& mu.innerTrack()->hitPattern().trackerLayersWithMeasurement() > 5 
	&& mu.globalTrack()->hitPattern().numberOfValidMuonHits() > 0 
	&& mu.innerTrack()->hitPattern().numberOfValidPixelHits() > 0 
	&& mu.numberOfMatchedStations() > 1 )
      {
	if(passDB && passPt && passEta) 
	  {
	    histContainer_["muonchreliso"]->Fill(relchIso);
	    histContainer_["muonchiso"]->Fill(mu.chargedHadronIso());
	    histContainer_["muonneuthadiso"]->Fill(mu.neutralHadronIso());
	    histContainer_["muonphotoniso"]->Fill(mu.photonIso());
	    histContainer_["muonpuchiso"]->Fill(mu.puChargedHadronIso());
	  }
	if(passIso && passPt && passEta) 
	  {
	    histContainer_["muondb"]->Fill(mu.dB());
	    histContainer_["muondz"]->Fill(dz);
	  }
	if(passDB && passIso)
	  {
	    if( passEta ) histContainer_["muonpt"]->Fill(mu.pt()); //N-1 plot
	    if( passPt )        histContainer_["muoneta"]->Fill(fabs(mu.eta()));
	    if( passPt && passEta )
	      {
		selectedMuons.push_back( &mu );
		leptonpt = mu.pt();
		leptonphi = mu.phi();
		
		//save the selected lepton
		ev_.l_id=13;
		ev_.l_charge=mu.charge();
		ev_.l_pt=mu.pt();
		ev_.l_eta=mu.eta();
		ev_.l_phi=mu.phi();
		ev_.l_chargedHadronIso=mu.chargedHadronIso();
                ev_.l_neutralHadronIso=mu.neutralHadronIso();
                ev_.l_photonIso=mu.photonIso();
                ev_.l_puChargedHadronIso=mu.puChargedHadronIso();
	      }
	    else if(passVetoPt && passEta)
	      {
		vetoMuons.push_back( &mu );
	      }
	  }
      }
  }
  histContainer_["nselmuons"]->Fill(selectedMuons.size());

  //
  // ELECTRONS
  // cf. https://twiki.cern.ch/twiki/bin/view/CMS/EgammaCutBasedIdentification  
  //
  edm::Handle<pat::ElectronCollection> electrons;
  iEvent.getByToken(electronToken_, electrons);
  histContainer_["nrecoelectrons"]->Fill(electrons->size());
  std::vector<const pat::Electron *> selectedElectrons,vetoElectrons;
  for (const pat::Electron &el : *electrons) {        	

    //kinematics cuts
    bool passPt(el.pt()>30);
    bool passVetoPt(el.pt()>20);
    bool passEta(fabs(el.eta()) < 2.5 && (fabs(el.superCluster()->eta()) < 1.4442 || fabs(el.superCluster()->eta()) > 1.5660));

    //use a cut based id
    bool passVetoId( EgammaCutBasedEleId::PassWP( EgammaCutBasedEleId::VETO, el.isEB(), el.pt(), el.eta(),
                                                  el.deltaEtaSuperClusterTrackAtVtx(), 
						  el.deltaPhiSuperClusterTrackAtVtx(),
						  el.sigmaIetaIeta(),
						  el.hadronicOverEm(),
						  (1.0/el.ecalEnergy() - el.eSuperClusterOverP()/el.ecalEnergy()),
						  fabs(el.gsfTrack()->dxy(primVtx.position())),
						  fabs(el.gsfTrack()->dz(primVtx.position())),
                                                  0., 0., 0., 
						  !(el.passConversionVeto()), 
						  el.gsfTrack()->trackerExpectedHitsInner().numberOfLostHits(),
						  rho) );
    bool passTightId( EgammaCutBasedEleId::PassWP( EgammaCutBasedEleId::TIGHT, el.isEB(), el.pt(), el.eta(),
						   el.deltaEtaSuperClusterTrackAtVtx(), 
						   el.deltaPhiSuperClusterTrackAtVtx(),
						   el.sigmaIetaIeta(),
						   el.hadronicOverEm(),
						   (1.0/el.ecalEnergy() - el.eSuperClusterOverP()/el.ecalEnergy()),
						   fabs(el.gsfTrack()->dxy(primVtx.position())),
						   fabs(el.gsfTrack()->dz(primVtx.position())),
						   0., 0., 0., 
						   !(el.passConversionVeto()), 
						   el.gsfTrack()->trackerExpectedHitsInner().numberOfLostHits(),
						   rho) );
    
    //isolation
    float relchIso((el.chargedHadronIso())/el.pt());
    bool passIso( relchIso<0.05 );
    bool passVetoIso( relchIso<0.15 );
    
    if( el.gsfTrack()->trackerExpectedHitsInner().numberOfHits() <= 0 
	&& el.dB() < 0.02 
	&& el.passConversionVeto() == true 
	)
      {
	
	if(passPt && passEta && passTightId)
	  {
	    histContainer_["electronchreliso"]->Fill(relchIso);
	    histContainer_["electronchiso"]->Fill(el.chargedHadronIso());
	    histContainer_["electronneuthadiso"]->Fill(el.neutralHadronIso());
	    histContainer_["electronphotoniso"]->Fill(el.photonIso());
	    histContainer_["electronpuchiso"]->Fill(el.puChargedHadronIso());
	  }
	
	if(passTightId && passIso) 
	  {
	    if(passEta) histContainer_["electronpt"]->Fill(el.pt());    //N-1 plot
	    if(passPt) histContainer_["electroneta"]->Fill(fabs(el.eta()));  //N-1 plot
	  }

	if(passPt && passEta && passTightId && passIso)
	  {
	    selectedElectrons.push_back(&el);
	    leptonpt = el.pt();
	    leptonphi = el.phi();
	    
	    //save the selected lepton
	    ev_.l_id=11;
	    ev_.l_charge=el.charge();
	    ev_.l_pt=el.pt();
	    ev_.l_eta=el.eta();
	    ev_.l_phi=el.phi();
	    ev_.l_chargedHadronIso=el.chargedHadronIso();
	    ev_.l_neutralHadronIso=el.neutralHadronIso();
	    ev_.l_photonIso=el.photonIso();
	    ev_.l_puChargedHadronIso=el.puChargedHadronIso();
	  }
	else if(passVetoPt && passEta && passVetoId && passVetoIso)
	  {
	    vetoElectrons.push_back(&el);
	  }
      }
  }
  histContainer_["nselelectrons"]->Fill(selectedElectrons.size());
  
  //require only 1 tight lepton in the event
  int nSelectedLeptons(selectedElectrons.size()+selectedMuons.size());
  if(nSelectedLeptons>1 || nSelectedLeptons==0) return;
  histContainer_["cutflow"]->Fill(1);
  if(selectedElectrons.size()==1)  histContainer_["ecutflow"]->Fill(1);
  if(selectedMuons.size()==1)      histContainer_["mucutflow"]->Fill(1);

  //require no other leptons in the event
  int nVetoLeptons(vetoElectrons.size()+vetoMuons.size());
  if(nVetoLeptons>0) return;
  histContainer_["cutflow"]->Fill(2);
  if(selectedElectrons.size()==1)  histContainer_["ecutflow"]->Fill(2);
  if(selectedMuons.size()==1)      histContainer_["mucutflow"]->Fill(2);  

  //
  // JETS
  //
  uint32_t nCSVMtags(0);
  edm::Handle<pat::JetCollection> jets;
  iEvent.getByToken(jetToken_, jets);
  std::vector<const pat::Jet *> selectedJets;
  histContainer_["nrecojets"]->Fill(jets->size());
  ev_.nj=0;
  for (const pat::Jet &j : *jets) {
	  
    float dR2lepton= selectedMuons.size()==1 ? 
      deltaR(j,*(selectedMuons[0])) :
      deltaR(j,*(selectedElectrons[0]));
    float rawEnergy(j.energy()*j.jecFactor("Uncorrected"));
    if ( j.numberOfDaughters() > 1 
	 && (j.neutralHadronEnergy() + j.HFHadronEnergy())/rawEnergy < 0.99 
	 && j.neutralEmEnergyFraction() < 0.99 
	 && (j.chargedEmEnergyFraction() < 0.99 || fabs(j.eta()) >= 2.4)
	 && (j.chargedHadronEnergyFraction() > 0. || fabs(j.eta()) >= 2.4) 
	 && (j.chargedMultiplicity() > 0 || fabs(j.eta()) >= 2.4)
	 && dR2lepton>0.4)
      {
	if(fabs(j.eta()) < 2.5) histContainer_["jetpt"]->Fill(j.pt());    //N-1 plots
	if(j.pt() > 30)         histContainer_["jeteta"]->Fill(fabs(j.eta())); 
	if( fabs(j.eta()) < 2.5 && j.pt() > 30)
	  {
	    selectedJets.push_back( &j );
	    float csv=j.bDiscriminator("combinedSecondaryVertexBJetTags");
	    if(csv>0.679) nCSVMtags++;
	    histContainer_["jetcsv"]->Fill(csv);
	    histContainer_["jetpileupid"]->Fill(j.userFloat("pileupJetId:fullDiscriminant"));
	    float svtxmass=j.userFloat("vtxMass");
	    histContainer_["jetsecvtxmass"]->Fill(svtxmass);
	    int vtxNtracks=j.userFloat("vtxNtracks");
	    histContainer_["jetvtxNtracks"]->Fill(vtxNtracks);
	    float vtx3DVal=j.userFloat("vtx3DVal");
	    histContainer_["jetvtx3DVal"]->Fill(vtx3DVal);
	    float vtx3DSig=j.userFloat("vtx3DSig");
	    histContainer_["jetvtx3DSig"]->Fill(vtx3DSig);

	    ev_.j_pt[ev_.nj]=j.pt();
	    ev_.j_eta[ev_.nj]=j.eta();
	    ev_.j_phi[ev_.nj]=j.phi();
	    ev_.j_csv[ev_.nj]=csv;
	    ev_.j_vtxmass[ev_.nj]=svtxmass;
	    ev_.j_vtxNtracks[ev_.nj]=vtxNtracks;
	    ev_.j_vtx3DVal[ev_.nj]=vtx3DVal;
	    ev_.j_vtx3DSig[ev_.nj]=vtx3DSig;
	    ev_.j_puid[ev_.nj]=j.userFloat("pileupJetId:fullDiscriminant");
	    ev_.j_flav[ev_.nj]=j.partonFlavour();
	    const reco::Candidate *genParton = j.genParton();
	    ev_.j_pid[ev_.nj]=genParton ? genParton->pdgId() : 0;
	    ev_.nj++;

	  }
      }
  }
  histContainer_["nseljets"]->Fill(selectedJets.size());


  //
  // MET
  //
  edm::Handle<pat::METCollection> mets;
  iEvent.getByToken(metToken_, mets);
  float metpt = mets->at(0).pt();
  float metphi = mets->at(0).phi();
  float dphi_met_lepton = deltaPhi(leptonphi, metphi); // use the function to restrict to the 0,pi range
  float mt=sqrt(2*leptonpt*metpt*(1-cos(dphi_met_lepton)));

  //save to ttree
  ev_.met_pt=metpt;
  ev_.met_phi=metphi;
  ev_.mt=mt;

  //
  // FINAL SELECTION PLOTS
  //
  if(selectedJets.size()>=2)
    {
      histContainer_["cutflow"]->Fill(3);
      if(selectedElectrons.size()==1)  histContainer_["ecutflow"]->Fill(3);
      if(selectedMuons.size()==1)      histContainer_["mucutflow"]->Fill(3);  
      tree_->Fill();
    }
  if(selectedJets.size()>=3) 
    {
      histContainer_["cutflow"]->Fill(4);
      if(selectedElectrons.size()==1)  histContainer_["ecutflow"]->Fill(4);
      if(selectedMuons.size()==1)      histContainer_["mucutflow"]->Fill(4);  
      histContainer_["mt_3"]->Fill(mt);
    }
  if(selectedJets.size()>=4) 
    {
      histContainer_["mt_4"]->Fill(mt);
      histContainer_["cutflow"]->Fill(5);
      if(selectedElectrons.size()==1)  histContainer_["ecutflow"]->Fill(5);
      if(selectedMuons.size()==1)      histContainer_["mucutflow"]->Fill(5);  
      histContainer_["ncsvmjets"]->Fill(nCSVMtags);
      histContainer_["nvertices"]->Fill(vertices->size());
      if(nCSVMtags>=1) 
	{
	  histContainer_["cutflow"]->Fill(6);
	  if(selectedElectrons.size()==1)  histContainer_["ecutflow"]->Fill(6);
	  if(selectedMuons.size()==1)      histContainer_["mucutflow"]->Fill(6);  
	  histContainer_["mt_41b"] ->Fill(mt);
	}
      if(nCSVMtags>=2) 
	{
	  histContainer_["cutflow"]->Fill(7);
	  if(selectedElectrons.size()==1)  histContainer_["ecutflow"]->Fill(7);
	  if(selectedMuons.size()==1)      histContainer_["mucutflow"]->Fill(7);  
	  histContainer_["mt_42b"] ->Fill(mt);
	  histContainer_["metpt"]->Fill(metpt);
	  histContainer_["metphi"] ->Fill(metphi);
	  histContainer_["dphimetlepton"]->Fill(dphi_met_lepton );
	}
    }
//  float pfpt=0.0;
//int id = 0.0;
    edm::Handle<pat::PackedCandidateCollection> pfs;
    iEvent.getByToken(pfToken_, pfs);
    // now loop on pf candidates
    for (unsigned int i = 0, n = pfs->size(); i < n; ++i) {
        const pat::PackedCandidate &pf = (*pfs)[i];
        if (pf.fromPV() > 0 && pf.charge()!= 0) 

histContainer_["pfpt"]->Fill(pf.pt());
histContainer_["pfpx"]->Fill(pf.px());
histContainer_["pfpy"]->Fill(pf.py());
histContainer_["pfpz"]->Fill(pf.pz());
histContainer_["sum_pt"]->Fill(sqrt(pow(pf.px(),2)+pow(pf.py(),2)+pow(pf.pz(),2)));
}
}
/*
if (pf.pdgId() == 13){
		pfpt = pf.pt();}
		cout<<"ID PF Candidates: "<<pf.pdgId()<<endl; //}
		cout<<"Pt of PF Candidates: "<<pfpt<<endl; 
		cout<<"All is well: "<<endl;}*/
//   }
//}
//


// ------------ method called once each job just before starting event loop  ------------
void 
MiniAnalyzer::beginJob()
{
  edm::Service<TFileService> fs;
  histContainer_["cutflow"]   = fs->make<TH1F>("cutflow",    ";Selection cut;Events", 8, 0., 8.); 
  histContainer_["ecutflow"]  = fs->make<TH1F>("ecutflow",   ";Selection cut;Events", 8, 0., 8.); 
  histContainer_["mucutflow"] = fs->make<TH1F>("mucutflow",  ";Selection cut;Events", 8, 0., 8.); 
  TString steps[]={"reco","=1 good lepton","=0 loose leptons","#geq2 jets","#geq3 jets","#geq4 jets","#geq1 b-tag","#geq2 b-tags"};
  for(size_t i=0; i<sizeof(steps)/sizeof(TString); i++) 
    {
      histContainer_["cutflow"]   -> GetXaxis()->SetBinLabel(i+1,steps[i]);
      histContainer_["ecutflow"]  -> GetXaxis()->SetBinLabel(i+1,steps[i]);
      histContainer_["mucutflow"] -> GetXaxis()->SetBinLabel(i+1,steps[i]);
    }

  histContainer_["nvertices"] = fs->make<TH1F>("nvertices",    ";# vertices;Events", 100, 0., 100.); 

  histContainer_["nrecomuons"]     = fs->make<TH1F>("nrecomuons",      ";# reconstructed muons; Events",             10, 0., 10.);
  histContainer_["muonpt"]         = fs->make<TH1F>("muonpt",          ";Transverse momentum [GeV];# muons",         100, 0., 300.);
  histContainer_["muoneta"]        = fs->make<TH1F>("muoneta",         ";Pseudo-rapidity;#muons ",                   100, 0, 3.);
  histContainer_["muondb"]         = fs->make<TH1F>("muondb",          ";d_{0} [cm];# muons",         100, 0.,0.3);
  histContainer_["muondz"]         = fs->make<TH1F>("muondz",          ";|d_{z}| [cm];# muons",       100, 0.,0.6);
  histContainer_["muonchreliso"]   = fs->make<TH1F>("muonchreliso",       ";Relative charged hadron isolation;# muons ",  100, 0, 0.1);
  histContainer_["muonchiso"]      = fs->make<TH1F>("muonchiso",       ";Charged hadron isolation [GeV];# muons ",  100, 0, 10);
  histContainer_["muonpuchiso"]    = fs->make<TH1F>("muonpuchiso",     ";Pileup charged hadron isolation [GeV];# muons ",  100, 0, 10);
  histContainer_["muonneuthadiso"] = fs->make<TH1F>("muonneuthadiso",  ";Neutral hadron isolation [GeV];# muons ",  100, 0, 10);
  histContainer_["muonphotoniso"]  = fs->make<TH1F>("muonphotoniso",   ";Photon isolation [GeV];# muons ",  100, 0, 10);
  histContainer_["nselmuons"]      = fs->make<TH1F>("nselmuons",       ";# reconstructed muons;Events",              5, 0.,5.);

  histContainer_["nrecoelectrons"]    = fs->make<TH1F>("nrecoelectrons",  ";# reconstructed electrons;Events",         5, 0., 5.);
  histContainer_["electronpt"]        = fs->make<TH1F>("electronpt",      ";Transverse momentum [GeV];# electrons",    100, 0., 300.);
  histContainer_["electroneta"]       = fs->make<TH1F>("electroneta",     ";Pseudo-rapidity;#electrons",               100, 0., 3.);
  histContainer_["electronchreliso"]   = fs->make<TH1F>("electronchreliso",       ";Relative charged hadron isolation;# electrons ",  100, 0, 0.1);
  histContainer_["electronchiso"]      = fs->make<TH1F>("electronchiso",       ";Charged hadron isolation [GeV];# electrons ",  100, 0, 10);
  histContainer_["electronpuchiso"]    = fs->make<TH1F>("electronpuchiso",     ";Pileup charged hadron isolation [GeV];# electrons ",  100, 0, 10);
  histContainer_["electronneuthadiso"] = fs->make<TH1F>("electronneuthadiso",  ";Neutral hadron isolation [GeV];# electrons ",  100, 0, 10);
  histContainer_["electronphotoniso"]  = fs->make<TH1F>("electronphotoniso",   ";Photon isolation [GeV];# electrons ",  100, 0, 10);
  histContainer_["nselelectrons"]     = fs->make<TH1F>("nselelectrons",   ";# selected electrons;Events", 5, 0., 5.);

  histContainer_["nrecojets"]      = fs->make<TH1F>("nrecojets",   ";#reconstructed jets;Events", 50, 0., 50.);
  histContainer_["jetpt"]          = fs->make<TH1F>("jetpt",       ";Transverse momentum [GeV];# jets", 100, 0., 300.);
  histContainer_["jeteta"]         = fs->make<TH1F>("jeteta",      ";Pseudo-rapidity;# jets", 100, 0., 3.);
  histContainer_["jetcsv"]         = fs->make<TH1F>("jetcsv",      ";Combined secondary vertes;# jets", 100, -1.2, 1.2);
  histContainer_["jetpileupid"]    = fs->make<TH1F>("jetpileupid", ";Pileup jet id;#jets", 100, -1.2, 1.2);
  histContainer_["jetsecvtxmass"]  = fs->make<TH1F>("jetvtxMass", ";Secondary vertex mass [GeV];#jets", 100, 0., 6.);
  histContainer_["jetvtxNtracks"]  = fs->make<TH1F>("jetvtxNtracks", ";Vertex Tracks;#jets", 6, 0., 6.);  
  histContainer_["jetvtx3DVal"]  = fs->make<TH1F>("jetvtx3DVal", ";vtx3DVal [cm];#jets", 100, -5., 5.);
  histContainer_["jetvtx3DSig"]  = fs->make<TH1F>("jetvtx3DSig", ";vtx3DSig;#jets", 100, 0., 5.);
  histContainer_["nseljets"]       = fs->make<TH1F>("nseljets",    ";#selected jets;Events", 6, 3., 10.);
  histContainer_["ncsvmjets"]      = fs->make<TH1F>("ncsvmjets",    ";b-tagged jets (CSVM);Events", 10, 0., 10.);
  
  histContainer_["metpt"] = fs->make<TH1F>("metpt",    ";Missing transverse energy [GeV];Events", 100, 0., 300.);
  histContainer_["metphi"] = fs->make<TH1F>("metphi",    ";Missing transverse energy #phi [rad];Events", 50, -3.2, 3.2);
  histContainer_["dphimetlepton"] = fs->make<TH1F>("dphimetlepton",    ";#Delta#phi(MET,lepton) [rad];Events", 50, -3.2, 3.2);
 
  histContainer_["mt_3"] = fs->make<TH1F>("mt_3",    ";Transverse mass [GeV]; Events", 100, 0., 200.);
  histContainer_["mt_4"] = fs->make<TH1F>("mt_4",    ";Transverse mass [GeV]; Events", 100, 0., 200.);
  histContainer_["mt_41b"] = fs->make<TH1F>("mt_41b",    ";Transverse mass [GeV]; Events", 100, 0., 200.);
  histContainer_["mt_42b"] = fs->make<TH1F>("mt_42b",    ";Transverse mass [GeV]; Events", 100, 0., 200.);

  histContainer_["pfpt"]         = fs->make<TH1F>("pfpt",          ";Transverse momentum [GeV];# pfpacked",         100, 0., 300.);
  histContainer_["pfpx"]         = fs->make<TH1F>("pfpx",          ";X_Component_Transverse momentum [GeV];# pfpacked",         100, 0., 300.);
  histContainer_["pfpy"]         = fs->make<TH1F>("pfpy",          ";Y_Component_Transverse momentum [GeV];# pfpacked",         100, 0., 300.);
  histContainer_["pfpz"]         = fs->make<TH1F>("pfpz",          ";Z_Component_Transverse momentum [GeV];# pfpacked",         100, 0., 300.);
  histContainer_["sum_pt"]         = fs->make<TH1F>("sum_pt",          ";SUM_Transverse momentum [GeV];# pfpacked",         100, 0., 300.);
  //instruct ROOT to compute the uncertainty from the square root of weights
  //http://root.cern.ch/root/html/TH1.html#TH1:Sumw2
//  for(std::unordered_map<std::string,TH1F*>::iterator it=histContainer_.begin();   it!=histContainer_.end();   it++) it->second->Sumw2();
//  for(std::unordered_map<std::string,TH2F*>::iterator it=histContainer2d_.begin(); it!=histContainer2d_.end(); it++) it->second->Sumw2();

  //create a tree for the selected events
  tree_ = fs->make<TTree>("AnaTree", "AnaTree");
  createMiniEventTree(tree_,ev_);
}

// ------------ method called once each job just after ending the event loop  ------------
void 
MiniAnalyzer::endJob() 
{
}

// ------------ method called when starting to processes a run  ------------
/*
  void 
  MiniAnalyzer::beginRun(edm::Run const&, edm::EventSetup const&)
  {
  }
*/

// ------------ method called when ending the processing of a run  ------------
/*
  void 
  MiniAnalyzer::endRun(edm::Run const&, edm::EventSetup const&)
  {
  }
*/

// ------------ method called when starting to processes a luminosity block  ------------
/*
  void 
  MiniAnalyzer::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
  {
  }
*/

// ------------ method called when ending the processing of a luminosity block  ------------
/*
  void 
  MiniAnalyzer::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
  {
  }
*/

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
MiniAnalyzer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(MiniAnalyzer);
