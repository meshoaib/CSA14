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
 
#include <TTree.h>
#include <TClonesArray.h>
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

  //TH1F* fHistnew_Histo;

  std::unordered_map<std::string,TH1F*> histContainer_;
  std::unordered_map<std::string,TH2F*> histContainer2d_; 

  //TH1F* fHistnew_Histo;

  bool muon_selection, electron_veto, jet_selection;

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
  metToken_(consumes<pat::METCollection>(iConfig.getParameter<edm::InputTag>("mets")))
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
  //   using namespace edm;

  muon_selection = false;
  electron_veto = false;
  jet_selection = false; 

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
  
  edm::Handle< double > rhoH;
  iEvent.getByToken(rhoToken_,rhoH);
  float rho=*rhoH;
		    
  histContainer_["cutflow"]->Fill(0);


  //
  // MUONS
  // cf. https://twiki.cern.ch/twiki/bin/view/CMSPublic/SWGuideMuonId  
  //
  float mupt(0),muphi(0);
  edm::Handle<pat::MuonCollection> muons;
  iEvent.getByToken(muonToken_, muons);
  histContainer_["nrecomuons"]->Fill(muons->size());
  std::vector<const pat::Muon *> selectedMuons;        
  for (const pat::Muon &mu : *muons) { 
    if( mu.isPFMuon() 
	&& mu.isGlobalMuon() 
	&& mu.normChi2() < 10 
	&& mu.innerTrack()->hitPattern().trackerLayersWithMeasurement() > 5 
	&& mu.globalTrack()->hitPattern().numberOfValidMuonHits() > 0 
	&& mu.dB() < 0.2 
	&& mu.innerTrack()->hitPattern().numberOfValidPixelHits() > 0 
	&& mu.numberOfMatchedStations() > 1 
	&& (mu.chargedHadronIso()+max(0.,mu.neutralHadronIso()+mu.photonIso()-0.50*mu.puChargedHadronIso()))/mu.pt() < 0.12 )
      {
	
	if( fabs(mu.eta()) < 2.1 ) histContainer_["muonpt"]->Fill(mu.pt()); //N-1 plot
	if( mu.pt() > 20 )        histContainer_["muoneta"]->Fill(fabs(mu.eta()));
	if( mu.pt() > 20 && fabs(mu.eta()) < 2.1 )
	  {
	    selectedMuons.push_back( &mu );
	    mupt = mu.pt();
	    muphi = mu.phi();
	  }
      }
  }
  histContainer_["nselmuons"]->Fill(selectedMuons.size());

  //require at least one muon
  if(selectedMuons.size()==1) histContainer_["cutflow"]->Fill(1);
  else return;

  //
  // ELECTRONS
  // cf. https://twiki.cern.ch/twiki/bin/view/CMS/EgammaCutBasedIdentification  
  //
  edm::Handle<pat::ElectronCollection> electrons;
  iEvent.getByToken(electronToken_, electrons);
  histContainer_["nrecoelectrons"]->Fill(electrons->size());
  std::vector<const pat::Electron *> selectedElectrons;
  for (const pat::Electron &el : *electrons) {        	

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
    
    if( passVetoId 
	&& fabs(el.superCluster()->eta()) > 1.4442 && fabs(el.superCluster()->eta()) < 1.5660 
	&& el.gsfTrack()->trackerExpectedHitsInner().numberOfHits() <= 0 
	&& el.dB() < 0.02 
	&& el.passConversionVeto() == true 
	&& (el.chargedHadronIso()+max(0.,el.neutralHadronIso()+el.photonIso()-0.50*el.puChargedHadronIso()))/el.pt() < 0.1 )
      {
	if(fabs(el.eta()) < 2.5) histContainer_["electronpt"]->Fill(el.pt());    //N-1 plot
	if(el.pt() > 20)        histContainer_["electroneta"]->Fill(fabs(el.eta()));  //N-1 plot
	if(el.pt() > 20 && fabs(el.eta()) < 2.5)
	  {
	    selectedElectrons.push_back(&el);
	  }
      }
  }
  histContainer_["nselelectrons"]->Fill(selectedElectrons.size());
  
  //require at least one electron
  if(selectedElectrons.size()==0) histContainer_["cutflow"]->Fill(2);
  else return;

  //
  // JETS
  //
  uint32_t nCSVMtags(0);
  edm::Handle<pat::JetCollection> jets;
  iEvent.getByToken(jetToken_, jets);
  std::vector<const pat::Jet *> selectedJets;
  histContainer_["nrecojets"]->Fill(jets->size());
  for (const pat::Jet &j : *jets) {
	  
    float dR2muon=deltaR(j,*(selectedMuons[0]));
    float rawEnergy(j.energy()*j.jecFactor("Uncorrected"));
    if ( j.numberOfDaughters() > 1 
	 && (j.neutralHadronEnergy() + j.HFHadronEnergy())/rawEnergy < 0.99 
	 && j.neutralEmEnergyFraction() < 0.99 
	 && (j.chargedEmEnergyFraction() < 0.99 || fabs(j.eta()) >= 2.4)
	 && (j.chargedHadronEnergyFraction() > 0. || fabs(j.eta()) >= 2.4) 
	 && (j.chargedMultiplicity() > 0 || fabs(j.eta()) >= 2.4)
	 && dR2muon>0.5)
      {
	if(fabs(j.eta()) < 2.5) histContainer_["jetpt"]->Fill(j.pt());    //N-1 plots
	if(j.pt() > 30)        histContainer_["jeteta"]->Fill(fabs(j.eta())); 
	if( fabs(j.eta()) < 2.5 && j.pt() > 30)
	  {
	    selectedJets.push_back( &j );
	    float csv=j.bDiscriminator("combinedSecondaryVertexBJetTags");
	    if(csv>0.679) nCSVMtags++;
	    histContainer_["jetcsv"]->Fill(csv);
	    histContainer_["jetpileupid"]->Fill(j.userFloat("pileupJetId:fullDiscriminant"));
	    float svtxmass=j.userFloat("vtxMass");
	    histContainer_["jetsecvtxmass"]->Fill(svtxmass);
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
  float dphi_met_mu = deltaPhi(muphi, metphi); // use the function to restrict to the 0,pi range
  float mt=sqrt(2*mupt*metpt*(1-cos(dphi_met_mu)));

  //
  // FINAL SELECTION PLOTS
  //
  if(selectedJets.size()>=3) 
    {
      histContainer_["cutflow"]->Fill(3);
      histContainer_["mt_3"]->Fill(mt);
    }
  if(selectedJets.size()>=4) 
    {
      histContainer_["mt_4"]->Fill(mt);
      histContainer_["cutflow"]->Fill(4);
      histContainer_["ncsvmjets"]->Fill(nCSVMtags);
      histContainer_["nvertices"]->Fill(vertices->size());
      if(nCSVMtags>=1) 
	{
	  histContainer_["cutflow"]->Fill(5);
	  histContainer_["mt_41b"] ->Fill(mt);
	}
      if(nCSVMtags>=2) 
	{
	  histContainer_["cutflow"]->Fill(6);
	  histContainer_["mt_42b"] ->Fill(mt);
	  histContainer_["metpt"]->Fill(metpt);
	  histContainer_["metphi"] ->Fill(metphi);
	  histContainer_["dphimetmuon"]->Fill(dphi_met_mu );
	}
    }
  
  
}


// ------------ method called once each job just before starting event loop  ------------
void 
MiniAnalyzer::beginJob()
{
  edm::Service<TFileService> fs;
  histContainer_["cutflow"] = fs->make<TH1F>("cutflow",    ";Selection cut;Events", 7, 0., 7.); 
  histContainer_["cutflow"] -> GetXaxis()->SetBinLabel(1,"reco");
  histContainer_["cutflow"] -> GetXaxis()->SetBinLabel(2,"=1 #mu");
  histContainer_["cutflow"] -> GetXaxis()->SetBinLabel(3,"=0 e");
  histContainer_["cutflow"] -> GetXaxis()->SetBinLabel(4,"#geq3 jets");
  histContainer_["cutflow"] -> GetXaxis()->SetBinLabel(5,"#geq4 jets");
  histContainer_["cutflow"] -> GetXaxis()->SetBinLabel(6,"#geq1 b-tag");
  histContainer_["cutflow"] -> GetXaxis()->SetBinLabel(7,"#geq2 b-tags");

  histContainer_["nvertices"] = fs->make<TH1F>("nvertices",    ";# vertices;Events", 100, 0., 100.); 

  histContainer_["nrecomuons"]  = fs->make<TH1F>("nrecomuons",   ";# reconstructed muons; Events",             10, 0., 10.);
  histContainer_["muonpt"]      = fs->make<TH1F>("muonpt",       ";Transverse momentum [GeV];# muons",         100, 0., 300.);
  histContainer_["muoneta"]     = fs->make<TH1F>("muoneta",      ";Pseudo-rapidity;#muons ",                   100, 0, 3.);
  histContainer_["nselmuons"]   = fs->make<TH1F>("nselmuons",    ";# reconstructed muons;Events",              5, 0.,5.);

  histContainer_["nrecoelectrons"]    = fs->make<TH1F>("nrecoelectrons",  ";# reconstructed electrons;Events",         5, 0., 5.);
  histContainer_["electronpt"]        = fs->make<TH1F>("electronpt",      ";Transverse momentum [GeV];# electrons",    100, 0., 300.);
  histContainer_["electroneta"]       = fs->make<TH1F>("electroneta",     ";Pseudo-rapidity;#electrons",               100, 0., 3.);
  histContainer_["nselelectrons"]     = fs->make<TH1F>("nselelectrons",   ";# selected electrons;Events", 5, 0., 5.);

  histContainer_["nrecojets"]      = fs->make<TH1F>("nrecojets",   ";#reconstructed jets;Events", 50, 0., 50.);
  histContainer_["jetpt"]          = fs->make<TH1F>("jetpt",       ";Transverse momentum [GeV];# jets", 100, 0., 300.);
  histContainer_["jeteta"]         = fs->make<TH1F>("jeteta",      ";Pseudo-rapidity;# jets", 100, 0., 3.);
  histContainer_["jetcsv"]         = fs->make<TH1F>("jetcsv",      ";Combined secondary vertes;# jets", 100, -1.2, 1.2);
  histContainer_["jetpileupid"]    = fs->make<TH1F>("jetpileupid", ";Pileup jet id;#jets", 100, -1.2, 1.2);
  histContainer_["jetsecvtxmass"]  = fs->make<TH1F>("jetvtxMass", ";Secondary vertex mass [GeV];#jets", 100, 0., 6.);
  histContainer_["nseljets"]       = fs->make<TH1F>("nseljets",    ";#selected jets;Events", 6, 3., 10.);
  histContainer_["ncsvmjets"]      = fs->make<TH1F>("ncsvmjets",    ";b-tagged jets (CSVM);Events", 10, 0., 10.);
  
  histContainer_["metpt"] = fs->make<TH1F>("metpt",    ";Missing transverse energy [GeV];Events", 100, 0., 300.);
  histContainer_["metphi"] = fs->make<TH1F>("metphi",    ";Missing transverse energy #phi [rad];Events", 50, -3.2, 3.2);
  histContainer_["dphimetmuon"] = fs->make<TH1F>("dphimetmuon",    ";#Delta#phi(MET,#mu) [rad];Events", 50, -3.2, 3.2);
 
  histContainer_["mt_3"] = fs->make<TH1F>("mt_3",    ";Transverse mass [GeV]; Events", 100, 0., 200.);
  histContainer_["mt_4"] = fs->make<TH1F>("mt_4",    ";Transverse mass [GeV]; Events", 100, 0., 200.);
  histContainer_["mt_41b"] = fs->make<TH1F>("mt_41b",    ";Transverse mass [GeV]; Events", 100, 0., 200.);
  histContainer_["mt_42b"] = fs->make<TH1F>("mt_42b",    ";Transverse mass [GeV]; Events", 100, 0., 200.);

  //instruct ROOT to compute the uncertainty from the square root of weights
  //http://root.cern.ch/root/html/TH1.html#TH1:Sumw2
  for(std::unordered_map<std::string,TH1F*>::iterator it=histContainer_.begin();   it!=histContainer_.end();   it++) it->second->Sumw2();
  for(std::unordered_map<std::string,TH2F*>::iterator it=histContainer2d_.begin(); it!=histContainer2d_.end(); it++) it->second->Sumw2();
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
