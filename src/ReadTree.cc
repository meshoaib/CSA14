#include <TFile.h>
#include <TROOT.h>
#include <TH1.h>
#include "UserCode/TopAnalysis/interface/MiniEvent.h"
 
#include <vector>
#include <iostream>
#include <algorithm>

bool sortBySignificance(std::pair<int,float> a,std::pair<int,float> b)
{
  return (a.second>b.second);
}

void ReadTree(TString filename,TString output)
{
  gROOT->Reset();

  TH1F *cutflow = new TH1F("cutflow",";Cut;Events" ,5,0.,5.);
  cutflow->GetXaxis()->SetBinLabel(1,"preselected");
  cutflow->GetXaxis()->SetBinLabel(2,"#geq 3j");
  cutflow->GetXaxis()->SetBinLabel(3,"#geq 4j");
  cutflow->GetXaxis()->SetBinLabel(4,"#geq 1b-tag");
  cutflow->GetXaxis()->SetBinLabel(5,"#geq 2b-tags");

  TH1F *csvbjetcutflow = new TH1F("csvbjetcutflow",";Cut;Events" ,6,0.,6.);
  csvbjetcutflow->GetXaxis()->SetBinLabel(1,"3j,=0b");
  csvbjetcutflow->GetXaxis()->SetBinLabel(2,"3j,=1b");
  csvbjetcutflow->GetXaxis()->SetBinLabel(3,"3j,#geq2b");
  csvbjetcutflow->GetXaxis()->SetBinLabel(4,"4j,=0b");
  csvbjetcutflow->GetXaxis()->SetBinLabel(5,"4j,=1b");
  csvbjetcutflow->GetXaxis()->SetBinLabel(6,"4j,#geq2b");
  TH1F *ssvbjetcutflow = (TH1F *)csvbjetcutflow->Clone("ssvbjetcutflow");

  TH1F *disvtx_3j_leading     = new TH1F("disvtx_3j_leading",";lxy;Events" ,100,0.,10.);
  TH1F *disvtx_3j_nextleading = (TH1F *)disvtx_3j_leading->Clone("disvtx_4j_nextleading");
  TH1F *disvtx_4j_leading     = (TH1F *)disvtx_3j_leading->Clone("disvtx_4j_leading");
  TH1F *disvtx_4j_nextleading = (TH1F *)disvtx_3j_leading->Clone("disvtx_4j_nextleading");

  TH1F *lxyz_sig_3j_leading     = new TH1F("lxyz_sig_3j_leading",";lxyz_sig;Events" ,100,0.,10.);
  TH1F *lxyz_sig_3j_nextleading = (TH1F *)lxyz_sig_3j_leading->Clone("lxyz_sig_3j_nextleading");
  TH1F *lxyz_sig_4j_leading     = (TH1F *)lxyz_sig_3j_leading->Clone("lxyz_sig_4j_leading");
  TH1F *lxyz_sig_4j_nextleading = (TH1F *)lxyz_sig_3j_leading->Clone("lxyz_sig_4j_nextleading");

  TH1F *vertexmass_3j_leading     = new TH1F("vertexmass_3j_leading",";vertexmass;Events" ,100,0.,10.);
  TH1F *vertexmass_3j_nextleading = (TH1F *)vertexmass_3j_leading->Clone("vertexmass_3j_nextleading");
  TH1F *vertexmass_4j_leading     = (TH1F *)vertexmass_3j_leading->Clone("vertexmass_4j_leading");
  TH1F *vertexmass_4j_nextleading = (TH1F *)vertexmass_3j_leading->Clone("vertexmass_3j_nextleading");

  TH1F *jetpt_3j_leading     = new TH1F("jetpt_3j_leading",";pt;Events" ,100,0.,300.);
  TH1F *jetpt_3j_nextleading = (TH1F *)jetpt_3j_leading->Clone("jetpt_3j_nextleading");
  TH1F *jetpt_4j_leading     = (TH1F *)jetpt_3j_leading->Clone("jetpt_4j_leading");
  TH1F *jetpt_4j_nextleading = (TH1F *)jetpt_3j_leading->Clone("jetpt_4j_nextleading");

  TH1F *jeteta_3j_leading     = new TH1F("jeteta_3j_leading",";eta;Events" ,100,0.,3.);
  TH1F *jeteta_3j_nextleading = (TH1F *) jeteta_3j_leading->Clone("jeteta_3j_nextleading");
  TH1F *jeteta_4j_leading     = (TH1F *) jeteta_3j_leading->Clone("jeteta_4j_leading");
  TH1F *jeteta_4j_nextleading = (TH1F *) jeteta_3j_leading->Clone("jeteta_4j_nextleading");
  

 
  //read tree from file
  MiniEvent_t ev;
  TFile *f = TFile::Open(filename);

  //get the original number of events in the dataset
  TH1F *origCutFlow=(TH1F *)f->Get("demo/cuflow");
  if(origCutFlow) {
    Float_t origEvents=origCutFlow->GetBinContent(1);
    cutflow->SetBinContent(1,origEvents);
  }

  //get the tree
  TTree *t = (TTree*)f->Get("demo/AnaTree");
  attachToMiniEventTree(t,ev);

  //fill histograms, loop over all entries
  Int_t nentries = (Int_t)t->GetEntriesFast();
  for (Int_t i=0;i<nentries;i++)
    {
      t->GetEntry(i);

      //select jets
      uint32_t nJets(0), nCSVMtags(0), nSSVMtags(0);
      std::vector< std::pair<int,float> > vlxyz_sig;
      for (int k=0; k<ev.nj;k++)
	{
	  //check pt and eta of this jet
	  float pt  = ev.j_pt[k];
	  float eta = ev.j_eta[k];
	  float csv = ev.j_csv[k];
	  float lxyz=ev.j_vtx3DVal[k];
	  float lxyz_sig= ev.j_vtx3DSig[k];

	  if (pt > 30 && abs(eta) < 2.5)        
	    {
	      nJets++;
	      if (csv>0.679) nCSVMtags++;
              if(lxyz>0) nSSVMtags++;
	      {
		vlxyz_sig.push_back( std::pair<int,float>(k,lxyz_sig) );
	      }
	    }
	}

      std::sort(vlxyz_sig.begin(),vlxyz_sig.end(),sortBySignificance);
	  
      
      //fill histos for CSVM
      if(nJets>=3)                  cutflow->Fill(1);
      if(nJets>=4)                  cutflow->Fill(2);
      if(nJets>=4 && nCSVMtags >=1) cutflow->Fill(3);
      if(nJets>=4 && nCSVMtags >=2) cutflow->Fill(4);
      if(nJets==3 && nCSVMtags ==0) csvbjetcutflow->Fill(0);	
      if(nJets==3 && nCSVMtags ==1) csvbjetcutflow->Fill(1);	
      if(nJets==3 && nCSVMtags >=2) csvbjetcutflow->Fill(2);	
      if(nJets>=4 && nCSVMtags ==0) csvbjetcutflow->Fill(3);
      if(nJets>=4 && nCSVMtags ==1) csvbjetcutflow->Fill(4);
      if(nJets>=4 && nCSVMtags >=2) csvbjetcutflow->Fill(5);
      
      //fill histos for SSV
      if(nJets==3 && nSSVMtags ==0) ssvbjetcutflow->Fill(0);
      if(nJets==3 && nSSVMtags ==1) ssvbjetcutflow->Fill(1);
      if(nJets==3 && nSSVMtags >=2) ssvbjetcutflow->Fill(2);
      if(nJets>=4 && nSSVMtags ==0) ssvbjetcutflow->Fill(3);
      if(nJets>=4 && nSSVMtags ==1) ssvbjetcutflow->Fill(4);
      if(nJets>=4 && nSSVMtags >=2) ssvbjetcutflow->Fill(5);
                                                                                                
      for(size_t v=0; v<vlxyz_sig.size(); v++) 
	{
	  if(v>1) break;

	  int k=vlxyz_sig[v].first;
	  float pt  = ev.j_pt[k];
	  float eta = ev.j_eta[k];
	  float lxyz=ev.j_vtx3DVal[k];
	  float lxyz_sig= ev.j_vtx3DSig[k];
	  float vtxmass = ev.j_vtxmass[k]; 

	  //most significantly displaced vertex
	  if(v==0) {
	    if(nJets >=3){
	      disvtx_3j_leading->Fill(lxyz);
	      lxyz_sig_3j_leading->Fill(lxyz_sig);
	      vertexmass_3j_leading->Fill(vtxmass);
	      jetpt_3j_leading->Fill(pt);
	      jeteta_3j_leading->Fill(fabs(eta));}

	    if(nJets >=4){
	      disvtx_4j_leading->Fill(lxyz);
	      lxyz_sig_4j_leading->Fill(lxyz_sig);
	      vertexmass_4j_leading->Fill(vtxmass);
	      jetpt_4j_leading->Fill(pt);
	      jeteta_4j_leading->Fill(fabs(eta));}
	  }

	  //second most significantly displaced vertex
	  if(v==1) {
	    if(nJets >=3){
	      disvtx_3j_nextleading->Fill(lxyz);
	      lxyz_sig_3j_nextleading->Fill(lxyz_sig);
	      vertexmass_3j_nextleading->Fill(vtxmass);
	      jetpt_3j_nextleading->Fill(pt);
	      jeteta_3j_nextleading->Fill(fabs(eta));}

	    if(nJets >=4){
	      disvtx_4j_nextleading->Fill(lxyz);
	      lxyz_sig_4j_nextleading->Fill(lxyz_sig);
	      vertexmass_4j_nextleading->Fill(vtxmass);
	      jetpt_4j_nextleading->Fill(pt);
	      jeteta_4j_nextleading->Fill(fabs(eta));}
	  }	
	}
    }
    

    
  
  //close file
  
  f->Close();
  
  //open output file
  TFile *fOut=TFile::Open(output+"/"+filename,"RECREATE");
  cutflow->Write();
  csvbjetcutflow->Write();
  ssvbjetcutflow->Write();

  disvtx_3j_leading->Write();
  lxyz_sig_3j_leading->Write();
  vertexmass_3j_leading->Write();
  jetpt_3j_leading->Write();
  jeteta_3j_leading->Write();
 
  disvtx_4j_leading->Write();
  lxyz_sig_4j_leading->Write();
  vertexmass_4j_leading->Write();
  jetpt_4j_leading->Write();
  jeteta_4j_leading->Write();

  disvtx_3j_nextleading->Write();
  lxyz_sig_3j_leading->Write();
  vertexmass_3j_nextleading->Write();
  jetpt_3j_nextleading->Write();
  jeteta_3j_nextleading->Write();

  disvtx_4j_nextleading->Write();
  lxyz_sig_4j_nextleading->Write();
  vertexmass_4j_nextleading->Write();
  jetpt_4j_nextleading->Write();
  jeteta_4j_nextleading->Write();
 
 
  fOut->Close();
}
