#include <TFile.h>
#include <TROOT.h>
#include <TH1.h>
#include "UserCode/TopAnalysis/interface/MiniEvent.h"
 
#include <vector>
#include <iostream>
#include <algorithm>

bool sortBySignificance(std::pair<int,std::pair<float,float> > a,
			std::pair<int,std::pair<float,float> > b)
{
  if( a.second.first>0 || b.second.first>0 ) return (a.second.first>b.second.first);
  return (a.second.second>b.second.second);
}

void ReadTree(TString filename,TString output)
{
  gROOT->Reset();

  TH1F *cutflow = new TH1F("cutflow",";Cut;Events" ,6,0.,6.);
  cutflow->GetXaxis()->SetBinLabel(1,"preselected");
  cutflow->GetXaxis()->SetBinLabel(2,"#geq 2j");
  cutflow->GetXaxis()->SetBinLabel(3,"#geq 3j");
  cutflow->GetXaxis()->SetBinLabel(4,"#geq 4j");
  cutflow->GetXaxis()->SetBinLabel(5,"#geq 1b-tag");
  cutflow->GetXaxis()->SetBinLabel(6,"#geq 2b-tags");

  TH1F *csvbjetcutflow = new TH1F("csvbjetcutflow",";Cut;Events" ,9,0.,9.);
  csvbjetcutflow->GetXaxis()->SetBinLabel(1,"2j,=0b");
  csvbjetcutflow->GetXaxis()->SetBinLabel(2,"2j,=1b");
  csvbjetcutflow->GetXaxis()->SetBinLabel(3,"2j,#geq2b");
  csvbjetcutflow->GetXaxis()->SetBinLabel(4,"3j,=0b");
  csvbjetcutflow->GetXaxis()->SetBinLabel(5,"3j,=1b");
  csvbjetcutflow->GetXaxis()->SetBinLabel(6,"3j,#geq2b");
  csvbjetcutflow->GetXaxis()->SetBinLabel(7,"4j,=0b");
  csvbjetcutflow->GetXaxis()->SetBinLabel(8,"4j,=1b");
  csvbjetcutflow->GetXaxis()->SetBinLabel(9,"4j,#geq2b");
  TH1F *ssvbjetcutflow = (TH1F *)csvbjetcutflow->Clone("ssvbjetcutflow");

  TH1F *disvtx_2j_leading     = new TH1F("disvtx_2j_leading",";lxy;Events" ,25,0.,10.);
  TH1F *disvtx_2j_nextleading = (TH1F *)disvtx_2j_leading->Clone("disvtx_2j_nextleading");
  TH1F *disvtx_3j_leading     = (TH1F *)disvtx_2j_leading->Clone("disvtx_3j_leading");
  TH1F *disvtx_3j_nextleading = (TH1F *)disvtx_2j_leading->Clone("disvtx_3j_nextleading");
  TH1F *disvtx_4j_leading     = (TH1F *)disvtx_2j_leading->Clone("disvtx_4j_leading");
  TH1F *disvtx_4j_nextleading = (TH1F *)disvtx_2j_leading->Clone("disvtx_4j_nextleading");

  TH1F *lxyz_sig_2j_leading     = new TH1F("lxyz_sig_2j_leading",";lxyz_sig;Events" ,25,0.,20.);
  TH1F *lxyz_sig_2j_nextleading = (TH1F *)lxyz_sig_2j_leading->Clone("lxyz_sig_2j_nextleading");
  TH1F *lxyz_sig_3j_leading     = (TH1F *)lxyz_sig_2j_leading->Clone("lxyz_sig_3j_leading");
  TH1F *lxyz_sig_3j_nextleading = (TH1F *)lxyz_sig_2j_leading->Clone("lxyz_sig_3j_nextleading");
  TH1F *lxyz_sig_4j_leading     = (TH1F *)lxyz_sig_2j_leading->Clone("lxyz_sig_4j_leading");
  TH1F *lxyz_sig_4j_nextleading = (TH1F *)lxyz_sig_2j_leading->Clone("lxyz_sig_4j_nextleading");

  TH1F *vertexmass_2j_leading     = new TH1F("vertexmass_2j_leading",";vertexmass;Events" ,25,0.,6.);
  TH1F *vertexmass_2j_nextleading = (TH1F *)vertexmass_2j_leading->Clone("vertexmass_2j_nextleading");
  TH1F *vertexmass_3j_leading     = (TH1F *)vertexmass_2j_leading->Clone("vertexmass_3j_leading");
  TH1F *vertexmass_3j_nextleading = (TH1F *)vertexmass_2j_leading->Clone("vertexmass_3j_nextleading");
  TH1F *vertexmass_4j_leading     = (TH1F *)vertexmass_2j_leading->Clone("vertexmass_4j_leading");
  TH1F *vertexmass_4j_nextleading = (TH1F *)vertexmass_2j_leading->Clone("vertexmass_4j_nextleading");

  TH1F *jetpt_2j_leading     = new TH1F("jetpt_2j_leading",";pt;Events" ,25,0.,300.);
  TH1F *jetpt_2j_nextleading = (TH1F *)jetpt_2j_leading->Clone("jetpt_2j_nextleading");
  TH1F *jetpt_3j_leading     = (TH1F *)jetpt_2j_leading->Clone("jetpt_3j_leading");
  TH1F *jetpt_3j_nextleading = (TH1F *)jetpt_2j_leading->Clone("jetpt_3j_nextleading");
  TH1F *jetpt_4j_leading     = (TH1F *)jetpt_2j_leading->Clone("jetpt_4j_leading");
  TH1F *jetpt_4j_nextleading = (TH1F *)jetpt_2j_leading->Clone("jetpt_4j_nextleading");

  TH1F *jeteta_2j_leading     = new TH1F("jeteta_2j_leading",";eta;Events" ,25,0.,3.);
  TH1F *jeteta_2j_nextleading = (TH1F *)jeteta_2j_leading->Clone("jeteta_2j_nextleading");
  TH1F *jeteta_3j_leading     = (TH1F *)jeteta_2j_leading->Clone("jeteta_3j_leading");
  TH1F *jeteta_3j_nextleading = (TH1F *)jeteta_2j_leading->Clone("jeteta_3j_nextleading");
  TH1F *jeteta_4j_leading     = (TH1F *)jeteta_2j_leading->Clone("jeteta_4j_leading");
  TH1F *jeteta_4j_nextleading = (TH1F *)jeteta_2j_leading->Clone("jeteta_4j_nextleading");
  

 
  //read tree from file
  MiniEvent_t ev;
  TFile *f = TFile::Open(filename);

  //get the original number of events in the dataset
  TH1F *origCutFlow=(TH1F *)f->Get("demo/cutflow");
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
      bool useChPt(ev.j_chpt);
//      int minVtx, maxVtx;
      int Nvertices(ev.nvtx); 	
      std::vector< std::pair<int,std::pair<float,float> > > vlxyz_sig;
//      std::vector< std::pair<int,std::pair<int,float> > > vlxyz_sig;
      for (int k=0; k<ev.nj;k++)
	{
	  //check pt and eta of this jet
	  float pt  = ev.j_pt[k];
	  float chPt  = ev.j_chpt[k];
	  float eta = ev.j_eta[k];
//	  int Nvertices = ev.nvtx;
	  float csv = ev.j_csv[k];
	  float lxyz=ev.j_vtx3DVal[k];
	  float lxyz_sig= ev.j_vtx3DSig[k];

	  if (useChPt) {
	   if(chPt > 15 && abs(eta) < 2.5 ) 
 	   {
              nJets++;
              if (csv>0.679) nCSVMtags++;
              if(lxyz>0) nSSVMtags++;
              {
                std::pair<float,float> jetKinematics(lxyz_sig,pt);
                vlxyz_sig.push_back( std::pair<int,std::pair<float,float> >(k,jetKinematics) );
              }
            }
            else
	     if (pt > 30 && abs(eta) < 2.5) 
	    {
	      nJets++;
	      if (csv>0.679) nCSVMtags++;
              if(lxyz>0) nSSVMtags++;
	      {
		std::pair<float,float> jetKinematics(lxyz_sig,pt);
		vlxyz_sig.push_back( std::pair<int,std::pair<float,float> >(k,jetKinematics) );
	      }
	    }
		}
	}

      std::sort(vlxyz_sig.begin(),vlxyz_sig.end(),sortBySignificance);
	  
      
      //fill histos for CSVM
      if(Nvertices >= 15 && Nvertices <= 30) 
      if(nJets>=2)                  cutflow->Fill(1);
      if(nJets>=3)                  cutflow->Fill(2);
      if(nJets>=4)                  cutflow->Fill(3);
      if(nJets>=4 && nCSVMtags >=1) cutflow->Fill(4);
      if(nJets>=4 && nCSVMtags >=2) cutflow->Fill(5);
      if(nJets==2 && nCSVMtags ==0) csvbjetcutflow->Fill(0);	
      if(nJets==2 && nCSVMtags ==1) csvbjetcutflow->Fill(1);	
      if(nJets==2 && nCSVMtags ==2) csvbjetcutflow->Fill(2);	
      if(nJets==3 && nCSVMtags ==0) csvbjetcutflow->Fill(3);	
      if(nJets==3 && nCSVMtags ==1) csvbjetcutflow->Fill(4);	
      if(nJets==3 && nCSVMtags >=2) csvbjetcutflow->Fill(5);	
      if(nJets>=4 && nCSVMtags ==0) csvbjetcutflow->Fill(6);
      if(nJets>=4 && nCSVMtags ==1) csvbjetcutflow->Fill(7);
      if(nJets>=4 && nCSVMtags >=2) csvbjetcutflow->Fill(8);
      
      //fill histos for SSV
      if(nJets==2 && nSSVMtags ==0) ssvbjetcutflow->Fill(0);
      if(nJets==2 && nSSVMtags ==1) ssvbjetcutflow->Fill(1);
      if(nJets==2 && nSSVMtags ==2) ssvbjetcutflow->Fill(2);
      if(nJets==3 && nSSVMtags ==0) ssvbjetcutflow->Fill(3);
      if(nJets==3 && nSSVMtags ==1) ssvbjetcutflow->Fill(4);
      if(nJets==3 && nSSVMtags >=2) ssvbjetcutflow->Fill(5);
      if(nJets>=4 && nSSVMtags ==0) ssvbjetcutflow->Fill(6);
      if(nJets>=4 && nSSVMtags ==1) ssvbjetcutflow->Fill(7);
      if(nJets>=4 && nSSVMtags >=2) ssvbjetcutflow->Fill(8);

	else continue;                                                                                                

      for(size_t v=0; v<vlxyz_sig.size(); v++) 
	{
	  if(v>1) break;

	  int k=vlxyz_sig[v].first;
	  float pt  = ev.j_pt[k];
//	  float pt  = ev.j_chpt[k];
	  float eta = ev.j_eta[k];
	  float lxyz=ev.j_vtx3DVal[k];
	  float lxyz_sig= ev.j_vtx3DSig[k];
	  float vtxmass = ev.j_vtxmass[k]; 

	  //most significantly displaced vertex
	  if(v==0) {
	    if(nJets >=2){
	      if(lxyz>0){
		disvtx_2j_leading->Fill(lxyz);
		lxyz_sig_2j_leading->Fill(lxyz_sig);
		vertexmass_2j_leading->Fill(vtxmass);
	      }
	      jetpt_2j_leading->Fill(pt);
	      jeteta_2j_leading->Fill(fabs(eta));}
	
	     if(nJets >=3){
              if(lxyz>0){
                disvtx_3j_leading->Fill(lxyz);
                lxyz_sig_3j_leading->Fill(lxyz_sig);
                vertexmass_3j_leading->Fill(vtxmass);
              }
              jetpt_3j_leading->Fill(pt);
              jeteta_3j_leading->Fill(fabs(eta));}

	    if(nJets >=4){
	      if(lxyz>0){
		disvtx_4j_leading->Fill(lxyz);
		lxyz_sig_4j_leading->Fill(lxyz_sig);
		vertexmass_4j_leading->Fill(vtxmass);
	      }
	      jetpt_4j_leading->Fill(pt);
	      jeteta_4j_leading->Fill(fabs(eta));}
	  }

	  //second most significantly displaced vertex
	  if(v==1) {
	    if(nJets >=2){
	      if(lxyz>0){
		disvtx_2j_nextleading->Fill(lxyz);
		lxyz_sig_2j_nextleading->Fill(lxyz_sig);
		vertexmass_2j_nextleading->Fill(vtxmass);
	      }
	      jetpt_2j_nextleading->Fill(pt);
	      jeteta_2j_nextleading->Fill(fabs(eta));}

	    if(nJets >=3){
              if(lxyz>0){
                disvtx_3j_nextleading->Fill(lxyz);
                lxyz_sig_3j_nextleading->Fill(lxyz_sig);
                vertexmass_3j_nextleading->Fill(vtxmass);
              }
              jetpt_3j_nextleading->Fill(pt);
              jeteta_3j_nextleading->Fill(fabs(eta));}

	    if(nJets >=4){
	      if(lxyz>0){
		disvtx_4j_nextleading->Fill(lxyz);
		lxyz_sig_4j_nextleading->Fill(lxyz_sig);
		vertexmass_4j_nextleading->Fill(vtxmass);
	      }
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

  disvtx_2j_leading->Write();
  lxyz_sig_2j_leading->Write();
  vertexmass_2j_leading->Write();
  jetpt_2j_leading->Write();
  jeteta_2j_leading->Write();

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

  disvtx_2j_nextleading->Write();
  lxyz_sig_2j_leading->Write();
  vertexmass_2j_nextleading->Write();
  jetpt_2j_nextleading->Write();
  jeteta_2j_nextleading->Write();

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


void RunOverSamples(TString output)
{
  TString files[]={
    "miniAOD_DY.root",
    "miniAOD_MC_v2.root",
    "miniAOD_QCD_1000_1400.root",
    "miniAOD_QCD_120_170.root",
    "miniAOD_QCD_1400_1800.root",
    "miniAOD_QCD_170_300.root",
    "miniAOD_QCD_1800_2400.root",
    "miniAOD_QCD_2400_3200.root",
    "miniAOD_QCD_300_470.root",
    "miniAOD_QCD_3200.root",
    "miniAOD_QCD_470_600.root",
    "miniAOD_QCD_600_800.root",
    "miniAOD_QCD_800_1000.root",
    "miniAOD_QCD_80_120.root",
    "miniAOD_TBarToLeptons_t_channel.root",
    "miniAOD_TToLeptons_s_channel.root",
    "miniAOD_T_tW_channel.root",
    "miniAOD_Tbar_tW_channel.root",
    "miniAOD_wjets.root"
  };
  for(size_t i=0; i<sizeof(files)/sizeof(TString); i++)
    {
      ReadTree(files[i],output);
    }
  

}
