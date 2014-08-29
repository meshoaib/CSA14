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
  csvbjetcutflow->GetXaxis()->SetBinLabel(1,"#equiv 0b-tag");
  csvbjetcutflow->GetXaxis()->SetBinLabel(2,"#equiv 1b-tag");
  csvbjetcutflow->GetXaxis()->SetBinLabel(3,"#geq 2b-tags");
  csvbjetcutflow->GetXaxis()->SetBinLabel(4,"#equiv 0b-tags");
  csvbjetcutflow->GetXaxis()->SetBinLabel(5,"#equiv 1b-tags");
  csvbjetcutflow->GetXaxis()->SetBinLabel(6,"#geq 2b-tags");

  TH1F *ssvbjetcutflow = new TH1F("ssvbjetcutflow",";Cut;Events" ,6,0.,6.);
  ssvbjetcutflow->GetXaxis()->SetBinLabel(1,"#equiv 0b-tag");
  ssvbjetcutflow->GetXaxis()->SetBinLabel(2,"#equiv 1b-tag");
  ssvbjetcutflow->GetXaxis()->SetBinLabel(3,"#geq 2b-tags");
  ssvbjetcutflow->GetXaxis()->SetBinLabel(4,"#equiv 0b-tags");
  ssvbjetcutflow->GetXaxis()->SetBinLabel(5,"#equiv 1b-tags");
  ssvbjetcutflow->GetXaxis()->SetBinLabel(6,"#geq 2b-tags");

  TH1F *disvtx_3j = new TH1F("disvtx_3j",";lxy;Events" ,100,0.,10.);
  TH1F *lxyz_sig_3j = new TH1F("lxyz_sig_3j",";lxyz_sig;Events" ,100,0.,10.);
  TH1F *vertexmass_3j = new TH1F("vertexmass_3j",";vertexmass;Events" ,100,0.,10.);
  TH1F *jetpt_3j = new TH1F("jetpt_3j",";pt;Events" ,100,0.,300.);
  TH1F *jeteta_3j = new TH1F("jeteta_3j",";eta;Events" ,100,0.,3.);

  TH1F *disvtx_4j = new TH1F("disvtx_4j",";lxy;Events" ,100,0.,10.);
  TH1F *lxyz_sig_4j = new TH1F("lxyz_sig_4j",";lxyz_sig;Events" ,100,0.,10.);
  TH1F *vertexmass_4j = new TH1F("vertexmass_4j",";vertexmass;Events" ,100,0.,10.);
  TH1F *jetpt_4j = new TH1F("jetpt_4j",";pt;Events" ,100,0.,300.);
  TH1F *jeteta_4j = new TH1F("jeteta_4j",";eta;Events" ,100,0.,3.);

  TH1F *disvtx_3j_leading = new TH1F("disvtx_3j_leading",";lxy;Events" ,100,0.,10.);
  TH1F *lxyz_sig_3j_leading = new TH1F("lxyz_sig_3j_leading",";lxyz_sig;Events" ,100,0.,10.);
  TH1F *vertexmass_3j_leading = new TH1F("vertexmass_3j_leading",";vertexmass;Events" ,100,0.,10.);
  TH1F *jetpt_3j_leading= new TH1F("jetpt_3j_leading",";pt;Events" ,100,0.,300.);
  TH1F *jeteta_3j_leading= new TH1F("jeteta_3j_leading",";eta;Events" ,100,0.,3.);

  TH1F *disvtx_4j_leading = new TH1F("disvtx_4j_leading",";lxy;Events" ,100,0.,10.);
  TH1F *lxyz_sig_4j_leading = new TH1F("lxyz_sig_4j_leading",";lxyz_sig;Events" ,100,0.,10.);
  TH1F *vertexmass_4j_leading = new TH1F("vertexmass_4j_leading",";vertexmass;Events" ,100,0.,10.);
  TH1F *jetpt_4j_leading= new TH1F("jetpt_4j_leading",";pt;Events" ,100,0.,300.);
  TH1F *jeteta_4j_leading= new TH1F("jeteta_4j_leading",";eta;Events" ,100,0.,3.);

  TH1F *disvtx_3j_nextleading = new TH1F("disvtx_3j_nextleading",";lxy;Events" ,100,0.,10.);
  TH1F *lxyz_sig_3j_nextleading = new TH1F("lxyz_sig_3j_nextleading",";lxyz_sig;Events" ,100,0.,10.);
  TH1F *vertexmass_3j_nextleading = new TH1F("vertexmass_3j_nextleading",";vertexmass;Events" ,100,0.,10.);
  TH1F *jetpt_3j_nextleading= new TH1F("jetpt_3j_nextleading",";pt;Events" ,100,0.,300.);
  TH1F *jeteta_3j_nextleading= new TH1F("jeteta_3j_nextleading",";eta;Events" ,100,0.,3.);

  TH1F *disvtx_4j_nextleading  = new TH1F("disvtx_4j_nextleading",";lxy;Events" ,100,0.,10.);
  TH1F *lxyz_sig_4j_nextleading = new TH1F("lxyz_sig_4j_nextleading",";lxyz_sig;Events" ,100,0.,10.);
  TH1F *vertexmass_4j_nextleading = new TH1F("vertexmass_4j_nextleading",";vertexmass;Events" ,100,0.,10.);
  TH1F *jetpt_4j_nextleading= new TH1F("jetpt_4j_nextleading",";pt;Events" ,100,0.,300.);
  TH1F *jeteta_4j_nextleading= new TH1F("jeteta_4j_nextleading",";eta;Events" ,100,0.,3.);

 
  //read tree from file
  MiniEvent_t ev;
  TFile *f = TFile::Open(filename);
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
      				    cutflow->Fill(0);
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

		if(v==0) {
                         if(nJets >=3){
                          disvtx_3j->Fill(lxyz);
                          lxyz_sig_3j->Fill(lxyz_sig);
                          vertexmass_3j->Fill(vtxmass);
                          jetpt_3j->Fill(pt);
                          jeteta_3j->Fill(fabs(eta));}

                          if(nJets >=4){
                          disvtx_4j->Fill(lxyz);
                          lxyz_sig_4j->Fill(lxyz_sig);
                          vertexmass_4j->Fill(vtxmass);
                          jetpt_4j->Fill(pt);
                          jeteta_4j->Fill(fabs(eta));}

                         }


		if(v==1) {
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

		 if(v==2) {
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

  disvtx_3j->Write();
  lxyz_sig_3j->Write();
  vertexmass_3j->Write();
  jetpt_3j->Write();
  jeteta_3j->Write();

  disvtx_4j->Write();
  lxyz_sig_4j->Write();
  vertexmass_4j->Write();
  jetpt_4j->Write();
  jeteta_4j->Write();

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

